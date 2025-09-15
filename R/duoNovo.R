#' Identify de novo variants from duo (proband - single parent) long read sequencing
#'
#' The `duoNovo` function uses phased variant calls from duo long-read sequencing to classify candidate variants (heterozygous in the proband; absent in the parent).
#' Variants are classifed as de novo, present on the haplotype inherited from the non-sequenced parent, or uncertain variants.
#'
#' @param LRS_phased_vcf_file_path File path to the vcf containing phased variant calls from long-read sequencing of the duo.
#' @param depth_cutoff A numeric value specifying the minimum sequencing depth for variants to be included in the analysis. The same cutoff applies to both proband and parent variants, for both LRS and SRS (if used).
#' @param GQ_cutoff A numeric value specifying the minimum GQ (genotype quality) for variants to be included in the analysis. The same cutoff applies to both proband and parent variants, for both LRS and SRS (if used).
#' @param proband_column_identifier A character corresponding to an identifier for the proband column in the metadata matrices. Should be the same for both the LRS vcf and (if used) the SRS vcf.
#' @param PS_width_cutoff A numeric value specifying the minimum width of the haplotype blocks included in the analysis.
#' @param boundary_cutoff A numeric value indicating the minimum distance from a haplotype block boundary (either start or end coordinate) for candidate variants to be analyzed.
#' @param IBD_distance_cutoff A numeric value specifying the minimum Hamming distance required for inferring
#'   the non-IBD status of a haplotype pair.
#' @param non_IBD_distance_cutoff A numeric value specifying the maximum Hamming distance allowed to infer
#'   the IBD status of a haplotype pair.
#' @param candidate_variants_concordant_with_SRS Logical value specifying if candidate variants should be concordant with short-read sequencing (default is `FALSE`).
#' @param SRS_vcf_file_path File path to the vcf containing variant calls from short-read sequencing of the duo.
#' @param test_reference_allele Logical value specifying if positions where the proband is heterozygous and the parent is homozygous for the variant allele (not the reference) should be tested (default is `FALSE`).
#' @param candidate_variant_coordinates Vector of coordinates for specific variants of interest (of the form c(chr1:1000, chr2:2000)).
#' @param problematic_regions BED file with coordinates of problematic regions (e.g. as defined by Genome-in-a-Bottle)
#' @param output_vcf_path File path for output vcf.
#' @param compress_output Logical value specifying if output_vcf should be compressed and indexed. appended .bgz to filename (default is `TRUE`).
#'
#' @return A `GRanges` containing additional columns for the variant classifications as well as supporting information.
#' @details The function ultimately works by detecting identical by descent haplotype blocks, to determine whether each candidate variant of interest is de novo, using the genotype of only one parent. If requested, concordance with short-read sequencing can be checked.
#'
#' @importFrom matrixStats colMins rowMaxs
#' @importFrom S4Vectors queryHits subjectHits DataFrame metadata append
#' @importFrom SummarizedExperiment rowRanges colData 
#' @import GenomicRanges
#' @import IRanges
#' @importFrom rtracklayer import
#' @import VariantAnnotation
#' @importFrom utils packageVersion
#' @export
#'
#' @examples
#' # duoNovo(LRS_phased_vcf_file_path = my_LRS_file_path, depth_cutoff = 20, GQ_cutoff = 30,
#' #          proband_column_identifier = my_proband_identifier,
#' #          PS_width_cutoff = 10000, boundary_cutoff = 2000, distance_cutoff = 40,
#' #          candidate_variants_concordant_with_SRS = FALSE,
#' #          SRS_vcf_file_path = my_SRS_file_path)
duoNovo <- function(LRS_phased_vcf_file_path, depth_cutoff = 20, GQ_cutoff = 30,
                    proband_column_identifier,
                    PS_width_cutoff = 10000, boundary_cutoff = 2000, 
                    IBD_distance_cutoff = 0, non_IBD_distance_cutoff = 40,
                    candidate_variants_concordant_with_SRS = FALSE, SRS_vcf_file_path = NULL, 
                    test_reference_allele = FALSE, 
                    candidate_variant_coordinates = NULL, problematic_regions = NULL,
                    output_vcf_path = NULL, compress_output = TRUE) {
  if (!file.exists(LRS_phased_vcf_file_path)) {
    stop("The LRS VCF file path does not exist: ", LRS_phased_vcf_file_path)
  }
  if (candidate_variants_concordant_with_SRS == TRUE && !file.exists(SRS_vcf_file_path)) {
    stop("The SRS VCF file path does not exist: ", SRS_vcf_file_path)
  }
  if (IBD_distance_cutoff >= non_IBD_distance_cutoff) {
    stop("The Hamming distance cutoff for IBD haplotypes exceeds that for non-IBD haplotypes")
  }

  message("Importing LRS VCF file...")
  vcf <- readVcf(LRS_phased_vcf_file_path)
  vcf_granges <- SummarizedExperiment::rowRanges(vcf)
  vcf_metadata <- geno(vcf)

  proband_column <- grep(proband_column_identifier, samples(header(vcf)))
  if (length(proband_column) == 0) {
    stop("Proband column identifier not found: ", proband_column_identifier)
  }

  if (dim(vcf)[2] != 2) {
    stop("VCF does not contain exactly 2 samples.  Number of samples: ", length(samples(header(vcf))))
  }

  if (!"GT" %in% names(vcf_metadata)) {
    stop("The 'GT' (genotype) field is missing in the VCF metadata from LRS.")
  }
  if (!"DP" %in% names(vcf_metadata)) {
    stop("The 'DP' (depth) field is missing in the VCF metadata from LRS.")
  }
  if (!"GQ" %in% names(vcf_metadata)) {
    stop("The 'GQ' (genotype quality) field is missing in the VCF metadata from LRS.")
  }
  if (!"PS" %in% names(vcf_metadata)) {
    stop("The 'PS' (phasing set) field is missing in the VCF metadata from LRS.")
  }

  vcf_granges$phasing1 <- vcf_metadata$GT[, proband_column]
  vcf_granges$phasing2 <- vcf_metadata$GT[, -proband_column]

  vcf_granges$depth1 <- vcf_metadata$DP[, proband_column]
  vcf_granges$depth2 <- vcf_metadata$DP[, -proband_column]

  vcf_granges$GQ1 <- vcf_metadata$GQ[, proband_column]
  vcf_granges$GQ2 <- vcf_metadata$GQ[, -proband_column]

  #the following extracts the phasing set to which each phased variant is assigned to in the proband and in the parent
  #and ensures the phasing set is represented by a unique identifier in the columns PS1 (for proband) and PS2 (for parent)
  vcf_granges$PS1 <- vcf_metadata$PS[, proband_column]
  has_PS1 <- !is.na(vcf_granges$PS1)
  vcf_granges$PS1[has_PS1] <- paste0(seqnames(vcf_granges[has_PS1]), "_", vcf_granges$PS1[has_PS1])
  vcf_granges$PS2 <- vcf_metadata$PS[, -proband_column]
  has_PS2 <- !is.na(vcf_granges$PS2)
  vcf_granges$PS2[has_PS2] <- paste0(seqnames(vcf_granges[has_PS2]), "_", vcf_granges$PS2[has_PS2])
  
  QC_fail_variants <- GRanges() # Initialize empty granges to store all variants that failed QC
  
  message("Reconstructing haplotypes...")
  hap_granges <- getHaplotypes(vcf_granges)

  low_depth_indices <- which(hap_granges$depth1 < depth_cutoff | hap_granges$depth2 < depth_cutoff)
  low_GQ_indices <- which(hap_granges$GQ1 < GQ_cutoff | hap_granges$GQ2 < GQ_cutoff)

  low_depth_or_GQ <- union(low_depth_indices, low_GQ_indices)
  low_depth_and_GQ <- intersect(low_depth_indices, low_GQ_indices)
  low_depth_only <- low_depth_indices[-which(low_depth_indices %in% low_GQ_indices)]
  low_GQ_only <- low_GQ_indices[-which(low_GQ_indices %in% low_depth_indices)]
  
  if (length(low_depth_only) > 0){
    hap_granges_low_depth <- hap_granges[low_depth_only]
    hap_granges_low_depth$QC_fail_step <- "low_depth"
  } else {
    hap_granges_low_depth <- GRanges()
  }
  if (length(low_GQ_only) > 0){
    hap_granges_low_GQ <- hap_granges[low_GQ_only]
    hap_granges_low_GQ$QC_fail_step <- "low_GQ"
  } else {
    hap_granges_low_GQ <- GRanges()
  }
  if (length(low_depth_and_GQ) > 0){
    hap_granges_low_depth_GQ <- hap_granges[low_depth_and_GQ]
    hap_granges_low_depth_GQ$QC_fail_step <- "low_depth_and_GQ"
  } else {
    hap_granges_low_depth_GQ<- GRanges()
  }
  QC_fail_variants <- c(hap_granges_low_GQ, hap_granges_low_depth, hap_granges_low_depth_GQ, QC_fail_variants)

  if (length(low_depth_or_GQ) == length(hap_granges)) {
    stop("No variants called from LRS pass depth and GQ thresholds.")
  }

  # If candidate variant coordinates are provided, filter hap_granges accordingly
  if (!is.null(candidate_variant_coordinates)) {
    split_coords <- strsplit(candidate_variant_coordinates, ":")
    seqnames <- sapply(split_coords, `[[`, 1)
    ranges <- sapply(split_coords, `[[`, 2)
    
    # Extract chromosome names and start/end coordinates
    ranges_split <- strsplit(ranges, "-")
    starts <- as.numeric(sapply(ranges_split, `[[`, 1))
    ends <- sapply(ranges_split, function(x) if (length(x) > 1) as.numeric(x[2]) else as.numeric(x[1]))
    
    # Create GRanges
    candidate_variant_granges <- GRanges(
      seqnames = seqnames,
      ranges = IRanges(
        start = starts,
        end = ends
      )
    )
    candidate_variant_indices <- unique(queryHits(findOverlaps(hap_granges, candidate_variant_granges)))
    ranges_to_subset <- hap_granges[candidate_variant_indices]
    
    candidate_variant_indices_left <- which(ranges_to_subset$phasing1 == "1|0" & ranges_to_subset$phasing2 == "0/0")  
    candidate_variant_indices_right <- which(ranges_to_subset$phasing1 == "0|1" & ranges_to_subset$phasing2 == "0/0")
    if (test_reference_allele == TRUE){
      candidate_variant_indices_left_ref <- which(ranges_to_subset$phasing1 == "0|1" & ranges_to_subset$phasing2 == "1/1")  
      candidate_variant_indices_right_ref <- which(ranges_to_subset$phasing1 == "1|0" & ranges_to_subset$phasing2 == "1/1")
      candidate_variant_indices_left <- c(candidate_variant_indices_left, candidate_variant_indices_left_ref)
      candidate_variant_indices_right <- c(candidate_variant_indices_right, candidate_variant_indices_right_ref)
    }
  } else { # Otherwise, identify candidate de novo variants directly from genotypes
    ranges_to_subset <- hap_granges
    
    candidate_variant_indices_left <- which(ranges_to_subset$phasing1 == "1|0" & ranges_to_subset$phasing2 == "0/0")  
    candidate_variant_indices_right <- which(ranges_to_subset$phasing1 == "0|1" & ranges_to_subset$phasing2 == "0/0")
    if (test_reference_allele == TRUE){
      candidate_variant_indices_left_ref <- which(ranges_to_subset$phasing1 == "0|1" & ranges_to_subset$phasing2 == "1/1")  
      candidate_variant_indices_right_ref <- which(ranges_to_subset$phasing1 == "1|0" & ranges_to_subset$phasing2 == "1/1")
      candidate_variant_indices_left <- c(candidate_variant_indices_left, candidate_variant_indices_left_ref)
      candidate_variant_indices_right <- c(candidate_variant_indices_right, candidate_variant_indices_right_ref)
    }
  }
  if (length(candidate_variant_indices_left) > 0){
    candidate_variant_granges_left <- ranges_to_subset[candidate_variant_indices_left]
    overlaps <- findOverlaps(QC_fail_variants, candidate_variant_granges_left)
    if (length(overlaps) > 0){
      QC_fail_variants_left <- QC_fail_variants[unique(queryHits(overlaps))]
      candidate_variant_granges_left <- candidate_variant_granges_left[-unique(subjectHits(overlaps))]
    } else {
      QC_fail_variants_left <- GRanges()
    }
  } else {
    candidate_variant_granges_left <- GRanges()
    QC_fail_variants_left <- GRanges()
  }
  if (length(candidate_variant_indices_right) > 0){
    candidate_variant_granges_right <- ranges_to_subset[candidate_variant_indices_right]
    overlaps <- findOverlaps(QC_fail_variants, candidate_variant_granges_right)
    if (length(overlaps) > 0){
      QC_fail_variants_right <- QC_fail_variants[unique(queryHits(overlaps))]
      candidate_variant_granges_right <- candidate_variant_granges_right[-unique(subjectHits(overlaps))]
    } else {
      QC_fail_variants_right <- GRanges()
    }  
  } else {
    candidate_variant_granges_right <- GRanges()
    QC_fail_variants_right <- GRanges()
  }  

  if (length(low_depth_or_GQ) > 0){
    hap_granges <- hap_granges[-low_depth_or_GQ]
  }
  hap_boundary_coordinates <- getHaplotypeBlockCoordinates(hap_granges)
  
  if (length(candidate_variant_granges_left) == 0 & length(candidate_variant_granges_right) == 0) {
    warning("No candidate variants passed QC.")
    return(c(QC_fail_variants_left, QC_fail_variants_right))
  }
  
  # Optional: Restrict to candidate variants concordant with short-read sequencing
  if (candidate_variants_concordant_with_SRS == TRUE) {
    message("Importing SRS VCF file...")
    vcf_SR <- readVcf(SRS_vcf_file_path)
    vcf_granges_SR <- SummarizedExperiment::rowRanges(vcf_SR)
    vcf_metadata_SR <- geno(vcf_SR)

    proband_column <- grep(proband_column_identifier, samples(header(vcf_SR)))
    if (length(proband_column) == 0) {
      stop("Proband column identifier not found in SRS vcf file: ", proband_column_identifier)
    }

    if (!"GT" %in% names(vcf_metadata_SR)) {
      stop("The 'GT' (genotype) field is missing in the VCF metadata from SRS.")
    }
    if (!"DP" %in% names(vcf_metadata_SR)) {
      stop("The 'DP' (depth) field is missing in the VCF metadata from SRS.")
    }
    if (!"GQ" %in% names(vcf_metadata_SR)) {
      stop("The 'GQ' (genotype quality) field is missing in the VCF metadata from SRS.")
    }

    vcf_granges_SR$gt1 <- vcf_metadata_SR$GT[, proband_column]
    vcf_granges_SR$gt2 <- vcf_metadata_SR$GT[, -proband_column]

    vcf_granges_SR$depth1 <- vcf_metadata_SR$DP[, proband_column]
    vcf_granges_SR$depth2 <- vcf_metadata_SR$DP[, -proband_column]

    vcf_granges_SR$GQ1 <- vcf_metadata_SR$GQ[, proband_column]
    vcf_granges_SR$GQ2 <- vcf_metadata_SR$GQ[, -proband_column]

    pass_depth_GQ_SR <- which(vcf_granges_SR$depth1 >= depth_cutoff & vcf_granges_SR$GQ1 >= GQ_cutoff &
                                vcf_granges_SR$depth2 >= depth_cutoff & vcf_granges_SR$GQ2 >= GQ_cutoff)
    if (length(pass_depth_GQ_SR) == 0) {
      stop("No variants called from SRS pass depth and GQ thresholds.")
    }
    vcf_granges_filtered_SR <- vcf_granges_SR[pass_depth_GQ_SR]
    candidate_indices_SR <- which(vcf_granges_filtered_SR$gt1 %in% c("1/0", "0/1") & 
                                    vcf_granges_filtered_SR$gt2 == "0/0")
    if (test_reference_allele == TRUE){
      candidate_indices_SR <- c(candidate_indices_SR, which(vcf_granges_filtered_SR$gt1 %in% c("1/0", "0/1") & 
                                      vcf_granges_filtered_SR$gt2 == "1/1"))
    }
    if (length(candidate_indices_SR) == 0) {
      stop("No variants called from SRS pass depth and GQ thresholds.")
    } 
    vcf_granges_filtered_SR <- vcf_granges_filtered_SR[candidate_indices_SR]
    
    QC_fail_variants_SR_left <- GRanges()
    QC_fail_variants_SR_right <- GRanges()
    
    SR_indices_ref <- which(vcf_granges_filtered_SR$gt2 == "0/0")
    SR_indices_alt <- which(vcf_granges_filtered_SR$gt2 == "1/1")
    
    if (length(candidate_variant_granges_left) > 0){
        concordant_genotype_overlaps_left_alt <- findOverlaps(
          candidate_variant_granges_left[which(candidate_variant_granges_left$phasing2 == "0/0")], 
          vcf_granges_filtered_SR[SR_indices_ref])
        
        concordant_genotype_overlaps_left_ref <- findOverlaps(
          candidate_variant_granges_left[which(candidate_variant_granges_left$phasing2 == "1/1")], 
          vcf_granges_filtered_SR[SR_indices_alt])
        
        queryHits_indices <- c(unique(queryHits(concordant_genotype_overlaps_left_alt)), 
                               unique(queryHits(concordant_genotype_overlaps_left_ref)))
        
        if (length(queryHits_indices) > 0){
          candidate_variant_granges_left <- candidate_variant_granges_left[queryHits_indices]
          
          QC_fail_variants_SR_left <- c(QC_fail_variants_SR_left, 
                                        candidate_variant_granges_left[-queryHits_indices])
          if (length(QC_fail_variants_SR_left) > 0){
            QC_fail_variants_SR_left$QC_fail_step <- "failed_SRS_concordance_QC"
          }
      } else {
        candidate_variant_granges_left <- GRanges()
        QC_fail_variants_SR_left <- c(QC_fail_variants_SR_left, 
                                      candidate_variant_granges_left)
        QC_fail_variants_SR_left$QC_fail_step <- "failed_SRS_concordance_QC"
      }
    }
    if (length(candidate_variant_granges_right) > 0){
      concordant_genotype_overlaps_right_alt <- findOverlaps(
        candidate_variant_granges_right[which(candidate_variant_granges_right$phasing2 == "0/0")], 
        vcf_granges_filtered_SR[SR_indices_ref])
      
      concordant_genotype_overlaps_right_ref <- findOverlaps(
        candidate_variant_granges_right[which(candidate_variant_granges_right$phasing2 == "1/1")], 
        vcf_granges_filtered_SR[SR_indices_alt])
      
      queryHits_indices <- c(unique(queryHits(concordant_genotype_overlaps_right_alt)), 
                             unique(queryHits(concordant_genotype_overlaps_right_ref)))
      
      if (length(queryHits_indices) > 0){
          candidate_variant_granges_right <- candidate_variant_granges_right[queryHits_indices]
          
          QC_fail_variants_SR_right <- c(QC_fail_variants_SR_right, 
                                        candidate_variant_granges_right[-queryHits_indices])
          if (length(QC_fail_variants_SR_right) > 0){
            QC_fail_variants_SR_right$QC_fail_step <- "failed_SRS_concordance_QC"
          }
      } else {
        candidate_variant_granges_right <- GRanges()
        QC_fail_variants_SR_right <- c(QC_fail_variants_SR_right, 
                                      candidate_variant_granges_right)
        QC_fail_variants_SR_right$QC_fail_step <- "failed_SRS_concordance_QC"
      }      
    }
    
    QC_fail_variants_left <- c(QC_fail_variants_left, QC_fail_variants_SR_left)
    QC_fail_variants_right <- c(QC_fail_variants_right, QC_fail_variants_SR_right)
    
    if (length(candidate_variant_granges_left) == 0 & length(candidate_variant_granges_right) == 0) {
      warning("No candidate variants passed QC.")
      return(c(QC_fail_variants_left, QC_fail_variants_right))
    }
  }
  
  # Finding overlaps and analyzing haplotypes
  message("Classifying variants...")
  if (length(candidate_variant_granges_left) > 0){
    classifications_left <- classifyVariants(candidate_variant_granges_left, phasing_orientation = "left", 
                                             haplotype_granges = hap_granges, 
                                             haplotype_boundary_coordinate_granges = hap_boundary_coordinates, 
                                             boundary_cutoff = boundary_cutoff, 
                                             non_IBD_distance_cutoff = non_IBD_distance_cutoff,
                                             IBD_distance_cutoff = IBD_distance_cutoff,  
                                             PS_width_cutoff = PS_width_cutoff, 
                                             QC_fail_variant_granges = QC_fail_variants_left)
  } else {
    if (length(QC_fail_variants_left) > 0){
      QC_fail_variants_left$duoNovo_classification <- "failed_QC"
      QC_fail_variants_left$supporting_hamming_distance <- NA
      QC_fail_variants_left$supporting_counts_het_hom <- NA
      QC_fail_variants_left$supporting_counts_het_het <- NA
      QC_fail_variants_left$supporting_counts_hom_het <- NA
      classifications_left <- QC_fail_variants_left
    } else {
      classifications_left <- GRanges()
    }
  }
  if (length(candidate_variant_granges_right) > 0){
    classifications_right <- classifyVariants(candidate_variant_granges_right, phasing_orientation = "right", 
                                              haplotype_granges = hap_granges, 
                                              haplotype_boundary_coordinate_granges = hap_boundary_coordinates, 
                                              boundary_cutoff = boundary_cutoff, 
                                              non_IBD_distance_cutoff = non_IBD_distance_cutoff,
                                              IBD_distance_cutoff = IBD_distance_cutoff, 
                                              PS_width_cutoff = PS_width_cutoff,
                                              QC_fail_variant_granges = QC_fail_variants_right)
  } else {
    if (length(QC_fail_variants_right) > 0){
      QC_fail_variants_right$duoNovo_classification <- "failed_QC"
      QC_fail_variants_right$supporting_hamming_distance <- NA
      QC_fail_variants_right$supporting_counts_het_hom <- NA
      QC_fail_variants_right$supporting_counts_het_het <- NA
      QC_fail_variants_right$supporting_counts_hom_het <- NA
      classifications_right <- QC_fail_variants_right
    } else {
      classifications_right <- GRanges()
    }
  }
  duo_novo_classifications <- c(classifications_left, classifications_right)
  output <- duo_novo_classifications
  output_sorted <- sort(output)
  
  #flag de novo variants that are clustered in the same phasing set -- these are likely false positives
  #first obtain phasing sets containing multiple variants classified as de novo
  output_sorted$n_de_novo_left_orientation_same_PS <- NA
  output_sorted$n_de_novo_right_orientation_same_PS <- NA
  de_novo_indices_left <- which(
    output_sorted$duoNovo_classification == "de_novo" & 
                                  (
                                    (output_sorted$phasing1 == "1|0" & output_sorted$phasing2 == "0/0") | 
                                      (output_sorted$phasing1 == "0|1" & output_sorted$phasing2 == "1/1")
                                    )
    )
  de_novo_indices_right <- which(
    output_sorted$duoNovo_classification == "de_novo" & 
      (
        (output_sorted$phasing1 == "0|1" & output_sorted$phasing2 == "0/0") | 
          (output_sorted$phasing1 == "1|0" & output_sorted$phasing2 == "1/1")
      )
  )
  if (length(de_novo_indices_left) > 0){
    de_novo_left <- output_sorted[de_novo_indices_left]
    dn_overlaps_left <- findOverlaps(de_novo_left, hap_boundary_coordinates)
    dn_by_ps_left <- split(queryHits(dn_overlaps_left), subjectHits(dn_overlaps_left))
    n_dn_in_ps_left <- lengths(dn_by_ps_left)
    indices_left <- unlist(dn_by_ps_left)
    counts_left <- rep(n_dn_in_ps_left, n_dn_in_ps_left)
    output_sorted$n_de_novo_left_orientation_same_PS[de_novo_indices_left[indices_left]] <- counts_left
  }
  if (length(de_novo_indices_right) > 0){
    de_novo_right <- output_sorted[de_novo_indices_right]
    dn_overlaps_right <- findOverlaps(de_novo_right, hap_boundary_coordinates)
    dn_by_ps_right <- split(queryHits(dn_overlaps_right), subjectHits(dn_overlaps_right))
    n_dn_in_ps_right <- lengths(dn_by_ps_right)
    indices_right <- unlist(dn_by_ps_right)
    counts_right <- rep(n_dn_in_ps_right, n_dn_in_ps_right)
    output_sorted$n_de_novo_right_orientation_same_PS[de_novo_indices_right[indices_right]] <- counts_right
  }
  multi_denovo_haplotype_indices <- which(output_sorted$n_de_novo_left_orientation_same_PS > 1 | 
                                               output_sorted$n_de_novo_right_orientation_same_PS > 1)
  if (length(multi_denovo_haplotype_indices) > 0){
    output_sorted$duoNovo_classification[multi_denovo_haplotype_indices] <- "on_multi_denovo_haplotype"
  }
  multi_denovo_mat <- as.matrix(
    mcols(output_sorted)[, c("n_de_novo_left_orientation_same_PS",
                  "n_de_novo_right_orientation_same_PS")]
  )
  if (all(is.na(multi_denovo_mat)) == FALSE){
    max_count <- rowMaxs(multi_denovo_mat, na.rm = TRUE)
    max_count[is.infinite(max_count)] <- NA
    output_sorted$n_de_novo_same_orientation_same_PS <- max_count
  } else {
    output_sorted$n_de_novo_same_orientation_same_PS <- NA
  }
  mcols(output_sorted)$n_de_novo_left_orientation_same_PS  <- NULL
  mcols(output_sorted)$n_de_novo_right_orientation_same_PS <- NULL
  
  if (!is.null(problematic_regions)){
    problematic_regions_bed <- rtracklayer::import(problematic_regions, format = "BED")
    problematic_region_overlap_indices <- unique(queryHits(findOverlaps(output_sorted, problematic_regions_bed)))
    if (length(problematic_region_overlap_indices) > 0){
      output_sorted$QC_fail_step[problematic_region_overlap_indices] <- paste0("classified_", 
                                            output_sorted$duoNovo_classification[problematic_region_overlap_indices], 
                                                                               "_in_problematic_region")
      output_sorted$duoNovo_classification[problematic_region_overlap_indices] <- "failed_QC"
    }
  }
  output_sorted$tested_allele <- 1
  if (test_reference_allele == TRUE){
    output_sorted$tested_allele[which(output_sorted$phasing1 %in% c("0|1", "1|0") & 
                                        output_sorted$phasing2 == "1/1")] <- 0
  }
  
  if (!is.null(output_vcf_path)){
    message("Writing classified variants into VCF...")
    # Add each new INFO field to the header
    vcf_header <- header(vcf)
    header_metadata <- meta(vcf_header)
    
    # Create new metadata entries in a similar format to existing entries
    # Adding the custom metadata lines to the header
    duoNovo_version <- "unknown"
    if (requireNamespace("duoNovo", quietly = TRUE)) {
       duoNovo_version <-packageVersion("duoNovo")
    }


    description_values <- c(
      paste0( ifelse(is.null(duoNovo_version), "NA", paste0(duoNovo_version, collapse="." ) )),
      paste0( ifelse(is.null(LRS_phased_vcf_file_path), "NA", LRS_phased_vcf_file_path)),
      paste0( ifelse(is.null(GQ_cutoff), "NA", GQ_cutoff)),
      paste0( ifelse(is.null(depth_cutoff), "NA", depth_cutoff)),
      paste0( ifelse(is.null(proband_column_identifier), "NA", proband_column_identifier)),
      paste0( ifelse(is.null(PS_width_cutoff), "NA", PS_width_cutoff)),
      paste0( ifelse(is.null(boundary_cutoff), "NA", boundary_cutoff)),
      paste0( ifelse(is.null(non_IBD_distance_cutoff), "NA", non_IBD_distance_cutoff)),
      paste0( ifelse(is.null(IBD_distance_cutoff), "NA", IBD_distance_cutoff)),
      paste0( ifelse(is.null(candidate_variants_concordant_with_SRS), "NA", candidate_variants_concordant_with_SRS)),
      paste0( ifelse(is.null(test_reference_allele), "NA", test_reference_allele)),
      paste0( ifelse(is.null(SRS_vcf_file_path), "NA", SRS_vcf_file_path)),
      paste0( ifelse(is.null(candidate_variant_coordinates), "NA", candidate_variant_coordinates)),
      paste0( ifelse(is.null(problematic_regions), "NA", problematic_regions)),
      paste0( ifelse(is.null(output_vcf_path), "NA", output_vcf_path)),
      paste0( ifelse(is.null(compress_output), "NA", compress_output))
    )
    
    additional_metadata <- DataFrame(
      Value  = description_values,
      row.names = c("Version",
                    "LRS_phased_input_vcf",
                    "minGQ", 
                    "minDepth",
                    "proband_column_identifier",
                    "PS_width_cutoff",
                    "boundary_cutoff",
                    "IBD_distance_cutoff",
                    "non_IBD_distance_cutoff",
                    "candidate_variants_concordant_with_SRS",
                    "test_reference_allele",
                    "SRS_vcf_file_path",
                    "candidate_variant_coordinates",
                    "problematic_regions",
                    "output_vcf_path",
                    "compress_output")
    )
    
    # Combine the metadata with the new entries
    combined_metadata <- S4Vectors::append(header_metadata, list(duoNovoPARAM=additional_metadata))
    
    # Update the header
    meta(vcf_header) <- combined_metadata
    
    new_info_fields <- DataFrame(
      row.names = c("phasing_proband", "phasing_parent", "depth_proband", "depth_parent",
                    "GQ_proband", "GQ_parent", "duoNovo_classification",
                    "supporting_hamming_distance", "supporting_counts_het_hom",
                    "supporting_counts_het_het", "supporting_counts_hom_het",
                    "QC_fail_step", "tested_allele", "n_de_novo_same_orientation_same_PS"),
      Number = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
      Type = c("String", "String", "Integer", "Integer", "Integer", "Integer", "String",
               "Integer", "Integer", "Integer", "Integer", "String", "Integer", "Integer"),
      Description = c("Phasing for proband", "Phasing for parent", "Depth for proband",
                      "Depth for parent", "Genotype quality for proband",
                      "Genotype quality for parent", "DuoNovo classification",
                      "Supporting Hamming distance", "Supporting counts (het-hom)",
                      "Supporting counts (het-het)", "Supporting counts (hom-het)",
                      "QC fail step (NA for variants that passed QC)", 
                      "Allele tested for de novo status (1 for ALT, 0 for REF)", 
                      "Total number of de novo variants in same orientation and same phasing set (NA for non-de novo variants)")
    )
    
    # Add each new INFO field to the header
    info(vcf_header) <- rbind(info(vcf_header), new_info_fields)

    fixed_fields <- fixed(vcf) 
    rownames(fixed_fields) <- rownames(info(vcf))
    fixed_fields <- fixed_fields[names(output_sorted), ]
    
    info_new <- DataFrame(
      phasing_proband = output_sorted$phasing1,
      phasing_parent = output_sorted$phasing2,
      depth_proband = output_sorted$depth1,
      depth_parent = output_sorted$depth2,
      GQ_proband = output_sorted$GQ1,
      GQ_parent = output_sorted$GQ2,
      duoNovo_classification = output_sorted$duoNovo_classification,
      supporting_hamming_distance = output_sorted$supporting_hamming_distance,
      supporting_counts_het_hom = output_sorted$supporting_counts_het_hom,
      supporting_counts_het_het = output_sorted$supporting_counts_het_het,
      supporting_counts_hom_het = output_sorted$supporting_counts_hom_het,
      QC_fail_step = output_sorted$QC_fail_step, 
      tested_allele = output_sorted$tested_allele,
      n_de_novo_same_orientation_same_PS = output_sorted$n_de_novo_same_orientation_same_PS
    )
    info <- cbind(info(vcf)[names(output_sorted), ], info_new)
    
    sample_info <- colData(vcf)
    geno_data <- geno(vcf)
    geno_data <- lapply(geno_data, function(mat) {
      if (length(dim(mat)) == 2) {
        # 2D: subset rows and retain all columns
        mat[names(output_sorted), , drop = FALSE]
      } else if (length(dim(mat)) == 3) {
        # 3D: subset rows and retain all slices and columns
        mat[names(output_sorted), , , drop = FALSE]
      }
    })
    
    vcf_out <- VCF(
      rowRanges = output_sorted, 
      fixed = fixed_fields,
      colData = sample_info,       # Retain the original sample information
      info = info,                 # Add the new INFO metadata fields
      geno = geno_data             # Add the geno_data after subsetting for the rows present in our output
    )
    
    S4Vectors::metadata(vcf_out)$header <- vcf_header
    writeVcf(vcf_out, output_vcf_path, index = compress_output)  
    }
  
  ###Also give a verbose summary of results
  # Generate the table of counts for each classification
  classification_counts <- table(output_sorted$duoNovo_classification)
  
  # Define the expected classification categories
  expected_classes <- c("de_novo", "on_other_parent_haplotype", "uncertain", 
                        "failed_QC", "on_multi_denovo_haplotype")
  
  # Add missing classes if they are not present in the table (to ensure all classes are covered)
  for (class in expected_classes) {
    if (!class %in% names(classification_counts)) {
      classification_counts[class] <- 0
    }
  }
  
  # Sort the counts in the expected order
  classification_counts <- classification_counts[expected_classes]
  
  # Create a verbose summary
  verbose_summary <- paste0(
    "Variant Classification Summary:\n",
    "--------------------------------\n",
    "Number of de novo variants: ", classification_counts["de_novo"], "\n",
    "Number of variants present on haplotype inherited from non-sequenced parent: ", classification_counts["on_other_parent_haplotype"], "\n",
    "Number of uncertain variants: ", classification_counts["uncertain"], "\n",
    "Number of variants on multi-denovo haplotypes: ", classification_counts["on_multi_denovo_haplotype"], "\n",
    "Number of variants that failed QC: ", classification_counts["failed_QC"], "\n",
    "--------------------------------"
  )
  # Print the verbose summary
  cat(verbose_summary)
  return(output_sorted)
}
