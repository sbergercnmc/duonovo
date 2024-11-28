#' Identify de novo variants from duo (proband - single parent) long read sequencing
#'
#' The `duoNovo` function uses phased variant calls from duo long-read sequencing to classify candidate variants (heterozygous in the proband; absent in the parent).
#' Variants are classifed as de novo, present on the haplotype inherited from the non-sequenced parent, or uncertain variants.
#'
#' @param LRS_phased_vcf_file_path File path to the vcf containing phased variant calls from long-read sequencing of the duo.
#' @param depth_cutoff A numeric value specifying the minimum sequencing depth for variants to be included in the analysis. The same cutoff applies to both proband and parent variants, for both LRS and SRS (if used).
#' @param GQ_cutoff A numeric value specifying the minimum GQ (genotype quality) for variants to be included in the analysis. The same cutoff applies to both proband and parent variants, for both LRS and SRS (if used).
#' @param proband_column_identifier A character corresponding to an identifier for the proband column in the metadata matrices. Should be the same for both the LRS vcf and (if used) the SRS vcf.
#' @param PS_width_cutoff A numeric value specifying the minimum width for phasing sets to be included in the analysis.
#' @param boundary_cutoff A numeric value indicating the minimum distance from a haplotype block boundary (either start or end coordinate) for candidate variants to be analyzed.
#' @param distance_cutoff A numeric value specifying the minimum hamming distance cutoff to determine that a proband-parent haplotype block are not identical by descent.
#' @param candidate_variants_concordant_with_SRS Logical value specifying if candidate variants should be concordant with short-read sequencing (default is `TRUE` or `FALSE`).
#' @param SRS_vcf_file_path File path to the vcf containing variant calls from short-read sequencing of the duo.
#' @param reference Reference genome name (e.g. hg38) used in the vcfs.
#' @param candidate_variant_coordinates List of coordinates for specific variants of interest.
#' @param output_vcf_path File path for output vcf.
#'
#' @return A `GRangesList` containing three elements: `de_novo`, `not_de_novo`, and `uncertain`.
#' @details The function ultimately works by detecting identical by descent haplotype blocks, to determine whether each candidate variant of interest is de novo, using the genotype of only one parent. If requested, concordance with short-read sequencing can be checked.
#'
#' @importFrom matrixStats colMins 
#' @importFrom S4Vectors queryHits subjectHits DataFrame
#' @importFrom SummarizedExperiment rowRanges
#' @import GenomicRanges
#' @import IRanges
#' @import VariantAnnotation
#' @export
#'
#' @examples
#' # duoNovo(LRS_phased_vcf_file_path = my_LRS_file_path,
#' #          proband_phasing = "1|0", proband_column_identifier = "-0$",
#' #          PS_width_cutoff = 10000, boundary_cutoff = 10000, distance_cutoff = 2000,
#' #          candidate_variants_concordant_with_SRS = TRUE,
#' #          SRS_vcf_file_path = my_SRS_file_path, reference = "hg38")
duoNovo <- function(LRS_phased_vcf_file_path, depth_cutoff = 20, GQ_cutoff = 30,
                    proband_column_identifier,
                    PS_width_cutoff = 10000, boundary_cutoff = 2000, distance_cutoff = 40,
                    candidate_variants_concordant_with_SRS = TRUE,
                    SRS_vcf_file_path, reference = "hg38", 
                    candidate_variant_coordinates = NULL, 
                    output_vcf_path = NULL) {
  if (!file.exists(LRS_phased_vcf_file_path)) {
    stop("The LRS VCF file path does not exist: ", LRS_phased_vcf_file_path)
  }
  if (candidate_variants_concordant_with_SRS == TRUE && !file.exists(SRS_vcf_file_path)) {
    stop("The SRS VCF file path does not exist: ", SRS_vcf_file_path)
  }

  message("Importing LRS VCF file...")
  vcf <- readVcf(LRS_phased_vcf_file_path, genome = reference)
  vcf_granges <- SummarizedExperiment::rowRanges(vcf)
  vcf_metadata <- geno(vcf)

  proband_column <- grep(proband_column_identifier, samples(header(vcf)))
  if (length(proband_column) == 0) {
    stop("Proband column identifier not found: ", proband_column_identifier)
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

  vcf_granges$PS1 <- vcf_metadata$PS[, proband_column]
  vcf_granges$PS2 <- vcf_metadata$PS[, -proband_column]

  QC_fail_variants <- GRanges() # Initialize empty granges to store all variants that failed QC
  
  message("Reconstructing haplotypes...")
  hap_granges <- getHaplotypes(vcf_granges)
  pass_depth_GQ <- which(hap_granges$depth1 >= depth_cutoff & hap_granges$GQ1 >= GQ_cutoff &
                           hap_granges$depth2 >= depth_cutoff & hap_granges$GQ2 >= GQ_cutoff)
  
  low_depth_indices <- which(hap_granges$depth1 < depth_cutoff | hap_granges$depth2 < depth_cutoff)
  hap_granges_low_depth <- hap_granges[low_depth_indices]
  hap_granges_low_depth$QC_fail_step <- "low_depth"
  
  low_GQ_indices <- which(hap_granges$GQ1 < GQ_cutoff | hap_granges$GQ2 < GQ_cutoff)
  hap_granges_low_GQ <- hap_granges[low_GQ_indices]
  hap_granges_low_GQ$QC_fail_step <- "low_GQ"
  
  hap_granges <- hap_granges[pass_depth_GQ]
  hap_boundary_coordinates <- getHaplotypeBlockCoordinates(hap_granges)

  # If candidate variant coordinates are provided, filter hap_granges accordingly
  if (!is.null(candidate_variant_coordinates)) {
    candidate_variant_granges <- GRanges(
      seqnames = sapply(strsplit(candidate_variant_coordinates, ":"), `[[`, 1),
      ranges = IRanges(
        start = as.numeric(sapply(strsplit(sapply(strsplit(candidate_variant_coordinates, ":"), `[[`, 2), "-"), `[[`, 1)),
        end = as.numeric(sapply(strsplit(sapply(strsplit(candidate_variant_coordinates, ":"), `[[`, 2), "-"), `[[`, 2))
      ))
    QC_fail_variants <- c(QC_fail_variants, hap_granges_low_depth[
      unique(queryHits(findOverlaps(hap_granges_low_depth, candidate_variant_granges)))])
    QC_fail_variants <- c(QC_fail_variants, hap_granges_low_GQ[
      unique(queryHits(findOverlaps(hap_granges_low_GQ, candidate_variant_granges)))])
    
    candidate_variant_indices <- unique(queryHits(findOverlaps(hap_granges, candidate_variant_granges)))
    ranges_to_subset <- hap_granges[candidate_variant_indices]
    candidate_variant_indices_left <- which(ranges_to_subset$phasing1 == "1|0" & ranges_to_subset$phasing2 == "0/0")  
    candidate_variant_indices_right <- which(ranges_to_subset$phasing1 == "0|1" & ranges_to_subset$phasing2 == "0/0")
  } else { # Otherwise, identify candidate de novo variants directly from genotypes
    QC_fail_variants <- c(QC_fail_variants, hap_granges_low_depth)
    QC_fail_variants <- c(QC_fail_variants, hap_granges_low_GQ)
    
    ranges_to_subset <- hap_granges
    candidate_variant_indices_left <- which(ranges_to_subset$phasing1 == "1|0" & ranges_to_subset$phasing2 == "0/0")  
    candidate_variant_indices_right <- which(ranges_to_subset$phasing1 == "0|1" & ranges_to_subset$phasing2 == "0/0")
  }
  if (length(candidate_variant_indices_left) > 0){
    candidate_variant_granges_left <- ranges_to_subset[candidate_variant_indices_left]
  } else {
    candidate_variant_granges_left <- GRanges()
  }
  if (length(candidate_variant_indices_right) > 0){
    candidate_variant_granges_right <- ranges_to_subset[candidate_variant_indices_right]
  } else {
    candidate_variant_granges_right <- GRanges()
  }  
  
  # Optional: Restrict to candidate variants concordant with short-read sequencing
  if (candidate_variants_concordant_with_SRS == TRUE) {
    message("Importing SRS VCF file...")
    vcf_SR <- readVcf(SRS_vcf_file_path, genome = reference)
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

    vcf_granges_filtered_SR <- vcf_granges_SR[which(vcf_granges_SR$depth1 >= depth_cutoff & vcf_granges_SR$GQ1 >= GQ_cutoff &
                                                      vcf_granges_SR$depth2 >= depth_cutoff & vcf_granges_SR$GQ2 >= GQ_cutoff)]

    vcf_granges_filtered_SR <- vcf_granges_filtered_SR[
      which(vcf_granges_filtered_SR$gt1 %in% c("1/0", "0/1") & vcf_granges_filtered_SR$gt2 == "0/0")]
    
    candidate_variant_granges_left <- candidate_variant_granges_left[
      unique(queryHits(findOverlaps(candidate_variant_granges_left, vcf_granges_filtered_SR)))]
    candidate_variant_granges_right <- candidate_variant_granges_right[
      unique(queryHits(findOverlaps(candidate_variant_granges_left, vcf_granges_filtered_SR)))]
  }

  # Finding overlaps and analyzing haplotypes
  message("Classifying variants...")
  if (length(candidate_variant_granges_left) > 0){
    classifications_left <- classifyVariants(candidate_variant_granges_left, phasing_orientation = "left", 
                                             haplotype_granges = hap_granges, 
                                             haplotype_boundary_coordinate_granges = hap_boundary_coordinates, 
                                             boundary_cutoff = boundary_cutoff, distance_cutoff = distance_cutoff, 
                                             PS_width_cutoff = PS_width_cutoff, 
                                             QC_fail_variant_granges = QC_fail_variants)
  } else {
    classifications_left <- GRanges()
  }
  if (length(candidate_variant_granges_right) > 0){
    classifications_right <- classifyVariants(candidate_variant_granges_right, phasing_orientation = "right", 
                                              haplotype_granges = hap_granges, 
                                              haplotype_boundary_coordinate_granges = hap_boundary_coordinates, 
                                              boundary_cutoff = boundary_cutoff, distance_cutoff = distance_cutoff, 
                                              PS_width_cutoff = PS_width_cutoff,
                                              QC_fail_variant_granges = QC_fail_variants)
  } else {
    classifications_right <- GRanges()
  }
  duo_novo_classifications <- c(classifications_left, classifications_right)
  output <- duo_novo_classifications
  output_sorted <- sort(output)
  
  if (!is.null(output_vcf_path)){
    message("Writing classified variants to VCF file...")
    info <- DataFrame(
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
      QC_fail_step = output_sorted$QC_fail_step
    )
    
    # Create a VCF header with the meta-information
    vcf_header <- VCFHeader()
    
    # Add INFO fields to the VCF header
    info_fields <- data.frame(
      ID = c("phasing_proband", "phasing_parent", "depth_proband", "depth_parent",
             "GQ_proband", "GQ_parent", "duoNovo_classification",
             "supporting_hamming_distance", "supporting_counts_het_hom",
             "supporting_counts_het_het", "supporting_counts_hom_het", "QC_fail_step"),
      Number = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
      Type = c("String", "String", "Integer", "Integer", "Integer", "Integer", "String",
               "Integer", "Integer", "Integer", "Integer", "String"),
      Description = c("Phasing of the proband",
                      "Phasing of the parent",
                      "Read depth for the proband",
                      "Read depth for the parent",
                      "Genotype Quality (GQ) for the proband",
                      "Genotype Quality (GQ) for the parent",
                      "Classification result by duoNovo ('de_novo', 'on_other_parent_haplotype', 'uncertain', 'failed_QC')",
                      "Hamming distance to support the classification",
                      "Count of heterozygous proband and homozygous parent variant pairs supporting the classification",
                      "Count of heterozygous proband and heterozygous parent variant pairs supporting the classification",
                      "Count of homozygous proband and heterozygous parent variant pairs supporting the classification",
                      "Reason for QC failure (NA for variants that passed QC)")
    )
    
    for (i in 1:nrow(info_fields)) {
      info(vcf_header) <- DataFrameRow(
        rownames = info_fields$ID[i],
        Number = info_fields$Number[i],
        Type = info_fields$Type[i],
        Description = info_fields$Description[i]
      )
    }
    
    # Create the VCF object with the row ranges and info
    vcf_out <- VCF(
      rowRanges = output_sorted,
      info = info,
      header = vcf_header  # Add the header to the VCF object
    )
    
    # Write the VCF to a file, with index
    writeVcf(vcf_out, output_vcf_path, index = TRUE)
  }
  return(output_sorted)
  
  ###Also give a verbose summary of results
  # Generate the table of counts for each classification
  classification_counts <- table(output_sorted$duoNovo_classification)
  
  # Define the expected classification categories
  expected_classes <- c("de_novo", "on_other_parent_haplotype", "uncertain", "failed_QC")
  
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
    "Number of variants that failed QC: ", classification_counts["failed_QC"], "\n",
    "--------------------------------"
  )
  # Print the verbose summary
  cat(verbose_summary)
}
