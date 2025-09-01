duoNovoSib <- function(LRS_phased_vcf_file_path, depth_cutoff = 20, GQ_cutoff = 30,
                       candidate_variant_coordinates = NULL, candidate_variants_duoNovo_output = NULL,
                       problematic_regions = NULL,
                       PS_width_cutoff = 10000, boundary_cutoff = 2000, distance_cutoff = 40,
                       output_vcf_path = NULL, compress_output = TRUE) {
  
  if (!file.exists(LRS_phased_vcf_file_path)) {
    stop("The LRS VCF file path does not exist: ", LRS_phased_vcf_file_path)
  }
  
  message("Importing LRS VCF file...")
  vcf <- readVcf(LRS_phased_vcf_file_path)
  vcf_granges <- SummarizedExperiment::rowRanges(vcf)
  vcf_metadata <- geno(vcf)
  
  proband_column <- 1
  sibling_column <- 2
  mother_column <- 3
  
  if (length(proband_column) == 0) {
    stop("Proband column identifier not found: ", proband_column_identifier)
  }
  
  if (dim(vcf)[2] != 3) {
    stop("VCF does not contain exactly 3 samples.  Number of samples: ", length(samples(header(vcf))))
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
  vcf_granges$phasing2 <- vcf_metadata$GT[, sibling_column]
  vcf_granges$phasing3 <- vcf_metadata$GT[, mother_column]
  
  vcf_granges$depth1 <- vcf_metadata$DP[, proband_column]
  vcf_granges$depth2 <- vcf_metadata$DP[, sibling_column]
  vcf_granges$depth3 <- vcf_metadata$DP[, mother_column]
  
  vcf_granges$GQ1 <- vcf_metadata$GQ[, proband_column]
  vcf_granges$GQ2 <- vcf_metadata$GQ[, sibling_column]
  vcf_granges$GQ3 <- vcf_metadata$GQ[, mother_column]
  
  #the following extracts the phasing set to which each phased variant is assigned to in the proband and in the parent
  #and ensures the phasing set is represented by a unique identifier in the columns PS1 (for proband) and PS2 (for parent)
  vcf_granges$PS1 <- vcf_metadata$PS[, proband_column]
  has_PS1 <- !is.na(vcf_granges$PS1)
  vcf_granges$PS1[has_PS1] <- paste0(seqnames(vcf_granges[has_PS1]), "_", vcf_granges$PS1[has_PS1])
  
  vcf_granges$PS2 <- vcf_metadata$PS[, sibling_column]
  has_PS2 <- !is.na(vcf_granges$PS2)
  vcf_granges$PS2[has_PS2] <- paste0(seqnames(vcf_granges[has_PS2]), "_", vcf_granges$PS2[has_PS2])
  
  vcf_granges$PS3 <- vcf_metadata$PS[, mother_column]
  has_PS3 <- !is.na(vcf_granges$PS3)
  vcf_granges$PS3[has_PS3] <- paste0(seqnames(vcf_granges[has_PS3]), "_", vcf_granges$PS3[has_PS3])
  
  QC_fail_variants <- GRanges() # Initialize empty granges to store all variants that failed QC
  
  message("Reconstructing haplotypes...")
  hap_granges <- getHaplotypesForTrio(vcf_granges)
  
  low_depth_indices <- which(hap_granges$depth1 < depth_cutoff | 
                               hap_granges$depth2 < depth_cutoff | hap_granges$depth3 < depth_cutoff)
  low_GQ_indices <- which(hap_granges$GQ1 < GQ_cutoff | hap_granges$GQ2 < GQ_cutoff | hap_granges$GQ3 < GQ_cutoff)
  
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
    
    candidate_variant_indices_left <- which(ranges_to_subset$phasing1 == "1|0" & 
                                              ranges_to_subset$phasing2 == "0/0" & 
                                              ranges_to_subset$phasing3 == "0/0")  
    candidate_variant_indices_right <- which(ranges_to_subset$phasing1 == "0|1" & 
                                               ranges_to_subset$phasing2 == "0/0" & 
                                               ranges_to_subset$phasing3 == "0/0")
    
  } else if (!is.null(candidate_variants_duoNovo_output)){
    candidate_variant_indices <- which(names(hap_granges) %in% names(candidate_variants_duoNovo_output))
    ranges_to_subset <- hap_granges[candidate_variant_indices] #this ensures we test same variants as previously classified
    
    candidate_variant_indices_left <- which(ranges_to_subset$phasing1 == "1|0" & 
                                              ranges_to_subset$phasing2 == "0/0" & 
                                              ranges_to_subset$phasing3 == "0/0")  
    candidate_variant_indices_right <- which(ranges_to_subset$phasing1 == "0|1" & 
                                               ranges_to_subset$phasing2 == "0/0" & 
                                               ranges_to_subset$phasing3 == "0/0")
    
  } else { # Otherwise, identify candidate de novo variants directly from genotypes
    ranges_to_subset <- hap_granges
    
    candidate_variant_indices_left <- which(ranges_to_subset$phasing1 == "1|0" & 
                                              ranges_to_subset$phasing2 == "0/0" & 
                                              ranges_to_subset$phasing3 == "0/0")  
    candidate_variant_indices_right <- which(ranges_to_subset$phasing1 == "0|1" & 
                                               ranges_to_subset$phasing2 == "0/0" & 
                                               ranges_to_subset$phasing3 == "0/0")
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
  hap_boundary_coordinates <- getHaplotypeBlockCoordinatesForTrio(hap_granges)
  
  if (length(candidate_variant_granges_left) == 0 & length(candidate_variant_granges_right) == 0) {
    warning("No candidate variants passed QC.")
    return(c(QC_fail_variants_left, QC_fail_variants_right))
  }
  
  # Finding overlaps and analyzing haplotypes
  message("Classifying variants...")
  if (length(candidate_variant_granges_left) > 0){
    classifications_left <- classifyVariantsTrio(candidate_variant_granges_left, phasing_orientation = "left", 
                                                 haplotype_granges = hap_granges, 
                                                 haplotype_boundary_coordinate_granges = hap_boundary_coordinates, 
                                                 boundary_cutoff = boundary_cutoff, distance_cutoff = distance_cutoff, 
                                                 PS_width_cutoff = PS_width_cutoff, 
                                                 QC_fail_variant_granges = QC_fail_variants_left)
  } else {
    classifications_left <- GRanges()
  }
  if (length(candidate_variant_granges_right) > 0){
    classifications_right <- classifyVariantsTrio(candidate_variant_granges_right, phasing_orientation = "right", 
                                                  haplotype_granges = hap_granges, 
                                                  haplotype_boundary_coordinate_granges = hap_boundary_coordinates, 
                                                  boundary_cutoff = boundary_cutoff, distance_cutoff = distance_cutoff, 
                                                  PS_width_cutoff = PS_width_cutoff,
                                                  QC_fail_variant_granges = QC_fail_variants_right)
  } else {
    classifications_right <- GRanges()
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
      (output_sorted$phasing1 == "1|0" & output_sorted$phasing2 == "0/0" & 
         output_sorted$phasing3 == "0/0"
      )
  )
  de_novo_indices_right <- which(
    output_sorted$duoNovo_classification == "de_novo" & 
      (
        output_sorted$phasing1 == "0|1" & output_sorted$phasing2 == "0/0" &
          output_sorted$phasing3 == "0/0"
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
  max_count <- rowMaxs(multi_denovo_mat, na.rm = TRUE)
  max_count[is.infinite(max_count)] <- NA
  
  output_sorted$n_de_novo_same_orientation_same_PS <- max_count
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
      paste0( ifelse(is.null(PS_width_cutoff), "NA", PS_width_cutoff)),
      paste0( ifelse(is.null(boundary_cutoff), "NA", boundary_cutoff)),
      paste0( ifelse(is.null(distance_cutoff), "NA", distance_cutoff)),
      paste0( ifelse(is.null(candidate_variant_coordinates), "NA", candidate_variant_coordinates)),
      paste0( ifelse(is.null(candidate_variants_duoNovo_output), "NA", candidate_variants_duoNovo_output)),
      paste0( ifelse(is.null(output_vcf_path), "NA", output_vcf_path)),
      paste0( ifelse(is.null(compress_output), "NA", compress_output))
    )
    
    additional_metadata <- DataFrame(
      Value  = description_values,
      row.names = c("Version",
                    "LRS_phased_input_vcf",
                    "minGQ", 
                    "minDepth",
                    "PS_width_cutoff",
                    "boundary_cutoff",
                    "distance_cutoff",
                    "candidate_variant_coordinates",
                    "candidate_variants_from_duoNovo_output",
                    "output_vcf_path",
                    "compress_output")
    )
    
    # Combine the metadata with the new entries
    combined_metadata <- S4Vectors::append(header_metadata, list(duoNovoPARAM=additional_metadata))
    
    # Update the header
    meta(vcf_header) <- combined_metadata
    
    new_info_fields <- DataFrame(
      row.names = c("phasing_proband", "phasing_sibling", "phasing_mother",
                    "depth_proband", "depth_sibling", "depth_mother",
                    "GQ_proband", "GQ_sibling", "GQ_mother",
                    "duoNovo_classification", "parent_of_origin",
                    "QC_fail_step", "n_de_novo_left_orientation_same_PS", "n_de_novo_right_orientation_same_PS"),
      Number = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
      Type = c("String", "String", "String", "Integer", "Integer", "Integer", 
               "Integer", "Integer", "Integer", "String", "String", "String",
               "Integer", "Integer"),
      Description = c("Phasing for proband", "Phasing for sibling", "Phasing for mother",
                      "Depth for proband", "Depth for sibling", "Depth for mother",
                      "Genotype quality for proband", "Genotype quality for sibling", 
                      "Genotype quality for mother",
                      "DuoNovo classification",
                      "Parent in whose haplotype the de novo variant arose (NA for non-de novo variants)",
                      "QC fail step (NA for variants that passed QC)", 
                      "Total number of de novo variants in left orientation and same phasing set (NA for non-de novo variants or de novo variants in right orientation)", 
                      "Total number of de novo variants in right orientation and same phasing set (NA for non-de novo variants or de novo variants in left orientation)")
    )
    
    # Add each new INFO field to the header
    info(vcf_header) <- rbind(info(vcf_header), new_info_fields)
    
    fixed_fields <- fixed(vcf) 
    rownames(fixed_fields) <- rownames(info(vcf))
    fixed_fields <- fixed_fields[names(output_sorted), ]
    
    info_new <- DataFrame(
      phasing_proband = output_sorted$phasing1,
      phasing_sibling = output_sorted$phasing2,
      phasing_mother = output_sorted$phasing3,
      depth_proband = output_sorted$depth1,
      depth_sibling = output_sorted$depth2,
      depth_mother = output_sorted$depth3,
      GQ_proband = output_sorted$GQ1,
      GQ_sibling = output_sorted$GQ2,
      GQ_mother = output_sorted$GQ3,
      duoNovo_classification = output_sorted$duoNovo_classification,
      parent_of_origin = output_sorted$parent_origin,
      QC_fail_step = output_sorted$QC_fail_step, 
      n_de_novo_left_orientation_same_PS = output_sorted$n_de_novo_left_orientation_same_PS,
      n_de_novo_right_orientation_same_PS = output_sorted$n_de_novo_right_orientation_same_PS
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
  return(output_sorted)
}


###
###
classifyVariantsTrio <- function(candidate_variant_granges, phasing_orientation = c("left", "right"),
                                 haplotype_granges, haplotype_boundary_coordinate_granges,
                                 boundary_cutoff, distance_cutoff, PS_width_cutoff,
                                 QC_fail_variant_granges){
  
  ###QC steps
  QC_fail_variants <- QC_fail_variant_granges
  
  #first obtain variants that do not overlap any of the haplotype blocks
  hap_overlap_indices <- unique(queryHits(findOverlaps(candidate_variant_granges,
                                                       haplotype_boundary_coordinate_granges)))
  if (length(hap_overlap_indices) == 0) {
    QC_fail_variants_no_hap_overlap <- candidate_variant_granges
    QC_fail_variants_no_hap_overlap$QC_fail_step <- "no_haplotype_block_overlap"
    QC_fail_variants <- c(QC_fail_variants, QC_fail_variants_no_hap_overlap)
    QC_fail_variants$duoNovo_classification <- "failed_QC"
    QC_fail_variants$parent_origin <- NA
    warning(paste0("No candidate variants of ", phasing_orientation, " phasing orientation passed QC."))
    return(QC_fail_variants)
  }
  if (length(hap_overlap_indices) < length(candidate_variant_granges)){
    QC_fail_variants_no_hap_overlap <- candidate_variant_granges[-hap_overlap_indices]
    QC_fail_variants_no_hap_overlap$QC_fail_step <- "no_haplotype_block_overlap"
    QC_fail_variants <- c(QC_fail_variants, QC_fail_variants_no_hap_overlap)
    candidate_variant_granges <- candidate_variant_granges[hap_overlap_indices]
  }
  
  #now obtain those that do not overlap any after filtering out small haplotype blocks
  haplotype_boundary_coordinate_granges <- haplotype_boundary_coordinate_granges[
    which(width(haplotype_boundary_coordinate_granges) > PS_width_cutoff)]
  
  hap_overlap_indices <- unique(queryHits(findOverlaps(candidate_variant_granges,
                                                       haplotype_boundary_coordinate_granges)))
  if (length(hap_overlap_indices) == 0) {
    QC_fail_variants_no_hap_overlap <- candidate_variant_granges
    QC_fail_variants_no_hap_overlap$QC_fail_step <- "in_small_haplotype_block"
    QC_fail_variants <- c(QC_fail_variants, QC_fail_variants_no_hap_overlap)
    QC_fail_variants$duoNovo_classification <- "failed_QC"
    QC_fail_variants$parent_origin <- NA
    warning(paste0("No candidate variants of ", phasing_orientation, " phasing orientation passed QC."))
    return(QC_fail_variants)
  }
  if (length(hap_overlap_indices) < length(candidate_variant_granges)){
    QC_fail_variants_no_hap_overlap <- candidate_variant_granges[-hap_overlap_indices]
    QC_fail_variants_no_hap_overlap$QC_fail_step <- "in_small_haplotype_block"
    QC_fail_variants <- c(QC_fail_variants, QC_fail_variants_no_hap_overlap)
    candidate_variant_granges <- candidate_variant_granges[hap_overlap_indices]
  }
  
  #now obtain those that fall too close within boundaries of haplotype blocks
  overlaps <- findOverlaps(candidate_variant_granges, haplotype_boundary_coordinate_granges - boundary_cutoff)
  if (length(overlaps) == 0) {
    QC_fail_variants_boundary_overlap <- candidate_variant_granges
    QC_fail_variants_boundary_overlap$QC_fail_step <- "in_haplotype_block_boundary"
    QC_fail_variants <- c(QC_fail_variants, QC_fail_variants_boundary_overlap)
    QC_fail_variants$duoNovo_classification <- "failed_QC"
    QC_fail_variants$parent_origin <- NA
    warning(paste0("No candidate variants of ", phasing_orientation, " phasing orientation passed QC."))
    return(QC_fail_variants)
  }
  no_boundary_overlap_indices <- unique(queryHits(overlaps))
  if (length(no_boundary_overlap_indices) < length(candidate_variant_granges)){
    QC_fail_variants_boundary_overlap <- candidate_variant_granges[-no_boundary_overlap_indices]
    QC_fail_variants_boundary_overlap$QC_fail_step <- "in_haplotype_block_boundary"
    QC_fail_variants <- c(QC_fail_variants, QC_fail_variants_boundary_overlap)
  }
  ###QC steps end here
  
  ###now proceed to variant classification
  overlapping_indices <- split(queryHits(overlaps), subjectHits(overlaps))
  
  haplotypes <- vector("list", length(overlapping_indices))
  names(haplotypes) <- names(overlapping_indices)
  indices <- as.numeric(names(overlapping_indices))
  
  haplotype_granges_no_denovo <- haplotype_granges[-unique(queryHits(findOverlaps(haplotype_granges,
                                                                                  candidate_variant_granges)))]
  
  #counts_het_hom <- rep(NA, length(indices))
  #counts_het_het <- rep(NA, length(indices))
  #counts_hom_het <- rep(NA, length(indices))
  het <- c("0/1", "1/0", "0|1", "1|0", "0|2", "2|0", "1|2", "2|1",
           "0|3", "3|0", "1|3", "3|1", "2|3", "3|2", "0|4", "4|0", "1|4", "4|1", "2|4", "4|2", "3|4", "4|3")
  hom <- c("0/0", "1/1", "2/2", "3/3", "4/4")
  
  selected_granges <- haplotype_boundary_coordinate_granges[indices]
  overlap_results <- findOverlaps(haplotype_granges_no_denovo, selected_granges)
  
  # Directly use subjectHits to access haplotype blocks with variants
  with_variants <- unique(subjectHits(overlap_results))
  
  # Precompute query hits split by subject hits
  overlapping_map <- split(queryHits(overlap_results), subjectHits(overlap_results))
  
  hap11_all <- haplotype_granges_no_denovo$hap11 #first haplotype proband
  hap12_all <- haplotype_granges_no_denovo$hap12 #first haplotype sibling
  hap13_all <- haplotype_granges_no_denovo$hap13 #first haplotype mother
  hap21_all <- haplotype_granges_no_denovo$hap21 #second haplotype proband
  hap22_all <- haplotype_granges_no_denovo$hap22 #second haplotype sibling
  hap23_all <- haplotype_granges_no_denovo$hap23 #second haplotype mother
  
  #is_het1_all <- haplotype_granges_no_denovo$phasing1 %in% het
  #is_het2_all <- haplotype_granges_no_denovo$phasing2 %in% het
  #is_het3_all <- haplotype_granges_no_denovo$phasing3 %in% het
  #is_hom1_all <- haplotype_granges_no_denovo$phasing1 %in% hom
  #is_hom2_all <- haplotype_granges_no_denovo$phasing2 %in% hom
  #is_hom3_all <- haplotype_granges_no_denovo$phasing3 %in% hom
  
  for (i in with_variants) {
    variant_indices <- overlapping_map[[as.character(i)]]
    
    haplotypes[[i]] <- cbind(
      hap11 = hap11_all[variant_indices],
      hap12 = hap12_all[variant_indices],
      hap13 = hap13_all[variant_indices],
      hap21 = hap21_all[variant_indices], 
      hap22 = hap22_all[variant_indices],
      hap23 = hap23_all[variant_indices]
    )
    
    # Compute logical indices once and reuse them
    #counts_het_hom1[i] <- sum(is_het1_all[variant_indices] & is_hom2_all[variant_indices], na.rm = TRUE)
    #counts_het_het1[i] <- sum(is_het1_all[variant_indices] & is_het2_all[variant_indices], na.rm = TRUE)
    #counts_hom_het1[i] <- sum(is_hom1_all[variant_indices] & is_het2_all[variant_indices], na.rm = TRUE)
    
    #counts_het_hom2[i] <- sum(is_het1_all[variant_indices] & is_hom3_all[variant_indices], na.rm = TRUE)
    #counts_het_het2[i] <- sum(is_het1_all[variant_indices] & is_het3_all[variant_indices], na.rm = TRUE)
    #counts_hom_het2[i] <- sum(is_hom1_all[variant_indices] & is_het3_all[variant_indices], na.rm = TRUE)
  }
  
  hamming_distance_mat <- sapply(haplotypes, function(xx)
    c(sum(xx[, "hap11"] != xx[, "hap12"], na.rm = TRUE),
      sum(xx[, "hap11"] != xx[, "hap22"], na.rm = TRUE),
      sum(xx[, "hap21"] != xx[, "hap12"], na.rm = TRUE),
      sum(xx[, "hap21"] != xx[, "hap22"], na.rm = TRUE), 
      sum(xx[, "hap11"] != xx[, "hap13"], na.rm = TRUE),
      sum(xx[, "hap11"] != xx[, "hap23"], na.rm = TRUE),
      sum(xx[, "hap21"] != xx[, "hap13"], na.rm = TRUE),
      sum(xx[, "hap21"] != xx[, "hap23"], na.rm = TRUE))
  )
  
  hamming_distance_mins_1vs1_proband_sib <- colMins(hamming_distance_mat[1, , drop = FALSE])
  hamming_distance_mins_1vs2_proband_sib <- colMins(hamming_distance_mat[2, , drop = FALSE])
  hamming_distance_mins_2vs1_proband_sib <- colMins(hamming_distance_mat[3, , drop = FALSE])
  hamming_distance_mins_2vs2_proband_sib <- colMins(hamming_distance_mat[4, , drop = FALSE])
  
  hamming_distance_mins_1vs1_proband_mother <- colMins(hamming_distance_mat[5, , drop = FALSE])
  hamming_distance_mins_1vs2_proband_mother <- colMins(hamming_distance_mat[6, , drop = FALSE])
  hamming_distance_mins_2vs1_proband_mother <- colMins(hamming_distance_mat[7, , drop = FALSE])
  hamming_distance_mins_2vs2_proband_mother <- colMins(hamming_distance_mat[8, , drop = FALSE])
  
  #hamming_distance_mins_hap1_proband_sib <- colMins(hamming_distance_mat[1:2, , drop = FALSE])
  #hamming_distance_mins_hap2_proband_sib <- colMins(hamming_distance_mat[3:4, , drop = FALSE])
  #hamming_distance_mins_hap1_proband_mother <- colMins(hamming_distance_mat[5:6, , drop = FALSE])
  #hamming_distance_mins_hap2_proband_mother <- colMins(hamming_distance_mat[7:8, , drop = FALSE])
  
  all_columns <- 1:dim(hamming_distance_mat)[2]
  
  similar_11_sib <- which(hamming_distance_mins_1vs1_proband_sib == 0)
  similar_11_mother <- which(hamming_distance_mins_1vs1_proband_mother == 0)
  similar_12_sib <- which(hamming_distance_mins_1vs2_proband_sib == 0)
  similar_12_mother <- which(hamming_distance_mins_1vs2_proband_mother == 0)
  similar_21_sib <- which(hamming_distance_mins_2vs1_proband_sib == 0)
  similar_21_mother <- which(hamming_distance_mins_2vs1_proband_mother == 0)
  similar_22_sib <- which(hamming_distance_mins_2vs2_proband_sib == 0)
  similar_22_mother <- which(hamming_distance_mins_2vs2_proband_mother == 0)
  
  dissimilar_11_sib <- which(hamming_distance_mins_1vs1_proband_sib > distance_cutoff)
  dissimilar_11_mother <- which(hamming_distance_mins_1vs1_proband_mother > distance_cutoff)
  dissimilar_12_sib <- which(hamming_distance_mins_1vs2_proband_sib > distance_cutoff)
  dissimilar_12_mother <- which(hamming_distance_mins_1vs2_proband_mother > distance_cutoff)
  dissimilar_21_sib <- which(hamming_distance_mins_2vs1_proband_sib > distance_cutoff)
  dissimilar_21_mother <- which(hamming_distance_mins_2vs1_proband_mother > distance_cutoff)
  dissimilar_22_sib <- which(hamming_distance_mins_2vs2_proband_sib > distance_cutoff)
  dissimilar_22_mother <- which(hamming_distance_mins_2vs2_proband_mother > distance_cutoff)
  
  not_perfect_1vs1_sib    <- which(hamming_distance_mins_1vs1_proband_sib    > 0)
  not_perfect_1vs2_sib    <- which(hamming_distance_mins_1vs2_proband_sib    > 0)
  not_perfect_2vs1_sib    <- which(hamming_distance_mins_2vs1_proband_sib    > 0)
  not_perfect_2vs2_sib    <- which(hamming_distance_mins_2vs2_proband_sib    > 0)
  not_perfect_1vs1_mother <- which(hamming_distance_mins_1vs1_proband_mother > 0)
  not_perfect_1vs2_mother <- which(hamming_distance_mins_1vs2_proband_mother > 0)
  not_perfect_2vs1_mother <- which(hamming_distance_mins_2vs1_proband_mother > 0)
  not_perfect_2vs2_mother <- which(hamming_distance_mins_2vs2_proband_mother > 0)
  
  
  ## ===================================================================
  ##  clean-inheritance patterns
  ## ===================================================================
  
  ###----------Case 1: haplotype 1 from sib (hap1), haplotype 2 from mom
  clean_inheritance_hap1vs1_case1 <- Reduce(intersect, list(
    similar_11_sib,
    dissimilar_11_mother, dissimilar_12_mother,
    not_perfect_1vs2_sib,
    dissimilar_21_sib,  dissimilar_22_sib,
    similar_21_mother,
    not_perfect_2vs2_mother
  ))
  
  clean_inheritance_hap1vs1_case2 <- Reduce(intersect, list(
    similar_11_sib,
    dissimilar_11_mother, dissimilar_12_mother,
    not_perfect_1vs2_sib,
    dissimilar_21_sib,  dissimilar_22_sib,
    similar_22_mother,
    not_perfect_2vs1_mother
  ))
  
  ###----------Case 2: haplotype 1 from sib (hap2), haplotype 2 from mom
  clean_inheritance_hap1vs2_case1 <- Reduce(intersect, list(
    similar_12_sib,
    dissimilar_11_mother, dissimilar_12_mother,
    not_perfect_1vs1_sib,
    dissimilar_21_sib,  dissimilar_22_sib,
    similar_21_mother,
    not_perfect_2vs2_mother
  ))
  
  clean_inheritance_hap1vs2_case2 <- Reduce(intersect, list(
    similar_12_sib,
    dissimilar_11_mother, dissimilar_12_mother,
    not_perfect_1vs1_sib,
    dissimilar_21_sib,  dissimilar_22_sib,
    similar_22_mother,
    not_perfect_2vs1_mother
  ))
  
  ###----------Case 3: haplotype 2 from sib (hap1), haplotype 1 from mom
  clean_inheritance_hap2vs1_case1 <- Reduce(intersect, list(
    similar_21_sib,
    dissimilar_21_mother, dissimilar_22_mother,
    not_perfect_2vs2_sib,
    dissimilar_11_sib,  dissimilar_12_sib,
    similar_11_mother,
    not_perfect_1vs2_mother
  ))
  
  clean_inheritance_hap2vs1_case2 <- Reduce(intersect, list(
    similar_21_sib,
    dissimilar_21_mother, dissimilar_22_mother,
    not_perfect_2vs2_sib,
    dissimilar_11_sib,  dissimilar_12_sib,
    similar_12_mother,
    not_perfect_1vs1_mother
  ))
  
  ###----------Case 4: haplotype 2 from sib (hap2), haplotype 1 from mom
  clean_inheritance_hap2vs2_case1 <- Reduce(intersect, list(
    similar_22_sib,
    dissimilar_21_mother, dissimilar_22_mother,
    not_perfect_2vs1_sib,
    dissimilar_11_sib,  dissimilar_12_sib,
    similar_11_mother,
    not_perfect_1vs2_mother
  ))
  
  clean_inheritance_hap2vs2_case2 <- Reduce(intersect, list(
    similar_22_sib,
    dissimilar_21_mother, dissimilar_22_mother,
    not_perfect_2vs1_sib,
    dissimilar_11_sib,  dissimilar_12_sib,
    similar_12_mother,
    not_perfect_1vs1_mother
  ))
  
  ## -------------------------------------------------------------------
  ##  Collapse the two “case” vectors per paternal/maternal combination
  ## -------------------------------------------------------------------
  clean_inheritance_hap1vs1 <- union(clean_inheritance_hap1vs1_case1,
                                     clean_inheritance_hap1vs1_case2)
  
  clean_inheritance_hap1vs2 <- union(clean_inheritance_hap1vs2_case1,
                                     clean_inheritance_hap1vs2_case2)
  
  clean_inheritance_hap2vs1 <- union(clean_inheritance_hap2vs1_case1,
                                     clean_inheritance_hap2vs1_case2)
  
  clean_inheritance_hap2vs2 <- union(clean_inheritance_hap2vs2_case1,
                                     clean_inheritance_hap2vs2_case2)
  
  ### combine all
  clean_inheritance_all <- Reduce(
    union,
    list(clean_inheritance_hap1vs1,
         clean_inheritance_hap1vs2,
         clean_inheritance_hap2vs1,
         clean_inheritance_hap2vs2)
  )
  uncertain_inheritance <- all_columns
  if (length(clean_inheritance_all) > 0) {
    uncertain_inheritance <- all_columns[-clean_inheritance_all]
  }
  
  if (length(clean_inheritance_hap1vs1) > 0){
    hap11_variants_by_hap_block <- overlapping_indices[clean_inheritance_hap1vs1]
    hap11_inherited <- unlist(hap11_variants_by_hap_block)
  } else {
    hap11_inherited <- NULL
  }
  
  if (length(clean_inheritance_hap1vs2) > 0){
    hap12_variants_by_hap_block <- overlapping_indices[clean_inheritance_hap1vs2]
    hap12_inherited <- unlist(hap12_variants_by_hap_block)
  } else {
    hap12_inherited <- NULL
  }
  
  if (length(clean_inheritance_hap2vs1) > 0){
    hap21_variants_by_hap_block <- overlapping_indices[clean_inheritance_hap2vs1]
    hap21_inherited <- unlist(hap21_variants_by_hap_block)
  } else {
    hap21_inherited <- NULL
  }
  
  if (length(clean_inheritance_hap2vs2) > 0){
    hap22_variants_by_hap_block <- overlapping_indices[clean_inheritance_hap2vs2]
    hap22_inherited <- unlist(hap22_variants_by_hap_block)
  } else {
    hap22_inherited <- NULL
  }
  
  if (length(uncertain_inheritance) > 0){
    uncertain <- unlist(overlapping_indices[uncertain_inheritance])
  } else {
    uncertain <- NULL
  }
  
  if (phasing_orientation == "left") {
    if (!is.null(hap11_inherited)){
      de_novo11 <- candidate_variant_granges[hap11_inherited]
      de_novo11$parent_origin <- "father"
    } else {
      de_novo11 <- GRanges()
    }
    
    if (!is.null(hap12_inherited)){
      de_novo12 <- candidate_variant_granges[hap12_inherited]
      de_novo12$parent_origin <- "father"
    } else {
      de_novo12 <- GRanges()
    }
    
    if (!is.null(hap21_inherited)){
      de_novo21 <- candidate_variant_granges[hap21_inherited]
      de_novo21$parent_origin <- "mother"
    } else {
      de_novo21 <- GRanges()
    }
    
    if (!is.null(hap22_inherited)){
      de_novo22 <- candidate_variant_granges[hap22_inherited]
      de_novo22$parent_origin <- "mother"
    } else {
      de_novo22 <- GRanges()
    }
    de_novo <- c(de_novo11, de_novo12, de_novo21, de_novo22)
    de_novo$duoNovo_classification <- "de_novo"
    
    if (!is.null(uncertain)){
      uncertain <- candidate_variant_granges[uncertain]
      uncertain$duoNovo_classification <- "uncertain"
      uncertain$parent_origin <- NA
    } else {
      uncertain <- GRanges()
    }
    
    output <- c(de_novo, uncertain)
  } else if (phasing_orientation == "right") {
    if (!is.null(hap21_inherited)){
      de_novo21 <- candidate_variant_granges[hap21_inherited]
      de_novo21$parent_origin <- "father"
    } else {
      de_novo21 <- GRanges()
    }
    
    if (!is.null(hap22_inherited)){
      de_novo22 <- candidate_variant_granges[hap22_inherited]
      de_novo22$parent_origin <- "father"
    } else {
      de_novo22 <- GRanges()
    }
    
    if (!is.null(hap11_inherited)){
      de_novo11 <- candidate_variant_granges[hap11_inherited]
      de_novo11$parent_origin <- "mother"
    } else {
      de_novo11 <- GRanges()
    }
    
    if (!is.null(hap12_inherited)){
      de_novo12 <- candidate_variant_granges[hap12_inherited]
      de_novo12$parent_origin <- "mother"
    } else {
      de_novo12 <- GRanges()
    }
    de_novo <- c(de_novo11, de_novo12, de_novo21, de_novo22)
    de_novo$duoNovo_classification <- "de_novo"
    
    if (!is.null(uncertain)){
      uncertain <- candidate_variant_granges[uncertain]
      uncertain$duoNovo_classification <- "uncertain"
      uncertain$parent_origin <- NA
    } else {
      uncertain <- GRanges()
    }
    
    output <- c(de_novo, uncertain)
  }
  
  if (length(QC_fail_variants) > 0){
    QC_fail_variants$duoNovo_classification <- "failed_QC"
    QC_fail_variants$parent_origin <- NA
  }
  
  output$QC_fail_step <- NA
  combined_output <- c(output, QC_fail_variants)
  combined_output
}

###
###
getHaplotypesForTrio <- function(phased_variant_granges_trio) {
  granges <- phased_variant_granges_trio
  granges$hap11 <- NA
  granges$hap21 <- NA
  
  ### Extracting haplotypes for the proband
  granges$hap11[grep("0\\|", granges$phasing1)] <- 0
  granges$hap11[grep("1\\|", granges$phasing1)] <- 1
  granges$hap11[grep("2\\|", granges$phasing1)] <- 2
  granges$hap11[grep("3\\|", granges$phasing1)] <- 3
  granges$hap11[grep("4\\|", granges$phasing1)] <- 4
  granges$hap11[grep("0/0", granges$phasing1)] <- 0
  granges$hap11[grep("1/1", granges$phasing1)] <- 1
  
  granges$hap21[grep("\\|0", granges$phasing1)] <- 0
  granges$hap21[grep("\\|1", granges$phasing1)] <- 1
  granges$hap21[grep("\\|2", granges$phasing1)] <- 2
  granges$hap21[grep("\\|3", granges$phasing1)] <- 3
  granges$hap21[grep("\\|4", granges$phasing1)] <- 4
  granges$hap21[grep("0/0", granges$phasing1)] <- 0
  granges$hap21[grep("1/1", granges$phasing1)] <- 1
  
  ### Extracting haplotypes for the sibling
  granges$hap12 <- NA
  granges$hap22 <- NA
  
  granges$hap12[grep("0\\|", granges$phasing2)] <- 0
  granges$hap12[grep("1\\|", granges$phasing2)] <- 1
  granges$hap12[grep("2\\|", granges$phasing2)] <- 2
  granges$hap12[grep("3\\|", granges$phasing2)] <- 3
  granges$hap12[grep("4\\|", granges$phasing2)] <- 4
  granges$hap12[grep("0/0", granges$phasing2)] <- 0
  granges$hap12[grep("1/1", granges$phasing2)] <- 1
  
  granges$hap22[grep("\\|0", granges$phasing2)] <- 0
  granges$hap22[grep("\\|1", granges$phasing2)] <- 1
  granges$hap22[grep("\\|2", granges$phasing2)] <- 2
  granges$hap22[grep("\\|3", granges$phasing2)] <- 3
  granges$hap22[grep("\\|4", granges$phasing2)] <- 4
  granges$hap22[grep("0/0", granges$phasing2)] <- 0
  granges$hap22[grep("1/1", granges$phasing2)] <- 1
  
  ### Extracting haplotypes for the mother
  granges$hap13 <- NA
  granges$hap23 <- NA
  
  granges$hap13[grep("0\\|", granges$phasing3)] <- 0
  granges$hap13[grep("1\\|", granges$phasing3)] <- 1
  granges$hap13[grep("2\\|", granges$phasing3)] <- 2
  granges$hap13[grep("3\\|", granges$phasing3)] <- 3
  granges$hap13[grep("4\\|", granges$phasing3)] <- 4
  granges$hap13[grep("0/0", granges$phasing3)] <- 0
  granges$hap13[grep("1/1", granges$phasing3)] <- 1
  
  granges$hap23[grep("\\|0", granges$phasing3)] <- 0
  granges$hap23[grep("\\|1", granges$phasing3)] <- 1
  granges$hap23[grep("\\|2", granges$phasing3)] <- 2
  granges$hap23[grep("\\|3", granges$phasing3)] <- 3
  granges$hap23[grep("\\|4", granges$phasing3)] <- 4
  granges$hap23[grep("0/0", granges$phasing3)] <- 0
  granges$hap23[grep("1/1", granges$phasing3)] <- 1
  
  granges
}



getHaplotypeBlockCoordinatesForTrio <- function(haplotype_granges) {
  hap_granges <- haplotype_granges
  
  # Split haplotype ranges by PS1 and calculate coordinates
  granges_by_ps1 <- split(hap_granges, hap_granges$PS1)
  start_coords <- start(granges_by_ps1)
  start_coords <- sapply(start_coords, min)
  end_coords <- end(granges_by_ps1)
  end_coords <- sapply(end_coords, max)
  all_chrs <- seqnames(granges_by_ps1)  
  all_chrs <- unlist(unique(all_chrs)) #Assuming all sequences are from the same chromosome in each list element
  
  ps1_boundaries <- GRanges(
    seqnames = all_chrs, ranges = IRanges(start = start_coords, end = end_coords),
    PS1 = names(granges_by_ps1)
  )
  
  # Split haplotype ranges by PS2 and calculate coordinates
  granges_by_ps2 <- split(hap_granges, hap_granges$PS2)
  
  start_coords <- start(granges_by_ps2)
  start_coords <- sapply(start_coords, min)
  end_coords <- end(granges_by_ps2)
  end_coords <- sapply(end_coords, max)
  all_chrs <- seqnames(granges_by_ps2)  
  all_chrs <- unlist(unique(all_chrs)) #Assuming all sequences are from the same chromosome in each list element
  
  ps2_boundaries <- GRanges(
    seqnames = all_chrs, ranges = IRanges(start = start_coords, end = end_coords),
    PS2 = names(granges_by_ps2)
  )
  
  # Split haplotype ranges by PS1 and calculate coordinates
  granges_by_ps3 <- split(hap_granges, hap_granges$PS3)
  start_coords <- start(granges_by_ps3)
  start_coords <- sapply(start_coords, min)
  end_coords <- end(granges_by_ps3)
  end_coords <- sapply(end_coords, max)
  all_chrs <- seqnames(granges_by_ps3)  
  all_chrs <- unlist(unique(all_chrs)) #Assuming all sequences are from the same chromosome in each list element
  
  ps3_boundaries <- GRanges(
    seqnames = all_chrs, ranges = IRanges(start = start_coords, end = end_coords),
    PS3 = names(granges_by_ps3)
  )
  
  
  ## ------------ pairwise intersection: PS1 ∩ PS2 ------------------
  ov12      <- findOverlaps(ps1_boundaries, ps2_boundaries, type = "any")
  inter12   <- pintersect(ps1_boundaries[queryHits(ov12)],
                          ps2_boundaries[subjectHits(ov12)],
                          drop.nohit.ranges = TRUE)
  
  ## ------------ add the third sample: (PS1 ∩ PS2) ∩ PS3 ----------
  ov123     <- findOverlaps(inter12, ps3_boundaries, type = "any")
  inter123  <- pintersect(inter12[queryHits(ov123)],
                          ps3_boundaries[subjectHits(ov123)],
                          drop.nohit.ranges = TRUE)
  
  ## ------------ assemble one row per 3-way shared segment ---------
  combined_PS <- inter123                              
  mcols(combined_PS) <- DataFrame(
    PS1 = ps1_boundaries$PS1[
      queryHits(ov12)[queryHits(ov123)]           ],
    PS2 = ps2_boundaries$PS2[
      subjectHits(ov12)[queryHits(ov123)]         ],
    PS3 = ps3_boundaries$PS3[
      subjectHits(ov123)                          ]
  )
  
  valid_ranges <- start(combined_PS) <= end(combined_PS)
  combined_PS <- combined_PS[valid_ranges]
  combined_PS
}

library(VariantAnnotation)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: falsepos_duoNovo_sib.R <run_directory>")
}
## set directory
current_dir <- args[1]
setwd(current_dir)


### run duoNovo (trio version) on the joint trio vcf (proband-sibling-mother and proband-sibling-father; samples have to be ordered this way in the vcf) 
message("getting de novo variants from PSM trio...")
surrogate_trio_vcf_filepath_psm <- list.files(pattern = "PSM\\.vcf.gz$")
sibling_trio_output_psm <- duoNovoSib(surrogate_trio_vcf_filepath_psm)

message("getting de novo variants from PSF trio...")
surrogate_trio_vcf_filepath_psf <- list.files(pattern = "PSF\\.vcf.gz$")
sibling_trio_output_psf <- duoNovoSib(surrogate_trio_vcf_filepath_psf)

# PSM
denovo_indices_sibling_trio_psm <- which(sibling_trio_output_psm$duoNovo_classification == "de_novo")
if (length(denovo_indices_sibling_trio_psm) > 0){
  denovo_sibling_psm <- sibling_trio_output_psm[denovo_indices_sibling_trio_psm]
  
  ### can optionally exclude de novos that are present in gnomAD
  # duoNovo_ranges$gnomad41_genome_AF <- as.numeric(unlist(duoNovo_ranges$gnomad41_genome_AF))
  # overlapping_indices <- unique(queryHits(findOverlaps(denovo_sibling, 
  #                                         duoNovo_ranges[which(duoNovo_ranges$gnomad41_genome_AF > 0)])))
  # if (length(overlapping_indices) > 0){
  # denovo_sibling <- denovo_sibling[-overlapping_indices]
  #}
  ###
} else {
  denovo_sibling_psm <- GRanges()
}

# PSF
denovo_indices_sibling_trio_psf <- which(sibling_trio_output_psf$duoNovo_classification == "de_novo")
if (length(denovo_indices_sibling_trio_psf) > 0){
  denovo_sibling_psf <- sibling_trio_output_psf[denovo_indices_sibling_trio_psf]
  
} else {
  denovo_sibling_psf <- GRanges()
}

### use full trio genotypes to assess false pos rate
### use duoNovo output from the mother-proband duo to separate candidate variants based on classification
duoNovo_output_filepath_pm <- list.files(pattern = "PM\\.duonovo\\.annovar\\.addedParent\\.dnm2.rda$")
load(file = duoNovo_output_filepath_pm)
dn_granges <- dn_granges_pm

dn_granges <- dn_granges[which(dn_granges$parentValidation_depth >= 20 & 
                                 dn_granges$parentValidation_GQ >= 30)]
dn_granges <- dn_granges[!grepl("\\.", dn_granges$parentValidation_gt)]
overlaps <- findOverlaps(dn_granges, denovo_sibling_psm)
if (length(queryHits(overlaps)) > 0){
  variants_to_assess <- dn_granges[unique(queryHits(overlaps))] 
  
  assessed <- length(variants_to_assess)
  false_dn <- length(grep("1", variants_to_assess$parentValidation_gt))
  rate <- false_dn/assessed
  output_vec <- c(false_dn, assessed, rate)
  names(output_vec) <- c("false_dn", "assessed", "false_pos_rate")
} else {
  output_vec <- NA
}
output_vec_psm <- output_vec

### now same from father-proband duo 
duoNovo_output_filepath_pf <- list.files(pattern = "PF\\.duonovo\\.annovar\\.addedParent\\.dnm2.rda$")
load(file = duoNovo_output_filepath_pf)
dn_granges <- dn_granges_pf

dn_granges <- dn_granges[which(dn_granges$parentValidation_depth >= 20 & 
                                 dn_granges$parentValidation_GQ >= 30)]
dn_granges <- dn_granges[!grepl("\\.", dn_granges$parentValidation_gt)]
overlaps <- findOverlaps(dn_granges, denovo_sibling_psf)
if (length(queryHits(overlaps)) > 0){
  variants_to_assess <- dn_granges[unique(queryHits(overlaps))] 
  
  assessed <- length(variants_to_assess)
  false_dn <- length(grep("1", variants_to_assess$parentValidation_gt))
  rate <- false_dn/assessed
  output_vec <- c(false_dn, assessed, rate)
  names(output_vec) <- c("false_dn", "assessed", "false_pos_rate")
} else {
  output_vec <- NA
}
output_vec_psf <- output_vec


sample_id <- sub("\\..*$", "", duoNovo_output_filepath_pm)
save(output_vec_psm, output_vec_psf, 
     file = paste0(sample_id, "_sibling_trio_falsepos_rate.rda"))







