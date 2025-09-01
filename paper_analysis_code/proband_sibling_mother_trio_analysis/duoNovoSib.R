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




