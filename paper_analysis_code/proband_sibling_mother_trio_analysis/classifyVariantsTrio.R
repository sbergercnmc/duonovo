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









