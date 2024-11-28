#' Classify Candidate Variants Based on Haplotype Analysis
#'
#' This function classifies candidate variants as de novo or inherited based on the comparison
#' of phased haplotype blocks between the proband and the sequenced parent.
#'
#' @param candidate_variant_granges A \code{GRanges} object containing the candidate variants to be classified.
#' @param phasing_orientation A character value indicating the phasing orientation ("left" or "right").
#' @param haplotype_granges A \code{GRanges} object containing haplotype information of the sequenced parent.
#' @param haplotype_boundary_coordinate_granges A \code{GRanges} object containing the haplotype block boundaries.
#' @param boundary_cutoff A numeric value indicating the distance from the boundary of a haplotype block
#'   within which candidate variants are considered.
#' @param distance_cutoff A numeric value specifying the minimum Hamming distance required for supporting
#'   the classification of variants.
#'
#' @return A \code{GRanges} object containing classified candidate variants. The classified variants include
#'   the following columns:
#'   \itemize{
#'     \item \code{duoNovo_classification}: Classification status ("de novo", "on other parent haplotype", or "uncertain").
#'     \item \code{hamming_distance_other_parent_hap}: The minimum Hamming distance for the other parent's haplotype.
#'     \item \code{supporting_counts_het_hom}, \code{supporting_counts_het_het}, \code{supporting_counts_hom_het}:
#'       Counts representing variant pairs used to calculate Hamming distance.
#'   }
#' @examples
#' # classifyVariants(candidate_variant_granges = candidate_variants,
#' #                phasing_orientation = "left",
#' #                 haplotype_granges = haplotypes,
#' #                haplotype_boundary_coordinate_granges = haplotype_boundaries,
#' #                  boundary_cutoff = 2000,
#' #                  distance_cutoff = 40)
#'
#' @export
classifyVariants <- function(candidate_variant_granges, phasing_orientation = c("left", "right"),
                             haplotype_granges, haplotype_boundary_coordinate_granges,
                             boundary_cutoff, distance_cutoff){

  ###QC steps
  #first obtain variants that do not overlap any of the haplotype blocks
  hap_overlap_indices <- unique(queryHits(findOverlaps(candidate_variant_granges,
                                                       haplotype_boundary_coordinate_granges)))
  QC_fail_variants_no_hap_overlap <- candidate_variant_granges[-hap_overlap_indices]
  QC_fail_variants_no_hap_overlap$QC_fail_step <- "no_haplotype_block_overlap"
  QC_fail_variants <- c(QC_fail_variants, QC_fail_variants_no_hap_overlap)
  candidate_variant_granges <- candidate_variant_granges[hap_overlap_indices]

  #now obtain those that do not overlap any after filtering out small haplotype blocks
  haplotype_boundary_coordinate_granges <- haplotype_boundary_coordinate_granges[
    which(width(haplotype_boundary_coordinate_granges) > PS_width_cutoff)]

  hap_overlap_indices <- unique(queryHits(findOverlaps(candidate_variant_granges,
                                                       haplotype_boundary_coordinate_granges)))
  QC_fail_variants_no_hap_overlap <- candidate_variant_granges[-hap_overlap_indices]
  QC_fail_variants_no_hap_overlap$QC_fail_step <- "in_small_haplotype_block"
  QC_fail_variants <- c(QC_fail_variants, QC_fail_variants_no_hap_overlap)
  candidate_variant_granges <- candidate_variant_granges[hap_overlap_indices]

  #now obtain those that fall too close within boundaries of haplotype blocks
  overlaps <- findOverlaps(candidate_variant_granges, haplotype_boundary_coordinate_granges - boundary_cutoff)
  no_boundary_overlap_indices <- unique(queryHits(overlaps))
  QC_fail_variants_boundary_overlap <- candidate_variant_granges[-no_boundary_overlap_indices]
  QC_fail_variants_boundary_overlap$QC_fail_step <- "in_haplotype_block_boundary"
  QC_fail_variants <- c(QC_fail_variants, QC_fail_variants_boundary_overlap)
  ###QC steps end here

  ###now proceed to variant classification
  overlapping_indices <- split(queryHits(overlaps), subjectHits(overlaps))
  if (length(overlapping_indices) == 0) {
    warning("Candidate variants cannot be classified because they all fall outside haplotype block boundaries after QC.")
    return(candidate_variant_granges)
  }

  haplotypes <- vector("list", length(overlapping_indices))
  names(haplotypes) <- names(overlapping_indices)

  indices <- as.numeric(names(overlapping_indices))
  PS1_ids <- haplotype_boundary_coordinate_granges$PS1[indices]
  PS2_ids <- haplotype_boundary_coordinate_granges$PS2[indices]

  haplotype_granges_no_denovo <- haplotype_granges[-unique(queryHits(findOverlaps(haplotype_granges,
                                                                                  candidate_variant_granges)))]

  counts_het_hom <- rep(NA, length(seq_along(indices)))
  counts_het_het <- rep(NA, length(seq_along(indices)))
  counts_hom_het <- rep(NA, length(seq_along(indices)))
  het <- c("0/1", "1/0", "0|1", "1|0", "0|2", "2|0", "1|2", "2|1",
           "0|3", "3|0", "1|3", "3|1", "2|3", "3|2", "0|4", "4|0", "1|4", "4|1", "2|4", "4|2", "3|4", "4|3")
  hom <- c("0/0", "1/1", "2/2", "3/3", "4/4")

  for(i in seq_along(indices)) {
    variant_granges <- haplotype_granges_no_denovo[
      unique(queryHits(findOverlaps(haplotype_granges_no_denovo, haplotype_boundary_coordinate_granges[
        which(haplotype_boundary_coordinate_granges$PS1 == PS1_ids[i] &
                haplotype_boundary_coordinate_granges$PS2 == PS2_ids[i])])))]

    hap11 <- variant_granges$hap11
    hap12 <- variant_granges$hap12
    hap21 <- variant_granges$hap21
    hap22 <- variant_granges$hap22

    hap_mat <- matrix(c(hap11, hap12, hap21, hap22), nrow = length(hap11), ncol = 4, byrow = FALSE)
    colnames(hap_mat) <- c("hap11", "hap12", "hap21", "hap22")
    haplotypes[[i]] <- hap_mat

    # Count the types of variant pair comparisons for each haplotype pair
    counts_het_hom[i] <- sum((variant_granges$phasing1 %in% het) &
                                 (variant_granges$phasing2 %in% hom), na.rm = TRUE)
    counts_het_het[i] <- sum((variant_granges$phasing1 %in% het) &
                                 (variant_granges$phasing2 %in% het), na.rm = TRUE)

    counts_hom_het[i] <- sum((variant_granges$phasing1 %in% hom) &
                                 (variant_granges$phasing2 %in% het), na.rm = TRUE)
  }

  hamming_distance_mat <- sapply(haplotypes, function(xx)
    c(sum(xx[, "hap11"] != xx[, "hap12"], na.rm = TRUE),
      sum(xx[, "hap11"] != xx[, "hap22"], na.rm = TRUE),
      sum(xx[, "hap21"] != xx[, "hap12"], na.rm = TRUE),
      sum(xx[, "hap21"] != xx[, "hap22"], na.rm = TRUE))
  )

  hamming_distance_mins <- colMins(hamming_distance_mat)

  hamming_distance_mins_1vs1 <- colMins(hamming_distance_mat[1, , drop = FALSE])
  hamming_distance_mins_1vs2 <- colMins(hamming_distance_mat[2, , drop = FALSE])
  hamming_distance_mins_2vs1 <- colMins(hamming_distance_mat[3, , drop = FALSE])
  hamming_distance_mins_2vs2 <- colMins(hamming_distance_mat[4, , drop = FALSE])

  hamming_distance_mins_hap1 <- colMins(hamming_distance_mat[1:2, ])
  hamming_distance_mins_hap2 <- colMins(hamming_distance_mat[3:4, ])

  clean_inheritance_hap1vs1 <- which(hamming_distance_mins_1vs1 == 0 & hamming_distance_mins_1vs2 > 0 &
                                        hamming_distance_mins_hap2 > distance_cutoff)
  clean_inheritance_hap1vs2 <- which(hamming_distance_mins_1vs1 > 0 & hamming_distance_mins_1vs2 == 0 &
                                       hamming_distance_mins_hap2 > distance_cutoff)
  clean_inheritance_hap2vs1 <- which(hamming_distance_mins_2vs1 == 0 & hamming_distance_mins_2vs2 > 0 &
                                       hamming_distance_mins_hap1 > distance_cutoff)
  clean_inheritance_hap2vs2 <- which(hamming_distance_mins_2vs1 > 0 & hamming_distance_mins_2vs2 == 0 &
                                       hamming_distance_mins_hap1 > distance_cutoff)
  uncertain_inheritance <- which(hamming_distance_mins > 0)

  ###4 different cases corresponding to 4 different clean inheritance pairs
  #proband haplotype 1 vs parent haplotype 1
  if (length(clean_inheritance_hap1vs1) > 0){
    hap11_variants_by_hap_block <- overlapping_indices[clean_inheritance_hap1vs1]
    hap11_inherited <- unlist(hap11_variants_by_hap_block)
    hap11_hamming_distance_other_parent_hap <- rep(hamming_distance_mins_hap2[clean_inheritance_hap1vs1],
                                                   lengths(hap11_variants_by_hap_block))
    hap11_supporting_counts_het_hom <- rep(counts_het_hom[clean_inheritance_hap1vs1], lengths(hap11_variants_by_hap_block))
    hap11_supporting_counts_het_het <- rep(counts_het_het[clean_inheritance_hap1vs1], lengths(hap11_variants_by_hap_block))
    hap11_supporting_counts_hom_het <- rep(counts_hom_het[clean_inheritance_hap1vs1], lengths(hap11_variants_by_hap_block))
  } else {
    hap11_inherited <- NULL
  }

  #proband haplotype 1 vs parent haplotype 2
  if (length(clean_inheritance_hap1vs2) > 0){
    hap12_variants_by_hap_block <- overlapping_indices[clean_inheritance_hap1vs2]
    hap12_inherited <- unlist(hap12_variants_by_hap_block)
    hap12_hamming_distance_other_parent_hap <- rep(hamming_distance_mins_hap2[clean_inheritance_hap1vs2],
                                                   lengths(hap12_variants_by_hap_block))
    hap12_supporting_counts_het_hom <- rep(counts_het_hom[clean_inheritance_hap1vs2], lengths(hap12_variants_by_hap_block))
    hap12_supporting_counts_het_het <- rep(counts_het_het[clean_inheritance_hap1vs2], lengths(hap12_variants_by_hap_block))
    hap12_supporting_counts_hom_het <- rep(counts_hom_het[clean_inheritance_hap1vs2], lengths(hap12_variants_by_hap_block))
  } else {
    hap12_inherited <- NULL
  }

  #proband haplotype 2 vs parent haplotype 1
  if (length(clean_inheritance_hap2vs1) > 0){
    hap21_variants_by_hap_block <- overlapping_indices[clean_inheritance_hap2vs1]
    hap21_inherited <- unlist(hap21_variants_by_hap_block)
    hap21_hamming_distance_other_parent_hap <- rep(hamming_distance_mins_hap1[clean_inheritance_hap2vs1],
                                                   lengths(hap21_variants_by_hap_block))
    hap21_supporting_counts_het_hom <- rep(counts_het_hom[clean_inheritance_hap2vs1], lengths(hap21_variants_by_hap_block))
    hap21_supporting_counts_het_het <- rep(counts_het_het[clean_inheritance_hap2vs1], lengths(hap21_variants_by_hap_block))
    hap21_supporting_counts_hom_het <- rep(counts_hom_het[clean_inheritance_hap2vs1], lengths(hap21_variants_by_hap_block))
  } else {
    hap21_inherited <- NULL
  }

  #proband haplotype 2 vs parent haplotype 2
  if (length(clean_inheritance_hap2vs2) > 0){
    hap22_variants_by_hap_block <- overlapping_indices[clean_inheritance_hap2vs2]
    hap22_inherited <- unlist(hap22_variants_by_hap_block)
    hap22_hamming_distance_other_parent_hap <- rep(hamming_distance_mins_hap1[clean_inheritance_hap2vs2],
                                                   lengths(hap22_variants_by_hap_block))
    hap22_supporting_counts_het_hom <- rep(counts_het_hom[clean_inheritance_hap2vs2], lengths(hap22_variants_by_hap_block))
    hap22_supporting_counts_het_het <- rep(counts_het_het[clean_inheritance_hap2vs2], lengths(hap22_variants_by_hap_block))
    hap22_supporting_counts_hom_het <- rep(counts_hom_het[clean_inheritance_hap2vs2], lengths(hap22_variants_by_hap_block))
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
      de_novo11$duoNovo_classification <- "de_novo"
      #add columns for supporting evidence of de novo classification
      de_novo11$supporting_hamming_distance <- hap11_hamming_distance_other_parent_hap
      de_novo11$supporting_counts_het_hom <- hap11_supporting_counts_het_hom
      de_novo11$supporting_counts_het_het <- hap11_supporting_counts_het_het
      de_novo11$supporting_counts_hom_het <- hap11_supporting_counts_hom_het
    } else {
      de_novo11 <- GRanges()
    }

    if (!is.null(hap12_inherited)){
      de_novo12 <- candidate_variant_granges[hap12_inherited]
      de_novo12$duoNovo_classification <- "de_novo"
      #add columns for supporting evidence of de novo classification
      de_novo12$supporting_hamming_distance <- hap12_hamming_distance_other_parent_hap
      de_novo12$supporting_counts_het_hom <- hap12_supporting_counts_het_hom
      de_novo12$supporting_counts_het_het <- hap12_supporting_counts_het_het
      de_novo12$supporting_counts_hom_het <- hap12_supporting_counts_hom_het
    } else {
      de_novo12 <- GRanges()
    }
    de_novo <- c(de_novo11, de_novo12)

    if (!is.null(hap21_inherited)){
      not_de_novo21 <- candidate_variant_granges[hap21_inherited]
      not_de_novo21$duoNovo_classification <- "on_other_parent_haplotype"
      #add columns for supporting evidence of de novo classification
      not_de_novo21$supporting_hamming_distance <- hap21_hamming_distance_other_parent_hap
      not_de_novo21$supporting_counts_het_hom <- hap21_supporting_counts_het_hom
      not_de_novo21$supporting_counts_het_het <- hap21_supporting_counts_het_het
      not_de_novo21$supporting_counts_hom_het <- hap21_supporting_counts_hom_het
    } else {
      not_de_novo21 <- GRanges()
    }

    if (!is.null(hap22_inherited)){
      not_de_novo22 <- candidate_variant_granges[hap22_inherited]
      not_de_novo22$duoNovo_classification <- "on_other_parent_haplotype"
      #add columns for supporting evidence of de novo classification
      not_de_novo22$supporting_hamming_distance <- hap22_hamming_distance_other_parent_hap
      not_de_novo22$supporting_counts_het_hom <- hap22_supporting_counts_het_hom
      not_de_novo22$supporting_counts_het_het <- hap22_supporting_counts_het_het
      not_de_novo22$supporting_counts_hom_het <- hap22_supporting_counts_hom_het
    } else {
      not_de_novo22 <- GRanges()
    }
    not_de_novo <- c(not_de_novo21, not_de_novo22)

    if (!is.null(uncertain)){
      uncertain <- candidate_variant_granges[uncertain]
      uncertain$duoNovo_classification <- "uncertain"
      uncertain$supporting_hamming_distance <- NA
      uncertain$supporting_counts_het_hom <- NA
      uncertain$supporting_counts_het_het <- NA
      uncertain$supporting_counts_hom_het <- NA
    } else {
      uncertain <- GRanges()
    }

    output <- c(de_novo, not_de_novo, uncertain)
  } else if (phasing_orientation == "right") {
    if (!is.null(hap21_inherited)){
      de_novo21 <- candidate_variant_granges[hap21_inherited]
      de_novo21$duoNovo_classification <- "de_novo"
      #add columns for supporting evidence of de novo classification
      de_novo21$supporting_hamming_distance <- hap21_hamming_distance_other_parent_hap
      de_novo21$supporting_counts_het_hom <- hap21_supporting_counts_het_hom
      de_novo21$supporting_counts_het_het <- hap21_supporting_counts_het_het
      de_novo21$supporting_counts_hom_het <- hap21_supporting_counts_hom_het
    } else {
      de_novo21 <- GRanges()
    }

    if (!is.null(hap22_inherited)){
      de_novo22 <- candidate_variant_granges[hap22_inherited]
      de_novo22$duoNovo_classification <- "de_novo"
      #add columns for supporting evidence of de novo classification
      de_novo22$supporting_hamming_distance <- hap22_hamming_distance_other_parent_hap
      de_novo22$supporting_counts_het_hom <- hap22_supporting_counts_het_hom
      de_novo22$supporting_counts_het_het <- hap22_supporting_counts_het_het
      de_novo22$supporting_counts_hom_het <- hap22_supporting_counts_hom_het
    } else {
      de_novo22 <- GRanges()
    }
    de_novo <- c(de_novo21, de_novo22)

    if (!is.null(hap11_inherited)){
      not_de_novo11 <- candidate_variant_granges[hap11_inherited]
      not_de_novo11$duoNovo_classification <- "on_other_parent_haplotype"
      #add columns for supporting evidence of de novo classification
      not_de_novo11$supporting_hamming_distance <- hap11_hamming_distance_other_parent_hap
      not_de_novo11$supporting_counts_het_hom <- hap11_supporting_counts_het_hom
      not_de_novo11$supporting_counts_het_het <- hap11_supporting_counts_het_het
      not_de_novo11$supporting_counts_hom_het <- hap11_supporting_counts_hom_het
    } else {
      not_de_novo11 <- GRanges()
    }

    if (!is.null(hap12_inherited)){
      not_de_novo12 <- candidate_variant_granges[hap12_inherited]
      not_de_novo12$duoNovo_classification <- "on_other_parent_haplotype"
      #add columns for supporting evidence of de novo classification
      not_de_novo12$supporting_hamming_distance <- hap12_hamming_distance_other_parent_hap
      not_de_novo12$supporting_counts_het_hom <- hap12_supporting_counts_het_hom
      not_de_novo12$supporting_counts_het_het <- hap12_supporting_counts_het_het
      not_de_novo12$supporting_counts_hom_het <- hap12_supporting_counts_hom_het
    } else {
      not_de_novo12 <- GRanges()
    }
    not_de_novo <- c(not_de_novo11, not_de_novo12)

    if (!is.null(uncertain)){
      uncertain <- candidate_variant_granges[uncertain]
      uncertain$duoNovo_classification <- "uncertain"
      uncertain$supporting_hamming_distance <- NA
      uncertain$supporting_counts_het_hom <- NA
      uncertain$supporting_counts_het_het <- NA
      uncertain$supporting_counts_hom_het <- NA
    } else {
      uncertain <- GRanges()
    }

    output <- c(de_novo, not_de_novo, uncertain)
  }

  QC_fail_variants$duoNovo_classification <- "failed_QC"
  QC_fail_variants$supporting_hamming_distance <- NA
  QC_fail_variants$supporting_counts_het_hom <- NA
  QC_fail_variants$supporting_counts_het_het <- NA
  QC_fail_variants$supporting_counts_hom_het <- NA

  output$QC_fail_step <- NA
  combined_output <- c(output, QC_fail_variants)
  combined_output
}


