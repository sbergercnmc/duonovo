#' Classify Candidate Variants Based on Haplotype Analysis
#'
#' This function classifies candidate variants as de novo or inherited based on the comparison
#' of phased haplotype blocks between the proband and the sequenced parent.
#'
#' @param candidate_variant_granges A \code{GRanges} object containing the candidate variants to be classified.
#' @param phasing_orientation A character value indicating the phasing orientation ("left" or "right").
#'   Default is c("left", "right"), specifying that both orientations can be used for classification.
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
#' # Example usage:
#' classifyVariants(candidate_variant_granges = candidate_variants,
#'                  phasing_orientation = "left",
#'                  haplotype_granges = haplotypes,
#'                  haplotype_boundary_coordinate_granges = haplotype_boundaries,
#'                  boundary_cutoff = 2000,
#'                  distance_cutoff = 40)
#'
#' @export
classifyVariants <- function(candidate_variant_granges, phasing_orientation = c("left", "right"),
                             haplotype_granges, haplotype_boundary_coordinate_granges,
                             boundary_cutoff, distance_cutoff){

  overlaps <- findOverlaps(candidate_variant_granges, haplotype_boundary_coordinate_granges - boundary_cutoff)
  overlapping_indices <- split(queryHits(overlaps), subjectHits(overlaps))
  if (length(overlapping_indices) == 0) {
    warning("Candidate variants cannot be classified because they fall outside haplotype boundaries.")
    return(NULL) # or return an empty GRangesList with appropriate names
  }

  haplotypes <- vector("list", length(overlapping_indices))
  names(haplotypes) <- names(overlapping_indices)

  indices <- as.numeric(names(overlapping_indices))
  PS1_ids <- haplotype_boundary_coordinate_granges$PS1[indices]
  PS2_ids <- haplotype_boundary_coordinate_granges$PS2[indices]

  haplotype_granges_no_denovo <- haplotype_granges[-unique(queryHits(findOverlaps(haplotype_granges,
                                                                                  candidate_variant_granges)))]

  counts_het_hom <- list()
  counts_het_het <- list()
  counts_hom_het <- list()
  het <- c("1|0", "0|1")
  hom <- c("0/0", "1/1")

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
    counts_het_hom[[i]] <- c(
      sum(hap_mat[, "hap11"] %in% het & hap_mat[, "hap12"] %in% hom, na.rm = TRUE),
      sum(hap_mat[, "hap11"] %in% het & hap_mat[, "hap22"] %in% hom, na.rm = TRUE),
      sum(hap_mat[, "hap21"] %in% het & hap_mat[, "hap12"] %in% hom, na.rm = TRUE),
      sum(hap_mat[, "hap21"] %in% het & hap_mat[, "hap22"] %in% hom, na.rm = TRUE)
    )

    counts_het_het[[i]] <- c(
      sum(hap_mat[, "hap11"] %in% het & hap_mat[, "hap12"] %in% het, na.rm = TRUE),
      sum(hap_mat[, "hap11"] %in% het & hap_mat[, "hap22"] %in% het, na.rm = TRUE),
      sum(hap_mat[, "hap21"] %in% het & hap_mat[, "hap12"] %in% het, na.rm = TRUE),
      sum(hap_mat[, "hap21"] %in% het & hap_mat[, "hap22"] %in% het, na.rm = TRUE)
    )

    counts_hom_het[[i]] <- c(
      sum(hap_mat[, "hap11"] %in% hom & hap_mat[, "hap12"] %in% het, na.rm = TRUE),
      sum(hap_mat[, "hap11"] %in% hom & hap_mat[, "hap22"] %in% het, na.rm = TRUE),
      sum(hap_mat[, "hap21"] %in% hom & hap_mat[, "hap12"] %in% het, na.rm = TRUE),
      sum(hap_mat[, "hap21"] %in% hom & hap_mat[, "hap22"] %in% het, na.rm = TRUE)
    )
  }

  hamming_distance_mat <- sapply(haplotypes, function(xx)
    c(sum(xx[, "hap11"] != xx[, "hap12"], na.rm = TRUE),
      sum(xx[, "hap11"] != xx[, "hap22"], na.rm = TRUE),
      sum(xx[, "hap21"] != xx[, "hap12"], na.rm = TRUE),
      sum(xx[, "hap21"] != xx[, "hap22"], na.rm = TRUE))
  )

  hamming_distance_mins <- colMins(hamming_distance_mat)

  hamming_distance_mins_hap1 <- colMins(hamming_distance_mat[, 1:2])
  hamming_distance_mins_hap2 <- colMins(hamming_distance_mat[, 3:4])

  clean_inheritance_hap1 <- which(hamming_distance_mins_hap1 == 0 & hamming_distance_mins_hap2 > distance_cutoff)
  clean_inheritance_hap2 <- which(hamming_distance_mins_hap2 == 0 & hamming_distance_mins_hap1 > distance_cutoff)

  hap1_variants_by_hap_block <- overlapping_indices[clean_inheritance_hap1]
  hap1_inherited <- unlist(hap1_variants_by_hap_block)
  hap1_hamming_distance_other_parent_hap <- rep(hamming_distance_mins_hap2[clean_inheritance_hap1],
                                                lengths(hap1_variants_by_hap_block))
  hap1_supporting_counts_het_hom <- rep(counts_het_hom[clean_inheritance_hap1], lengths(hap1_variants_by_hap_block))
  hap1_supporting_counts_het_het <- rep(counts_het_het[clean_inheritance_hap1], lengths(hap1_variants_by_hap_block))
  hap1_supporting_counts_hom_het <- rep(counts_hom_het[clean_inheritance_hap1], lengths(hap1_variants_by_hap_block))


  hap2_variants_by_hap_block <- overlapping_indices[clean_inheritance_hap2]
  hap2_inherited <- unlist(hap2_variants_by_hap_block)
  hap2_hamming_distance_other_parent_hap <- rep(hamming_distance_mins_hap1[clean_inheritance_hap2],
                                                lengths(hap2_variants_by_hap_block))
  hap2_supporting_counts_het_hom <- rep(counts_het_hom[clean_inheritance_hap2], lengths(hap2_variants_by_hap_block))
  hap2_supporting_counts_het_het <- rep(counts_het_het[clean_inheritance_hap2], lengths(hap2_variants_by_hap_block))
  hap2_supporting_counts_hom_het <- rep(counts_hom_het[clean_inheritance_hap2], lengths(hap2_variants_by_hap_block))


  uncertain <- unlist(overlapping_indices[which(hamming_distance_mins > 0)])

  if (phasing_orientation == "left") {
    de_novo <- candidate_variant_granges[hap1_inherited]
    de_novo$duoNovo_classification <- "de novo"
    de_novo$hamming_distance_other_parent_hap <- hap1_hamming_distance_other_parent_hap

    de_novo$supporting_counts_het_hom <- hap1_supporting_counts_het_hom
    de_novo$supporting_counts_het_het <- hap1_supporting_counts_het_het
    de_novo$supporting_counts_hom_het <- hap1_supporting_counts_hom_het


    not_de_novo <- candidate_variant_granges[hap2_inherited]
    not_de_novo$duoNovo_classification <- "on other parent haplotype"

    uncertain <- candidate_variant_granges[uncertain]
    uncertain$duoNovo_classification <- "uncertain"

    output <- c(de_novo, not_de_novo, uncertain)
  } else if (phasing_orientation == "right") {
    de_novo <- candidate_variant_granges[hap2_inherited]
    de_novo$duoNovo_classification <- "de novo"
    de_novo$hamming_distance_other_parent_hap <- hap2_hamming_distance_other_parent_hap

    de_novo$supporting_counts_het_hom <- hap2_supporting_counts_het_hom
    de_novo$supporting_counts_het_het <- hap2_supporting_counts_het_het
    de_novo$supporting_counts_hom_het <- hap2_supporting_counts_hom_het

    not_de_novo <- candidate_variant_granges[hap1_inherited]
    not_de_novo$duoNovo_classification <- "on other parent haplotype"

    uncertain <- candidate_variant_granges[uncertain]
    uncertain$duoNovo_classification <- "uncertain"

    output <- c(de_novo, not_de_novo, uncertain)
  }
  output_sorted <- sort(output)

  message("Writing classified variants to VCF file...")
  # Create a DataFrame for info columns, adding the classifications
  info <- DataFrame(duoNovo_classification = output_sorted$duoNovo_classification,
                    output_sorted$hamming_distance_other_parent_hap,
                    output_sorted$supporting_counts_het_hom,
                    output_sorted$supporting_counts_het_het,
                    output_sorted$supporting_counts_hom_het)

  # Create the VCF object and save to the current directory
  vcf_out <- VCF(rowRanges = output_sorted, info = info)
  writeVcf(vcf_out, output_vcf_name)

  return(output_sorted)
}


