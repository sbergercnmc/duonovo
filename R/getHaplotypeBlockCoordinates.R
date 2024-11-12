#' Get haplotype block boundaries
#'
#' The `getHaplotypeBlockCoordinates` function takes as input a `GRanges` object containing haplotype information, and identifies genomic segments where proband variants are assigned to the same phasing set, and parental variants are assigned to the same phasing set.
#' The function combines phasing set boundaries (`PS1` and `PS2`) from the proband and parent phasing (using HiPhase), creating a unified `GRanges` object.
#'
#' @param haplotype_granges A `GRanges` object containing haplotype information with `PS1` and `PS2` columns representing the phasing sets of the proband and parent, respectively.
#' @return A `GRanges` object with combined coordinates of overlapping phasing sets, including `PS1` and `PS2` identifiers.
#' @details The function splits haplotype ranges by `PS1` and `PS2` identifiers, calculates start and end coordinates for each group,
#' and then finds overlapping phasing sets between `PS1` and `PS2`. The resulting overlapping coordinates are stored in a new `GRanges` object.
#'
#' @importFrom S4Vectors queryHits subjectHits
#' @import GenomicRanges
#' @import IRanges
#' @import VariantAnnotation
#' @export
#'
#' @examples
#' # Assuming `haplotype_granges` is defined and contains `PS1` and `PS2` columns:
#' # combined_coordinates <- getPhasingSetCoordinates(haplotype_granges)
getHaplotypeBlockCoordinates <- function(haplotype_granges) {
  hap_granges <- haplotype_granges

  # Split haplotype ranges by PS1 and calculate coordinates
  granges_by_ps1 <- split(hap_granges, hap_granges$PS1)
  all_chrs <- sapply(granges_by_ps1, function(xx) unique(seqnames(xx)))
  indices_to_use <- which(lengths(all_chrs) == 1)
  start_coords <- sapply(granges_by_ps1[indices_to_use], function(xx) min(start(xx)))
  end_coords <- sapply(granges_by_ps1[indices_to_use], function(xx) max(end(xx)))
  all_chrs <- unlist(all_chrs[indices_to_use])

  ps1_boundaries <- GRanges(
    seqnames = all_chrs, ranges = IRanges(start = start_coords, end = end_coords),
    PS1 = names(granges_by_ps1)[indices_to_use]
  )

  # Split haplotype ranges by PS2 and calculate coordinates
  granges_by_ps2 <- split(hap_granges, hap_granges$PS2)
  start_coords <- sapply(granges_by_ps2, function(xx) min(start(xx)))
  end_coords <- sapply(granges_by_ps2, function(xx) max(end(xx)))
  all_chrs <- sapply(granges_by_ps2, function(xx) unique(seqnames(xx)))

  ps2_boundaries <- GRanges(
    seqnames = all_chrs, ranges = IRanges(start = start_coords, end = end_coords),
    PS2 = names(granges_by_ps2)
  )

  # Find overlaps between ps1_boundaries and ps2_boundaries
  overlaps <- findOverlaps(ps1_boundaries, ps2_boundaries)

  # Create a new GRanges object to store the combined information
  combined_PS <- GRanges(
    seqnames = seqnames(ps1_boundaries)[queryHits(overlaps)],
    ranges = IRanges(
      start = pmax(start(ps1_boundaries)[queryHits(overlaps)], start(ps2_boundaries)[subjectHits(overlaps)]),
      end = pmin(end(ps1_boundaries)[queryHits(overlaps)], end(ps2_boundaries)[subjectHits(overlaps)])
    ),
    PS1 = ps1_boundaries$PS1[queryHits(overlaps)],
    PS2 = ps2_boundaries$PS2[subjectHits(overlaps)]
  )

  # Remove any ranges that have become invalid (start > end) due to the pmax/pmin operation
  valid_ranges <- start(combined_PS) <= end(combined_PS)
  combined_PS <- combined_PS[valid_ranges]
  combined_PS
}
