#' Extract haplotypes from phased variant calls
#'
#' The `getHaplotypes` function extracts haplotype information from a `GRanges` object containing phased variant calls for both proband and parent.
#' Each of the two haplotypes of an individual is represented in a separate column.
#'
#' @param phased_variant_granges_duo A `GRanges` object containing phased variant calls for both the proband and parent. Generated with the VariantAnnotation package from a phased vcf as input. We also recommend filtering based on sequencing depth and genotype quality prior to extracting haplotypes.
#' @return A `GRanges` object with additional columns (`hap11`, `hap21`, `hap12`, `hap22`) indicating phased alleles for both individuals.
#' @details The function identifies the haplotype data for the proband and parent from the `phasing1` and `phasing2` columns.
#' The phased alleles are categorized for each variant using regular expressions to determine the phase (e.g., `0|1`, `1|0`, etc.).
#'
#' @import GenomicRanges
#' @import IRanges
#' @import VariantAnnotation
#' @export
#'
#' @examples
#' # Assuming `phased_variant_granges_duo` is defined and contains the appropriate phasing columns:
#' # phased_haplotypes <- getHaplotypes(phased_variant_granges_duo)
getHaplotypes <- function(phased_variant_granges_duo) {
  granges <- phased_variant_granges_duo
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

  ### Extracting haplotypes for the parent
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

  granges
}
