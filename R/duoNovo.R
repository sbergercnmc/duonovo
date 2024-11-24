#' Identify de novo variants from duo (proband - single parent) long read sequencing
#'
#' The `duoNovo` function uses phased variant calls from duo long-read sequencing to classify candidate variants (heterozygous in the proband; absent in the parent).
#' Variants are classifed as de novo, present on the haplotype inherited from the non-sequenced parent, or uncertain variants.
#'
#' @param LRS_phased_vcf_file_path File path to the vcf containing phased variant calls from long-read sequencing of the duo.
#' @param depth_cutoff A numeric value specifying the minimum sequencing depth for variants to be included in the analysis. The same cutoff applies to both proband and parent variants, for both LRS and SRS (if used).
#' @param GQ_cutoff A numeric value specifying the minimum GQ (genotype quality) for variants to be included in the analysis. The same cutoff applies to both proband and parent variants, for both LRS and SRS (if used).
#' @param proband_phasing A character string indicating the phasing of the proband, either "1|0" or "0|1".
#' @param proband_column_identifier A character corresponding to an identifier for the proband column in the metadata matrices. Should be the same for both the LRS vcf and (if used) the SRS vcf.
#' @param PS_width_cutoff A numeric value specifying the minimum width for phasing sets to be included in the analysis.
#' @param boundary_cutoff A numeric value indicating the minimum distance from a haplotype block boundary (either start or end coordinate) for candidate variants to be analyzed.
#' @param distance_cutoff A numeric value specifying the minimum hamming distance cutoff to determine that a proband-parent haplotype block are not identical by descent.
#' @param candidate_variants_concordant_with_SRS Logical value specifying if candidate variants should be concordant with short-read sequencing (default is `TRUE` or `FALSE`).
#' @param SRS_vcf_file_path File path to the vcf containing variant calls from short-read sequencing of the duo.
#'
#' @return A `GRangesList` containing three elements: `de_novo`, `not_de_novo`, and `uncertain`.
#' @details The function ultimately works by detecting identical by descent haplotype blocks, to determine whether each candidate variant of interest is de novo, using the genotype of only one parent. If requested, concordance with short-read sequencing can be checked.
#'
#' @importFrom matrixStats colMins
#' @importFrom S4Vectors queryHits subjectHits
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
#' #          SRS_vcf_file_path = my_SRS_file_path)
duoNovo <- function(LRS_phased_vcf_file_path, depth_cutoff = 20, GQ_cutoff = 30,
                    proband_phasing, proband_column_identifier,
                    PS_width_cutoff = 10000, boundary_cutoff = 2000, distance_cutoff = 40,
                    candidate_variants_concordant_with_SRS = TRUE,
                    SRS_vcf_file_path) {
  if (!file.exists(LRS_phased_vcf_file_path)) {
    stop("The LRS VCF file path does not exist: ", LRS_phased_vcf_file_path)
  }
  if (candidate_variants_concordant_with_SRS == TRUE && !file.exists(SRS_vcf_file_path)) {
    stop("The SRS VCF file path does not exist: ", SRS_vcf_file_path)
  }

  message("Importing LRS VCF file...")
  vcf <- readVcf(LRS_phased_vcf_file_path, genome = "hg38")
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

  vcf_granges_filtered <- vcf_granges[which(vcf_granges$depth1 >= depth_cutoff & vcf_granges$GQ1 >= GQ_cutoff &
                                              vcf_granges$depth2 >= depth_cutoff & vcf_granges$GQ2 >= GQ_cutoff)]

  message("Reconstructing haplotypes...")
  hap_granges <- getHaplotypes(vcf_granges_filtered)
  combined_PS <- getHaplotypeBlockCoordinates(hap_granges)
  combined_PS <- combined_PS[which(width(combined_PS) > PS_width_cutoff)]

  # Identifying potential de novo variants
  potential_de_novo_indices_allele1 <- which(hap_granges$phasing1 == proband_phasing & hap_granges$phasing2 == "0/0")
  hap_granges1 <- hap_granges[potential_de_novo_indices_allele1]

  # Optional: Restrict to candidate variants concordant with short-read sequencing
  if (candidate_variants_concordant_with_SRS == TRUE) {
    message("Importing SRS VCF file...")
    vcf_SR <- readVcf(SRS_vcf_file_path, genome = "hg38")
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
    hap_granges1 <- hap_granges1[unique(queryHits(findOverlaps(hap_granges1, vcf_granges_filtered_SR)))]
  }

  # Finding overlaps and analyzing haplotypes
  message("Classifying variants...")

  overlaps <- findOverlaps(hap_granges1, combined_PS - boundary_cutoff)
  overlapping_indices <- split(queryHits(overlaps), subjectHits(overlaps))

  haplotypes <- vector("list", length(overlapping_indices))
  names(haplotypes) <- names(overlapping_indices)

  indices <- as.numeric(names(overlapping_indices))
  PS1_ids <- combined_PS$PS1[indices]
  PS2_ids <- combined_PS$PS2[indices]

  hap_granges_no_denovo <- hap_granges[-potential_de_novo_indices_allele1]
  for(i in seq_along(indices)) {
    variant_granges <- hap_granges_no_denovo[
      unique(queryHits(findOverlaps(hap_granges_no_denovo, combined_PS[
        which(combined_PS$PS1 == PS1_ids[i] & combined_PS$PS2 == PS2_ids[i])])))]

    hap11 <- variant_granges$hap11
    hap12 <- variant_granges$hap12
    hap21 <- variant_granges$hap21
    hap22 <- variant_granges$hap22

    hap_mat <- matrix(c(hap11, hap12, hap21, hap22), nrow = length(hap11), ncol = 4, byrow = FALSE)
    colnames(hap_mat) <- c("hap11", "hap12", "hap21", "hap22")
    haplotypes[[i]] <- hap_mat
  }

  hamming_distance_mat <- sapply(haplotypes, function(xx)
    c(sum(xx[, "hap11"] != xx[, "hap12"], na.rm = TRUE),
      sum(xx[, "hap11"] != xx[, "hap22"], na.rm = TRUE),
      sum(xx[, "hap21"] != xx[, "hap12"], na.rm = TRUE),
      sum(xx[, "hap21"] != xx[, "hap22"], na.rm = TRUE))
  )

  hamming_distance_mins <- colMins(hamming_distance_mat)

  inheritance_stringent_hap1 <- apply(hamming_distance_mat, 2,
                                      function(xx) (xx[1] == 0 | xx[2] == 0) &
                                        xx[3] > distance_cutoff & xx[4] > distance_cutoff)

  inheritance_stringent_hap2 <- apply(hamming_distance_mat, 2,
                                      function(xx) (xx[3] == 0 | xx[4] == 0) &
                                        xx[1] > distance_cutoff & xx[2] > distance_cutoff)

  hap1_inherited <- unlist(overlapping_indices[which(inheritance_stringent_hap1 == TRUE)])
  hap2_inherited <- unlist(overlapping_indices[which(inheritance_stringent_hap2 == TRUE)])
  uncertain <- unlist(overlapping_indices[which(hamming_distance_mins > 0)])

  if(proband_phasing == "1|0") {
    output_list <- GRangesList(de_novo = hap_granges1[hap1_inherited], not_de_novo = hap_granges1[hap2_inherited],
                               uncertain = hap_granges1[uncertain])
  } else if (proband_phasing == "0|1") {
    output_list <- GRangesList(de_novo = hap_granges1[hap2_inherited], not_de_novo = hap_granges1[hap1_inherited],
                               uncertain = hap_granges1[uncertain])
  }
  output_list
}
