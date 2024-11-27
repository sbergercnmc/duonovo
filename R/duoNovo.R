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
                    output_vcf_path) {
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

  vcf_granges_filtered <- vcf_granges[which(vcf_granges$depth1 >= depth_cutoff & vcf_granges$GQ1 >= GQ_cutoff &
                                              vcf_granges$depth2 >= depth_cutoff & vcf_granges$GQ2 >= GQ_cutoff)]

  message("Reconstructing haplotypes...")
  hap_granges <- getHaplotypes(vcf_granges_filtered)
  hap_boundary_coordinates <- getHaplotypeBlockCoordinates(hap_granges)
  hap_boundary_coordinates <- hap_boundary_coordinates[which(width(hap_boundary_coordinates) > PS_width_cutoff)]

  # If candidate variant coordinates are provided, filter hap_granges accordingly
  if (!is.null(candidate_variant_coordinates)) {
    candidate_variant_granges <- GRanges(
      seqnames = sapply(strsplit(candidate_variant_coordinates, ":"), `[[`, 1),
      ranges = IRanges(
        start = as.numeric(sapply(strsplit(sapply(strsplit(candidate_variant_coordinates, ":"), `[[`, 2), "-"), `[[`, 1)),
        end = as.numeric(sapply(strsplit(sapply(strsplit(candidate_variant_coordinates, ":"), `[[`, 2), "-"), `[[`, 2))
      ))
    
    candidate_variant_indices <- unique(queryHits(findOverlaps(hap_granges, candidate_variant_granges)))
    candidate_variant_granges_combined <- hap_granges[candidate_variant_indices]
    candidate_variant_granges_left <- candidate_variant_granges_combined[
      which(candidate_variant_granges_combined$phasing1 == "1|0")]
    candidate_variant_granges_right <- candidate_variant_granges_combined[
      which(candidate_variant_granges_combined$phasing1 == "0|1")]
    
  } else { # Otherwise, identify candidate de novo variants directly from genotypes
    candidate_variant_indices_left <- which(hap_granges$phasing1 == "1|0" & hap_granges$phasing2 == "0/0")  
    candidate_variant_indices_right <- which(hap_granges$phasing1 == "0|1" & hap_granges$phasing2 == "0/0")
    candidate_variant_granges_left <- hap_granges[candidate_variant_indices_left]
    candidate_variant_granges_right <- hap_granges[candidate_variant_indices_right]
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
  classifications_left <- classifyVariants(candidate_variant_granges_left, phasing_orientation = "left", 
                                           haplotype_granges = hap_granges, 
                                           haplotype_boundary_coordinate_granges = hap_boundary_coordinates, 
                                           boundary_cutoff = boundary_cutoff, distance_cutoff = distance_cutoff)
  classifications_right <- classifyVariants(candidate_variant_granges_right, phasing_orientation = "right", 
                                           haplotype_granges = hap_granges, 
                                           haplotype_boundary_coordinate_granges = hap_boundary_coordinates, 
                                           boundary_cutoff = boundary_cutoff, distance_cutoff = distance_cutoff)
  
  duo_novo_classifications <- c(classifications_left, classifications_right)
  output <- duo_novo_classifications
  output_sorted <- sort(output)
  
  message("Writing classified variants to VCF file...")
  # Create a DataFrame for info columns, adding the classifications
  info <- DataFrame(
    duoNovo_classification = output_sorted$duoNovo_classification,
    hamming_distance_other_parent_hap = output_sorted$hamming_distance_other_parent_hap,
    supporting_counts_het_hom = output_sorted$supporting_counts_het_hom,
    supporting_counts_het_het = output_sorted$supporting_counts_het_het,
    supporting_counts_hom_het = output_sorted$supporting_counts_hom_het
  )
  
  # Create the VCF object and save it to the current directory
  vcf_out <- VCF(rowRanges = output_sorted, info = info)
  writeVcf(vcf_out, output_vcf_path)
  
  return(output_sorted)
}
