% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/duoNovo.R
\name{duoNovo}
\alias{duoNovo}
\title{Identify de novo variants from duo (proband - single parent) long read sequencing}
\usage{
duoNovo(
  LRS_phased_vcf_file_path,
  depth_cutoff = 20,
  GQ_cutoff = 30,
  proband_column_identifier,
  PS_width_cutoff = 10000,
  boundary_cutoff = 2000,
  distance_cutoff = 40,
  candidate_variants_concordant_with_SRS = FALSE,
  SRS_vcf_file_path = NULL,
  test_reference_allele = FALSE,
  candidate_variant_coordinates = NULL,
  problematic_regions = NULL,
  output_vcf_path = NULL,
  compress_output = TRUE
)
}
\arguments{
\item{LRS_phased_vcf_file_path}{File path to the vcf containing phased variant calls from long-read sequencing of the duo.}

\item{depth_cutoff}{A numeric value specifying the minimum sequencing depth for variants to be included in the analysis. The same cutoff applies to both proband and parent variants, for both LRS and SRS (if used).}

\item{GQ_cutoff}{A numeric value specifying the minimum GQ (genotype quality) for variants to be included in the analysis. The same cutoff applies to both proband and parent variants, for both LRS and SRS (if used).}

\item{proband_column_identifier}{A character corresponding to an identifier for the proband column in the metadata matrices. Should be the same for both the LRS vcf and (if used) the SRS vcf.}

\item{PS_width_cutoff}{A numeric value specifying the minimum width of the haplotype blocks included in the analysis.}

\item{boundary_cutoff}{A numeric value indicating the minimum distance from a haplotype block boundary (either start or end coordinate) for candidate variants to be analyzed.}

\item{distance_cutoff}{A numeric value specifying the minimum hamming distance cutoff to determine that a proband-parent haplotype block are not identical by descent.}

\item{candidate_variants_concordant_with_SRS}{Logical value specifying if candidate variants should be concordant with short-read sequencing (default is \code{FALSE}).}

\item{SRS_vcf_file_path}{File path to the vcf containing variant calls from short-read sequencing of the duo.}

\item{test_reference_allele}{Logical value specifying if positions where the proband is heterozygous and the parent is homozygous for the variant allele (not the reference) should be tested (default is \code{FALSE}).}

\item{candidate_variant_coordinates}{Vector of coordinates for specific variants of interest (of the form c(chr1:1000, chr2:2000)).}

\item{problematic_regions}{BED file with coordinates of problematic regions (e.g. as defined by Genome-in-a-Bottle)}

\item{output_vcf_path}{File path for output vcf.}

\item{compress_output}{Logical value specifying if output_vcf should be compressed and indexed. appended .bgz to filename (default is \code{TRUE}).}
}
\value{
A \code{GRanges} containing additional columns for the variant classifications as well as supporting information.
}
\description{
The \code{duoNovo} function uses phased variant calls from duo long-read sequencing to classify candidate variants (heterozygous in the proband; absent in the parent).
Variants are classifed as de novo, present on the haplotype inherited from the non-sequenced parent, or uncertain variants.
}
\details{
The function ultimately works by detecting identical by descent haplotype blocks, to determine whether each candidate variant of interest is de novo, using the genotype of only one parent. If requested, concordance with short-read sequencing can be checked.
}
\examples{
# duoNovo(LRS_phased_vcf_file_path = my_LRS_file_path, depth_cutoff = 20, GQ_cutoff = 30,
#          proband_column_identifier = my_proband_identifier,
#          PS_width_cutoff = 10000, boundary_cutoff = 2000, distance_cutoff = 40,
#          candidate_variants_concordant_with_SRS = FALSE,
#          SRS_vcf_file_path = my_SRS_file_path)
}
