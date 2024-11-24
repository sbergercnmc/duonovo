
<!-- README.md is generated from README.Rmd. Please edit that file -->

# duoNovo

<!-- badges: start -->
<!-- badges: end -->

duoNovo is an R package that identifies de novo variants from single
parent-proband duos, that is, without having to sequence both biological
parents. duoNovo uses phased variant calls generated via long-read
sequencing followed by variant calling and read-backed phasing (for
example, using HiPhase), in order to determine whether candidate
variants of interest (heterozygous in the proband; absent in the parent)
are present on haplotypes that are identical by descent between the
sequenced parent and the proband. If that is the case, duoNovo
classifies these candidate variants as de novo, whereas if that is not
the case it infers that they have arisen on the haplotype inherited from
the non-sequenced parent (thus ascertaining their de novo status is not
possible).

## Installation

duoNovo can be installed as follows:

``` r
# Step 1: Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("VariantAnnotation", "GenomicRanges", "IRanges", "S4Vectors"))

# Step 2: Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Or alternatively using remotes
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

# Step 3: Install duoNovo from GitHub
devtools::install_github("sbergercnmc/duoNovo", dependencies = TRUE)

# Or alternatively
remotes::install_github("sbergercnmc/duoNovo", dependencies = TRUE)
```

## Run duoNovo

duoNovo is simple to use. All it requires is a VCF containing phased
variant calls. Optionally, a VCF containing variant calls from
short-read sequencing can be used, in case one wants to only classify
candidate variants whose genotype calls are concordant with the genotype
calls from short-read sequencing. Some additional parameters which can
optionally be adjusted are described below. The following is a simple
example of how to run duoNovo.

### Prepare Input VCF from Long-Read Sequencing Data

The first step is to prepare the input VCF file from long-read
sequencing data. Below are steps to generate a suitable VCF:

1.  **Call Variants**: Generate gVCF files separately for the proband
    and the parent. A variant caller such as **DeepVariant** can be used
    for this step.

2.  **Joint Calling**: Run the generated gVCF files through **GLnexus**
    to produce a joint-called VCF file for the duo.

3.  **Phase Variants**: Use **HiPhase** (available at [HiPhase
    GitHub](https://github.com/PacificBiosciences/HiPhase)) to phase the
    variants. This step will annotate each phased variant with a phasing
    set ID, which is required for reconstruction of haplotypes with
    `duoNovo`.

Optionally, a VCF from short-read sequencing data can also be prepared
using steps 1 and 2 above.

### Classify candidate variants with duoNovo

The following shows how to run duoNovo, assuming we have a VCF
“duo_proband_father.longread.hiphase.vcf.gz” from long-read sequencing,
and a VCF “duo_proband_father.shortRead.vcf.gz” from short-read
sequencing.

``` r
library(duoNovo)
duoNovo_results <- duoNovo(
  LRS_phased_vcf_file_path = "duo_proband_father.longread.hiphase.vcf.gz", 
  depth_cutoff = 20, 
  GQ_cutoff = 30,
  proband_phasing = "1|0", 
  proband_column_identifier = "-0$",
  PS_width_cutoff = 10000, 
  boundary_cutoff = 2000, 
  distance_cutoff = 40,
  candidate_variants_concordant_with_SRS = TRUE,
  SRS_vcf_file_path = "duo_proband_father.shortRead.vcf.gz",
  reference = "hg38"
)
```

duoNovo returns a `GRangesList` containing three elements: `de_novo`,
`not_de_novo`, and `uncertain`. Each of these elements corresponds to a
`GRanges` containing the variants that have been classified into the
corresponding category.

Below is a detailed description of each argument of `duoNovo()`:

- **LRS_phased_vcf_file_path**: File path to the VCF containing phased
  variant calls from long-read sequencing of the duo.
- **depth_cutoff**: A numeric value specifying the minimum sequencing
  depth for variants to be included in the analysis. The same cutoff
  applies to both proband and parent variants, for both LRS and SRS (if
  used).
- **GQ_cutoff**: A numeric value specifying the minimum GQ (genotype
  quality) for variants to be included in the analysis. The same cutoff
  applies to both proband and parent variants, for both LRS and SRS (if
  used).
- **proband_phasing**: A character string indicating the phasing of the
  proband, either “1\|0” or “0\|1”.
- **proband_column_identifier**: A character corresponding to an
  identifier for the proband column in the metadata matrices. Should be
  the same for both the LRS VCF and (if used) the SRS VCF.
- **PS_width_cutoff**: A numeric value specifying the minimum width for
  phasing sets to be included in the analysis.
- **boundary_cutoff**: A numeric value indicating the minimum distance
  from a haplotype block boundary (either start or end coordinate) for
  candidate variants to be analyzed.
- **distance_cutoff**: A numeric value specifying the minimum Hamming
  distance cutoff to determine that a proband-parent haplotype block is
  not identical by descent.
- **candidate_variants_concordant_with_SRS**: Logical value specifying
  if candidate variants should be concordant with short-read sequencing
  (default is `TRUE` or `FALSE`).
- **SRS_vcf_file_path**: File path to the VCF containing variant calls
  from short-read sequencing of the duo.
- **reference**: name of reference genome used, e.g. hg38
