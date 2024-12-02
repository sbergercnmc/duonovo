
<!-- README.md is generated from README.Rmd. Please edit that file -->

# *duoNovo*: Identification of *de novo* variants from single parent-proband duos

<!-- badges: start -->
<!-- badges: end -->

*duoNovo* is an R package that identifies *de novo* variants from single
parent-proband duos, that is, without requiring sequencing of both
biological parents. *duoNovo* uses phased variant calls generated via
long-read sequencing followed by variant calling and read-backed phasing
(for example, using HiPhase), in order to determine whether candidate
variants of interest (heterozygous in the proband; absent in the parent)
are present on haplotypes that are identical by descent between the
sequenced parent and the proband. If that is the case, *duoNovo*
classifies these candidate variants as *de novo*, whereas if that is not
the case it infers that they have arisen on the haplotype inherited from
the non-sequenced parent (thus ascertaining their *de novo* status is
not possible).

## Installation

*duoNovo* can be installed as follows:

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

# Or alternatively install remotes
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

# Step 3: Install duoNovo from GitHub using devtools
devtools::install_github("sbergercnmc/duoNovo", dependencies = TRUE)

# Or alternatively using remotes
remotes::install_github("sbergercnmc/duoNovo", dependencies = TRUE)
```

## Run *duoNovo*

*duoNovo* is simple to use. All it requires is a VCF containing phased
variant calls. Optionally, a VCF containing variant calls from
short-read sequencing can be used, in case one wants to only classify
candidate variants whose genotype calls are concordant with the genotype
calls from short-read sequencing. Some additional parameters which can
optionally be adjusted are described below. The following is a simple
example of how to run *duoNovo*.

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

### Classify candidate variants with *duoNovo*

The following shows how to run duoNovo, assuming we have a VCF
“duo_proband_father.longread.hiphase.vcf.gz” from long-read sequencing.

``` r
library(duoNovo)
duoNovo_results <- duoNovo(
  LRS_phased_vcf_file_path = "duo_proband_father.longread.hiphase.vcf.gz", 
  proband_column_identifier = "-0$"
  )
```

*duoNovo* returns a `GRanges` containing all candidate variants tested
for *de novo* status. This `GRanges` is a subset of the `rowRanges`
corresponding to the input LRS VCF, but contains additional metadata
columns. A column named `duoNovo_classification` provides the
classification of each variant as *de novo*, present on the haplotype
inherited from the non-sequenced parent, or uncertain. Additional
columns provide further pertinent information, such as the minimum of
the two Hamming distances supporting the non-IBD status of the proband
haplotype not containing the variant of interest and the two parental
haplotypes, as well as the counts of the different proband-parent
genotype classes supporting the classification. For variants that were
not classified because they failed QC, a column named `QC_fail_step`
provides the specific QC step they failed.

The user can optionally provide a vector of coordinates of variants of
interest. In this case, only these variants will be tested for *de novo*
status, provided they are appropriate candidates (that is, they are
heterozygous in the proband and absent in the parent), and they pass QC.
If specific variants of interest are not provided, then *duoNovo* tests
all variants heterozygous in the proband and absent in the parent that
pass QC.

Optionally, the output is written into a vcf file.

Below is a description of each argument that `duoNovo()` accepts:

- **LRS_phased_vcf_file_path**: File path to the VCF containing phased
  variant calls from long-read sequencing of the duo.
- **depth_cutoff**: A numeric value specifying the minimum sequencing
  depth for variants to be included in the analysis. The same cutoff
  applies to both proband and parent variants, for both LRS and SRS (if
  used). The default is 20.
- **GQ_cutoff**: A numeric value specifying the minimum GQ (genotype
  quality) for variants to be included in the analysis. The same cutoff
  applies to both proband and parent variants, for both LRS and SRS (if
  used). The default is 30.
- **proband_column_identifier**: A character corresponding to an
  identifier for the proband column in the metadata matrices of the VCF.
  Should be the same for both the LRS VCF and (if used) the SRS VCF.
- **PS_width_cutoff**: A numeric value specifying the minimum width of
  haplotype blocks used in the analysis. The default is 10000.
- **boundary_cutoff**: A numeric value indicating the minimum distance
  from a haplotype block boundary (either start or end coordinate) for
  candidate variants to be analyzed. The default is 2000.
- **distance_cutoff**: A numeric value specifying the minimum Hamming
  distance cutoff to determine that a proband-parent haplotype block is
  not identical by descent. The default is 40.
- **candidate_variants_concordant_with_SRS**: Logical value specifying
  if candidate variant genotype calls should be concordant with
  short-read sequencing (the default is `FALSE`).
- **SRS_vcf_file_path**: File path to the VCF containing variant calls
  from short-read sequencing of the duo.
- **test_reference_allele**: Logical value specifying if positions where
  the proband is heterozygous and the parent is homozygous for the
  variant allele should also be tested (the default is `FALSE`).
- **reference**: Name of reference genome (e.g. hg38) used by the vcf.
- **candidate_variant_coordinates**: A vector of coordinates
  (e.g. c(chr1:1000, chr2:2000)) of candidate variants of interest.
- **output_vcf_path**: File path for the output VCF file.
- **compress_output**: Logical value specifying whether or not to
  compress the output VCF file. The default is `TRUE`
