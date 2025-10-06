
<!-- README.md is generated from README.Rmd. Please edit that file -->

# *duoNovo*: Identification of *de novo* variants from single parent-proband duos

<!-- badges: start -->
<!-- badges: end -->

*duoNovo* is an R package that identifies *de novo* variants from
parent-proband duos, that is, without requiring sequencing of both
biological parents. *duoNovo* uses phased variant calls generated via
long-read sequencing followed by read-backed phasing (for example, using
HiPhase), in order to determine whether candidate variants of interest
(heterozygous in the proband; absent in the parent) are present on
haplotypes that are identical by descent between the sequenced parent
and the proband. If that is the case, *duoNovo* classifies these
candidate variants as *de novo*, whereas if that is not the case it
infers that they have arisen on the haplotype inherited from the missing
parent (thus ascertaining their *de novo* status is not possible).

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

# Step 4: Optionally, to use the duoNovo command line interface, install argparse
if (!requireNamespace("argparse", quietly = TRUE)) {
    install.packages("argparse")
}
```

## Run *duoNovo*

*duoNovo* is simple to use. All it requires is a joint-called VCF for
the parent-proband duo with phased variant calls. The following shows
the basic steps of running *duoNovo*.

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

(Optionally, a VCF from short-read sequencing data can also be prepared
using steps 1 and 2 above, in case one only wants to classify variants
whose genotype is concordant between long-read and short-read
sequencing).

#### Preprocessing Script

The preprocess folder includes an optional shell script to automate the
processing steps above for the long-read sequencing duos.

    Usage: prep_duo_input.sh PROBAND_GVCF PROBAND_BAM PARENT_GVCF PARENT_BAM REFERENCE_FASTA OUTPUT_VCF [THREADS:default nproc]
    Requires: bcftools, glnexus_cli, hiphase
    Inputs:
    PROBAND_GVCF: g.vcf file from the proband (can be g.vcf or g.vcf.gz)
    PROBAND_BAM: bam (or cram) file from the proband
    PARENT_GVCF: g.vcf file from the parent (can be g.vcg or g.vcf.gz)
    PARENT_BAM: bam (or cram) file from the proband
    REFERENCE_FASTA: reference fasta file for the aligned sequences
    OUTPUT_VCF: vcf file to be created, a joint called (glnexus DeepVariant_unfiltered) duo vcf with both samples with phasing added from hiphase
    optional THREADS: number of thread to use for glnexus and hiphase. Default is is nproc - 88

### Classify candidate variants with *duoNovo*

The following shows how to run duoNovo within R, assuming we have a VCF
“duo_proband_father.longread.hiphase.vcf.gz” from long-read sequencing.

``` r
library(duoNovo)
duoNovo_results <- duoNovo(
  LRS_phased_vcf_file_path = "duo_proband_father.longread.hiphase.vcf.gz", 
  proband_column_identifier = "-0$"
  )
```

*duoNovo* returns a `GRanges` containing all candidate variants
identified from the joint-called duo vcf. This `GRanges` is a subset of
the `rowRanges` corresponding to the input LRS VCF, but contains
additional metadata columns. A column named `duoNovo_classification`
provides the classification of each variant as *de novo*, present on the
haplotype inherited from the missing parent, or uncertain. See below for
a more detailed description of the output.

Optionally, the output is written into a vcf.

The user can optionally provide a vector of coordinates of variants of
interest. In this case, only these variants will be tested for *de novo*
status, provided they are appropriate candidates (that is, they are
heterozygous in the proband and absent in the parent), and pass QC
filters. If specific variants of interest are not provided, then
*duoNovo* tests all variants heterozygous in the proband and absent in
the parent that pass QC.

Below is a description of each argument that `duoNovo()` accepts:

- **LRS_phased_vcf_file_path**: Required input. Path to the VCF
  containing phased variant calls from long-read sequencing of the duo.
- **depth_cutoff**: A numeric value specifying the minimum sequencing
  depth for variants to be included in the analysis. The same cutoff
  applies to both proband and parent variants, for both LRS and SRS (if
  used). The default is 20.
- **GQ_cutoff**: A numeric value specifying the minimum GQ (PHRED
  quality) for variants to be included in the analysis. The same cutoff
  applies to both proband and parent variants, for both LRS and SRS (if
  used). The default is 30.
- **proband_column_identifier**: Required input. A character
  corresponding to an identifier for the proband column in the VCF.
  Should be the same for both the LRS VCF and (if used) the SRS VCF.
- **PS_width_cutoff**: A numeric value specifying the minimum width of
  haplotype blocks used in the analysis. The default is 10000.
- **boundary_cutoff**: A numeric value indicating the minimum distance
  from a haplotype block boundary (either start or end coordinate) for
  candidate variants to be evaluated. The default is 2000.
- **IBD_distance_cutoff**: A numeric value specifying the maximum
  Hamming distance cutoff to infer that a proband-parent haplotype block
  is identical by descent. The default is 0.
- **non_IBD_distance_cutoff**: A numeric value specifying the minimum
  Hamming distance cutoff to infer that a proband-parent haplotype block
  is not identical by descent. The default is 40.
- **candidate_variants_concordant_with_SRS**: Logical value specifying
  if candidate variant genotype calls should be concordant with
  short-read sequencing (the default is `FALSE`).
- **SRS_vcf_file_path**: Path to the VCF containing variant calls from
  short-read sequencing of the duo.
- **test_reference_allele**: Logical value specifying if positions where
  the proband is heterozygous and the parent is homozygous for the
  variant allele should also be tested (the default is `FALSE`).
- **candidate_variant_coordinates**: A vector of 1-based coordinates
  (e.g. `c("chr1:1000", "chr2:2000")`) of candidate variants of
  interest.
- **problematic_regions**: Bed file with regions to exclude from
  candidates and data (not supplied by default).
- **output_vcf_path**: Path to write an output VCF file. Optional. If
  used, by default it appends \_duoNovo to input vcf file path.
- **compress_output**: Logical value specifying whether to compress the
  output VCF file, append .bgz to file name, and tabix index. The
  default is `TRUE`

### Run *duoNovo* from the Command line

The *duoNovo* command line interface provides a tool to run *duoNovo*
directly from the command line. After installing the ‘*duoNovo*’ and the
‘argparse’ packages, the first step is to clone the `duoNovo_cli.R`
script from the cli directory. The command line script can then be
invoked using Rscript from the command line and has command line
parameters to set each of the `duoNovo()` arguments.

Below is the output of `Rscript duoNovo_cli.R -h`.

    usage: cli/duoNovo_cli.R [-h] [-v] [-q] [-d number] [-g number] -f FILE -p
                             PROBAND_SAMPLE_ID [-w number] [-b number] [-c number]
                             [-t] [-n coordinates] [-s FILE] [-o FILE] [-z]

    options:
      -h, --help            show this help message and exit
      -v, --verbose         Print extra output/parameter values [default TRUE]
      -q, --quietly         Print little output
      -d number, --depth_cutoff number
                            Numeric value specifying the minimum sequencing depth for variants to be included in the analysis [default 20]
      -g number, --GQ_cutoff number
                            Numeric value specifying the minimum GQ (PHRED quality) for variants to be included in the analysis [default 30]
      -f FILE, --LRS_phased_vcf_file_path FILE
                            Path to joint called phased duo vcf [REQUIRED]
      -p PROBAND_SAMPLE_ID, --proband_id PROBAND_SAMPLE_ID
                            Identifier for the proband column in the metadata matrices of the VCF [REQUIRED]
      -w number, --PS_width_cutoff number
                            Numeric value specifying the minimum width of haplotype blocks used in the analysis. [default
                            10000]
      -b number, --boundary_cutoff number
                            Numeric value indicating the minimum distance from a
                            haplotype block boundary (either start or end
                            coordinate) for candidate variants to be evaluated.
                            [default 2000]
      -i number, --IBD_distance_cutoff number
                            A numeric value specifying the maximum Hamming
                            distance cutoff to infer that a proband-parent
                            haplotype block is identical by descent. [default 0]                      
      -c number, --non_IBD_distance_cutoff number
                            A numeric value specifying the minimum Hamming
                            distance cutoff to infer that a proband-parent
                            haplotype block is not identical by descent. [default 40]
      -t, --test_reference_allele
                            Logical value specifying if positions where the proband is heterozygous and the parent is homozygous for the variant allele should also be tested [default FALSE]
      -n coordinates, --candidate_variant_coordinates coordinates
                            Vector of 1-based coordinates of candidate variants, e.g. chr1:12345-12345,chr2:65430-65430.
                            [Optional]
      -x BED_FILE, --exclude_regions_bed_path BED_FILE
                            Path to bedfile for regions to exclude [Optional]
      -s FILE, --SRS_vcf_file_path FILE
                            Path to short read duo vcf [Optional]
      -o FILE, --output_vcf FILE
                            Path to file to write output vcf [Optional: Default
                            appends _duoNovo to input vcf]
      -z, --compress_output
                            Logical value specifying whether to compress output, append .bgz to filename, and tabix
                            index [default TRUE]

### The output of *duoNovo*

`duoNovo()` outputs a `GRanges` (or optionally, a vcf) containing all
candidate variants identified from the joint-called duo vcf (or the
specified list of candidate variants, if such a list was provided as
input). This `GRanges` is a subset of the `rowRanges` corresponding to
the input LRS VCF, but contains additional metadata columns. If an
output vcf is written, these columns are part of the vcf INFO. The
column of main interest is `duoNovo_classification`. It provides the
classification of each candidate variant as *de novo* vs present on the
haplotype inherited from the missing vs uncertain. Some of the other
metadata columns are self-explanatory, and provide information about the
phasing of the proband and the parent, the sequencing depth, and GQ.
Additional columns include:

- **QC_fail_step**, which describes the specific QC step that variants
  that didn’t pass QC failed. Its values are self-explanatory
  (e.g. `low_depth`), and are `NA` for variants that passed QC)
- **non_IBD_hamming_distance**, which corresponds to the minimum of the
  two Hamming distances between the proband haplotype block inferred to
  have been transmitted from the missing parent and the two parental
  haplotype blocks.
- **IBD_hamming_distance**, which corresponds to the Hamming distance
  between the proband-parent haplotype block inferred to be
  identical-by-descent.
- **supporting_counts_het_hom**, which corresponds to the number of
  positions that went into the Hamming distance calculation where the
  proband is heterozygous and the parent is homozygous (see the Methods
  section of the preprint for details).
- **supporting_counts_het_het**, which corresponds to the number of
  positions that went into the Hamming distance calculation where both
  the proband and the parent are heterozygous (see preprint Methods).
- **supporting_counts_hom_het**, which corresponds to the number of
  positions that went into the Hamming distance calculation where the
  proband is homozygous and the parent is heterozygous (see preprint
  Methods).
- **n_de_novo_same_orientation_same_PS**, which corresponds to the total
  number of *de novo* variants in the same orientation and same
  haplotype block (NA for non- *de novo* variants). If that number is
  greater than 1, then the variant is classified as
  `on_multi_denovo_haplotype`. In our evaluation, we have found that
  these variants often tend to be false positives.

### Interpretation of *duoNovo* results

Classifications generated by *duoNovo* should be evaluated according to
variant quality, sequence quality, and regional quality.

Moreover, as is true for any classifier, the prior probability of each
variant being *de novo* affects the predictive value of *duoNovo*’s
classifications. For instance, variants absent from gnomAD have a higher
prior probability of being *de novo* (and pathogenic) and thus the
positive predictive value of *duoNovo* for these variants is higher. In
our evaluations, *duoNovo* achieved near-perfect accuracy among these
variants.

### *duoNovo* performance

We systematically evaluated *duoNovo*’s performance on a cohort of 117
trios which we used to construct 234 duos (117 father-proband and 117
mother-proband duos). Details can be found in our preprint
(<https://www.medrxiv.org/content/10.1101/2025.02.24.25322424v1>).
