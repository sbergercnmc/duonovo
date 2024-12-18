
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

# Step 4: Optionally, to use the duoNovo command line interface, install argparse
if (!requireNamespace("argparse", quietly = TRUE)) {
    install.packages("argparse")
}
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
“duo\_proband\_father.longread.hiphase.vcf.gz” from long-read
sequencing.

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

  - **LRS\_phased\_vcf\_file\_path**: Required input. Path to the VCF
    containing phased variant calls from long-read sequencing of the
    duo.
  - **depth\_cutoff**: A numeric value specifying the minimum sequencing
    depth for variants to be included in the analysis. The same cutoff
    applies to both proband and parent variants, for both LRS and SRS
    (if used). The default is 20.
  - **GQ\_cutoff**: A numeric value specifying the minimum GQ (genotype
    quality) for variants to be included in the analysis. The same
    cutoff applies to both proband and parent variants, for both LRS and
    SRS (if used). The default is 30.
  - **proband\_column\_identifier**: Required input. A character
    corresponding to an identifier for the proband column in the
    metadata matrices of the VCF. Should be the same for both the LRS
    VCF and (if used) the SRS VCF.
  - **PS\_width\_cutoff**: A numeric value specifying the minimum width
    of haplotype blocks used in the analysis. The default is 10000.
  - **boundary\_cutoff**: A numeric value indicating the minimum
    distance from a haplotype block boundary (either start or end
    coordinate) for candidate variants to be analyzed. The default is
    2000.
  - **distance\_cutoff**: A numeric value specifying the minimum Hamming
    distance cutoff to determine that a proband-parent haplotype block
    is not identical by descent. The default is 40.
  - **candidate\_variants\_concordant\_with\_SRS**: Logical value
    specifying if candidate variant genotype calls should be concordant
    with short-read sequencing (the default is `FALSE`).
  - **SRS\_vcf\_file\_path**: Path to the VCF containing variant calls
    from short-read sequencing of the duo.
  - **test\_reference\_allele**: Logical value specifying if positions
    where the proband is heterozygous and the parent is homozygous for
    the variant allele should also be tested (the default is `FALSE`).
  - **candidate\_variant\_coordinates**: A vector of coordinates (e.g.
    `c("chr1:1000", "chr2:2000")`) of candidate variants of interest.
  - **output\_vcf\_path**: Path to write the output VCF file.
  - **compress\_output**: Logical value specifying whether or not to
    compress the output VCF file. The default is `TRUE`

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
                            Depth cutoff for variant evaluation [default 20]
      -g number, --GQ_cutoff number
                            GQ cutoff for variant evaluation [default 30]
      -f FILE, --LRS_phased_vcf_file_path FILE
                            Path to joint called phased duo vcf [REQUIRED]
      -p PROBAND_SAMPLE_ID, --proband_id PROBAND_SAMPLE_ID
                            VCF column heading for Proband [REQUIRED]
      -w number, --PS_width_cutoff number
                            A numeric value specifying the minimum width for
                            phasing sets to be included in the analysis. [default
                            10000]
      -b number, --boundary_cutoff number
                            A numeric value indicating the minimum distance from a
                            haplotype block boundary (either start or end
                            coordinate) for candidate variants to be analyzed.
                            [default 2000]
      -c number, --distance_cutoff number
                            A numeric value specifying the minimum Hamming
                            distance cutoff to determine that a proband-parent
                            haplotype block is not identical by descent. [default
                            40]
      -t, --test_reference_allele
                            Test for deNovo Reference reversions (parent is
                            hom_var, proband is het) [default FALSE]
      -n coordinates, --candidate_variant_coordinates coordinates
                            1-based list of chromosome ranges to evaluate for
                            variants, e.g. chr1:12345-12345,chr2:65430-65430.
                            [Optional]
      -s FILE, --SRS_vcf_file_path FILE
                            Path to short read duo vcf [Optional]
      -o FILE, --output_vcf FILE
                            Path to file to write output vcf [Optional: Default
                            appends _duoNovo to input vcf]
      -z, --compress_output
                            Compress output, append .bgz to filename, and tabix
                            index [default FALSE]
