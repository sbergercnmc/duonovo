#!/usr/bin/Rscript --vanilla

suppressPackageStartupMessages(library("duoNovo"))
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action="store_true",  default=TRUE,  help="Print extra output/parameter values [default TRUE]")
parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", help="Print little output")
parser$add_argument("-d", "--depth_cutoff", type="integer", default=20, 
                     help="Depth cutoff for variant evaluation [default %(default)d]",
                     metavar="number")
parser$add_argument("-g", "--GQ_cutoff", type="integer", default=30, 
                     help="GQ cutoff for variant evaluation [default %(default)d]",
                     metavar="number")
parser$add_argument("-f", "--LRS_phased_vcf_file_path", required=TRUE,  metavar="FILE", help = "Path to joint called phased duo vcf [REQUIRED]")
parser$add_argument("-p", "--proband_id",  metavar="PROBAND_SAMPLE_ID", required=TRUE, help = "VCF column heading for Proband [REQUIRED]")
parser$add_argument("-w", "--PS_width_cutoff", type="integer", default=10000, 
                     help="A numeric value specifying the minimum width for phasing sets to be included in the analysis. [default %(default)d]",
                     metavar="number")
parser$add_argument("-b", "--boundary_cutoff", type="integer", default=2000, 
                     help="A numeric value indicating the minimum distance from a haplotype block boundary (either start or end coordinate) for candidate variants to be analyzed. [default %(default)d]",
                     metavar="number")
parser$add_argument("-c", "--distance_cutoff", type="integer", default=40, 
                     help= "A numeric value specifying the minimum Hamming distance cutoff to determine that a proband-parent haplotype block is not identical by descent. [default %(default)d]",
                     metavar="number")
parser$add_argument("-t", "--test_reference_allele", action="store_true",  default=FALSE,  help="Test for deNovo Reference reversions (parent is hom_var, proband is het) [default FALSE]")
parser$add_argument("-n", "--candidate_variant_coordinates", metavar="coordinates", 
                     help= "1-based list of chromosome ranges to evaluate for variants, e.g.  chr1:12345-12345,chr2:65430-65430. [Optional]")
parser$add_argument("-x", "--exclude_regions_bed_path",  metavar="BED_FILE", help = "Path to bedfile for regions to exclude [Optional]")
parser$add_argument("-s", "--SRS_vcf_file_path",  metavar="FILE", help = "Path to short read duo vcf [Optional]")
parser$add_argument("-o", "--output_vcf",  metavar="FILE", help = "Path to file to write output vcf [Optional: Default appends _duoNovo to input vcf]")
parser$add_argument("-z", "--compress_output", action="store_true",  default=FALSE,  help="Compress output, append .bgz to filename, and tabix index [default FALSE]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

args$use_SRS <- FALSE
if (! is.null(args$SRS_vcf_file_path)){
   args$use_SRS <- TRUE
}

candidateCoords <- NULL
if (!is.null(args$candidate_variant_coordinates)){
     candidateCoords <- unlist(strsplit(args$candidate_variant_coordinates , '[,;]+'))
}

if (is.null(args$output_vcf) ) {
 args$output_vcf <- args$LRS_phased_vcf_file_path
 args$output_vcf <- sub('.vcf','_duoNovo.vcf',args$output_vcf )
}


#need to process output file

# print values of command line parameters to stderr if "quietly" wasn't requested
if ( args$verbose ) { 
    write(paste(names(args),args, sep=" --> "), stderr() )
}


duoNovo_results <- duoNovo(
  LRS_phased_vcf_file_path = args$LRS_phased_vcf_file_path, 
  depth_cutoff = args$depth_cutoff, 
  GQ_cutoff = args$GQ_cutoff,
  proband_column_identifier = args$proband_id,
  PS_width_cutoff = args$PS_width_cutoff, 
  boundary_cutoff = args$boundary_cutoff, 
  distance_cutoff = args$distance_cutoff,
  candidate_variants_concordant_with_SRS = args$use_SRS,
  SRS_vcf_file_path = args$SRS_vcf_file_path,
  test_reference_allele = args$test_reference_allele,
  candidate_variant_coordinates=candidateCoords,
  output_vcf_path=args$output_vcf,
  problematic_regions = args$exclude_regions_bed_path, 
  compress_output=args$compress_output
)
