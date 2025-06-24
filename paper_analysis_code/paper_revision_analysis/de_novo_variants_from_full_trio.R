getDeNovosFromTrio <- function(trio_vcf_filepath){
  if (!file.exists(trio_vcf_filepath)) {
    stop("The trio VCF file path does not exist: ", trio_vcf_filepath)
  }
  
  message("Reading trio VCF file...")
  trio_vcf <- readVcf(trio_vcf_filepath)
  if (length(samples(header(trio_vcf))) != 3) {
    stop("Trio VCF does not contain exactly 3 samples. Number of samples: ", length(samples(header(trio_vcf))))
  }
  proband_column <- 1
  if (!"GT" %in% names(geno(trio_vcf)) ) {
    stop("The 'GT' (genotype) field is missing in the VCF metadata.")
  }
  if (!"DNM" %in% names(geno(trio_vcf)) ) {
    stop("The 'DNM' field is missing in the VCF metadata.")
  }
  
  trio_ranges <- rowRanges(trio_vcf)
  vcf_info <- info(trio_vcf)
  columns_to_add <- c("cpg", "problematic_region")
  
  missing_info <- setdiff(columns_to_add, colnames(vcf_info))
  if (length(missing_info))
    stop("VCF INFO is missing fields: ", paste(missing_info, collapse=", "))
  
  mcols(trio_ranges) <- vcf_info[, columns_to_add]
  
  g <- geno(trio_vcf)
  trio_ranges$proband_gt <- g$GT[, proband_column]
  trio_ranges$proband_AD <- g$AD[, proband_column]
  trio_ranges$proband_dp <- g$DP[, proband_column]
  trio_ranges$proband_GQ <- g$GQ[, proband_column]
  
  trio_ranges$parent1_gt <- g$GT[, 2]
  trio_ranges$parent1_AD <- g$AD[, 2]
  trio_ranges$parent1_dp <- g$DP[, 2]
  trio_ranges$parent1_GQ <- g$GQ[, 2]
  
  trio_ranges$parent2_gt <- g$GT[, 3]
  trio_ranges$parent2_AD <- g$AD[, 3]
  trio_ranges$parent2_dp <- g$DP[, 3]
  trio_ranges$parent2_GQ <- g$GQ[, 3]
  
  trio_ranges$DNM <- g$DNM[, 1]
  
  de_novo_indices <- which(trio_ranges$DNM == 1)
  if (length(de_novo_indices) > 0){
    trio_de_novo <- trio_ranges[de_novo_indices]
  } else {
    trio_de_novo <- GRanges()
  }
  trio_de_novo
}

library(VariantAnnotation)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: process_duoNovo_output.R <run_directory>")
}

## set directory
current_dir <- args[1]
setwd(current_dir)


trio_vcf_filepath <- list.files(pattern = "PFM\\.annovar\\.dnm2\\.vcf\\.gz$")
if (length(trio_vcf_filepath) != 1L ) {
  stop("Couldn’t unambiguously detect trio vcf file")
}
trio_de_novo <- getDeNovosFromTrio(trio_vcf_filepath)

trio_de_novo_basename <- basename(trio_vcf_filepath)
output_filename <- paste0(sub("\\.vcf\\.gz$", "", trio_de_novo_basename), ".rda")
save(trio_de_novo, file = output_filename)

