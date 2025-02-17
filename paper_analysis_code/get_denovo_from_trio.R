library(VariantAnnotation)
setwd('/scratch/sberger/pmgrc_lr_data/duoNovoOutputs_ANNOVAR')
directories <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
dir_index <- as.integer(args[1])
setwd(directories[dir_index])

trio_vcf_filepath <- list.files(pattern = "TRIO_PFM")
duo_vcf_filepath <- list.files(pattern = "DUO_PM") #only needed to get the proband ID


trio_vcf <- readVcf(trio_vcf_filepath)
duo_vcf <- readVcf(duo_vcf_filepath)

if (length(samples(header(trio_vcf))) != 3) {
  stop("Trio VCF does not contain exactly 3 samples. Number of samples: ", length(samples(header(trio_vcf))))
}

proband_column_identifier <- meta(header(duo_vcf))$duoNovoPARAM['proband_column_identifier', ]
trio_proband_column <- grep(proband_column_identifier, samples(header(trio_vcf)))
trio_parent_columns <- setdiff(seq(3), trio_proband_column)

trio_ranges <- rowRanges(trio_vcf)
vcf_info <- info(trio_vcf)
mcols(trio_ranges) <- vcf_info[, c("CpG", "giab_problematic")]

trio_ranges$proband_gt <- geno(trio_vcf)$GT[, trio_proband_column]
trio_ranges$proband_vaf <- geno(trio_vcf)$VAF[, trio_proband_column]
trio_ranges$proband_dp <- geno(trio_vcf)$DP[, trio_proband_column]
trio_ranges$proband_gq <- geno(trio_vcf)$GQ[, trio_proband_column]

trio_ranges$parent1_gt <- geno(trio_vcf)$GT[, trio_parent_columns[1]]
trio_ranges$parent1_vaf <- geno(trio_vcf)$VAF[, trio_parent_columns[1]]
trio_ranges$parent1_dp <- geno(trio_vcf)$DP[, trio_parent_columns[1]]
trio_ranges$parent1_gq <- geno(trio_vcf)$GQ[, trio_parent_columns[1]]

trio_ranges$parent2_gt <- geno(trio_vcf)$GT[, trio_parent_columns[2]]
trio_ranges$parent2_vaf <- geno(trio_vcf)$VAF[, trio_parent_columns[2]]
trio_ranges$parent2_dp <- geno(trio_vcf)$DP[, trio_parent_columns[2]]
trio_ranges$parent2_gq <- geno(trio_vcf)$GQ[, trio_parent_columns[2]]

de_novo_indices <- which(trio_ranges$proband_gt %in% c("1|0", "0|1", "1/0", "0/1") & 
                           trio_ranges$parent1_gt == "0/0" & 
                           trio_ranges$parent2_gt == "0/0" & 
                           trio_ranges$proband_dp >= 20 & trio_ranges$parent1_dp >= 20 & 
                               trio_ranges$parent2_dp >= 20 & 
                           trio_ranges$proband_gq >= 30 & trio_ranges$parent1_gq >= 30 &
                               trio_ranges$parent2_gq >= 30)

if (length(de_novo_indices) == 0){
  trio_de_novo <- GenomicRanges()
} else {
  trio_de_novo <- trio_ranges[de_novo_indices]
}

output_dir <- '/home/lboukas/PMGRC_longread/duo/duoNovo_outputs/trio_de_novo'

trio_basename <- basename(trio_vcf_filepath)

output_filename <- paste0(sub("\\.vcf$", "", trio_basename), ".rda")
output_file <- file.path(output_dir, output_filename)
save(trio_de_novo, file = output_file)

