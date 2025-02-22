addOtherParentGenotype <- function(duoNovo_file_path, trio_vcf_file_path){
  if (!file.exists(duoNovo_file_path)) {
    stop("The duoNovo VCF file path does not exist: ", duoNovo_file_path )
  }
  if (!file.exists(trio_vcf_file_path)) {
    stop("The Trio VCF file path does not exist: ", trio_vcf_file_path )
  }
  
  message("Reading duoNovo VCF file...")
  dnvcf <- readVcf(duoNovo_file_path)
  if (is.null(meta(header(dnvcf))$duoNovoPARAM)){
    stop("DuoNovoHeader Parameters not found")
  }
  
  proband_column_identifier <- meta(header(dnvcf))$duoNovoPARAM['proband_column_identifier', ]
  proband_column <- grep(proband_column_identifier, samples(header(dnvcf)))
  if (length(proband_column) != 1) {
    stop("Proband column identifier not found unambiguously in duoNovo Vcf: ", proband_column_identifier)
  }
  
  proband_id <- samples(header(dnvcf))[proband_column]
  parent_id <- samples(header(dnvcf))[-proband_column]
  if (!"GT" %in% names(geno(dnvcf)) ) {
    stop("The 'GT' (genotype) field is missing in the VCF metadata from LRS.")
  }
  
  message("Reading trio VCF file...")
  trio_vcf <- readVcf(trio_vcf_file_path)
  if (length(samples(header(trio_vcf))) != 3) {
    stop("Trio VCF does not contain exactly 3 samples. Number of samples: ", length(samples(header(trio_vcf))))
  }
  proband_trio_column <- grep(proband_column_identifier, samples(header(trio_vcf)))
  if (length(proband_trio_column) != 1) {
    stop("Proband column identifier not found unambiguously in trio Vcf: ", proband_column_identifier)
  }
  parentTested_trio_column <- grep(parent_id, samples(header(trio_vcf)))
  if (length(parentTested_trio_column) != 1) {
    stop("Tested Parent column identifier not found in duoNovo Vcf: ", proband_column_identifier)
  }
  parentValidation_trio_column <- setdiff(seq(3), c(proband_trio_column, parentTested_trio_column))
  parentValidation_trio_gt <- geno(trio_vcf)$GT[, parentValidation_trio_column]
  parentValidation_trio_ad <- geno(trio_vcf)$AD[, parentValidation_trio_column]
  parentValidation_trio_dp <- geno(trio_vcf)$DP[, parentValidation_trio_column]
  parentValidation_trio_gq <- geno(trio_vcf)$GQ[, parentValidation_trio_column]
  
  dnvcf_info <- info(dnvcf)
  duoNovo_ranges <- rowRanges(dnvcf)
  columns_to_add <- c(
    "phasing_proband", "phasing_parent", "depth_proband", "depth_parent",
    "GQ_proband", "GQ_parent", "duoNovo_classification",
    "supporting_hamming_distance", "supporting_counts_het_hom", "supporting_counts_het_het",
    "supporting_counts_hom_het", "QC_fail_step", 
    "n_de_novo_left_orientation_same_PS", "n_de_novo_right_orientation_same_PS",
    "ANNOVAR_DATE",
    "Func.refGeneWithVer", "Gene.refGeneWithVer", "GeneDetail.refGeneWithVer",
    "ExonicFunc.refGeneWithVer", "AAChange.refGeneWithVer",
    "CpG",
    "giab_problematic",
    "gnomad41_genome_AF", "gnomad41_genome_AF_raw", "gnomad41_genome_AF_XX",
    "gnomad41_genome_AF_XY", "gnomad41_genome_AF_grpmax",
    "gnomad41_genome_faf95", "gnomad41_genome_faf99",
    "gnomad41_genome_fafmax_faf95_max", "gnomad41_genome_fafmax_faf99_max",
    "gnomad41_genome_AF_afr", "gnomad41_genome_AF_ami",
    "gnomad41_genome_AF_amr", "gnomad41_genome_AF_asj",
    "gnomad41_genome_AF_eas", "gnomad41_genome_AF_fin",
    "gnomad41_genome_AF_mid", "gnomad41_genome_AF_nfe",
    "gnomad41_genome_AF_remaining", "gnomad41_genome_AF_sas",
    "ALLELE_END")
  mcols(duoNovo_ranges) <- dnvcf_info[, columns_to_add]
  
  duoNovo_ranges$AD_proband <-  geno(dnvcf)$AD[, proband_column]
  duoNovo_ranges$AD_parent <-  geno(dnvcf)$AD[, -proband_column]
  
  indices_in_trio <- which(names(duoNovo_ranges) %in% names(parentValidation_trio_gt))
  duoNovo_ranges <- duoNovo_ranges[indices_in_trio]
  duoNovo_ranges$parentValidation_gt <- parentValidation_trio_gt[names(duoNovo_ranges)]
  duoNovo_ranges$parentValidation_ad <- parentValidation_trio_ad[names(duoNovo_ranges)]
  duoNovo_ranges$parentValidation_depth <- parentValidation_trio_dp[names(duoNovo_ranges)]
  duoNovo_ranges$parentValidation_GQ <- parentValidation_trio_gq[names(duoNovo_ranges)]
  duoNovo_ranges
}



library(VariantAnnotation)
setwd('/scratch/sberger/pmgrc_lr_data/duoNovoOutputs_ANNOVAR')
directories <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
dir_index <- as.integer(args[1])
setwd(directories[dir_index])

duo_vcf <- list.files(pattern = "DUO_PF|DUO_PM")
trio_vcf <- list.files(pattern = "TRIO_PFM")


dn_granges_pf <- addOtherParentGenotype(duo_vcf[1], trio_vcf)
dn_granges_pm <- addOtherParentGenotype(duo_vcf[2], trio_vcf)


output_dir <- '/home/lboukas/PMGRC_longread/duo/duoNovo_outputs'

duo_basename_pf <- basename(duo_vcf[1])
duo_basename_pm <- basename(duo_vcf[2])

output_filename_pf <- paste0(sub("\\.vcf$", "", duo_basename_pf), ".rda")
output_filename_pm <- paste0(sub("\\.vcf$", "", duo_basename_pm), ".rda")

output_file_pf <- file.path(output_dir, output_filename_pf)
output_file_pm <- file.path(output_dir, output_filename_pm)

save(dn_granges_pf, file = output_file_pf)
save(dn_granges_pm, file = output_file_pm)


