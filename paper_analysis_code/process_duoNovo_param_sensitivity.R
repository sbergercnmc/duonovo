addOtherParentGenotype <- function(duoNovo_file_path, trio_vcf){
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
  parentValidation_trio_dp <- geno(trio_vcf)$DP[, parentValidation_trio_column]
  parentValidation_trio_gq <- geno(trio_vcf)$GQ[, parentValidation_trio_column]
  giab_problematic <- info(trio_vcf)[, "giab_problematic"]
  giab_problematic <- unlist(giab_problematic)
  names(giab_problematic) <- rownames(info(trio_vcf))
  
  duoNovo_ranges <- rowRanges(dnvcf)
  dnvcf_info <- info(dnvcf)
  columns_to_add <- c(
    "phasing_proband", "phasing_parent", "depth_proband", "depth_parent",
    "GQ_proband", "GQ_parent", "duoNovo_classification",
    "supporting_hamming_distance", "supporting_counts_het_hom", "supporting_counts_het_het",
    "supporting_counts_hom_het", "QC_fail_step", 
    "n_de_novo_left_orientation_same_PS", "n_de_novo_right_orientation_same_PS")
  mcols(duoNovo_ranges) <- dnvcf_info[, columns_to_add]

  indices_in_trio <- which(names(duoNovo_ranges) %in% names(parentValidation_trio_gt))
  duoNovo_ranges <- duoNovo_ranges[indices_in_trio]
  duoNovo_ranges$parentValidation_gt <- parentValidation_trio_gt[names(duoNovo_ranges)]
  duoNovo_ranges$parentValidation_depth <- parentValidation_trio_dp[names(duoNovo_ranges)]
  duoNovo_ranges$parentValidation_GQ <- parentValidation_trio_gq[names(duoNovo_ranges)]
  duoNovo_ranges$giab_problematic <- giab_problematic[names(duoNovo_ranges)]
  duoNovo_ranges
}

getPPV <- function(duoNovo_granges){
  dn_granges <- duoNovo_granges
  dn_granges <- dn_granges[which(dn_granges$duoNovo_classification == "de_novo" & 
                                   dn_granges$GQ_proband >= 40 & dn_granges$GQ_parent >= 40)]
  dn_granges <- dn_granges[which(dn_granges$giab_problematic == ".")]
  dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0")]
  dn_granges <- dn_granges[which(dn_granges$n_de_novo_left_orientation_same_PS == 1 | 
                                   dn_granges$n_de_novo_right_orientation_same_PS == 1)]
  total <- length(dn_granges)
  assessed_indices <- which(dn_granges$parentValidation_depth >= 20 & 
                              dn_granges$parentValidation_GQ >= 40 & 
                              !grepl("\\.", dn_granges$parentValidation_gt))
  assessed <- length(assessed_indices)
  
  dn_granges <- dn_granges[assessed_indices]
  dn <- length(dn_granges) - length(grep("1", dn_granges$parentValidation_gt))
  c(dn, assessed, total, length(dn_granges))
}


library(VariantAnnotation)
setwd('/scratch/sberger/pmgrc_lr_data/duoNovoOutputs_ANNOVAR')
directories <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
dir_index <- as.integer(args[1])
setwd(directories[dir_index])
trio_vcf_file_path <- list.files(pattern = "TRIO_PFM")
trio_vcf <- readVcf(trio_vcf_file_path)

output_dir <- '/home/lboukas/PMGRC_longread/duo/duoNovo_outputs/parameter_sensitivity'

setwd('/scratch/sberger/pmgrc_lr_data/duoNovoOutputsPARAMS')
setwd(directories[dir_index])
subdirectories <- list.files(pattern = "(D20_G30_W10000_C.*_B2000)|(D20_G30_W10000_C40_B.*)|(D.*_G30_W10000_C40_B2000)|(D20_G.*_W10000_C40_B2000)|(D20_G30_W5000_C40_B2000)|(D20_G30_W15000_C40_B2000))")


ppv_df_pf <- matrix(NA, nrow = 4, ncol = length(subdirectories))
rownames(ppv_df_pf) <- c("dn", "assessed", "total", "all_candidates")
colnames(ppv_df_pf) <- subdirectories
ppv_df_pm <- ppv_df_pf
for (i in 1:length(subdirectories)){
  setwd(subdirectories[i])
  duo_vcf <- list.files(pattern = "(DUO_PF|DUO_PM).*bgz$")
  
  dn_granges_pf <- addOtherParentGenotype(duo_vcf[1], trio_vcf)
  ppv_df_pf[, i] <- getPPV(dn_granges_pf)
  
  dn_granges_pm <- addOtherParentGenotype(duo_vcf[2], trio_vcf)
  ppv_df_pm[, i] <- getPPV(dn_granges_pm)
  
  setwd('/scratch/sberger/pmgrc_lr_data/duoNovoOutputsPARAMS')
  setwd(directories[dir_index])
}

setwd(subdirectories[1])
duo_vcf <- list.files(pattern = "(DUO_PF|DUO_PM).*bgz$")
output_filename_pf <- paste0(sub(".*Duo_(.*)_hiphase.*", "\\1", 
                                 basename(duo_vcf[1])), ".rda")
output_filename_pm <- paste0(sub(".*Duo_(.*)_hiphase.*", "\\1", 
                                 basename(duo_vcf[2])), ".rda")

save(ppv_df_pf, file = file.path(output_dir, output_filename_pf))
save(ppv_df_pm, file = file.path(output_dir, output_filename_pm))
