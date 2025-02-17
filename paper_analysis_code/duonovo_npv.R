setwd("C:/Users/lboukas/duoNovo_results/PF")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

npv_pf_alt <- matrix(NA, nrow = 3, ncol = length(all_duonovo_outputs_pass_QC))
rownames(npv_pf_alt) <- c("NPV", "naive_NPV", "total")
colnames(npv_pf_alt) <- sub(".*PF_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pf
  dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0")]
  dn_granges <- dn_granges[-which(dn_granges$QC_fail_step %in% 
                                    c("low_depth", "low_GQ", "low_depth_and_GQ"))]
  dn_granges$giab_problematic <- unlist(dn_granges$giab_problematic)
  dn_granges <- dn_granges[which(dn_granges$giab_problematic == ".")]
  
  dn_granges <- dn_granges[!grepl("\\.", dn_granges$parentValidation_gt)]
  
  true_dn <- grep("1", dn_granges$parentValidation_gt, invert = TRUE)
  naive_npv <- 1 - length(true_dn)/length(dn_granges)
  
  classified_ndn_1 <- which(dn_granges$duoNovo_classification == "on_other_parent_haplotype")
  classified_ndn_2 <- which(dn_granges$duoNovo_classification == "de_novo" & 
                              dn_granges$n_de_novo_left_orientation_same_PS > 1)
  classified_ndn_3 <- which(dn_granges$duoNovo_classification == "de_novo" & 
                              dn_granges$n_de_novo_right_orientation_same_PS > 1)
  
  classified_ndn_granges <- dn_granges[c(classified_ndn_1)]
  true_dn <- grep("1", classified_ndn_granges$parentValidation_gt, invert = TRUE)
  duonovo_npv <- 1 - length(true_dn)/length(classified_ndn_granges)
  total <- length(classified_ndn_granges)
  
  npv_pf_alt[, i] <- c(duonovo_npv, naive_npv, total)
  npv_pf_alt
}

setwd("C:/Users/lboukas/duoNovo_results/PM")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

npv_pm_alt <- matrix(NA, nrow = 3, ncol = length(all_duonovo_outputs_pass_QC))
rownames(npv_pm_alt) <- c("NPV", "naive_NPV", "total")
colnames(npv_pm_alt) <- sub(".*PF_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pm
  dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0")]
  dn_granges <- dn_granges[-which(dn_granges$QC_fail_step %in% 
                                    c("low_depth", "low_GQ", "low_depth_and_GQ"))]
  dn_granges$giab_problematic <- unlist(dn_granges$giab_problematic)
  dn_granges <- dn_granges[which(dn_granges$giab_problematic == ".")]
  
  dn_granges <- dn_granges[!grepl("\\.", dn_granges$parentValidation_gt)]
  
  true_dn <- grep("1", dn_granges$parentValidation_gt, invert = TRUE)
  naive_npv <- 1 - length(true_dn)/length(dn_granges)
  
  classified_ndn_1 <- which(dn_granges$duoNovo_classification == "on_other_parent_haplotype")
  classified_ndn_2 <- which(dn_granges$duoNovo_classification == "de_novo" & 
                              dn_granges$n_de_novo_left_orientation_same_PS > 1)
  classified_ndn_3 <- which(dn_granges$duoNovo_classification == "de_novo" & 
                              dn_granges$n_de_novo_right_orientation_same_PS > 1)
  classified_ndn_granges <- dn_granges[c(classified_ndn_1)]
  
  true_dn <- grep("1", classified_ndn_granges$parentValidation_gt, invert = TRUE)
  duonovo_npv <- 1 - length(true_dn)/length(classified_ndn_granges)
  total <- length(classified_ndn_granges)
  
  npv_pm_alt[, i] <- c(duonovo_npv, naive_npv, total)
  npv_pm_alt
}



