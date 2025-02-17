###Father - proband
setwd("C:/Users/lboukas/duoNovo_results/PF")
all_duonovo_outputs <- list.files()

failed_qc_pf <- rep(NA, length(all_duonovo_outputs))
for (i in 1:length(all_duonovo_outputs)){
  load(file = all_duonovo_outputs[i])
  gr <- dn_granges_pf
  tab <- prop.table(table(gr$duoNovo_classification))
  failed_qc_pf[i] <- tab['failed_QC'] > 0.75
  failed_qc_pf
}


### Mother-proband
setwd("C:/Users/lboukas/duoNovo_results/PM")
all_duonovo_outputs <- list.files()

failed_qc_pm <- rep(NA, length(all_duonovo_outputs))
for (i in 1:length(all_duonovo_outputs)){
  load(file = all_duonovo_outputs[i])
  gr <- dn_granges_pm
  tab <- prop.table(table(gr$duoNovo_classification))
  failed_qc_pm[i] <- tab['failed_QC'] > 0.75
  failed_qc_pm
}

### PPV assessment
pass_QC <- failed_qc_pf == FALSE & failed_qc_pm == FALSE

setwd("C:/Users/lboukas/duoNovo_results/PF")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

ppv_pf_alt <- matrix(NA, nrow = 4, ncol = length(all_duonovo_outputs_pass_QC))
rownames(ppv_pf_alt) <- c("dn", "assessed", "total", "all_candidates")
colnames(ppv_pf_alt) <- sub(".*PF_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pf
  dn_granges <- dn_granges[which(dn_granges$duoNovo_classification == "de_novo" & 
                                   dn_granges$GQ_proband >= 40 & dn_granges$GQ_parent >= 40)]
  dn_granges$giab_problematic <- unlist(dn_granges$giab_problematic)
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
  #precision <- 1 - length(grep("1", dn_granges$parentValidation_gt))/length(dn_granges)
  dn <- length(dn_granges) - length(grep("1", dn_granges$parentValidation_gt))
  ppv_pf_alt[, i] <- c(dn, assessed, total, length(dn_granges_pf))
  ppv_pf_alt
}

##
setwd("C:/Users/lboukas/duoNovo_results/PM")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

ppv_pm_alt <- matrix(NA, nrow = 4, ncol = length(all_duonovo_outputs_pass_QC))
rownames(ppv_pm_alt) <- c("dn", "assessed", "total", "all_candidates")
colnames(ppv_pm_alt) <- sub(".*PM_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pm
  dn_granges <- dn_granges[which(dn_granges$duoNovo_classification == "de_novo" & 
                                   dn_granges$GQ_proband >= 40 & dn_granges$GQ_parent >= 40)]
  dn_granges$giab_problematic <- unlist(dn_granges$giab_problematic)
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
  #precision <- 1 - length(grep("1", dn_granges$parentValidation_gt))/length(dn_granges)
  
  dn <- length(dn_granges) - length(grep("1", dn_granges$parentValidation_gt))
  ppv_pm_alt[, i] <- c(dn, assessed, total, length(dn_granges_pm))
  ppv_pm_alt
}



