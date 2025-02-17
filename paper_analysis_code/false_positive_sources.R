########## Without exluding giab problematic regions
##########
setwd("C:/Users/lboukas/duoNovo_results/PF")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

ppv_pf_alt2 <- matrix(NA, nrow = 4, ncol = length(all_duonovo_outputs_pass_QC))
rownames(ppv_pf_alt2) <- c("dn", "assessed", "total", "all_candidates")
colnames(ppv_pf_alt2) <- sub(".*PF_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pf
  dn_granges <- dn_granges[which(dn_granges$duoNovo_classification == "de_novo" & 
                                   dn_granges$GQ_proband >= 40 & dn_granges$GQ_parent >= 40)]
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
  ppv_pf_alt2[, i] <- c(dn, assessed, total, length(dn_granges_pf))
  ppv_pf_alt2
}

##
setwd("C:/Users/lboukas/duoNovo_results/PM")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

ppv_pm_alt2 <- matrix(NA, nrow = 4, ncol = length(all_duonovo_outputs_pass_QC))
rownames(ppv_pm_alt2) <- c("dn", "assessed", "total", "all_candidates")
colnames(ppv_pm_alt2) <- sub(".*PM_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pm
  dn_granges <- dn_granges[which(dn_granges$duoNovo_classification == "de_novo" & 
                                   dn_granges$GQ_proband >= 40 & dn_granges$GQ_parent >= 40)]
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
  ppv_pm_alt2[, i] <- c(dn, assessed, total, length(dn_granges_pm))
  ppv_pm_alt2
}



############# Restricting to exonic/intronic regions only
#############
setwd("C:/Users/lboukas/duoNovo_results/PF")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

ppv_pf_alt3 <- matrix(NA, nrow = 4, ncol = length(all_duonovo_outputs_pass_QC))
rownames(ppv_pf_alt3) <- c("dn", "assessed", "total", "all_candidates")
colnames(ppv_pf_alt3) <- sub(".*PF_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)

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
  dn_granges$Func.refGeneWithVer <- unlist(dn_granges$Func.refGeneWithVer)
  dn_granges <- dn_granges[which(dn_granges$Func.refGeneWithVer %in% 
                                   c("exonic", "intronic", "ncRNA_exonic", "ncRNA_intronic", 
                                     "UTR3", "UTR5"))]
  
  total <- length(dn_granges)
  assessed_indices <- which(dn_granges$parentValidation_depth >= 20 & 
                              dn_granges$parentValidation_GQ >= 40 & 
                              !grepl("\\.", dn_granges$parentValidation_gt))
  assessed <- length(assessed_indices)
  
  dn_granges <- dn_granges[assessed_indices]
  #precision <- 1 - length(grep("1", dn_granges$parentValidation_gt))/length(dn_granges)
  
  dn <- length(dn_granges) - length(grep("1", dn_granges$parentValidation_gt))
  ppv_pf_alt3[, i] <- c(dn, assessed, total, length(dn_granges_pf))
  ppv_pf_alt3
}

##
setwd("C:/Users/lboukas/duoNovo_results/PM")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

ppv_pm_alt3 <- matrix(NA, nrow = 4, ncol = length(all_duonovo_outputs_pass_QC))
rownames(ppv_pm_alt3) <- c("dn", "assessed", "total", "all_candidates")
colnames(ppv_pm_alt3) <- sub(".*PM_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)

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
  dn_granges$Func.refGeneWithVer <- unlist(dn_granges$Func.refGeneWithVer)
  dn_granges <- dn_granges[which(dn_granges$Func.refGeneWithVer %in% 
                                   c("exonic", "intronic", "ncRNA_exonic", "ncRNA_intronic", 
                                     "UTR3", "UTR5"))]
  
  total <- length(dn_granges)
  assessed_indices <- which(dn_granges$parentValidation_depth >= 20 & 
                              dn_granges$parentValidation_GQ >= 40 & 
                              !grepl("\\.", dn_granges$parentValidation_gt))
  assessed <- length(assessed_indices)
  
  dn_granges <- dn_granges[assessed_indices]
  #precision <- 1 - length(grep("1", dn_granges$parentValidation_gt))/length(dn_granges)
  dn <- length(dn_granges) - length(grep("1", dn_granges$parentValidation_gt))
  ppv_pm_alt3[, i] <- c(dn, assessed, total, length(dn_granges_pm))
  ppv_pm_alt3
}


############### Without excluding clustered de novos
###############
setwd("C:/Users/lboukas/duoNovo_results/PF")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

ppv_pf_alt4 <- matrix(NA, nrow = 4, ncol = length(all_duonovo_outputs_pass_QC))
rownames(ppv_pf_alt4) <- c("dn", "assessed", "total", "all_candidates")
colnames(ppv_pf_alt4) <- sub(".*PF_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pf
  dn_granges <- dn_granges[which(dn_granges$duoNovo_classification == "de_novo" & 
                                   dn_granges$GQ_proband >= 40 & dn_granges$GQ_parent >= 40)]
  dn_granges$giab_problematic <- unlist(dn_granges$giab_problematic)
  dn_granges <- dn_granges[which(dn_granges$giab_problematic == ".")]
  dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0")]
  total <- length(dn_granges)
  assessed_indices <- which(dn_granges$parentValidation_depth >= 20 & 
                              dn_granges$parentValidation_GQ >= 40 & 
                              !grepl("\\.", dn_granges$parentValidation_gt))
  assessed <- length(assessed_indices)
  
  dn_granges <- dn_granges[assessed_indices]
  #precision <- 1 - length(grep("1", dn_granges$parentValidation_gt))/length(dn_granges)
  dn <- length(dn_granges) - length(grep("1", dn_granges$parentValidation_gt))
  ppv_pf_alt4[, i] <- c(dn, assessed, total, length(dn_granges_pf))
  ppv_pf_alt4
}

##
setwd("C:/Users/lboukas/duoNovo_results/PM")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

ppv_pm_alt4 <- matrix(NA, nrow = 4, ncol = length(all_duonovo_outputs_pass_QC))
rownames(ppv_pm_alt4) <- c("dn", "assessed", "total", "all_candidates")
colnames(ppv_pm_alt4) <- sub(".*PM_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pm
  dn_granges <- dn_granges[which(dn_granges$duoNovo_classification == "de_novo" & 
                                   dn_granges$GQ_proband >= 40 & dn_granges$GQ_parent >= 40)]
  dn_granges$giab_problematic <- unlist(dn_granges$giab_problematic)
  dn_granges <- dn_granges[which(dn_granges$giab_problematic == ".")]
  dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0")]

  total <- length(dn_granges)
  assessed_indices <- which(dn_granges$parentValidation_depth >= 20 & 
                              dn_granges$parentValidation_GQ >= 40 & 
                              !grepl("\\.", dn_granges$parentValidation_gt))
  assessed <- length(assessed_indices)
  
  dn_granges <- dn_granges[assessed_indices]
  #precision <- 1 - length(grep("1", dn_granges$parentValidation_gt))/length(dn_granges)
  dn <- length(dn_granges) - length(grep("1", dn_granges$parentValidation_gt))
  ppv_pm_alt4[, i] <- c(dn, assessed, total, length(dn_granges_pm))
  ppv_pm_alt4
}


getCollectivePPV <- function(ppv_df){
  dn_count_total <-  sum(ppv_df['dn', ])
  assessed_total <-  sum(ppv_df['assessed', ])
  dn_count_total/assessed_total
}

collective_ppv_pf <- getCollectivePPV(ppv_pf_alt) #from duonovo_ppv.R script
collective_ppv_pm <- getCollectivePPV(ppv_pm_alt)

collective_ppv_pf_with_giab_problematic <- getCollectivePPV(ppv_pf_alt2)
collective_ppv_pm_with_giab_problematic <- getCollectivePPV(ppv_pm_alt2)

collective_ppv_pf_genes_only <- getCollectivePPV(ppv_pf_alt3)
collective_ppv_pm_genes_only <- getCollectivePPV(ppv_pm_alt3)

collective_ppv_pf_with_clustered <- getCollectivePPV(ppv_pf_alt4)
collective_ppv_pm_with_clustered <- getCollectivePPV(ppv_pm_alt4)

