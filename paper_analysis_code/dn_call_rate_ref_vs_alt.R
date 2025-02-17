setwd("C:/Users/lboukas/duoNovo_results/PF")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

dn_call_rate_alt_pf <- matrix(NA, nrow = 3, ncol = length(all_duonovo_outputs_pass_QC))
colnames(dn_call_rate_alt_pf) <- sub(".*PF_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)
rownames(dn_call_rate_alt_pf) <- c("dn_count", "total_count", "rate")
  
dn_call_rate_ref_pf <- dn_call_rate_alt_pf

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pf
  
  dn_granges_alt <- dn_granges[which(dn_granges$phasing_parent == "0/0" & dn_granges$GQ_parent >= 40)]
  dn_granges_alt <- dn_granges_alt[-which(dn_granges_alt$duoNovo_classification == "failed_QC")]
  dn_granges_alt$giab_problematic <- unlist(dn_granges_alt$giab_problematic)
  dn_granges_alt <- dn_granges_alt[which(dn_granges_alt$giab_problematic == ".")]

  dn_alt_1 <- which(dn_granges_alt$duoNovo_classification == "de_novo" & 
                      dn_granges_alt$n_de_novo_left_orientation_same_PS == 1 & 
                      dn_granges_alt$GQ_proband >= 40)
  dn_alt_2 <- which(dn_granges_alt$duoNovo_classification == "de_novo" & 
                      dn_granges_alt$n_de_novo_right_orientation_same_PS == 1 & 
                      dn_granges_alt$GQ_proband >= 40)
  dn_call_rate_alt_pf[3, i] <- (length(dn_alt_1) + length(dn_alt_2))/length(dn_granges_alt)
  dn_call_rate_alt_pf[1, i] <- length(dn_alt_1) + length(dn_alt_2)
  dn_call_rate_alt_pf[2, i] <- length(dn_granges_alt)
  
  dn_granges_ref <- dn_granges[which(dn_granges$phasing_parent == "1/1" & dn_granges$GQ_parent >= 40)]
  dn_granges_ref <- dn_granges_ref[-which(dn_granges_ref$duoNovo_classification == "failed_QC")]
  dn_granges_ref$giab_problematic <- unlist(dn_granges_ref$giab_problematic)
  dn_granges_ref <- dn_granges_ref[which(dn_granges_ref$giab_problematic == ".")]
  
  dn_ref_1 <- which(dn_granges_ref$duoNovo_classification == "de_novo" & 
                      dn_granges_ref$n_de_novo_left_orientation_same_PS == 1 & 
                      dn_granges_ref$GQ_proband >= 40)
  dn_ref_2 <- which(dn_granges_ref$duoNovo_classification == "de_novo" & 
                      dn_granges_ref$n_de_novo_right_orientation_same_PS == 1 & 
                      dn_granges_ref$GQ_proband >= 40)
  dn_call_rate_ref_pf[3, i] <- (length(dn_ref_1) + length(dn_ref_2))/length(dn_granges_ref)
  dn_call_rate_ref_pf[1, i] <- length(dn_ref_1) + length(dn_ref_2)
  dn_call_rate_ref_pf[2, i] <- length(dn_granges_ref)
  
  dn_call_rate_alt_pf
  dn_call_rate_ref_pf
}

#######

setwd("C:/Users/lboukas/duoNovo_results/PM")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

dn_call_rate_alt_pm <- matrix(NA, nrow = 3, ncol = length(all_duonovo_outputs_pass_QC))
colnames(dn_call_rate_alt_pm) <- sub(".*PF_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)
rownames(dn_call_rate_alt_pm) <- c("dn_count", "total_count", "rate")

dn_call_rate_ref_pm <- dn_call_rate_alt_pm

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pm
  
  dn_granges_alt <- dn_granges[which(dn_granges$phasing_parent == "0/0" & dn_granges$GQ_parent >= 40)]
  dn_granges_alt <- dn_granges_alt[-which(dn_granges_alt$duoNovo_classification == "failed_QC")]
  dn_granges_alt$giab_problematic <- unlist(dn_granges_alt$giab_problematic)
  dn_granges_alt <- dn_granges_alt[which(dn_granges_alt$giab_problematic == ".")]
  
  dn_alt_1 <- which(dn_granges_alt$duoNovo_classification == "de_novo" & 
                      dn_granges_alt$n_de_novo_left_orientation_same_PS == 1 & 
                      dn_granges_alt$GQ_proband >= 40)
  dn_alt_2 <- which(dn_granges_alt$duoNovo_classification == "de_novo" & 
                      dn_granges_alt$n_de_novo_right_orientation_same_PS == 1 & 
                      dn_granges_alt$GQ_proband >= 40)
  dn_call_rate_alt_pm[3, i] <- (length(dn_alt_1) + length(dn_alt_2))/length(dn_granges_alt)
  dn_call_rate_alt_pm[1, i] <- length(dn_alt_1) + length(dn_alt_2)
  dn_call_rate_alt_pm[2, i] <- length(dn_granges_alt)
  
  dn_granges_ref <- dn_granges[which(dn_granges$phasing_parent == "1/1" & dn_granges$GQ_parent >= 40)]
  dn_granges_ref <- dn_granges_ref[-which(dn_granges_ref$duoNovo_classification == "failed_QC")]
  dn_granges_ref$giab_problematic <- unlist(dn_granges_ref$giab_problematic)
  dn_granges_ref <- dn_granges_ref[which(dn_granges_ref$giab_problematic == ".")]
  
  dn_ref_1 <- which(dn_granges_ref$duoNovo_classification == "de_novo" & 
                      dn_granges_ref$n_de_novo_left_orientation_same_PS == 1 & 
                      dn_granges_ref$GQ_proband >= 40)
  dn_ref_2 <- which(dn_granges_ref$duoNovo_classification == "de_novo" & 
                      dn_granges_ref$n_de_novo_right_orientation_same_PS == 1 & 
                      dn_granges_ref$GQ_proband >= 40)
  dn_call_rate_ref_pm[3, i] <- (length(dn_ref_1) + length(dn_ref_2))/length(dn_granges_ref)
  dn_call_rate_ref_pm[1, i] <- length(dn_ref_1) + length(dn_ref_2)
  dn_call_rate_ref_pm[2, i] <- length(dn_granges_ref)
  
  dn_call_rate_alt_pm
  dn_call_rate_ref_pm
}




