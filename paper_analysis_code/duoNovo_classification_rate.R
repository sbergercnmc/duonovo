setwd("C:/Users/lboukas/duoNovo_results/PF")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

dn_call_rate_alt_pf <- rep(NA, length(all_duonovo_outputs_pass_QC))
names(dn_call_rate_alt_pf) <- sub(".*PF_(.*?)_hiphase.*", "\\1", all_duonovo_outputs_pass_QC)


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
  ndn_alt <- which(dn_granges_alt$duoNovo_classification == "on_other_parent_haplotype")
  dn_call_rate_alt_pf[i] <- (length(dn_alt_1) + length(dn_alt_2) + length(ndn_alt))/length(dn_granges_alt)
  dn_call_rate_alt_pf
}

setwd("C:/Users/lboukas/duoNovo_results/PM")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

dn_call_rate_alt_pm <- rep(NA, length(all_duonovo_outputs_pass_QC))
names(dn_call_rate_alt_pm) <- sub(".*PM_(.*?)_hiphase.*", "\\1", all_duonovo_outputs_pass_QC)


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
  ndn_alt <- which(dn_granges_alt$duoNovo_classification == "on_other_parent_haplotype")
  dn_call_rate_alt_pm[i] <- (length(dn_alt_1) + length(dn_alt_2) + length(ndn_alt))/length(dn_granges_alt)
  dn_call_rate_alt_pm
}
