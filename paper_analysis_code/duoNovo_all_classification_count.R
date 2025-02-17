all_classified <- rep(NA, 40)

for (i in 1:40){
  setwd("C:/Users/lboukas/duoNovo_results/PF")
  all_duonovo_outputs <- list.files()
  all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]
  
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pf
  dn_granges <- dn_granges[which(dn_granges$duoNovo_classification %in% 
                                   c("de_novo", "on_other_parent_haplotype") & 
                                   dn_granges$GQ_proband >= 40 & dn_granges$GQ_parent >= 40)]
  dn_granges$giab_problematic <- unlist(dn_granges$giab_problematic)
  dn_granges <- dn_granges[which(dn_granges$giab_problematic == ".")]
  dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0")]
  clustered <- which(dn_granges$n_de_novo_left_orientation_same_PS > 1 | 
                       dn_granges$n_de_novo_right_orientation_same_PS > 1)
  if (length(clustered) > 0){
    dn_granges <- dn_granges[-clustered]
  }  
  dn_granges_pf_duo <- dn_granges
  
  setwd("C:/Users/lboukas/duoNovo_results/PM")
  all_duonovo_outputs <- list.files()
  all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]
  
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pm
  dn_granges <- dn_granges[which(dn_granges$duoNovo_classification %in% 
                                   c("de_novo", "on_other_parent_haplotype") & 
                                   dn_granges$GQ_proband >= 40 & dn_granges$GQ_parent >= 40)]
  dn_granges$giab_problematic <- unlist(dn_granges$giab_problematic)
  dn_granges <- dn_granges[which(dn_granges$giab_problematic == ".")]
  dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0")]
  clustered <- which(dn_granges$n_de_novo_left_orientation_same_PS > 1 | 
                       dn_granges$n_de_novo_right_orientation_same_PS > 1)
  if (length(clustered) > 0){
    dn_granges <- dn_granges[-clustered]
  }
  dn_granges_pm_duo <- dn_granges
  
  all_classified[i] <- length(unique(c(dn_granges_pf_duo, dn_granges_pm_duo)))
}

