setwd("C:/Users/lboukas/duoNovo_results/PF")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

false_pos_pf <- matrix(NA, nrow = 3, ncol = length(all_duonovo_outputs_pass_QC))
colnames(false_pos_pf) <- sub(".*PF_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)
rownames(false_pos_pf) <- c("dn", "ndn", "uncertain")

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pf
  dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0" & dn_granges$GQ_parent >= 40)]
  inherited_from_mother <- grep("1", dn_granges$parentValidation_gt)
  all_inherited <- dn_granges[inherited_from_mother]
  all_inherited <- all_inherited[-which(all_inherited$duoNovo_classification == "failed_QC")]
  all_inherited$giab_problematic <- unlist(all_inherited$giab_problematic)
  all_inherited <- all_inherited[which(all_inherited$giab_problematic == ".")]
  
  all_inherited <- all_inherited[which(all_inherited$parentValidation_GQ >= 40)]
  
  classified_dn <- which(all_inherited$duoNovo_classification == "de_novo" & 
                           all_inherited$GQ_proband >= 40 &
                         (all_inherited$n_de_novo_left_orientation_same_PS == 1 | 
                             all_inherited$n_de_novo_right_orientation_same_PS == 1))
  dn <- length(classified_dn)
  
  classified_ndn <- which(all_inherited$duoNovo_classification == "on_other_parent_haplotype")
  clustered_1 <- which(all_inherited$duoNovo_classification == "de_novo" & 
                         all_inherited$n_de_novo_left_orientation_same_PS > 1)
  clustered_2 <- which(all_inherited$duoNovo_classification == "de_novo" & 
                         all_inherited$n_de_novo_right_orientation_same_PS > 1)
  ndn <- length(classified_ndn)
  
  uncertain <- length(which(dn_granges$duoNovo_classification == "uncertain"))
  uncertain2 <- length(which(all_inherited$duoNovo_classification == "de_novo" & 
                        all_inherited$GQ_proband < 40 & 
                      (all_inherited$n_de_novo_left_orientation_same_PS == 1 | 
                          all_inherited$n_de_novo_right_orientation_same_PS == 1)))
  
  false_pos_pf[, i] <- c(dn, ndn, uncertain + uncertain2)
  false_pos_pf
}

###

setwd("C:/Users/lboukas/duoNovo_results/PM")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

false_pos_pm <- matrix(NA, nrow = 3, ncol = length(all_duonovo_outputs_pass_QC))
colnames(false_pos_pm) <- sub(".*PF_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)
rownames(false_pos_pm) <- c("dn", "ndn", "uncertain")

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pm
  dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0" & dn_granges$GQ_parent >= 40)]
  inherited_from_father <- grep("1", dn_granges$parentValidation_gt)
  all_inherited <- dn_granges[inherited_from_father]
  all_inherited <- all_inherited[-which(all_inherited$duoNovo_classification == "failed_QC")]
  all_inherited$giab_problematic <- unlist(all_inherited$giab_problematic)
  all_inherited <- all_inherited[which(all_inherited$giab_problematic == ".")]
  
  all_inherited <- all_inherited[which(all_inherited$parentValidation_GQ >= 40)]
  
  classified_dn <- which(all_inherited$duoNovo_classification == "de_novo" & 
                           all_inherited$GQ_proband >= 40 & all_inherited$GQ_parent >= 40 &
                           (all_inherited$n_de_novo_left_orientation_same_PS == 1 | 
                              all_inherited$n_de_novo_right_orientation_same_PS == 1))
  dn <- length(classified_dn)
  
  classified_ndn <- which(all_inherited$duoNovo_classification == "on_other_parent_haplotype")
  clustered_1 <- which(all_inherited$duoNovo_classification == "de_novo" & 
                              all_inherited$n_de_novo_left_orientation_same_PS > 1)
  clustered_2 <- which(all_inherited$duoNovo_classification == "de_novo" & 
                              all_inherited$n_de_novo_right_orientation_same_PS > 1)
  ndn <- length(classified_ndn)
  
  uncertain <- length(which(dn_granges$duoNovo_classification == "uncertain"))
  uncertain2 <- length(which(all_inherited$duoNovo_classification == "de_novo" & 
                               all_inherited$GQ_proband < 40 & 
                               (all_inherited$n_de_novo_left_orientation_same_PS == 1 | 
                                  all_inherited$n_de_novo_right_orientation_same_PS == 1)))  
  false_pos_pm[, i] <- c(dn, ndn, uncertain + uncertain2)
  false_pos_pm
}


aggregate_transmitted_f <- rowSums(false_pos_pf)
percentages_transmitted_f <- aggregate_transmitted_f/sum(aggregate_transmitted_f)
per_sample_false_pos_rate_f <- apply(false_pos_pf, 2, function(xx) xx[1]/sum(xx))

aggregate_transmitted_m <- rowSums(false_pos_pm)
percentages_transmitted_m <- aggregate_transmitted_m/sum(aggregate_transmitted_m)
per_sample_false_pos_rate_m <- apply(false_pos_pm, 2, function(xx) xx[1]/sum(xx))




