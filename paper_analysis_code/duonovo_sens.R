###all trios
setwd("C:/Users/lboukas/duoNovo_results/Trio")
all_trios <- list.files()
all_trios_pass_QC <- all_trios[pass_QC]


sens <- matrix(NA, nrow = 12, ncol = length(all_trios_pass_QC))
colnames(sens) <- sub(".*PF_([^_]+)_.*", "\\1", all_trios_pass_QC)
rownames(sens) <- c("dn_f", "ndn_f", "uncertain_f", "clustered_f", "failed_qc_f", "not_called_f", 
                    "dn_m", "ndn_m", "uncertain_m", "clustered_m", "failed_qc_m", "not_called_m")
setwd("C:/Users/lboukas/duoNovo_results/Trio")

for (i in 1:length(all_trios_pass_QC)){
  setwd("C:/Users/lboukas/duoNovo_results/Trio")
  load(file = all_trios_pass_QC[i])
  trio_de_novo <- trio_de_novo[which(trio_de_novo$proband_gq >= 40 & 
                                       trio_de_novo$parent1_gq >= 40 & 
                                       trio_de_novo$parent2_gq >= 40)]
  trio_de_novo$giab_problematic <- unlist(trio_de_novo$giab_problematic)
  trio_de_novo <- trio_de_novo[which(trio_de_novo$giab_problematic == ".")]
  
  setwd("C:/Users/lboukas/duoNovo_results/PF")
  all_duonovo_outputs <- list.files()
  all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pf

  overlaps_f <- unique(queryHits(findOverlaps(trio_de_novo, dn_granges)))
  if (length(overlaps_f) == length(trio_de_novo)){
    not_called_f <- 0
  } else {
    not_called_f <- length(trio_de_novo[-overlaps_f])
  }

  classified_dn <- which(dn_granges$duoNovo_classification == "de_novo" & 
                              (dn_granges$n_de_novo_left_orientation_same_PS == 1 | 
                                 dn_granges$n_de_novo_right_orientation_same_PS == 1))
  dn_f <- length(unique(queryHits(findOverlaps(trio_de_novo, dn_granges[classified_dn]))))
  
  classified_ndn <- which(dn_granges$duoNovo_classification == "on_other_parent_haplotype")
  ndn_f <- length(unique(queryHits(findOverlaps(
    trio_de_novo, dn_granges[classified_ndn]))))
  
  clustered_f_1 <- which(dn_granges$duoNovo_classification == "de_novo" & 
                              dn_granges$n_de_novo_left_orientation_same_PS > 1)
  clustered_f_2 <- which(dn_granges$duoNovo_classification == "de_novo" & 
                              dn_granges$n_de_novo_right_orientation_same_PS > 1)
  clustered_f <- length(unique(queryHits(findOverlaps(
    trio_de_novo, dn_granges[c(clustered_f_1, clustered_f_2)]))))

  uncertain_f <- length(unique(queryHits(findOverlaps(
    trio_de_novo, dn_granges[which(dn_granges$duoNovo_classification == "uncertain")]))))
  failed_qc_f <- length(unique(queryHits(findOverlaps(
    trio_de_novo, dn_granges[which(dn_granges$duoNovo_classification == "failed_QC")]))))

  
  setwd("C:/Users/lboukas/duoNovo_results/PM")
  all_duonovo_outputs <- list.files()
  all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pm
  
  overlaps_m <- unique(queryHits(findOverlaps(trio_de_novo, dn_granges)))
  if (length(overlaps_m) == length(trio_de_novo)){
    not_called_m <- 0
  } else {
    not_called_m <- length(trio_de_novo[-overlaps_m])
  }
  
  classified_dn <- which(dn_granges$duoNovo_classification == "de_novo" & 
                           (dn_granges$n_de_novo_left_orientation_same_PS == 1 | 
                              dn_granges$n_de_novo_right_orientation_same_PS == 1))
  dn_m <- length(unique(queryHits(findOverlaps(trio_de_novo, dn_granges[classified_dn]))))
  
  classified_ndn <- which(dn_granges$duoNovo_classification == "on_other_parent_haplotype")
  ndn_m <- length(unique(queryHits(findOverlaps(
    trio_de_novo, dn_granges[classified_ndn]))))
  
  clustered_m_1 <- which(dn_granges$duoNovo_classification == "de_novo" & 
                              dn_granges$n_de_novo_left_orientation_same_PS > 1)
  clustered_m_2 <- which(dn_granges$duoNovo_classification == "de_novo" & 
                              dn_granges$n_de_novo_right_orientation_same_PS > 1)
  clustered_m <- length(unique(queryHits(findOverlaps(
    trio_de_novo, dn_granges[c(clustered_m_1, clustered_m_2)]))))

  uncertain_m <- length(unique(queryHits(findOverlaps(
    trio_de_novo, dn_granges[which(dn_granges$duoNovo_classification == "uncertain")]))))
  failed_qc_m <- length(unique(queryHits(findOverlaps(
    trio_de_novo, dn_granges[which(dn_granges$duoNovo_classification == "failed_QC")]))))
  
  sens[, i] <- c(dn_f, ndn_f, uncertain_f, clustered_f, failed_qc_f, not_called_f, 
                 dn_m, ndn_m, uncertain_m, clustered_m, failed_qc_m, not_called_m)
  sens
}

aggregate <- rowSums(sens)
candidate_percentages_f <- aggregate[1:3]/sum(aggregate[1:3])
candidate_percentages_m <- aggregate[7:9]/sum(aggregate[7:9])


#####################
#####################
sens2 <- matrix(NA, nrow = 6, ncol = length(all_trios_pass_QC))
colnames(sens2) <- sub(".*PF_([^_]+)_.*", "\\1", all_trios_pass_QC)
rownames(sens2) <- c("dn_either", "ndn_both", "clustered_f", "clustered_m", "clustered_either", "total")
setwd("C:/Users/lboukas/duoNovo_results/Trio")

for (i in 1:length(all_trios_pass_QC)){
  setwd("C:/Users/lboukas/duoNovo_results/Trio")
  load(file = all_trios_pass_QC[i])
  trio_de_novo <- trio_de_novo[which(trio_de_novo$proband_gq >= 40 & 
                                       trio_de_novo$parent1_gq >= 40 & 
                                       trio_de_novo$parent2_gq >= 40)]
  trio_de_novo$giab_problematic <- unlist(trio_de_novo$giab_problematic)
  trio_de_novo <- trio_de_novo[which(trio_de_novo$giab_problematic == ".")]
  
  setwd("C:/Users/lboukas/duoNovo_results/PF")
  all_duonovo_outputs <- list.files()
  all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pf
  
  classified_dn <- which(dn_granges$duoNovo_classification == "de_novo" & 
                           (dn_granges$n_de_novo_left_orientation_same_PS == 1 | 
                              dn_granges$n_de_novo_right_orientation_same_PS == 1))
  trio_dn_dn_f <- trio_de_novo[unique(queryHits(findOverlaps(trio_de_novo, 
                                                             dn_granges[classified_dn])))]
  
  classified_ndn <- which(dn_granges$duoNovo_classification == "on_other_parent_haplotype")
  ndn_f_overlaps <- unique(queryHits(findOverlaps(
    trio_de_novo, dn_granges[classified_ndn])))
  
  clustered_f_1 <- which(dn_granges$duoNovo_classification == "de_novo" & 
                           dn_granges$n_de_novo_left_orientation_same_PS > 1)
  clustered_f_2 <- which(dn_granges$duoNovo_classification == "de_novo" & 
                           dn_granges$n_de_novo_right_orientation_same_PS > 1)
  clustered_f_indices <- unique(queryHits(findOverlaps(
    trio_de_novo, dn_granges[c(clustered_f_1, clustered_f_2)])))
  clustered_f <- length(clustered_f_indices)  

  setwd("C:/Users/lboukas/duoNovo_results/PM")
  all_duonovo_outputs <- list.files()
  all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pm
  
  classified_dn <- which(dn_granges$duoNovo_classification == "de_novo" & 
                           (dn_granges$n_de_novo_left_orientation_same_PS == 1 | 
                              dn_granges$n_de_novo_right_orientation_same_PS == 1))
  trio_dn_dn_m <- trio_de_novo[unique(queryHits(findOverlaps(trio_de_novo, 
                                                             dn_granges[classified_dn])))]
  
  classified_ndn <- which(dn_granges$duoNovo_classification == "on_other_parent_haplotype")
  ndn_m_overlaps <- unique(queryHits(findOverlaps(
    trio_de_novo, dn_granges[classified_ndn])))
  
  clustered_m_1 <- which(dn_granges$duoNovo_classification == "de_novo" & 
                           dn_granges$n_de_novo_left_orientation_same_PS > 1)
  clustered_m_2 <- which(dn_granges$duoNovo_classification == "de_novo" & 
                           dn_granges$n_de_novo_right_orientation_same_PS > 1)
  clustered_m_indices <- unique(queryHits(findOverlaps(
    trio_de_novo, dn_granges[c(clustered_m_1, clustered_m_2)])))
  clustered_m <- length(clustered_m_indices)
  
  
  dn_either <- length(unique(c(trio_dn_dn_f, trio_dn_dn_m)))
  ndn_both <- length(intersect(ndn_f_overlaps, ndn_m_overlaps))
  clustered_either <- length(unique(clustered_f_indices, clustered_m_indices))
  
  sens2[, i] <- c(dn_either, ndn_both, clustered_f, clustered_m, 
                  clustered_either, length(trio_de_novo))
  sens2
}

sensitivity <- apply(sens2, 2, function(xx) xx['dn_either']/(xx['total'] - xx['clustered_either']))
false_ndn <- apply(sens2, 2, function(xx) xx['ndn_both']/(xx['total'] - xx['clustered_either']))
other <- apply(sens2, 2, function(xx) 
  (xx['total'] - xx['clustered_either'] - xx['dn_either'] - xx['ndn_both'])/(xx['total'] - xx['clustered_either']))

n_denovo <- apply(sens2, 2, function(xx) xx['total'] - xx['clustered_either'])
plot(jitter(rep(1, 40), 2), n_denovo, xlim = c(0.8, 1.2), ylim = c(0, 157), bty = 'l', yaxt = 'n')
axis(2, at = c(0, 70, 140))

