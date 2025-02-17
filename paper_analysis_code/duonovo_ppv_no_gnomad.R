setwd("C:/Users/lboukas/duoNovo_results/PF")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

ppv_pf_alt_no_gnomad <- matrix(NA, nrow = 4, ncol = length(all_duonovo_outputs_pass_QC))
rownames(ppv_pf_alt_no_gnomad) <- c("dn", "assessed", "total", "all_candidates")
colnames(ppv_pf_alt_no_gnomad) <- sub(".*PF_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pf
  dn_granges$gnomad41_genome_AF <- as.numeric(unlist(dn_granges$gnomad41_genome_AF))
  dn_granges <- dn_granges[-which(dn_granges$gnomad41_genome_AF > 0)]
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
  dn <- length(dn_granges) - length(grep("1", dn_granges$parentValidation_gt))
  ppv_pf_alt_no_gnomad[, i] <- c(dn, assessed, total, length(dn_granges_pf))
  ppv_pf_alt_no_gnomad
}


setwd("C:/Users/lboukas/duoNovo_results/PM")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

ppv_pm_alt_no_gnomad <- matrix(NA, nrow = 4, ncol = length(all_duonovo_outputs_pass_QC))
rownames(ppv_pm_alt_no_gnomad) <- c("dn", "assessed", "total", "all_candidates")
colnames(ppv_pm_alt_no_gnomad) <- sub(".*PF_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_pm
  dn_granges$gnomad41_genome_AF <- as.numeric(unlist(dn_granges$gnomad41_genome_AF))
  dn_granges <- dn_granges[-which(dn_granges$gnomad41_genome_AF > 0)]
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
  dn <- length(dn_granges) - length(grep("1", dn_granges$parentValidation_gt))
  ppv_pm_alt_no_gnomad[, i] <- c(dn, assessed, total, length(dn_granges_pf))
  ppv_pm_alt_no_gnomad
}



pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_father_proband_ppv_no_gnomad.pdf", 
    height = 2.2, width = 2.2, pointsize = 8)

precision <- ppv_pf_alt_no_gnomad["dn", ]/ppv_pf_alt_no_gnomad["assessed", ]

plot(sort(precision, decreasing = TRUE), pch = 19, 
     cex = 0.5, bty = 'l', col = "dark orange", ylim = c(0, 1),
     yaxt = 'n', xaxt = 'n', ylab = "positive predictive value", xlab = "father-proband duos")
axis(2, at = c(0, 0.5, 1))
dev.off()

pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_father_proband_n_called_no_gnomad.pdf", 
    height = 1.5, width = 2.2, pointsize = 8)

number_called <- ppv_pf_alt_no_gnomad["assessed", ]

plot(number_called[order(precision, decreasing = TRUE)], pch = 19, type = 'h', lwd = 1,
     cex = 0.5, bty = 'l', col = "dark orange", ylim = c(0, max(number_called)), 
     yaxt = 'n', xaxt = 'n', ylab = "# de novos called", xlab = "father-proband duos")
axis(2, at = c(0, 12, 24))
dev.off()

pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_mother_proband_ppv_no_gnomad.pdf", 
    height = 2.2, width = 2.2, pointsize = 8)

precision <- ppv_pm_alt_no_gnomad["dn", ]/ppv_pm_alt_no_gnomad["assessed", ]

plot(sort(precision, decreasing = TRUE), pch = 19, 
     cex = 0.5, bty = 'l', col = "dark orange", ylim = c(0, 1),
     yaxt = 'n', xaxt = 'n', ylab = "positive predictive value", xlab = "father-proband duos")
axis(2, at = c(0, 0.5, 1))
dev.off()

pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_mother_proband_n_called_no_gnomad.pdf", 
    height = 1.5, width = 2.2, pointsize = 8)

number_called <- ppv_pm_alt_no_gnomad["assessed", ]

plot(number_called[order(precision, decreasing = TRUE)], pch = 19, type = 'h', lwd = 1,
     cex = 0.5, bty = 'l', col = "dark orange", ylim = c(0, max(number_called)), 
     yaxt = 'n', xaxt = 'n', ylab = "# de novos called", xlab = "father-proband duos")
axis(2, at = c(0, 4, 8))
dev.off()





