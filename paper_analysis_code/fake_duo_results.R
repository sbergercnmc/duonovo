setwd("C:/Users/lboukas/duoNovo_results/Parents")
all_duonovo_outputs <- list.files(pattern = 'FM')
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[-grep("UCI-064|UCI-084", all_duonovo_outputs)] #exclude problematic maternal samples

dn_call_rate_alt_fm <- rep(NA, length(all_duonovo_outputs_pass_QC))
names(dn_call_rate_alt_fm) <- sub(".*FM_(.*?)_hiphase.*", "\\1", all_duonovo_outputs_pass_QC)

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_fm
  
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
  dn_call_rate_alt_fm[i] <- (length(dn_alt_1) + length(dn_alt_2) + length(ndn_alt))/length(dn_granges_alt)
  dn_call_rate_alt_fm
}

####
####
all_duonovo_outputs <- list.files(pattern = 'MF')
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[-grep("UCI-064|UCI-084", all_duonovo_outputs)] #exclude problematic maternal samples

dn_call_rate_alt_mf <- rep(NA, length(all_duonovo_outputs_pass_QC))
names(dn_call_rate_alt_mf) <- sub(".*MF_(.*?)_hiphase.*", "\\1", all_duonovo_outputs_pass_QC)

for (i in 1:length(all_duonovo_outputs_pass_QC)){
  load(file = all_duonovo_outputs_pass_QC[i])
  dn_granges <- dn_granges_mf
  
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
  dn_call_rate_alt_mf[i] <- (length(dn_alt_1) + length(dn_alt_2) + length(ndn_alt))/length(dn_granges_alt)
  dn_call_rate_alt_mf
}


###
### Evaluate these results in the context of relatedness between parents
somalier <- read_tsv('C:\\Users\\lboukas\\somalier.pairs.tsv')
colnames(somalier)[1] <- "sample_a"
somalier$pair1 <- paste0(somalier$sample_a, "_", somalier$sample_b)
somalier$pair2 <- paste0(somalier$sample_b, "_", somalier$sample_a)


split_names <- strsplit(names(dn_call_rate_alt_fm), "_")
df_fm <- data.frame(
  father_ID = sapply(split_names, `[`, 1),
  mother_ID = sapply(split_names, `[`, 2),
  stringsAsFactors = FALSE
)

df_fm$relatedness <- unlist(sapply(1:dim(df_fm)[1], function(xx) 
  somalier$relatedness[which(somalier$pair1 == names(dn_call_rate_alt_fm)[xx] | 
                               somalier$pair2 == names(dn_call_rate_alt_fm)[xx])]))
df_fm$hom_concordance <- unlist(sapply(1:dim(df_fm)[1], function(xx) 
  somalier$hom_concordance[which(somalier$pair1 == names(dn_call_rate_alt_fm)[xx] | 
                               somalier$pair2 == names(dn_call_rate_alt_fm)[xx])]))
df_fm$IBS2 <- unlist(sapply(1:dim(df_fm)[1], function(xx) 
  somalier$ibs2[which(somalier$pair1 == names(dn_call_rate_alt_fm)[xx] | 
                                   somalier$pair2 == names(dn_call_rate_alt_fm)[xx])]))
df_fm$ibd_rate <- dn_call_rate_alt_fm

df_fm$call_rate_pf <- sapply(df_fm$father_ID, 
                             function(xx) dn_call_rate_alt_pf[grep(xx, names(dn_call_rate_alt_pf))])
df_fm$call_rate_pm <- sapply(df_fm$mother_ID, 
                             function(xx) dn_call_rate_alt_pm[grep(xx, names(dn_call_rate_alt_pm))])

pdf(file = "C:/Users/lboukas/duoNovo_figures/parent_relatedness.pdf", 
    height = 2.4, width = 1.8, pointsize = 8)
vec <- df_fm$relatedness
boxplot(vec, staplewex = 0, medlty = 0, ylim = c(-0.8, 0), 
        yaxt = 'n', bty = 'l', ylab = "relatedness")
axis(2, at = c(-0.8, -0.4, 0))
dev.off()



### Before the next plot, reorder names to match the same fake duo
reverse_name <- function(name) {
  parts <- strsplit(name, "_")[[1]]  # Split the name into parts
  reversed <- paste(rev(parts), collapse = "_")  # Reverse and concatenate
  return(reversed)
}
reversed_names <- sapply(names(dn_call_rate_alt_fm), reverse_name)
dn_call_rate_alt_mf_reordered <- dn_call_rate_alt_mf[reversed_names]


pdf(file = "C:/Users/lboukas/duoNovo_figures/fake_duos.pdf", 
    height = 2.4, width = 7, pointsize = 8)
par(mfrow = c(1, 2))
plot(1:38, 100*dn_call_rate_alt_fm,
     pch = 19, col = "deep pink", xlab = "fake duos", cex = 0.8,
     ylab = "% variants classified", 
     bty = 'l', xaxt = 'n', yaxt = 'n', xlim = c(0.8, 40), ylim = c(0, 70))
axis(2, at = c(0, 20, 40, 60))
abline(h = 100*median(c(dn_call_rate_alt_pf, dn_call_rate_alt_pm)), 
       lty = "longdash", col = rgb(0,0,0,0.4))
#abline(v = seq(1.5, 39.5, by = 1), lty = "longdash", col = rgb(0,0,0,0.4))

plot(1:38, 100*dn_call_rate_alt_mf_reordered, 
     pch = 19, col = "deep pink", xlab = "fake duos", cex = 0.8,
     ylab = "% variants classified", 
     bty = 'l', xaxt = 'n', yaxt = 'n', xlim = c(0.8, 40), ylim = c(0, 70))
axis(2, at = c(0, 20, 40, 60))
abline(h = 100*median(c(dn_call_rate_alt_pf, dn_call_rate_alt_pm)), 
       lty = "longdash", col = rgb(0,0,0,0.4))
dev.off()





