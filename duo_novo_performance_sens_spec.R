subdirs <- c("PMGRC-107-107-0", "PMGRC-323-323-0", "PMGRC-5-5-0",
             "PMGRC-191-191-0", "PMGRC-615-615-0",
             "PMGRC-291-291-0", "PMGRC-351-351-0", "PMGRC-645-645-0",
             "PMGRC-296-296-0", "PMGRC-372-372-0", "PMGRC-482-482-0", "PMGRC-661-661-0",
             "PMGRC-320-320-0", "PMGRC-417-417-0", "PMGRC-522-522-0", "PMGRC-770-770-0")

de_novo_count <- sapply(subdirs, function(xx) {
  file_path_de_novo <- paste0("C:\\Users\\lboukas\\ground_truth_de_novo_SR_", xx, ".rda")
  load(file = file_path_de_novo)
  keep <- which(trio_de_novo$child_depth >= 20 & trio_de_novo$father_depth >= 20 & 
                  trio_de_novo$mother_depth >= 20 & trio_de_novo$child_gq >= 30 & 
                  trio_de_novo$father_gq >= 30 & trio_de_novo$mother_gq >= 30)
  
  length(trio_de_novo[keep])
})

#the following is used to compute sensitivity
duo_novo_performance_vs_short_read_dn <- sapply(subdirs, function(xx) {
  file_path_duo_novo <- paste0("C:\\Users\\lboukas\\alt_father_output_candidate_de_novo_", xx, ".rda")
  load(file = file_path_duo_novo)
  duo_novo_f_dn <- c(de_novo_candidates1$de_novo, de_novo_candidates2$de_novo)
  duo_novo_f_ndn <- c(de_novo_candidates1$not_de_novo, de_novo_candidates2$not_de_novo)
  duo_novo_f_unc <- c(de_novo_candidates1$uncertain, de_novo_candidates2$uncertain)
  duo_novo_f_all <- c(duo_novo_f_dn, duo_novo_f_ndn, duo_novo_f_unc)
  
  file_path_duo_novo <- paste0("C:\\Users\\lboukas\\alt_mother_output_candidate_de_novo_", xx, ".rda")
  load(file = file_path_duo_novo)
  duo_novo_m_dn <- c(de_novo_candidates1$de_novo, de_novo_candidates2$de_novo)
  duo_novo_m_ndn <- c(de_novo_candidates1$not_de_novo, de_novo_candidates2$not_de_novo)
  duo_novo_m_unc <- c(de_novo_candidates1$uncertain, de_novo_candidates2$uncertain)
  duo_novo_m_all <- c(duo_novo_m_dn, duo_novo_m_ndn, duo_novo_m_unc)
  
  file_path_de_novo <- paste0("C:\\Users\\lboukas\\ground_truth_de_novo_SR_", xx, ".rda")
  load(file = file_path_de_novo)
  
  keep <- which(trio_de_novo$child_depth >= 20 & trio_de_novo$father_depth >= 20 & 
                  trio_de_novo$mother_depth >= 20 & trio_de_novo$child_gq >= 30 & 
                  trio_de_novo$father_gq >= 30 & trio_de_novo$mother_gq >= 30)
  
  trio_de_novo <- trio_de_novo[keep]
  
  overlaps_f_dn <- findOverlaps(trio_de_novo, duo_novo_f_dn)
  overlaps_f_ndn <- findOverlaps(trio_de_novo, duo_novo_f_ndn)
  overlaps_f_unc <- findOverlaps(trio_de_novo, duo_novo_f_unc)
  overlaps_f_all <- findOverlaps(trio_de_novo, duo_novo_f_all)
  
  overlaps_m_dn <- findOverlaps(trio_de_novo, duo_novo_m_dn)
  overlaps_m_ndn <- findOverlaps(trio_de_novo, duo_novo_m_ndn)
  overlaps_m_unc <- findOverlaps(trio_de_novo, duo_novo_m_unc)
  overlaps_m_all <- findOverlaps(trio_de_novo, duo_novo_m_all)
  
  dn_f <- trio_de_novo[unique(queryHits(overlaps_f_dn))]
  dn_m <- trio_de_novo[unique(queryHits(overlaps_m_dn))]
  dn_both <- intersect(names(dn_f), names(dn_m))
  
  ndn_f <- trio_de_novo[unique(queryHits(overlaps_f_ndn))]
  ndn_m <- trio_de_novo[unique(queryHits(overlaps_m_ndn))]
  ndn_both <- intersect(names(ndn_f), names(ndn_m))
  
  uncertain_f <- trio_de_novo[unique(queryHits(overlaps_f_unc))]
  uncertain_m <- trio_de_novo[unique(queryHits(overlaps_m_unc))]
  uncertain_both <- intersect(names(uncertain_f), names(uncertain_m))
  
  not_tested_f <- trio_de_novo[-unique(queryHits(overlaps_f_all))]
  not_tested_m <- trio_de_novo[-unique(queryHits(overlaps_m_all))]
  not_assessed_both <- intersect(names(not_tested_f), names(not_tested_m))
  
  uncertain_f_not_tested_m <- intersect(names(uncertain_f), names(not_tested_m))
  uncertain_m_not_tested_f <- intersect(names(uncertain_m), names(not_tested_f))
  
  uncertain_f_ndn_m <- intersect(names(uncertain_f), names(ndn_m))
  uncertain_m_ndn_f <- intersect(names(uncertain_m), names(ndn_f))
  
  ndn_f_not_tested_m <- intersect(names(ndn_f), names(not_tested_m))
  ndn_m_not_tested_f <- intersect(names(ndn_m), names(not_tested_f))
  
  c((length(unique(queryHits(overlaps_f_dn))) + length(unique(queryHits(overlaps_m_dn))))/length(trio_de_novo),
    (length(not_assessed_both) + length(uncertain_f_not_tested_m) + 
       length(uncertain_m_not_tested_f) + length(ndn_f_not_tested_m) + 
       length(ndn_m_not_tested_f))/length(trio_de_novo),
    (length(uncertain_both) + length(uncertain_f_ndn_m) + length(uncertain_m_ndn_f))/length(trio_de_novo),
    length(unique(queryHits(overlaps_f_dn)))/length(trio_de_novo),
    length(unique(queryHits(overlaps_m_dn)))/length(trio_de_novo),
    length(dn_both)/length(trio_de_novo), 
    length(ndn_both)/length(trio_de_novo), 
    length(not_assessed_both)/length(trio_de_novo),
    length(uncertain_both)/length(trio_de_novo),
    length(uncertain_f_not_tested_m)/length(trio_de_novo), 
    length(uncertain_m_not_tested_f)/length(trio_de_novo),
    length(uncertain_f_ndn_m)/length(trio_de_novo),
    length(uncertain_m_ndn_f)/length(trio_de_novo),
    length(ndn_f_not_tested_m)/length(trio_de_novo),
    length(ndn_m_not_tested_f)/length(trio_de_novo))
})
rownames(duo_novo_performance_vs_short_read_dn) <- c("de_novo_either", "not_tested", "uncertain",
                                                     "de_novo_f", "de_novo_m", 
                                                     "dn_in_both", "ndn_in_both",
                                                     "not_assessed_either", "uncertain_both",
                                                     "uncertain_f_not_tested_m", 
                                                     "uncertain_m_not_tested_f", 
                                                     "uncertain_f_ndn_m", "uncertain_m_ndn_f", 
                                                     "ndn_f_not_tested_m", "ndn_m_not_tested_f")




############# This computes the percentages after excluding variants not tested
############# in either duo because of QC 
duo_novo_performance_vs_short_read_dn2 <- sapply(subdirs, function(xx) {
  file_path_duo_novo <- paste0("C:\\Users\\lboukas\\alt_father_output_candidate_de_novo_", xx, ".rda")
  load(file = file_path_duo_novo)
  duo_novo_f_dn <- c(de_novo_candidates1$de_novo, de_novo_candidates2$de_novo)
  duo_novo_f_ndn <- c(de_novo_candidates1$not_de_novo, de_novo_candidates2$not_de_novo)
  duo_novo_f_unc <- c(de_novo_candidates1$uncertain, de_novo_candidates2$uncertain)
  duo_novo_f_all <- c(duo_novo_f_dn, duo_novo_f_ndn, duo_novo_f_unc)
  
  file_path_duo_novo <- paste0("C:\\Users\\lboukas\\alt_mother_output_candidate_de_novo_", xx, ".rda")
  load(file = file_path_duo_novo)
  duo_novo_m_dn <- c(de_novo_candidates1$de_novo, de_novo_candidates2$de_novo)
  duo_novo_m_ndn <- c(de_novo_candidates1$not_de_novo, de_novo_candidates2$not_de_novo)
  duo_novo_m_unc <- c(de_novo_candidates1$uncertain, de_novo_candidates2$uncertain)
  duo_novo_m_all <- c(duo_novo_m_dn, duo_novo_m_ndn, duo_novo_m_unc)
  
  file_path_de_novo <- paste0("C:\\Users\\lboukas\\ground_truth_de_novo_SR_", xx, ".rda")
  load(file = file_path_de_novo)
  
  overlaps_f_dn <- findOverlaps(trio_de_novo, duo_novo_f_dn)
  overlaps_f_ndn <- findOverlaps(trio_de_novo, duo_novo_f_ndn)
  overlaps_f_unc <- findOverlaps(trio_de_novo, duo_novo_f_unc)
  overlaps_f_all <- findOverlaps(trio_de_novo, duo_novo_f_all)
  
  overlaps_m_dn <- findOverlaps(trio_de_novo, duo_novo_m_dn)
  overlaps_m_ndn <- findOverlaps(trio_de_novo, duo_novo_m_ndn)
  overlaps_m_unc <- findOverlaps(trio_de_novo, duo_novo_m_unc)
  overlaps_m_all <- findOverlaps(trio_de_novo, duo_novo_m_all)
  
  dn_f <- trio_de_novo[unique(queryHits(overlaps_f_dn))]
  dn_m <- trio_de_novo[unique(queryHits(overlaps_m_dn))]
  dn_both <- intersect(names(dn_f), names(dn_m))
  
  ndn_f <- trio_de_novo[unique(queryHits(overlaps_f_ndn))]
  ndn_m <- trio_de_novo[unique(queryHits(overlaps_m_ndn))]
  ndn_both <- intersect(names(ndn_f), names(ndn_m))
  
  uncertain_f <- trio_de_novo[unique(queryHits(overlaps_f_unc))]
  uncertain_m <- trio_de_novo[unique(queryHits(overlaps_m_unc))]
  uncertain_both <- intersect(names(uncertain_f), names(uncertain_m))
  
  
  not_tested_f <- trio_de_novo[-unique(queryHits(overlaps_f_all))]
  not_tested_m <- trio_de_novo[-unique(queryHits(overlaps_m_all))]
  not_assessed_both <- intersect(names(not_tested_f), names(not_tested_m))
  
  uncertain_f_not_tested_m <- intersect(names(uncertain_f), names(not_tested_m))
  uncertain_m_not_tested_f <- intersect(names(uncertain_m), names(not_tested_f))
  
  uncertain_f_ndn_m <- intersect(names(uncertain_f), names(ndn_m))
  uncertain_m_ndn_f <- intersect(names(uncertain_m), names(ndn_f))
  
  ndn_f_not_tested_m <- intersect(names(ndn_f), names(not_tested_m))
  ndn_m_not_tested_f <- intersect(names(ndn_m), names(not_tested_f))
  
  trio_de_novo <- trio_de_novo[-which(names(trio_de_novo) %in% 
                                        c(not_assessed_both, uncertain_f_not_tested_m, 
                                          uncertain_m_not_tested_f, ndn_f_not_tested_m, 
                                          ndn_m_not_tested_f))]
  
  c((length(unique(queryHits(overlaps_f_dn))) + length(unique(queryHits(overlaps_m_dn))))/length(trio_de_novo),
    (length(uncertain_both) + length(uncertain_f_ndn_m) + length(uncertain_m_ndn_f))/length(trio_de_novo),
    length(unique(queryHits(overlaps_f_dn)))/length(trio_de_novo),
    length(unique(queryHits(overlaps_m_dn)))/length(trio_de_novo),
    length(dn_both)/length(trio_de_novo), 
    length(ndn_both)/length(trio_de_novo), 
    length(uncertain_both)/length(trio_de_novo),
    length(uncertain_f_ndn_m)/length(trio_de_novo),
    length(uncertain_m_ndn_f)/length(trio_de_novo))
})
rownames(duo_novo_performance_vs_short_read_dn2) <- c("de_novo_either", "uncertain",
                                                     "de_novo_f", "de_novo_m", 
                                                     "dn_in_both", "ndn_in_both",
                                                     "uncertain_both",
                                                     "uncertain_f_ndn_m", "uncertain_m_ndn_f")




#################################################
#################################################
duo_novo_performance_false_pos <- sapply(subdirs, function(xx) {
  file_path_duo_novo <- paste0("C:\\Users\\lboukas\\alt_father_output_candidate_de_novo_", xx, ".rda")
  load(file = file_path_duo_novo)
  duo_novo_f_dn <- c(de_novo_candidates1$de_novo, de_novo_candidates2$de_novo)
  duo_novo_f_ndn <- c(de_novo_candidates1$not_de_novo, de_novo_candidates2$not_de_novo)
  duo_novo_f_unc <- c(de_novo_candidates1$uncertain, de_novo_candidates2$uncertain)
  duo_novo_f_all <- c(duo_novo_f_dn, duo_novo_f_ndn, duo_novo_f_unc)
  
  suspicious_PS <- unique(duo_novo_f_dn$PS1[which(duplicated(duo_novo_f_dn$PS1))])
  if (length(suspicious_PS) > 0){
    duo_novo_f_dn <- duo_novo_f_dn[-which(duo_novo_f_dn$PS1 %in% suspicious_PS)]
  }
  
  file_path_duo_novo <- paste0("C:\\Users\\lboukas\\alt_mother_output_candidate_de_novo_", xx, ".rda")
  load(file = file_path_duo_novo)
  duo_novo_m_dn <- c(de_novo_candidates1$de_novo, de_novo_candidates2$de_novo)
  duo_novo_m_ndn <- c(de_novo_candidates1$not_de_novo, de_novo_candidates2$not_de_novo)
  duo_novo_m_unc <- c(de_novo_candidates1$uncertain, de_novo_candidates2$uncertain)
  duo_novo_m_all <- c(duo_novo_m_dn, duo_novo_m_ndn, duo_novo_m_unc)
  
  suspicious_PS <- unique(duo_novo_m_dn$PS1[which(duplicated(duo_novo_m_dn$PS1))])
  if (length(suspicious_PS) > 0){
    duo_novo_m_dn <- duo_novo_m_dn[-which(duo_novo_m_dn$PS1 %in% suspicious_PS)]
  }
  
  file_path_other_parent <- paste0("C:\\Users\\lboukas\\maternal_genotype_short_reads_", xx, ".rda")
  load(file = file_path_other_parent)
  maternal_gt <- other_parent_genotype
  
  file_path_other_parent <- paste0("C:\\Users\\lboukas\\paternal_genotype_short_reads_", xx, ".rda")
  load(file = file_path_other_parent)
  paternal_gt <- other_parent_genotype
  
  transmitted_variant_names_m <- names(maternal_gt)[
    which(maternal_gt %in% c("0/1", "1/0", "1/1"))]
  transmitted_variant_names_m <- transmitted_variant_names_m[
    which(transmitted_variant_names_m %in% names(duo_novo_f_all))]
  
  transmitted_variant_names_f <- names(paternal_gt)[
    which(paternal_gt %in% c("0/1", "1/0", "1/1"))]
  transmitted_variant_names_f <- transmitted_variant_names_f[
    which(transmitted_variant_names_f %in% names(duo_novo_m_all))]
  
  c(length(which(transmitted_variant_names_m %in% names(duo_novo_f_dn)))/length(transmitted_variant_names_m), 
    length(which(transmitted_variant_names_f %in% names(duo_novo_m_dn)))/length(transmitted_variant_names_f), 
    length(which(transmitted_variant_names_m %in% names(duo_novo_f_ndn)))/length(transmitted_variant_names_m), 
    length(which(transmitted_variant_names_f %in% names(duo_novo_m_ndn)))/length(transmitted_variant_names_f), 
    length(which(transmitted_variant_names_m %in% names(duo_novo_f_unc)))/length(transmitted_variant_names_m), 
    length(which(transmitted_variant_names_f %in% names(duo_novo_m_unc)))/length(transmitted_variant_names_f))
})
rownames(duo_novo_performance_false_pos) <- c("false_pos_f", "false_pos_m", 
                                              "true_neg_f", "true_neg_m", 
                                              "uncertain_f", "uncertain_m")




df <- as.data.frame(duo_novo_performance_false_pos)
df <- cbind(RowNames = rownames(df), df)
write_csv(df, file = "C:\\Users\\lboukas\\duo_novo_inherited_variants.csv")

df <- as.data.frame(duo_novo_performance_vs_short_read_dn2)
df <- cbind(RowNames = rownames(df), df)
write_csv(df, file = "C:\\Users\\lboukas\\duo_novo_SRdenovo_variants.csv")








performance_to_plot <- duo_novo_performance2[, -c(2, 10, 11)] #exclude these with too few or too many de novos from the full trio


pdf(file = "duo_novo_performance2.pdf", height = 3, width = 3, pointsize = 8)
plot(seq(1, 30, by = 2), performance_to_plot[1, ][order(performance_to_plot[1, ], decreasing = TRUE)], pch = 19, 
     cex = 0.75, bty = 'l', col = "dark orange", ylim = c(0, 1), xlab = "duos", xlim = c(0.8, 31.2),
     ylab = "% of true de novo variants", xaxt = 'n', yaxt = 'n')
points(seq(1.5, 30.5, by = 2), performance_to_plot[2, ][order(performance_to_plot[1, ], decreasing = TRUE)], pch = 19, 
     cex = 0.75, bty = 'l', col = "deep pink")
points(seq(2, 31, by = 2), performance_to_plot[3, ][order(performance_to_plot[1, ], decreasing = TRUE)], pch = 19, 
       cex = 0.75, bty = 'l', col = "forest green")
abline(v = seq(2.5, 29.5, by = 2), lty = "longdash")
axis(2, at = c(0, 0.5, 1))
legend <- legend("topright", legend = c("called de novo", "called not de novo", "uncertain call"), 
                 pch = 19, bty = 'n', cex = 0.84, col = c("dark orange", "deep pink", 
                                                          "forest green"))
dev.off()

#######################







###
sapply(subdirs, function(xx) {
  file_path_de_novo <- paste0("C:\\Users\\lboukas\\ground_truth_de_novo_SR_", xx, ".rda")
  load(file = file_path_de_novo)
  
  length(trio_de_novo)
})




