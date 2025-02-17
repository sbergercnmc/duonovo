library(MutationalPatterns)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)


###
###

setwd("C:/Users/lboukas/duoNovo_results/PF")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

sample_names <- sub(".*PF_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)
pf_dn_grl <- GRangesList(
  setNames(
    replicate(length(sample_names), GRanges(), simplify = FALSE),
    sample_names
))

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
  pf_dn_grl[[i]] <- dn_granges
}
genome(pf_dn_grl) <- "hg38"
chromosomes <- paste0('chr', c(1:22,'X', 'Y')) 
seqlevels(pf_dn_grl, pruning.mode = 'tidy') <- chromosomes

pf_dn_grl <- endoapply(pf_dn_grl, add_ref_alt)
pf_dn_grl_snvs <- get_mut_type(pf_dn_grl, type = "snv")
type_occurrences_pf <- mut_type_occurrences(pf_dn_grl_snvs, ref_genome)
type_occurrences_pf$`C>T` <- NULL
p1 <- plot_mutation_aggregate_percentages(type_occurrences_pf)



####
####
setwd("C:/Users/lboukas/duoNovo_results/PM")
all_duonovo_outputs <- list.files()
all_duonovo_outputs_pass_QC <- all_duonovo_outputs[pass_QC]

sample_names <- sub(".*PM_([^_]+)_.*", "\\1", all_duonovo_outputs_pass_QC)
pm_dn_grl <- GRangesList(
  setNames(
    replicate(length(sample_names), GRanges(), simplify = FALSE),
    sample_names
  ))

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
  pm_dn_grl[[i]] <- dn_granges
}
genome(pm_dn_grl) <- "hg38"
chromosomes <- paste0('chr', c(1:22,'X', 'Y')) 
seqlevels(pm_dn_grl, pruning.mode = 'tidy') <- chromosomes

pm_dn_grl <- endoapply(pm_dn_grl, add_ref_alt)
pm_dn_grl_snvs <- get_mut_type(pm_dn_grl, type = "snv")
type_occurrences_pm <- mut_type_occurrences(pm_dn_grl_snvs, ref_genome)
type_occurrences_pm$`C>T` <- NULL

p2 <- plot_mutation_aggregate_percentages(type_occurrences_pm)

library("gridExtra")
grid.arrange(p1, p2, ncol = 2, widths = c(3, 3))

pdf(file = "C:/Users/lboukas/duoNovo_figures/mutation_subtypes.pdf", 
    height = 3, width = 8, pointsize = 8)
grid.arrange(p1, p2, ncol = 2, widths = c(3, 3))
dev.off()


###
### Statistical test for difference in subtype proportions
###between father-proband duo SNVs and mother-proband duo SNVs
pf_all <- colSums(type_occurrences_pf, na.rm = TRUE)
pm_all <- colSums(type_occurrences_pm, na.rm = TRUE)
data_matrix <- rbind(pf_all, pm_all)
chisq.test(data_matrix)


###
### Indels
pf_dn_grl_indels <- get_mut_type(pf_dn_grl, type = "indel")
context_pf_indels <- get_indel_context(pf_dn_grl_indels, ref_genome)
type_occurrences_pf_indels <- count_indel_contexts(context_pf_indels)

pm_dn_grl_indels <- get_mut_type(pm_dn_grl, type = "indel")
type_occurrences_pm_indels <- mut_type_occurrences(pm_dn_grl_indels, ref_genome)
context_pm_indels <- get_indel_context(p_dn_grl_indels, ref_genome)
type_occurrences_pm_indels <- count_indel_contexts(context_pm_indels)













###
### Analysis with validated de novo variants only
pf_dn_grl_validated_only <- endoapply(pf_dn_grl, function(xx) {
  not_validated <- grep("1", xx$parentValidation_gt)
  if (length(not_validated) == 0){
    out <- xx
  }
  else {
    out <- xx[-not_validated]
  }
  out
})

pm_dn_grl_validated_only <- endoapply(pm_dn_grl, function(xx) {
  not_validated <- grep("1", xx$parentValidation_gt)
  if (length(not_validated) == 0){
    out <- xx
  }
  else {
    out <- xx[-not_validated]
  }
  out
})







