### Assess precision (pos predictive value)
getPrecision <- function(subdirectory_id, de_novo_candidates_filepath, 
                         other_parent_genotype_filepath, 
                         exclude = c(TRUE, FALSE), ranges_for_exclusion, 
                         include = c(TRUE, FALSE), ranges_for_inclusion,
                         exclude_mult_denovo_PS = c(TRUE, FALSE)){
  
    file_path_duo_novo <- paste0(de_novo_candidates_filepath, subdirectory_id, ".rda")
    load(file = file_path_duo_novo)
    duo_novo <- c(de_novo_candidates1$de_novo, de_novo_candidates2$de_novo)
    
    file_path_other_parent <- paste0(other_parent_genotype_filepath, subdirectory_id, ".rda")
    load(file = file_path_other_parent)
    other_parent_gt <- other_parent_genotype
  
  if (exclude == TRUE){
    seqlevels(duo_novo) <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)
    overlap_indices <- unique(queryHits(findOverlaps(duo_novo, ranges_for_exclusion)))
    if(length(overlap_indices) > 0){
      duo_novo <- duo_novo[-overlap_indices]
    }
  }
  if (include == TRUE){
    duo_novo <- duo_novo[unique(queryHits(findOverlaps(duo_novo, ranges_for_inclusion)))]
  }
  
  mult_denovo_PS <- unique(duo_novo$PS1[which(duplicated(duo_novo$PS1))])
  
  if (exclude_mult_denovo_PS == TRUE & length(mult_denovo_PS) > 0){
    duo_novo <- duo_novo[-which(duo_novo$PS1 %in% mult_denovo_PS)]
  }
  
  ###now compute precision
  dnv <- duo_novo[which(names(duo_novo) %in% names(other_parent_gt))]
  other_parent_gt_table <- table(other_parent_gt[names(dnv)])
  
  if ("0/0" %in% names(other_parent_gt_table)){
    output <- c(other_parent_gt_table["0/0"]/sum(other_parent_gt_table), 
                sum(other_parent_gt_table), length(duo_novo))
  } else {
    output <- c(0, sum(other_parent_gt_table), length(duo_novo))
  }
  output
}

duo_novo_performance_pos_f <- sapply(subdirs, function(xx) {
  getPrecision(xx, 
               de_novo_candidates_filepath = "C:\\Users\\lboukas\\alt_father_output_candidate_de_novo_",
               other_parent_genotype_filepath = "C:\\Users\\lboukas\\maternal_genotype_short_reads_",
               exclude = FALSE, include = FALSE, 
               exclude_mult_denovo_PS = TRUE)
})
rownames(duo_novo_performance_pos_f) <- c("precision", "total_assessed", "total")


duo_novo_performance_pos_f2 <- sapply(subdirs, function(xx) {
  getPrecision(xx, 
               de_novo_candidates_filepath = "C:\\Users\\lboukas\\alt_father_output_candidate_de_novo_",
               other_parent_genotype_filepath = "C:\\Users\\lboukas\\maternal_genotype_short_reads_",
               exclude = FALSE, include = TRUE, 
               ranges_for_inclusion = tx_all, 
               exclude_mult_denovo_PS = TRUE)
})
rownames(duo_novo_performance_pos_f2) <- c("precision", "total_assessed", "total")

duo_novo_performance_pos_f3 <- sapply(subdirs, function(xx) {
  getPrecision(xx, 
               de_novo_candidates_filepath = "C:\\Users\\lboukas\\alt_father_output_candidate_de_novo_",
               other_parent_genotype_filepath = "C:\\Users\\lboukas\\maternal_genotype_short_reads_",
               exclude = FALSE, include = TRUE, 
               ranges_for_inclusion = tx_all[which(tx_all$tx_biotype == "protein_coding")], 
               exclude_mult_denovo_PS = TRUE)
})
rownames(duo_novo_performance_pos_f3) <- c("precision", "total_assessed", "total")

duo_novo_performance_pos_f4 <- sapply(subdirs, function(xx) {
  getPrecision(xx, 
               de_novo_candidates_filepath = "C:\\Users\\lboukas\\alt_father_output_candidate_de_novo_",
               other_parent_genotype_filepath = "C:\\Users\\lboukas\\maternal_genotype_short_reads_",
               exclude = TRUE, include = FALSE, 
               ranges_for_exclusion = seg_dups_granges, 
               exclude_mult_denovo_PS = TRUE)
})
rownames(duo_novo_performance_pos_f4) <- c("precision", "total_assessed", "total")

duo_novo_performance_pos_f5 <- sapply(subdirs, function(xx) {
  getPrecision(xx, 
               de_novo_candidates_filepath = "C:\\Users\\lboukas\\alt_father_output_candidate_de_novo_",
               other_parent_genotype_filepath = "C:\\Users\\lboukas\\maternal_genotype_short_reads_",
               exclude = FALSE, include = FALSE, 
               exclude_mult_denovo_PS = FALSE)
})
rownames(duo_novo_performance_pos_f5) <- c("precision", "total_assessed", "total")

duo_novo_performance_pos_f6 <- sapply(subdirs, function(xx) {
  getPrecision(xx, 
               de_novo_candidates_filepath = "C:\\Users\\lboukas\\alt_father_output_candidate_de_novo_0cutoff_boundaries",
               other_parent_genotype_filepath = "C:\\Users\\lboukas\\maternal_genotype_short_reads_",
               exclude = FALSE, include = FALSE, 
               exclude_mult_denovo_PS = TRUE)
})
rownames(duo_novo_performance_pos_f6) <- c("precision", "total_assessed", "total")

duo_novo_performance_pos_f7 <- sapply(subdirs, function(xx) {
  getPrecision(xx, 
               de_novo_candidates_filepath = "C:\\Users\\lboukas\\alt_father_output_candidate_de_novo_1000cutoff_boundaries",
               other_parent_genotype_filepath = "C:\\Users\\lboukas\\maternal_genotype_short_reads_",
               exclude = FALSE, include = FALSE, 
               exclude_mult_denovo_PS = TRUE)
})
rownames(duo_novo_performance_pos_f7) <- c("precision", "total_assessed", "total")
plot(1- duo_novo_performance_pos_f["precision", ], 1 - duo_novo_performance_pos_f5["precision", ], 
     xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1)

###mother proband duos
duo_novo_performance_pos_m <- sapply(subdirs, function(xx) {
  getPrecision(xx, 
               de_novo_candidates_filepath = "C:\\Users\\lboukas\\alt_mother_output_candidate_de_novo_",
               other_parent_genotype_filepath = "C:\\Users\\lboukas\\paternal_genotype_short_reads_",
               exclude = FALSE, include = FALSE, 
               exclude_mult_denovo_PS = TRUE)
})
rownames(duo_novo_performance_pos_m) <- c("precision", "total_assessed", "total")


duo_novo_performance_pos_m2 <- sapply(subdirs, function(xx) {
  getPrecision(xx, 
               de_novo_candidates_filepath = "C:\\Users\\lboukas\\alt_mother_output_candidate_de_novo_",
               other_parent_genotype_filepath = "C:\\Users\\lboukas\\paternal_genotype_short_reads_",
               exclude = FALSE, include = TRUE, 
               ranges_for_inclusion = tx_all, 
               exclude_mult_denovo_PS = TRUE)
})
rownames(duo_novo_performance_pos_m2) <- c("precision", "total_assessed", "total")



duo_novo_performance_pos_m5 <- sapply(subdirs, function(xx) {
  getPrecision(xx, 
               de_novo_candidates_filepath = "C:\\Users\\lboukas\\alt_mother_output_candidate_de_novo_",
               other_parent_genotype_filepath = "C:\\Users\\lboukas\\paternal_genotype_short_reads_",
               exclude = FALSE, include = FALSE, 
               exclude_mult_denovo_PS = FALSE)
})
rownames(duo_novo_performance_pos_m5) <- c("precision", "total_assessed", "total")

duo_novo_performance_pos_m6 <- sapply(subdirs, function(xx) {
  getPrecision(xx, 
               de_novo_candidates_filepath = "C:\\Users\\lboukas\\alt_mother_output_candidate_de_novo_0cutoff_boundaries",
               other_parent_genotype_filepath = "C:\\Users\\lboukas\\paternal_genotype_short_reads_",
               exclude = FALSE, include = FALSE, 
               exclude_mult_denovo_PS = TRUE)
})
rownames(duo_novo_performance_pos_m6) <- c("precision", "total_assessed", "total")

duo_novo_performance_pos_m7 <- sapply(subdirs, function(xx) {
  getPrecision(xx, 
               de_novo_candidates_filepath = "C:\\Users\\lboukas\\alt_mother_output_candidate_de_novo_1000cutoff_boundaries",
               other_parent_genotype_filepath = "C:\\Users\\lboukas\\paternal_genotype_short_reads_",
               exclude = FALSE, include = FALSE, 
               exclude_mult_denovo_PS = TRUE)
})
rownames(duo_novo_performance_pos_m7) <- c("precision", "total_assessed", "total")

plot(1- duo_novo_performance_pos_m["precision", ], 1 - duo_novo_performance_pos_m5["precision", ], 
     xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1)




##########################
##########################
##########################
###Assess negative predictive value
duo_novo_performance_neg_f <- sapply(subdirs, function(xx) {
  file_path_other_parent <- paste0("C:\\Users\\lboukas\\maternal_genotype_short_reads_", xx, ".rda")
  load(file = file_path_other_parent)
  other_parent_gt <- other_parent_genotype
  
  file_path_duo_novo <- paste0("C:\\Users\\lboukas\\alt_father_output_candidate_de_novo_", xx, ".rda")
  load(file = file_path_duo_novo)
  
  duo_novo <- c(de_novo_candidates1$not_de_novo, de_novo_candidates2$not_de_novo)

  dnv <- duo_novo[which(names(duo_novo) %in% names(other_parent_gt))]
  other_parent_gt_table <- table(other_parent_gt[names(dnv)])
  
  duo_novo_baseline <- c(de_novo_candidates1$not_de_novo, de_novo_candidates2$not_de_novo, 
                         de_novo_candidates1$de_novo, de_novo_candidates2$de_novo, 
                         de_novo_candidates1$uncertain, de_novo_candidates2$uncertain)
  dnv_baseline <- duo_novo_baseline[which(names(duo_novo_baseline) %in% names(other_parent_gt))]
  other_parent_gt_table_baseline <- table(other_parent_gt[names(dnv_baseline)])
  if ("0/0" %in% names(other_parent_gt_table_baseline)){
    baseline_npv <- 1 - other_parent_gt_table_baseline["0/0"]/sum(other_parent_gt_table_baseline)
  } else {
    baseline_npv <- 1
  }
  
  if ("0/0" %in% names(other_parent_gt_table)){
    output <- c(1 - other_parent_gt_table["0/0"]/sum(other_parent_gt_table), 
                baseline_npv,
                sum(other_parent_gt_table), length(duo_novo))
  } else {
    output <- c(1, baseline_npv, 
                sum(other_parent_gt_table), length(duo_novo))
  }
  output
})

rownames(duo_novo_performance_neg_f) <- c("neg_pred_value", "baseline", "total_assessed", "total")









