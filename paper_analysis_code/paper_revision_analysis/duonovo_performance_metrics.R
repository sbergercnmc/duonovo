getDuoNovoPerformanceMetric <- function(duoNovo_granges_output_filepath, duo_type = c("PM", "PF"), 
                   include_gnomad = c(TRUE, FALSE), allele_type = c("ALT", "REF"),
                   metric = c("PPV", "NPV", "false positive rate", 
                              "de novo call rate", "n of total variants classified"), 
                   validation_GQ_cutoff = 30){
  load(file = duoNovo_granges_output_filepath)
  
  if (duo_type == "PF"){
    dn_granges <- dn_granges_pf
  } else if (duo_type == "PM"){
    dn_granges <- dn_granges_pm
  }
  
  if (metric == "PPV"){
    
    if (include_gnomad == FALSE){
      dn_granges$gnomad41_genome_AF <- as.numeric(unlist(dn_granges$gnomad41_genome_AF))
      dn_granges <- dn_granges[-which(dn_granges$gnomad41_genome_AF > 0)]
    }
    dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0" & dn_granges$GQ_parent >= validation_GQ_cutoff)]
    dn_granges <- dn_granges[!grepl("\\.", dn_granges$parentValidation_gt)]
    all_candidates <- length(dn_granges)
    
    classified_dn <- which(dn_granges$duoNovo_classification == "de_novo" & 
                             dn_granges$GQ_proband >= validation_GQ_cutoff)
    if (length(classified_dn) > 0){
      dn_granges <- dn_granges[classified_dn]
      total <- length(dn_granges)
      assessed_indices <- which(dn_granges$parentValidation_depth >= 20 & 
                                  dn_granges$parentValidation_GQ >= validation_GQ_cutoff)
      if (length(assessed_indices) > 0){
        assessed <- length(assessed_indices)
        dn_granges <- dn_granges[assessed_indices]
        dn <- length(dn_granges) - length(grep("1", dn_granges$parentValidation_gt))
      } else {
        dn <- NA
        assessed <- 0
      }
    } else {
      dn <- NA
      assessed <- NA
      total <- 0
    }
    ppv <- c(dn, assessed, total, all_candidates)
    names(ppv) <- c("dn", "assessed", "total", "all_candidates")
    out <- ppv
  } else if (metric == "NPV"){
    
    dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0" & dn_granges$GQ_parent >= validation_GQ_cutoff)]
    
    ## -- QC fails that naive approach would discard as well
    dn_granges <- dn_granges[!(dn_granges$QC_fail_step %in% 
                                      c("low_depth", "low_GQ", "low_depth_and_GQ"))] 
    problematic <- grep("problematic_region", dn_granges$QC_fail_step)
    dn_granges  <- dn_granges[-problematic]
    ##
    
    dn_granges <- dn_granges[which(dn_granges$parentValidation_depth >= 20 & 
                                     dn_granges$parentValidation_GQ >= validation_GQ_cutoff)]
    dn_granges <- dn_granges[!grepl("\\.", dn_granges$parentValidation_gt)]
    
    true_dn <- grep("1", dn_granges$parentValidation_gt, invert = TRUE)
    naive_npv <- 1 - length(true_dn)/length(dn_granges)
    
    classified_ndn <- which(dn_granges$duoNovo_classification == "on_other_parent_haplotype")
    classified_ndn_granges <- dn_granges[c(classified_ndn)]
    true_dn <- grep("1", classified_ndn_granges$parentValidation_gt, invert = TRUE)
    duonovo_npv <- 1 - length(true_dn)/length(classified_ndn_granges)

    npv <- c(duonovo_npv, naive_npv, length(classified_ndn_granges))
    names(npv) <- c("duonovo_NPV", "naive_NPV", "n_other_parent")
    out <- npv
  } else if (metric == "false positive rate"){
    
    dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0" & dn_granges$GQ_parent >= validation_GQ_cutoff)]
    inherited_from_missing_parent <- grep("1", dn_granges$parentValidation_gt)
    all_inherited <- dn_granges[inherited_from_missing_parent]
    all_inherited <- all_inherited[all_inherited$duoNovo_classification != "failed_QC"]
    all_inherited <- all_inherited[which(all_inherited$parentValidation_depth >= 20 & 
                                           all_inherited$parentValidation_GQ >= validation_GQ_cutoff)]
    
    classified_dn <- which(all_inherited$duoNovo_classification == "de_novo" & 
                             all_inherited$GQ_proband >= validation_GQ_cutoff)
    dn <- length(classified_dn)
    
    classified_ndn <- which(all_inherited$duoNovo_classification == "on_other_parent_haplotype")
    ndn <- length(classified_ndn)
    
    uncertain <- length(which(dn_granges$duoNovo_classification == "uncertain"))
    uncertain2 <- length(which(all_inherited$duoNovo_classification == "de_novo" & 
                                 all_inherited$GQ_proband < validation_GQ_cutoff))
    
    false_pos <- c(dn, ndn, uncertain + uncertain2)
    names(false_pos) <- c("dn", "ndn", "uncertain")
    out <- false_pos
  } else if (metric == "de novo call rate"){
    if (allele_type == "ALT"){
      
      dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0" & dn_granges$GQ_parent >= validation_GQ_cutoff)]
    } else if (allele_type == "REF"){
      
      dn_granges <- dn_granges[which(dn_granges$phasing_parent == "1/1" & dn_granges$GQ_parent >= validation_GQ_cutoff)]
    }
    dn_granges <- dn_granges[dn_granges$duoNovo_classification != "failed_QC"]
    classified_dn <- which(dn_granges$duoNovo_classification == "de_novo" & 
                             dn_granges$GQ_proband >= validation_GQ_cutoff)
    
    dn_call_rate <- c(length(classified_dn), length(dn_granges), 
                          length(classified_dn)/length(dn_granges))
    out <- dn_call_rate
  } else if (metric == "n of total variants classified"){
    dn_granges$problematic_region <- unlist(dn_granges$problematic_region)
    dn_granges <- dn_granges[which(dn_granges$problematic_region == FALSE)]
    dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0" & dn_granges$GQ_parent >= validation_GQ_cutoff)]
    dn_granges <- dn_granges[which(dn_granges$duoNovo_classification %in% 
                                     c("de_novo", "on_other_parent_haplotype") & 
                                     dn_granges$GQ_proband >= validation_GQ_cutoff)]
    out <- dn_granges
  }
  out
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: process_duoNovo_output.R <run_directory>")
}

## set directory
current_dir <- args[1]
setwd(current_dir)

duonovo_granges_output_filepaths <- list.files(pattern = "PF\\.duonovo\\.annovar\\.addedParent\\.dnm2.rda$|PM\\.duonovo\\.annovar\\.addedParent\\.dnm2.rda$")
duoNovo_output_filepath_pm <- grep("^.*\\.PM\\.", duonovo_granges_output_filepaths, value = TRUE)
duoNovo_output_filepath_pf <- grep("^.*\\.PF\\.", duonovo_granges_output_filepaths, value = TRUE)

if (length(duoNovo_output_filepath_pf) != 1L ||
    length(duoNovo_output_filepath_pm) != 1L) {
  stop("Couldn’t unambiguously detect PF vs PM GRanges files")
}

### --- father-proband duos
###
ppv_pf_alt <- getDuoNovoPerformanceMetric(duoNovo_output_filepath_pf, duo_type = "PF", 
                                          include_gnomad = TRUE, metric = "PPV")
ppv_pf_alt_no_gnomad <- getDuoNovoPerformanceMetric(duoNovo_output_filepath_pf, duo_type = "PF", 
                                                    include_gnomad = FALSE, metric = "PPV")

npv_pf <- getDuoNovoPerformanceMetric(duoNovo_output_filepath_pf, duo_type = "PF", metric = "NPV")

false_pos_pf <- getDuoNovoPerformanceMetric(duoNovo_output_filepath_pf, duo_type = "PF", 
                                            metric = "false positive rate")

dn_call_rate_alt_pf <- getDuoNovoPerformanceMetric(duoNovo_output_filepath_pf, duo_type = "PF", 
                                               allele_type = "ALT", metric = "de novo call rate")
dn_call_rate_ref_pf <- getDuoNovoPerformanceMetric(duoNovo_output_filepath_pf, duo_type = "PF", 
                                                   allele_type = "REF", metric = "de novo call rate")
variants_classified_pf <- getDuoNovoPerformanceMetric(duoNovo_output_filepath_pf, duo_type = "PF", 
                                                            metric = "n of total variants classified")

### --- mother-proband duos
###
ppv_pm_alt <- getDuoNovoPerformanceMetric(duoNovo_output_filepath_pm, duo_type = "PM", 
                                          include_gnomad = TRUE, metric = "PPV")
ppv_pm_alt_no_gnomad <- getDuoNovoPerformanceMetric(duoNovo_output_filepath_pm, duo_type = "PM", 
                                                    include_gnomad = FALSE, metric = "PPV")

npv_pm <- getDuoNovoPerformanceMetric(duoNovo_output_filepath_pm, duo_type = "PM", metric = "NPV")

false_pos_pm <- getDuoNovoPerformanceMetric(duoNovo_output_filepath_pm, duo_type = "PM", 
                                            metric = "false positive rate")

dn_call_rate_alt_pm <- getDuoNovoPerformanceMetric(duoNovo_output_filepath_pm, duo_type = "PM", 
                                                   allele_type = "ALT", metric = "de novo call rate")
dn_call_rate_ref_pm <- getDuoNovoPerformanceMetric(duoNovo_output_filepath_pm, duo_type = "PM", 
                                                   allele_type = "REF", metric = "de novo call rate")

variants_classified_pm <- getDuoNovoPerformanceMetric(duoNovo_output_filepath_pm, duo_type = "PM", 
                                                            metric = "n of total variants classified")

### both duos
total_variants_classified <- length(unique(c(variants_classified_pf, variants_classified_pm)))


### --- save results
###
sample_id <- sub("\\..*$", "", duoNovo_output_filepath_pf)
save(ppv_pf_alt, ppv_pm_alt, ppv_pf_alt_no_gnomad, ppv_pm_alt_no_gnomad, 
     file = paste0(sample_id, "_ppv_pm_and_pf.rda"))
save(npv_pf, npv_pm, 
     file = paste0(sample_id, "_npv_pm_and_pf.rda"))
save(false_pos_pf, false_pos_pm, 
     file = paste0(sample_id, "_false_pos_pm_and_pf.rda"))
save(dn_call_rate_alt_pf, dn_call_rate_ref_pf, dn_call_rate_alt_pm, dn_call_rate_ref_pm,
     file = paste0(sample_id, "_denovo_call_rate_alt_vs_ref.rda"))
save(total_variants_classified, file = paste0(sample_id, "_n_total_variants_classified.rda"))
