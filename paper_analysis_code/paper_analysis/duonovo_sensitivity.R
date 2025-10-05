getSensitivity <- function(trio_denovo_filepath, duoNovo_granges_output_filepath_pm, duoNovo_granges_output_filepath_pf, 
                           duo_type = c("PM", "PF", "both"), validation_GQ_cutoff = 30, 
                           allele_depth_filter = TRUE){
  load(file = trio_denovo_filepath)
  trio_de_novo$problematic_region <- unlist(trio_de_novo$problematic_region)
  trio_de_novo <- trio_de_novo[which(trio_de_novo$problematic_region == "." & 
                                       trio_de_novo$proband_GQ >= validation_GQ_cutoff & 
                                       trio_de_novo$parent1_GQ >= validation_GQ_cutoff & 
                                       trio_de_novo$parent2_GQ >= validation_GQ_cutoff & 
                                       trio_de_novo$proband_dp >= 20 & 
                                       trio_de_novo$parent1_dp >= 20 & 
                                       trio_de_novo$parent2_dp >= 20 & 
                                       trio_de_novo$proband_gt %in% c("0|1", "1|0") & 
                                       trio_de_novo$parent1_gt == "0/0" & 
                                       trio_de_novo$parent2_gt == "0/0")]
  if (length(trio_de_novo) == 0) {
    stop("No variant allele de novos from the trio vcf pass GQ and depth filters")
  }
  
  if (allele_depth_filter == TRUE){
    allele_depth_mat_1 <- matrix(
      unlist(trio_de_novo$parent1_AD, use.names = FALSE),
      ncol    = 2,      # two columns: ref, alt
      byrow   = TRUE
    )
    colnames(allele_depth_mat_1) <- c("AD_ref", "AD_alt")

    allele_depth_mat_2 <- matrix(
      unlist(trio_de_novo$parent2_AD, use.names = FALSE),
      ncol    = 2,      # two columns: ref, alt
      byrow   = TRUE
    )
    colnames(allele_depth_mat_2) <- c("AD_ref", "AD_alt")

    de_novos_to_keep <- which(allele_depth_mat_1[, 'AD_alt'] == 0 & allele_depth_mat_2[, 'AD_alt'] == 0)
    if (length(de_novos_to_keep) == 0) {
      stop("No variant allele de novos from the trio vcf pass allele depth filter")
    }
    trio_de_novo <- trio_de_novo[de_novos_to_keep]
  }
  total <- length(trio_de_novo)
  
  if (duo_type == "PF"){
    
    load(file = duoNovo_granges_output_filepath_pf)
    dn_granges <- dn_granges_pf
    
    classified_dn <- which(dn_granges$duoNovo_classification == "de_novo")
    classified_ndn <- which(dn_granges$duoNovo_classification == "on_other_parent_haplotype")
    classified_uncertain <- which(dn_granges$duoNovo_classification == "uncertain")
    if (length(classified_dn) > 0){
      dn <- length(queryHits(findOverlaps(trio_de_novo, dn_granges[classified_dn])))
    } else {
      dn <- 0
    }
    if (length(classified_ndn) > 0){
      ndn <- length(queryHits(findOverlaps(trio_de_novo, dn_granges[classified_ndn])))
    } else {
      ndn <- 0
    }
    if (length(classified_uncertain) > 0){
      uncertain <- length(queryHits(findOverlaps(trio_de_novo, dn_granges[classified_uncertain])))
    } else {
      uncertain <- 0
    }
    sens <- c(dn, ndn, uncertain, total)
    names(sens) <- c("de_novo", "on_other_parent_haplotype", "uncertain", "total")
    out <- sens
  } else if (duo_type == "PM"){
    
    load(file = duoNovo_granges_output_filepath_pm)
    dn_granges <- dn_granges_pm
    
    classified_dn <- which(dn_granges$duoNovo_classification == "de_novo")
    classified_ndn <- which(dn_granges$duoNovo_classification == "on_other_parent_haplotype")
    classified_uncertain <- which(dn_granges$duoNovo_classification == "uncertain")
    if (length(classified_dn) > 0){
      dn <- length(queryHits(findOverlaps(trio_de_novo, dn_granges[classified_dn])))
    } else {
      dn <- 0
    }
    if (length(classified_ndn) > 0){
      ndn <- length(queryHits(findOverlaps(trio_de_novo, dn_granges[classified_ndn])))
    } else {
      ndn <- 0
    }
    if (length(classified_uncertain) > 0){
      uncertain <- length(queryHits(findOverlaps(trio_de_novo, dn_granges[classified_uncertain])))
    } else {
      uncertain <- 0
    }
    sens <- c(dn, ndn, uncertain, total)
    names(sens) <- c("de_novo", "on_other_parent_haplotype", "uncertain", "total")
    out <- sens
  } else if (duo_type == "both"){
    
    load(file = duoNovo_granges_output_filepath_pf)
    load(file = duoNovo_granges_output_filepath_pm)
    
    dn_f <- unique(queryHits(findOverlaps(trio_de_novo, 
                                          dn_granges_pf[dn_granges_pf$duoNovo_classification == "de_novo"])))
    dn_m <- unique(queryHits(findOverlaps(trio_de_novo, 
                                          dn_granges_pm[dn_granges_pm$duoNovo_classification == "de_novo"])))
    dn_both <- intersect(dn_f, dn_m)
    dn_either <- union(dn_f, dn_m)
    
    ndn_f <- unique(queryHits(findOverlaps(trio_de_novo, 
                                           dn_granges_pf[dn_granges_pf$duoNovo_classification == "on_other_parent_haplotype"])))
    ndn_m <- unique(queryHits(findOverlaps(trio_de_novo, 
                                           dn_granges_pm[dn_granges_pm$duoNovo_classification == "on_other_parent_haplotype"])))
    ndn_both <- intersect(ndn_f, ndn_m)
    sens <- c(length(dn_either), length(dn_both), length(ndn_both), total)
    names(sens) <- c("de_novo_either", "de_novo_both", "non_denovo_both", "total")
    out <- sens
  }
  out
}

library(VariantAnnotation)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: duonovo_sensitivity.R <run_directory>")
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

trio_denovo_filepath <- list.files(pattern = "PFM\\.annovar\\.dnm2\\.rda$")
if (length(trio_denovo_filepath) != 1L ) {
  stop("Couldn’t unambiguously detect trio de novo GRanges file")
}

### --- father-proband duo
message("getting sensitivity from PF duo...")
sens_pf <- getSensitivity(trio_denovo_filepath, duoNovo_granges_output_filepath_pf = duoNovo_output_filepath_pf, 
                          duo_type = "PF")
message("getting sensitivity from PF duo without allele depth filter...")
sens_pf_no_ad_filter <- getSensitivity(trio_denovo_filepath, duoNovo_granges_output_filepath_pf = duoNovo_output_filepath_pf, 
                          duo_type = "PF", allele_depth_filter = FALSE)
### --- mother-proband duo
message("getting sensitivity from PM duo...")
sens_pm <- getSensitivity(trio_denovo_filepath, duoNovo_granges_output_filepath_pm = duoNovo_output_filepath_pm, 
                          duo_type = "PM")
message("getting sensitivity from PM duo without allele depth filter...")
sens_pm_no_ad_filter <- getSensitivity(trio_denovo_filepath, duoNovo_granges_output_filepath_pm = duoNovo_output_filepath_pm, 
                          duo_type = "PM", allele_depth_filter = FALSE)

### --- both duos
message("getting sensitivity when combining both duos...")
sens_both <- getSensitivity(trio_denovo_filepath, duoNovo_granges_output_filepath_pm = duoNovo_output_filepath_pm, 
                            duoNovo_granges_output_filepath_pf = duoNovo_output_filepath_pf, duo_type = "both")
message("getting sensitivity when combining both duos without allele depth filter...")
sens_both_no_ad_filter <- getSensitivity(trio_denovo_filepath, duoNovo_granges_output_filepath_pm = duoNovo_output_filepath_pm, 
                            duoNovo_granges_output_filepath_pf = duoNovo_output_filepath_pf, duo_type = "both", 
                            allele_depth_filter = FALSE)


### --- save results
###
sample_id <- sub("\\..*$", "", duoNovo_output_filepath_pf)
save(sens_pf, sens_pm, sens_both, sens_pf_no_ad_filter, sens_pm_no_ad_filter, sens_both_no_ad_filter,
     file = paste0(sample_id, "_sensitivity.rda"))



