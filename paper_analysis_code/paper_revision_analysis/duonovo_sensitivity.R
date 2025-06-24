getSensitivity <- function(trio_denovo_filepath, duoNovo_granges_output_filepath_pm, duoNovo_granges_output_filepath_pf, 
                           duo_type = c("PM", "PF", "both"), validation_GQ_cutoff = 30){
  load(file = trio_denovo_filepath)
  trio_de_novo$problematic_region <- as.logical(trio_de_novo$problematic_region)
  trio_de_novo <- trio_de_novo[which(trio_de_novo$problematic_region == FALSE & 
                                       trio_de_novo$proband_GQ >= validation_GQ_cutoff & 
                                       trio_de_novo$parent1_GQ >= validation_GQ_cutoff & 
                                       trio_de_novo$parent2_GQ >= validation_GQ_cutoff)]
  
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
sens_pf <- getSensitivity(trio_denovo_filepath, duoNovo_granges_output_filepath_pf = duoNovo_output_filepath_pf, 
                          duo_type = "PF")
### --- mother-proband duo
sens_pm <- getSensitivity(trio_denovo_filepath, duoNovo_granges_output_filepath_pm = duoNovo_output_filepath_pm, 
                          duo_type = "PM")
### --- both duos
sens_both <- getSensitivity(trio_denovo_filepath, duoNovo_output_filepath_pm, 
                            duoNovo_output_filepath_pf, duo_type = "both")


### --- save results
###
sample_id <- sub("\\..*$", "", duoNovo_output_filepath_pf)
save(sens_pf, sens_pm, sens_both,
     file = paste0(sample_id, "_sensitivity.rda"))



