source("duoNovoSib.R")
source("getHaplotypesAndBlockBoundariesForTrio.R")
source("classifyVariantsTrio.R")

library(VariantAnnotation)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: run_duoNovo_sib.R <run_directory>")
}
## set directory
current_dir <- args[1]
setwd(current_dir)

### use duoNovo output from the mother-proband duo to separate candidate variants based on classification
duoNovo_output_filepath_pm <- list.files(pattern = "PM\\.duonovo\\.annovar\\.addedParent\\.dnm2.rda$")
load(file = duonovo_granges_output_filepath_pm)
duoNovo_ranges <- dn_granges_pm

classified_denovo_indices <- which(duoNovo_ranges$duoNovo_classification == "de_novo")
if (length(classified_denovo_indices) > 0){
  classified_denovo <- duoNovo_ranges[classified_denovo_indices]
} else {
  classified_denovo <- GRanges()
}
classified_other_parent_indices <- which(duoNovo_ranges$duoNovo_classification == "on_other_parent_haplotype")
classified_uncertain_indices <- which(duoNovo_ranges$duoNovo_classification == "uncertain")
classified_other_parent_or_uncertain <- duoNovo_ranges[c(classified_uncertain_indices, 
                                                         classified_other_parent_indices)]


### Now run duoNovo (trio version) on the joint trio vcf (proband-sibling-mother; samples have to be ordered this way in the vcf) 
### to classify the candidate variants labeled as "on other parent haplotype" or uncertain
message("getting de novo variants from PSM trio...")
surrogate_trio_vcf_filepath <- list.files(pattern = "PSM\\.vcf.gz$")
sibling_trio_output <- duoNovoSib(surrogate_trio_vcf_filepath, 
                                  candidate_variants_duoNovo_output = classified_other_parent_or_uncertain)

# Obtain de novo variants from the output
denovo_indices_sibling_trio <- which(sibling_trio_output$duoNovo_classification == "de_novo")
if (length(denovo_indices_sibling_trio) > 0){
  denovo_sibling <- sibling_trio_output[denovo_indices_sibling_trio]
  
  ### can optionally exclude de novos that are present in gnomAD
  # duoNovo_ranges$gnomad41_genome_AF <- as.numeric(unlist(duoNovo_ranges$gnomad41_genome_AF))
  # overlapping_indices <- unique(queryHits(findOverlaps(denovo_sibling, 
  #                                         duoNovo_ranges[which(duoNovo_ranges$gnomad41_genome_AF > 0)])))
  # if (length(overlapping_indices) > 0){
  # denovo_sibling <- denovo_sibling[-overlapping_indices]
  #}
  ###
} else {
  denovo_sibling <- GRanges()
}

### now calculate sensitivity based on de novo variants identified from PM duo as well as PSM trio
###
getSensitivity <- function(trio_denovo_filepath, duoNovo_granges_output_filepath_pm, denovo_psm_trio, 
                           validation_GQ_cutoff = 30, 
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
  
  load(file = duoNovo_granges_output_filepath_pm)
  dn_granges <- dn_granges_pm
  
  classified_dn <- which(dn_granges$duoNovo_classification == "de_novo")
  if (length(classified_dn) > 0){
    dn_duo <- length(queryHits(findOverlaps(trio_de_novo, dn_granges[classified_dn])))
  } else {
    dn_duo <- 0
  }
  if (length(denovo_psm_trio) > 0){
    dn_psm_trio <- length(queryHits(findOverlaps(trio_de_novo, denovo_psm_trio)))
  }  else {
    dn_psm_trio <- 0
  }
  
  sens <- c(dn_duo + dn_psm_trio, total)
  names(sens) <- c("de_novo", "total")
  out <- sens
  out
}

trio_denovo_filepath <- list.files(pattern = "PFM\\.annovar\\.dnm2\\.rda$")
if (length(trio_denovo_filepath) != 1L ) {
  stop("Couldn’t unambiguously detect trio de novo GRanges file")
}

message("calculating updated sensitivity...")
sens <- getSensitivity(trio_denovo_filepath, duoNovo_granges_output_filepath_pm = duoNovo_output_filepath_pm, 
                       denovo_psm_trio = denovo_sibling)

sample_id <- sub("\\..*$", "", duoNovo_output_filepath_pm)
save(denovo_sibling, sens,
     file = paste0(sample_id, "_psm_de_novo_variants.rda"))





