source("duoNovoSib.R")
source("getHaplotypesAndBlockBoundariesForTrio.R")
source("classifyVariantsTrio.R")

library(VariantAnnotation)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: falsepos_duoNovo_sib.R <run_directory>")
}
## set directory
current_dir <- args[1]
setwd(current_dir)


### run duoNovo (trio version) on the joint trio vcf (proband-sibling-mother; samples have to be ordered this way in the vcf) 
surrogate_trio_vcf_filepath <- list.files(pattern = "PSM\\.vcf.gz$")
sibling_trio_output <- duoNovoSib(surrogate_trio_vcf_filepath)

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

### use full trio genotypes to assess false pos rate
### use duoNovo output from the mother-proband duo to separate candidate variants based on classification
duoNovo_output_filepath_pm <- list.files(pattern = "PM\\.duonovo\\.annovar\\.addedParent\\.dnm2.rda$")
load(file = duonovo_granges_output_filepath_pm)
gn_ranges <- dn_granges_pm

dn_granges <- dn_granges[which(dn_granges$parentValidation_depth >= 20 & 
                                 dn_granges$parentValidation_GQ >= 30)]
dn_granges <- dn_granges[!grepl("\\.", dn_granges$parentValidation_gt)]
overlaps <- findOverlaps(dn_granges, denovo_sibling)
if (length(queryHits(overlaps) > 0)){
  variants_to_assess <- dn_granges[unique(queryHits(overlaps))] 
  
  assessed <- length(variants_to_assess)
  false_dn <- length(grep("1", variants_to_assess$parentValidation_gt))
  output_vec <- c(false_dn, assessed, rate)
  names(output_vec) <- c("false_dn", "assessed", "false_pos_rate")
} else {
  output_vec <- NA
}

sample_id <- sub("\\..*$", "", duoNovo_output_filepath_pm)
save(output_vec, 
     file = paste0(sample_id, "_psm_falsepos_rate.rda"))







