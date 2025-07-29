source("duoNovoSib.R")
source("getHaplotypesAndBlockBoundariesForTrio.R")
source("classifyVariantsTrio.R")

library(VariantAnnotation)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: run_duoNovo_sib.R <run_directory>")
}

psm_trio_dir <- args[1]
duo_dir <- args[2]

### use duoNovo output from the mother-proband duo to separate candidate variants based on classification
setwd(duo_dir)
duonovo_granges_output_filepath <- list.files(pattern = "PM\\.duonovo\\.annovar\\.addedParent\\.dnm2.rda$")
load(file = duonovo_granges_output_filepath)
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
setwd(trio_dir)
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


