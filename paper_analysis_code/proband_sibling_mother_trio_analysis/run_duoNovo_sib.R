source("duoNovoSib.R")
source("getHaplotypesAndBlockBoundariesForTrio.R")
source("classifyVariantsTrio.R")

library(VariantAnnotation)
setwd('duoNovoOutputs_ANNOVAR')
directories <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE)

# retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)
dir_index <- as.integer(args[1]) #will be used later
setwd(directories[dir_index])

duoNovo_output_filepaths <- list.files( #output from process_duoNovo_output.R
  path       = "duoNovoOutputsWithOtherParentGeno/PM", #replace with correct directory if necessary
  pattern    = "duoNovo\\.annovar\\.rda$",
  full.names = TRUE
)
duoNovo_output_filepath <- duoNovo_output_filepaths[dir_index]

actual_trio_denovo_filepaths <- list.files( #output from get_denovo_from_trio.R
  path        = "trio_de_novo", #replace with correct directory if necessary
  pattern     = "^TRIO.*annovar\\.rda$",
  full.names  = TRUE
)
actual_trio_denovo_filepath <- actual_trio_denovo_filepaths[dir_index]

### use duoNovo output from the mother-proband duo to separate candidate variants based on classification
load(file = duoNovo_output_filepath) 
duoNovo_ranges <- dn_granges_pm

classified_denovo_indices <- which(duoNovo_ranges$duoNovo_classification == "de_novo" & 
                         (duoNovo_ranges$n_de_novo_left_orientation_same_PS == 1 | 
                            duoNovo_ranges$n_de_novo_right_orientation_same_PS == 1))
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
surrogate_trio_vcf_filepath <- list.files(pattern = "TRIO_PSM")
sibling_trio_output <- duoNovoSib(surrogate_trio_vcf_filepath, 
                                  candidate_variants_duoNovo_output = classified_other_parent_or_uncertain)

# Obtain de novo variants from the output
denovo_indices_sibling_trio <- which(sibling_trio_output$duoNovo_classification == "de_novo" & 
                                   (sibling_trio_output$n_de_novo_left_orientation_same_PS == 1 | 
                                      sibling_trio_output$n_de_novo_right_orientation_same_PS == 1))
if (length(denovo_indices_sibling_trio) > 0){
  denovo_sibling <- sibling_trio_output[denovo_indices_sibling_trio]
} else {
  denovo_sibling <- GRanges()
}

### Now import actual trio, and calculate what fraction of the ground truth de novo variants are detected:
### a) from mother-proband duo only; b) using the proband-sibling-mother trio
load(file = actual_trio_denovo_filepath) #name of imported file is "trio_de_novo"
trio_de_novo <- trio_de_novo[which(trio_de_novo$proband_gq >= 40 & 
                                     trio_de_novo$parent1_gq >= 40 & 
                                     trio_de_novo$parent2_gq >= 40)]
trio_de_novo$giab_problematic <- unlist(trio_de_novo$giab_problematic)
trio_de_novo <- trio_de_novo[which(trio_de_novo$giab_problematic == ".")]

sensitivity_duo_only <- length(unique(queryHits(findOverlaps(trio_de_novo, 
                                                             classified_denovo))))/ length(trio_de_novo)
combined_denovos <- c(classified_denovo, denovo_sibling)
sensitivity_duo_plus_sibling <- length(unique(queryHits(findOverlaps(trio_de_novo, 
                                                                     combined_denovos))))/ length(trio_de_novo)

### get false positive rate (calculate fraction of de novos detected from the proband-sibling-mother trio
### that are not de novo when we look at the actual father's genotype)
if (length(denovo_indices_sibling_trio) > 0){
  denovo_sibling_false <- unique(queryHits(findOverlaps(denovo_sibling, 
                                                        duoNovo_ranges[grep("1", duoNovo_ranges$parentValidation_gt)])))
  false_positive_rate <- length(denovo_sibling_false)/length(denovo_sibling)
  parental_origin <- table(denovo_sibling$parent_origin)
} else {
  false_positive_rate <- NA
}

### save results
save(sensitivity_duo_only, sensitivity_duo_plus_sibling, false_positive_rate, parental_origin, 
     file = "sensitivity_sibling.rda")


