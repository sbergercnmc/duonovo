###
library(VariantAnnotation)
library(matrixStats)
###
###Also needs the main function (defined in duoNovo_main_function.R) and the auxiliary functions
###

args <- commandArgs(trailingOnly = TRUE)
vcf_index <- as.integer(args[1])

setwd('/home/lboukas/PMGRC_longread/duo/vcfs/short_read_vcfs')
all_vcfs_short_read <- list.files(pattern = '*short_read_vcfs.rda')
load(file = all_vcfs_short_read[vcf_index])
v_sr <- vcfs[[1]] #replace 1 with 2 for mother-proband duo
v_granges_duo_SR <- getQualFilteredVariantGrangesSR(v_sr, depth_threshold = 20, GQ_threshold = 30)

setwd('/home/lboukas/PMGRC_longread/duo/vcfs')
all_vcfs <- list.files(pattern = '*vcfs.rda*')
all_vcfs <- all_vcfs[-grep("de_novo", all_vcfs)]


load(file = all_vcfs[vcf_index])
v <- vcfs[[1]] #replace 1 with 2 for mother-proband duo
v_granges_duo <- getQualFilteredVariantGranges(v, depth_threshold = 20, GQ_threshold = 30)
hap_granges <- getHaplotypes(v_granges_duo)
combined_PS <- getPhasingSetCoordinates(hap_granges)

de_novo_candidates1 <- duoNovo(hap_granges, combined_PS, "0|1", 10000, 2000, 
                                           distance_cutoff = 40, v_granges_duo_SR)
de_novo_candidates2 <- duoNovo(hap_granges, combined_PS, "1|0", 10000, 2000, 
                                           distance_cutoff = 40, v_granges_duo_SR)

file_name_id <- sub("_.*", "", all_vcfs[vcf_index])
save(list = c("de_novo_candidates1", "de_novo_candidates2"),
     file = paste0("de_novo_candidates/alt_father_output_candidate_de_novo_", file_name_id, ".rda"))


