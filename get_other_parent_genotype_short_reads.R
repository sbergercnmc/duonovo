getOtherParentGenotype <- function(vcf_path, genotype = c("paternal", "maternal")){
  load(file = vcf_path)
  if (genotype == "maternal"){
    vcf <- vcfs[[2]]
  } else if (genotype == "paternal"){
    vcf <- vcfs[[1]]
  }
  
  gt <- geno(vcf)$GT
  depth <- geno(vcf)$DP
  gq <- geno(vcf)$GQ
  
  
  proband_column <- grep("-0$", colnames(gt))
  good_quality_gt_indices <- which(depth[, -proband_column] >= 20 & gq[, -proband_column] >= 30)
  gt[good_quality_gt_indices, -proband_column]
}

library(VariantAnnotation)
library(matrixStats)
args <- commandArgs(trailingOnly = TRUE)
vcf_index <- as.integer(args[1])

setwd('/home/lboukas/PMGRC_longread/duo/vcfs/short_read_vcfs')

all_vcfs <- list.files(pattern = '*short_read_vcfs.rda')

other_parent_genotype <- getOtherParentGenotype(all_vcfs[vcf_index], genotype = "maternal") #run once for maternal, once for paternal

file_name_id <- sub("_.*", "", all_vcfs[vcf_index])
save(other_parent_genotype, file = paste0("maternal_genotype_short_reads_", file_name_id, ".rda"))

