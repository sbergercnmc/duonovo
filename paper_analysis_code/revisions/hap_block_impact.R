calculateHapBlockImpact <- function(duoNovo_granges_output_filepath, duo_type = c("PM", "PF")){
  load(file = duoNovo_granges_output_filepath)
  
  if (duo_type == "PF"){
    dn_granges <- dn_granges_pf
  } else if (duo_type == "PM"){
    dn_granges <- dn_granges_pm
  }
  
  all <- length(dn_granges)
  all_fail <- length(which(dn_granges$duoNovo_classification == "failed_QC"))
  
  
  hap_block_overlap_fail <- length(which(dn_granges$QC_fail_step == "no_haplotype_block_overlap"))
  small_hap_block_fail <- length(which(dn_granges$QC_fail_step == "in_small_haplotype_block"))
  boundary_fail <- length(which(dn_granges$QC_fail_step == "in_haplotype_block_boundary"))
  depth_fail <- length(which(dn_granges$QC_fail_step == "low_depth"))
  GQ_fail <- length(which(dn_granges$QC_fail_step == "low_GQ"))
  depth_and_GQ_fail <- length(which(dn_granges$QC_fail_step == "low_depth_and_GQ"))
  problematic_region_fail1 <- length(which(dn_granges$QC_fail_step %in% c("classified_failed_QC_in_problematic_region")))
  problematic_region_fail2 <- length(which(dn_granges$QC_fail_step %in% 
                                             c("classified_de_novo_in_problematic_region", 
                                               "classified_uncertain_in_problematic_region", 
                                               "classified_on_other_parent_haplotype_in_problematic_region")))
  
  qc_fail_lengths <- c(hap_block_overlap_fail, small_hap_block_fail, boundary_fail, 
                       depth_fail, GQ_fail, depth_and_GQ_fail, problematic_region_fail1, problematic_region_fail2,
                       all_fail, all)
  names(qc_fail_lengths) <- c("hap_block_overlap_qc_fail", "small_hap_block_fail", "boundary_fail",
                              "depth_fail", "GQ_fail", "depth_and_GQ_fail", 
                              "problematic_region_fail_QC", "problematic_region_classified",
                              "all_qc_fail", "total")
  out <- qc_fail_lengths
  out
}

library(VariantAnnotation)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: hap_block_impact.R <run_directory>")
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

### --- father-proband duo
hap_block_impact_pf <- calculateHapBlockImpact(duoNovo_output_filepath_pf, duo_type = "PF")
### --- mother-proband duo
hap_block_impact_pm <- calculateHapBlockImpact(duoNovo_output_filepath_pm, duo_type = "PM")

sample_id <- sub("\\..*$", "", duoNovo_output_filepath_pf)
save(hap_block_impact_pf, hap_block_impact_pm, file = paste0(sample_id, "_hap_block_impact.rda"))




