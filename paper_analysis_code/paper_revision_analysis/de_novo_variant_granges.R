getDeNovoVariantGRanges <- function(duoNovo_granges_output_filepath, duo_type = c("PM", "PF"), 
                                    filter_problematic_regions = TRUE, 
                                    exclude_clustered_denovos = TRUE,
                                    genomic_annotation = NULL, 
                                    validation_GQ_cutoff = 30){
  load(file = duoNovo_granges_output_filepath)
  
  if (duo_type == "PF"){
    dn_granges <- dn_granges_pf
  } else if (duo_type == "PM"){
    dn_granges <- dn_granges_pm
  }
  
  dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0" & dn_granges$GQ_parent >= validation_GQ_cutoff)]
  if (filter_problematic_regions == FALSE){
    indices <- grep("de_novo", dn_granges$QC_fail_step)
    if (length(indices) > 0){
      dn_granges$duoNovo_classification[indices] <- "de_novo"
    }
  }
  if (exclude_clustered_denovos == FALSE){
    indices <- which(dn_granges$duoNovo_classification == "on_multi_denovo_haplotype")
    if (length(indices) > 0){
      dn_granges$duoNovo_classification[indices] <- "de_novo"
    }
  }
  classified_dn <- which(dn_granges$duoNovo_classification == "de_novo" & 
                           dn_granges$GQ_proband >= validation_GQ_cutoff)
  if (length(classified_dn) > 0){
    dn_granges <- dn_granges[classified_dn]
  } else {
    dn_granges <- GRanges()
  }
  if (!is.null(genomic_annotation)){
    dn_granges$Func.refGeneWithVer <- unlist(dn_granges$Func.refGeneWithVer)
    variants_in_annotation <- which(dn_granges$Func.refGeneWithVer %in% genomic_annotation)
    if (length(variants_in_annotation) > 0){
      dn_granges <- dn_granges[variants_in_annotation]
    } else {
      dn_granges <- GRanges()
    }
  }
  dn_granges
}

library(VariantAnnotation)
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
message("getting variants classified de novo from PF duo...")
de_novos_pf <- getDeNovoVariantGRanges(duoNovo_output_filepath_pf, duo_type = "PF")
message("getting variants classified de novo from PF duo without excluding problematic regions...")
de_novos_pf_with_giab_problematic <- getDeNovoVariantGRanges(duoNovo_output_filepath_pf, duo_type = "PF", 
                                                             filter_problematic_regions = FALSE)
message("getting variants classified de novo from PF duo without excluding those on multi-denovo haplotype...")
de_novos_pf_with_clustered <- getDeNovoVariantGRanges(duoNovo_output_filepath_pf, duo_type = "PF", 
                                                             exclude_clustered_denovos = FALSE)
message("getting variants classified de novo from PF duo in exonic/intronic regions...")
de_novos_pf_exonic_intronic_only <- getDeNovoVariantGRanges(duoNovo_output_filepath_pf, duo_type = "PF", 
                                                            genomic_annotation = c("exonic", "intronic", "ncRNA_exonic", "ncRNA_intronic", 
                                                                                   "UTR3", "UTR5"))

### --- mother-proband duos
message("getting variants classified de novo from PM duo...")
de_novos_pm <- getDeNovoVariantGRanges(duoNovo_output_filepath_pm, duo_type = "PM")
message("getting variants classified de novo from PM duo without excluding problematic regions...")
de_novos_pm_with_giab_problematic <- getDeNovoVariantGRanges(duoNovo_output_filepath_pm, duo_type = "PM", 
                                                             filter_problematic_regions = FALSE)
message("getting variants classified de novo from PM duo without excluding those on multi-denovo haplotype...")
de_novos_pm_with_clustered <- getDeNovoVariantGRanges(duoNovo_output_filepath_pm, duo_type = "PM", 
                                                      exclude_clustered_denovos = FALSE)
message("getting variants classified de novo from PM duo in exonic/intronic regions...")
de_novos_pm_exonic_intronic_only <- getDeNovoVariantGRanges(duoNovo_output_filepath_pm, duo_type = "PM", 
                                                            genomic_annotation = c("exonic", "intronic", "ncRNA_exonic", "ncRNA_intronic", 
                                                                                   "UTR3", "UTR5"))


### save
sample_id <- sub("\\..*$", "", duoNovo_output_filepath_pf)
save(de_novos_pf, de_novos_pf_with_giab_problematic, de_novos_pf_with_clustered, de_novos_pf_exonic_intronic_only,
     de_novos_pm, de_novos_pm_with_giab_problematic, de_novos_pm_with_clustered, de_novos_pm_exonic_intronic_only,
     file = paste0(sample_id, "_de_novo_variant_granges.rda"))

