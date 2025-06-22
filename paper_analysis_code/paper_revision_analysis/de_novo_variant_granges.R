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
  if (filter_problematic_regions == TRUE){
    dn_granges$problematic_region <- as.logical(dn_granges$problematic_region)
    dn_granges <- dn_granges[which(dn_granges$problematic_region == FALSE)]
  }
  if (exclude_clustered_denovos == TRUE){
    clustered <- which(dn_granges$n_de_novo_left_orientation_same_PS > 1 | 
                         dn_granges$n_de_novo_right_orientation_same_PS > 1)
    if(length(clustered) > 0){
      dn_granges <- dn_granges[-clustered]
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
    variants_in_annotation <- which(dn_granges$Func.refGeneWithVer %in% genomic_annotation)
    if (length(variants_in_annotation) > 0){
      dn_granges <- dn_granges[variants_in_annotation]
    } else {
      dn_granges <- GRanges()
    }
  }
  dn_granges
}

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
de_novos_pf <- getDeNovoVariantGRanges(duoNovo_output_filepath_pf, duo_type = "PF")
de_novos_pf_with_giab_problematic <- getDeNovoVariantGRanges(duoNovo_output_filepath_pf, duo_type = "PF", 
                                                             filter_problematic_regions = FALSE)
de_novos_pf_with_clustered <- getDeNovoVariantGRanges(duoNovo_output_filepath_pf, duo_type = "PF", 
                                                             exclude_clustered_denovos = FALSE)
de_novos_pf_exonic_intronic_only <- getDeNovoVariantGRanges(duoNovo_output_filepath_pf, duo_type = "PF", 
                                                            genomic_annotation = c("exonic", "intronic", "ncRNA_exonic", "ncRNA_intronic", 
                                                                                   "UTR3", "UTR5"))

### --- mother-proband duos
de_novos_pm <- getDeNovoVariantGRanges(duoNovo_output_filepath_pm, duo_type = "PM")
de_novos_pm_with_giab_problematic <- getDeNovoVariantGRanges(duoNovo_output_filepath_pm, duo_type = "PM", 
                                                             filter_problematic_regions = FALSE)
de_novos_pm_with_clustered <- getDeNovoVariantGRanges(duoNovo_output_filepath_pm, duo_type = "PM", 
                                                      exclude_clustered_denovos = FALSE)
de_novos_pm_exonic_intronic_only <- getDeNovoVariantGRanges(duoNovo_output_filepath_pm, duo_type = "PM", 
                                                            genomic_annotation = c("exonic", "intronic", "ncRNA_exonic", "ncRNA_intronic", 
                                                                                   "UTR3", "UTR5"))


### save
sample_id <- sub("\\..*$", "", duoNovo_output_filepath_pf)
save(de_novos_pf, de_novos_pf_with_giab_problematic, de_novos_pf_with_clustered, de_novos_pf_exonic_intronic_only,
     de_novos_pm, de_novos_pm_with_giab_problematic, de_novos_pm_with_clustered, de_novos_pm_exonic_intronic_only,
     file = paste0(sample_id, "_de_novo_variant_granges.rda"))

