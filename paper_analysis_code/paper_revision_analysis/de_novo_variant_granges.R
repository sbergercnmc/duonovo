getDeNovoVariantGRanges <- function(duoNovo_granges_output_filepath, duo_type = c("PM", "PF"), 
                                    filter_problematic_regions = TRUE, 
                                    exclude_clustered_denovos = TRUE,
                                    genomic_context = NULL){
  load(file = duoNovo_granges_output_filepath)
  if (duo_type == "PF"){
    dn_granges <- dn_granges_pf
  } else if (duo_type == "PM"){
    dn_granges <- dn_granges_pm
  }
  dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0" & dn_granges$GQ_parent >= 40)]
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
                           dn_granges$GQ_proband >= 40)
  if (length(classified_dn) > 0){
    dn_granges <- dn_granges[classified_dn]
  } else {
    dn_granges <- GRanges()
  }
  if (!is.null(genomic_context)){
    variants_in_genomic_context <- which(dn_granges$Func.refGeneWithVer %in% genomic_context)
    if (length(variants_in_genomic_context) > 0){
      dn_granges <- dn_granges[variants_in_genomic_context]
    } else {
      dn_granges <- GRanges()
    }
  }
  dn_granges
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: duonovo_performance_metrics.R <trio_directories.txt> <index>")
}

dir_file   <- args[1]
dir_index  <- as.integer(args[2])

## read directory list
dirs <- trimws(readLines(dir_file))
if (dir_index < 1L || dir_index > length(dirs))
  stop("Index ", dir_index, " is out of range 1–", length(dirs))

current_dir <- dirs[dir_index]
setwd(current_dir)

duonovo_granges_output_filepaths <- list.files(pattern = "PF\\.duonovo\\.addedParent\\.dnm2.rda$|PM\\.duonovo\\.addedParent\\.dnm2.rda$")
duoNovo_output_filepath_pm <- grep("^.*\\.PM\\.", duonovo_granges_output_filepaths, value = TRUE)
duoNovo_output_filepath_pf <- grep("^.*\\.PF\\.", duonovo_granges_output_filepaths, value = TRUE)

if (length(duoNovo_output_filepath_pf) != 1L ||
    length(duoNovo_output_filepath_pm) != 1L) {
  stop("Couldn’t unambiguously detect PF vs PM GRanges files")
}

### --- father-proband duos
de_novos_pf <- getDeNovoVariantGRanges(duoNovo_output_filepath_pf, duo_type = "PF")

### --- mother-proband duos
de_novos_pm <- getDeNovoVariantGRanges(duoNovo_output_filepath_pm, duo_type = "PF")

### save
sample_id <- sub("\\..*$", "", duoNovo_output_filepath_pf)
save(de_novos_pf, de_novos_pm, 
     file = paste0(sample_id, "_de_novo_variant_granges.rda"))

