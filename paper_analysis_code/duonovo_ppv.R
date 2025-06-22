getPPV <- function(duoNovo_granges_output_filepath, duo_type = c("PM", "PF"), 
                   include_gnomad = c(TRUE, FALSE)){
  load(file = duoNovo_granges_output_filepath)
  if (duo_type == "PF"){
    dn_granges <- dn_granges_pf
  } else if (duo_type == "PM"){
    dn_granges <- dn_granges_pm
  }
  if (include_gnomad = FALSE){
    dn_granges$gnomad41_genome_AF <- as.numeric(unlist(dn_granges$gnomad41_genome_AF))
    dn_granges <- dn_granges[-which(dn_granges$gnomad41_genome_AF > 0)]
  }
  dn_granges <- dn_granges[which(dn_granges$duoNovo_classification == "de_novo" & 
                                   dn_granges$GQ_proband >= 40 & dn_granges$GQ_parent >= 40)]
  dn_granges$giab_problematic <- unlist(dn_granges$giab_problematic)
  dn_granges <- dn_granges[which(dn_granges$giab_problematic == ".")]
  dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0")]
  dn_granges <- dn_granges[which(dn_granges$n_de_novo_left_orientation_same_PS == 1 | 
                                   dn_granges$n_de_novo_right_orientation_same_PS == 1)]
  
  total <- length(dn_granges)
  assessed_indices <- which(dn_granges$parentValidation_depth >= 20 & 
                              dn_granges$parentValidation_GQ >= 40 & 
                              !grepl("\\.", dn_granges$parentValidation_gt))
  assessed <- length(assessed_indices)
  
  dn_granges <- dn_granges[assessed_indices]
  dn <- length(dn_granges) - length(grep("1", dn_granges$parentValidation_gt))
  ppv <- c(dn, assessed, total, length(dn_granges_pf))
  names(ppv) <- c("dn", "assessed", "total", "all_candidates")
  ppv
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: duonovo_ppv.R <dir_list.txt> <index>")
}

dir_file   <- args[1]
dir_index  <- as.integer(args[2])

## read directory list
dirs <- scan(dir_file, what = character(), sep = "\n", quiet = TRUE)
if (dir_index < 1L || dir_index > length(dirs))
  stop("Index ", dir_index, " is out of range 1–", length(dirs))

current_dir <- dirs[dir_index]
setwd(current_dir)

duoNovo_output_filepath_pm <- list.files(
  pattern = "^duo_pm.*\\.rda$",   # prefix + anything + .rda
  ignore.case = TRUE,             # make it case-insensitive if you wish
  full.names = FALSE)
duoNovo_output_filepath_pf <- list.files(
  pattern = "^duo_pf.*\\.rda$",   # prefix + anything + .rda
  ignore.case = TRUE,             # make it case-insensitive if you wish
  full.names = FALSE)


### --- father-proband duos
###
ppv_pf_alt <- getPPV(duoNovo_output_filepath_pf, duo_type = "PF", include_gnomad = TRUE)
ppv_pf_alt_no_gnomad <- getPPV(duoNovo_output_filepath_pf, duo_type = "PF", include_gnomad = FALSE)

### --- mother-proband duos
###
ppv_pm_alt <- getPPV(duoNovo_output_filepath_pm, duo_type = "PM", include_gnomad = TRUE)
ppv_pm_alt_no_gnomad <- getPPV(duoNovo_output_filepath_pm, duo_type = "PM", include_gnomad = FALSE)

dir_tag <- basename(current_dir)
save(ppv_pf_alt, ppv_pm_alt, ppv_pf_alt_no_gnomad, ppv_pm_alt_no_gnomad, 
     file = paste0(dir_tag, "_ppv_pm_and_pf.rda"))

