getFalsePositives <- function(duoNovo_granges_output_filepath, duo_type = c("PM", "PF")){
  load(file = duoNovo_granges_output_filepath)
  if (duo_type == "PF"){
    dn_granges <- dn_granges_pf
  } else if (duo_type == "PM"){
    dn_granges <- dn_granges_pm
  }
  dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0" & dn_granges$GQ_parent >= 40)]
  inherited_from_missing_parent <- grep("1", dn_granges$parentValidation_gt)
  all_inherited <- dn_granges[inherited_from_missing_parent]
  all_inherited <- all_inherited[-which(all_inherited$duoNovo_classification == "failed_QC")]
  all_inherited$giab_problematic <- unlist(all_inherited$giab_problematic)
  all_inherited <- all_inherited[which(all_inherited$giab_problematic == ".")]
  
  all_inherited <- all_inherited[which(all_inherited$parentValidation_GQ >= 40)]
  
  classified_dn <- which(all_inherited$duoNovo_classification == "de_novo" & 
                           all_inherited$GQ_proband >= 40 &
                           (all_inherited$n_de_novo_left_orientation_same_PS == 1 | 
                              all_inherited$n_de_novo_right_orientation_same_PS == 1))
  dn <- length(classified_dn)
  
  classified_ndn <- which(all_inherited$duoNovo_classification == "on_other_parent_haplotype")
  clustered_1 <- which(all_inherited$duoNovo_classification == "de_novo" & 
                         all_inherited$n_de_novo_left_orientation_same_PS > 1)
  clustered_2 <- which(all_inherited$duoNovo_classification == "de_novo" & 
                         all_inherited$n_de_novo_right_orientation_same_PS > 1)
  ndn <- length(classified_ndn)
  
  uncertain <- length(which(all_inherited$duoNovo_classification == "uncertain"))
  uncertain2 <- length(which(all_inherited$duoNovo_classification == "de_novo" & 
                               all_inherited$GQ_proband < 40 & 
                               (all_inherited$n_de_novo_left_orientation_same_PS == 1 | 
                                  all_inherited$n_de_novo_right_orientation_same_PS == 1)))
  
  false_pos <- c(dn, ndn, uncertain + uncertain2)
  names(false_pos) <- c("dn", "ndn", "uncertain")
  false_pos
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
false_pos_pf <- getFalsePositives(duoNovo_output_filepath_pf, duo_type = "PF")

### --- mother-proband duos
false_pos_pm <- getFalsePositives(duoNovo_output_filepath_pm, duo_type = "PM")

dir_tag <- basename(current_dir)
save(false_pos_pf, false_pos_pm, 
     file = paste0(dir_tag, "_false_pos_pm_and_pf.rda"))






