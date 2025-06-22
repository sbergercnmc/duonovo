getNPV <- function(duoNovo_granges_output_filepath, duo_type = c("PM", "PF")){
  load(file = duoNovo_granges_output_filepath)
  if (duo_type == "PF"){
    dn_granges <- dn_granges_pf
  } else if (duo_type == "PM"){
    dn_granges <- dn_granges_pm
  }
  dn_granges <- dn_granges[which(dn_granges$phasing_parent == "0/0")]
  dn_granges <- dn_granges[-which(dn_granges$QC_fail_step %in% 
                                    c("low_depth", "low_GQ", "low_depth_and_GQ"))] #this is because the naive approach would discard these as well
  dn_granges$giab_problematic <- unlist(dn_granges$giab_problematic)
  dn_granges <- dn_granges[which(dn_granges$giab_problematic == ".")]
  
  dn_granges <- dn_granges[!grepl("\\.", dn_granges$parentValidation_gt)]
  
  true_dn <- grep("1", dn_granges$parentValidation_gt, invert = TRUE)
  naive_npv <- 1 - length(true_dn)/length(dn_granges)
  
  classified_ndn_1 <- which(dn_granges$duoNovo_classification == "on_other_parent_haplotype")
  classified_ndn_2 <- which(dn_granges$duoNovo_classification == "de_novo" & 
                              dn_granges$n_de_novo_left_orientation_same_PS > 1)
  classified_ndn_3 <- which(dn_granges$duoNovo_classification == "de_novo" & 
                              dn_granges$n_de_novo_right_orientation_same_PS > 1)
  
  classified_ndn_granges <- dn_granges[c(classified_ndn_1)]
  true_dn <- grep("1", classified_ndn_granges$parentValidation_gt, invert = TRUE)
  duonovo_npv <- 1 - length(true_dn)/length(classified_ndn_granges)
  total <- length(classified_ndn_granges)
  
  npv <- c(duonovo_npv, naive_npv, total)
  names(npv) <- c("NPV", "naive_NPV", "total")
  npv
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
npv_pf <- getNPV(duoNovo_output_filepath_pf, duo_type = "PF")

### --- mother-proband duos
npv_pm <- getNPV(duoNovo_output_filepath_pm, duo_type = "PM")


dir_tag <- basename(current_dir)
save(npv_pf, npv_pm, 
     file = paste0(dir_tag, "_npv_pm_and_pf.rda"))


