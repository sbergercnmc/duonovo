args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: dn_call_rate_ref_vs_alt.R <dir_list.txt> <index>")
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

dn_call_rate_alt_pf <- rep(NA, 3)
names(dn_call_rate_alt_pf) <- c("dn_count", "total_count", "rate")
dn_call_rate_ref_pf <- dn_call_rate_alt_pf
dn_call_rate_alt_pm <- dn_call_rate_alt_pf
dn_call_rate_ref_pm <- dn_call_rate_alt_pf

### --- father-proband duos
###
  load(file = duoNovo_output_filepath_pf)
  dn_granges <- dn_granges_pf
  
  dn_granges_alt <- dn_granges[which(dn_granges$phasing_parent == "0/0" & dn_granges$GQ_parent >= 40)]
  dn_granges_alt <- dn_granges_alt[-which(dn_granges_alt$duoNovo_classification == "failed_QC")]
  dn_granges_alt$giab_problematic <- unlist(dn_granges_alt$giab_problematic)
  dn_granges_alt <- dn_granges_alt[which(dn_granges_alt$giab_problematic == ".")]

  dn_alt_1 <- which(dn_granges_alt$duoNovo_classification == "de_novo" & 
                      dn_granges_alt$n_de_novo_left_orientation_same_PS == 1 & 
                      dn_granges_alt$GQ_proband >= 40)
  dn_alt_2 <- which(dn_granges_alt$duoNovo_classification == "de_novo" & 
                      dn_granges_alt$n_de_novo_right_orientation_same_PS == 1 & 
                      dn_granges_alt$GQ_proband >= 40)
  dn_call_rate_alt_pf[3] <- (length(dn_alt_1) + length(dn_alt_2))/length(dn_granges_alt)
  dn_call_rate_alt_pf[1] <- length(dn_alt_1) + length(dn_alt_2)
  dn_call_rate_alt_pf[2] <- length(dn_granges_alt)
  
  dn_granges_ref <- dn_granges[which(dn_granges$phasing_parent == "1/1" & dn_granges$GQ_parent >= 40)]
  dn_granges_ref <- dn_granges_ref[-which(dn_granges_ref$duoNovo_classification == "failed_QC")]
  dn_granges_ref$giab_problematic <- unlist(dn_granges_ref$giab_problematic)
  dn_granges_ref <- dn_granges_ref[which(dn_granges_ref$giab_problematic == ".")]
  
  dn_ref_1 <- which(dn_granges_ref$duoNovo_classification == "de_novo" & 
                      dn_granges_ref$n_de_novo_left_orientation_same_PS == 1 & 
                      dn_granges_ref$GQ_proband >= 40)
  dn_ref_2 <- which(dn_granges_ref$duoNovo_classification == "de_novo" & 
                      dn_granges_ref$n_de_novo_right_orientation_same_PS == 1 & 
                      dn_granges_ref$GQ_proband >= 40)
  dn_call_rate_ref_pf[3] <- (length(dn_ref_1) + length(dn_ref_2))/length(dn_granges_ref)
  dn_call_rate_ref_pf[1] <- length(dn_ref_1) + length(dn_ref_2)
  dn_call_rate_ref_pf[2] <- length(dn_granges_ref)

### --- mother-proband duos
###
  load(file = duoNovo_output_filepath_pm)
  dn_granges <- dn_granges_pm
  
  dn_granges_alt <- dn_granges[which(dn_granges$phasing_parent == "0/0" & dn_granges$GQ_parent >= 40)]
  dn_granges_alt <- dn_granges_alt[-which(dn_granges_alt$duoNovo_classification == "failed_QC")]
  dn_granges_alt$giab_problematic <- unlist(dn_granges_alt$giab_problematic)
  dn_granges_alt <- dn_granges_alt[which(dn_granges_alt$giab_problematic == ".")]
  
  dn_alt_1 <- which(dn_granges_alt$duoNovo_classification == "de_novo" & 
                      dn_granges_alt$n_de_novo_left_orientation_same_PS == 1 & 
                      dn_granges_alt$GQ_proband >= 40)
  dn_alt_2 <- which(dn_granges_alt$duoNovo_classification == "de_novo" & 
                      dn_granges_alt$n_de_novo_right_orientation_same_PS == 1 & 
                      dn_granges_alt$GQ_proband >= 40)
  dn_call_rate_alt_pm[3] <- (length(dn_alt_1) + length(dn_alt_2))/length(dn_granges_alt)
  dn_call_rate_alt_pm[1] <- length(dn_alt_1) + length(dn_alt_2)
  dn_call_rate_alt_pm[2] <- length(dn_granges_alt)
  
  dn_granges_ref <- dn_granges[which(dn_granges$phasing_parent == "1/1" & dn_granges$GQ_parent >= 40)]
  dn_granges_ref <- dn_granges_ref[-which(dn_granges_ref$duoNovo_classification == "failed_QC")]
  dn_granges_ref$giab_problematic <- unlist(dn_granges_ref$giab_problematic)
  dn_granges_ref <- dn_granges_ref[which(dn_granges_ref$giab_problematic == ".")]
  
  dn_ref_1 <- which(dn_granges_ref$duoNovo_classification == "de_novo" & 
                      dn_granges_ref$n_de_novo_left_orientation_same_PS == 1 & 
                      dn_granges_ref$GQ_proband >= 40)
  dn_ref_2 <- which(dn_granges_ref$duoNovo_classification == "de_novo" & 
                      dn_granges_ref$n_de_novo_right_orientation_same_PS == 1 & 
                      dn_granges_ref$GQ_proband >= 40)
  dn_call_rate_ref_pm[3] <- (length(dn_ref_1) + length(dn_ref_2))/length(dn_granges_ref)
  dn_call_rate_ref_pm[1] <- length(dn_ref_1) + length(dn_ref_2)
  dn_call_rate_ref_pm[2] <- length(dn_granges_ref)

dir_tag <- basename(current_dir)
save(dn_call_rate_alt_pf, dn_call_rate_ref_pf, dn_call_rate_alt_pm, dn_call_rate_ref_pm,
     file = paste0(dir_tag, "_denovo_call_rate_alt_vs_ref.rda"))


