getDuoNovoOutputGRanges <- function(duoNovo_file_path){
  if (!file.exists(duoNovo_file_path)) {
    stop("The duoNovo VCF file path does not exist: ", duoNovo_file_path )
  }
  
  message("Reading duoNovo VCF file...")
  dnvcf <- readVcf(duoNovo_file_path)
  if (is.null(meta(header(dnvcf))$duoNovoPARAM)){
    stop("DuoNovoHeader Parameters not found")
  }
  
  if (!"GT" %in% names(geno(dnvcf)) ) {
    stop("The 'GT' (genotype) field is missing in the VCF metadata from LRS.")
  }
  
  proband_column_identifier <- meta(header(dnvcf))$duoNovoPARAM['proband_column_identifier', ]
  proband_column <- grep(proband_column_identifier, samples(header(dnvcf)))
  if (length(proband_column) != 1) {
    stop("Proband column identifier not found unambiguously in duoNovo Vcf: ", proband_column_identifier)
  }
  proband_id <- samples(header(dnvcf))[proband_column]
  
  missing_parent_column <- 3
  sequenced_parent_column <- setdiff(seq(3), c(proband_column, missing_parent_column))
  
  g <- geno(dnvcf)
  parentValidation_trio_gt <- g$GT[, missing_parent_column]
  parentValidation_trio_dp <- g$DP[, missing_parent_column]
  parentValidation_trio_gq <- g$GQ[, missing_parent_column]
  
  dnvcf_info <- info(dnvcf)
  duoNovo_ranges <- rowRanges(dnvcf)
  columns_to_add <- c(
    "phasing_proband", "phasing_parent", "depth_proband", "depth_parent",
    "GQ_proband", "GQ_parent", "duoNovo_classification",
    "supporting_hamming_distance", "supporting_counts_het_hom", "supporting_counts_het_het",
    "supporting_counts_hom_het", "QC_fail_step", 
    "n_de_novo_left_orientation_same_PS", "n_de_novo_right_orientation_same_PS",
    "ANNOVAR_DATE",
    "Func.refGeneWithVer", "Gene.refGeneWithVer", "GeneDetail.refGeneWithVer",
    "ExonicFunc.refGeneWithVer", "AAChange.refGeneWithVer",
    "CpG",
    "problematic_region",
    "gnomad41_genome_AF", "gnomad41_genome_AF_raw", "gnomad41_genome_AF_XX",
    "gnomad41_genome_AF_XY", "gnomad41_genome_AF_grpmax",
    "gnomad41_genome_faf95", "gnomad41_genome_faf99",
    "gnomad41_genome_fafmax_faf95_max", "gnomad41_genome_fafmax_faf99_max",
    "gnomad41_genome_AF_afr", "gnomad41_genome_AF_ami",
    "gnomad41_genome_AF_amr", "gnomad41_genome_AF_asj",
    "gnomad41_genome_AF_eas", "gnomad41_genome_AF_fin",
    "gnomad41_genome_AF_mid", "gnomad41_genome_AF_nfe",
    "gnomad41_genome_AF_remaining", "gnomad41_genome_AF_sas",
    "ALLELE_END")
  missing_info <- setdiff(columns_to_add, colnames(info(dnvcf)))
  if (length(missing_info))
    stop("VCF INFO is missing fields: ", paste(missing_info, collapse=", "))
  
  mcols(duoNovo_ranges) <- dnvcf_info[, columns_to_add]
  
  duoNovo_ranges$parentValidation_gt <- parentValidation_trio_gt
  duoNovo_ranges$parentValidation_depth <- parentValidation_trio_dp
  duoNovo_ranges$parentValidation_GQ <- parentValidation_trio_gq
  duoNovo_ranges
}


library(VariantAnnotation)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: process_duoNovo_output.R <trio_directories.txt> <index>")
}

dir_file   <- args[1]
dir_index  <- as.integer(args[2])

## read directory list
dirs <- trimws(readLines(dir_file))
if (dir_index < 1L || dir_index > length(dirs))
  stop("Index ", dir_index, " is out of range 1–", length(dirs))

current_dir <- dirs[dir_index]
setwd(current_dir)

duonovo_vcf_output_filepaths <- list.files(pattern = "PF\\.duonovo\\.addedParent\\.dnm2|PM\\.duonovo\\.addedParent\\.dnm2")
dn_granges_pf <- getDuoNovoOutputGRanges(duonovo_vcf_output_filepaths[1])
dn_granges_pm <- getDuoNovoOutputGRanges(duonovo_vcf_output_filepaths[2])

duo_basename_pf <- basename(duonovo_vcf_output_filepaths[1])
duo_basename_pm <- basename(duonovo_vcf_output_filepaths[2])
output_filename_pf <- paste0(sub("\\.vcf\\.gz$", "", duo_basename_pf), ".rda")
output_filename_pm <- paste0(sub("\\.vcf\\.gz$", "", duo_basename_pm), ".rda")

save(dn_granges_pf, file = output_filename_pf)
save(dn_granges_pm, file = output_filename_pm)


