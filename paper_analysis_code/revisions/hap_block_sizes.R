getHaplotypes <- function(phased_variant_granges_duo) {
  granges <- phased_variant_granges_duo
  granges$hap11 <- NA
  granges$hap21 <- NA
  
  ### Extracting haplotypes for the proband
  granges$hap11[grep("0\\|", granges$phasing1)] <- 0
  granges$hap11[grep("1\\|", granges$phasing1)] <- 1
  granges$hap11[grep("2\\|", granges$phasing1)] <- 2
  granges$hap11[grep("3\\|", granges$phasing1)] <- 3
  granges$hap11[grep("4\\|", granges$phasing1)] <- 4
  granges$hap11[grep("0/0", granges$phasing1)] <- 0
  granges$hap11[grep("1/1", granges$phasing1)] <- 1
  
  granges$hap21[grep("\\|0", granges$phasing1)] <- 0
  granges$hap21[grep("\\|1", granges$phasing1)] <- 1
  granges$hap21[grep("\\|2", granges$phasing1)] <- 2
  granges$hap21[grep("\\|3", granges$phasing1)] <- 3
  granges$hap21[grep("\\|4", granges$phasing1)] <- 4
  granges$hap21[grep("0/0", granges$phasing1)] <- 0
  granges$hap21[grep("1/1", granges$phasing1)] <- 1
  
  ### Extracting haplotypes for the parent
  granges$hap12 <- NA
  granges$hap22 <- NA
  
  granges$hap12[grep("0\\|", granges$phasing2)] <- 0
  granges$hap12[grep("1\\|", granges$phasing2)] <- 1
  granges$hap12[grep("2\\|", granges$phasing2)] <- 2
  granges$hap12[grep("3\\|", granges$phasing2)] <- 3
  granges$hap12[grep("4\\|", granges$phasing2)] <- 4
  granges$hap12[grep("0/0", granges$phasing2)] <- 0
  granges$hap12[grep("1/1", granges$phasing2)] <- 1
  
  granges$hap22[grep("\\|0", granges$phasing2)] <- 0
  granges$hap22[grep("\\|1", granges$phasing2)] <- 1
  granges$hap22[grep("\\|2", granges$phasing2)] <- 2
  granges$hap22[grep("\\|3", granges$phasing2)] <- 3
  granges$hap22[grep("\\|4", granges$phasing2)] <- 4
  granges$hap22[grep("0/0", granges$phasing2)] <- 0
  granges$hap22[grep("1/1", granges$phasing2)] <- 1
  
  granges
}

getHaplotypeBlockCoordinates <- function(haplotype_granges) {
  hap_granges <- haplotype_granges
  
  # Split haplotype ranges by PS1 and calculate coordinates
  granges_by_ps1 <- split(hap_granges, hap_granges$PS1)
  start_coords <- start(granges_by_ps1)
  start_coords <- sapply(start_coords, min)
  end_coords <- end(granges_by_ps1)
  end_coords <- sapply(end_coords, max)
  all_chrs <- seqnames(granges_by_ps1)  
  all_chrs <- unlist(unique(all_chrs)) #Assuming all sequences are from the same chromosome in each list element
  
  ps1_boundaries <- GRanges(
    seqnames = all_chrs, ranges = IRanges(start = start_coords, end = end_coords),
    PS1 = names(granges_by_ps1)
  )
  
  # Split haplotype ranges by PS2 and calculate coordinates
  granges_by_ps2 <- split(hap_granges, hap_granges$PS2)
  
  start_coords <- start(granges_by_ps2)
  start_coords <- sapply(start_coords, min)
  end_coords <- end(granges_by_ps2)
  end_coords <- sapply(end_coords, max)
  all_chrs <- seqnames(granges_by_ps2)  
  all_chrs <- unlist(unique(all_chrs)) #Assuming all sequences are from the same chromosome in each list element
  
  ps2_boundaries <- GRanges(
    seqnames = all_chrs, ranges = IRanges(start = start_coords, end = end_coords),
    PS2 = names(granges_by_ps2)
  )
  
  # Find overlaps between ps1_boundaries and ps2_boundaries
  overlaps <- findOverlaps(ps1_boundaries, ps2_boundaries)
  
  # Create a new GRanges object to store the combined information
  combined_PS <- GRanges(
    seqnames = seqnames(ps1_boundaries)[queryHits(overlaps)],
    ranges = IRanges(
      start = pmax(start(ps1_boundaries)[queryHits(overlaps)], start(ps2_boundaries)[subjectHits(overlaps)]),
      end = pmin(end(ps1_boundaries)[queryHits(overlaps)], end(ps2_boundaries)[subjectHits(overlaps)])
    ),
    PS1 = ps1_boundaries$PS1[queryHits(overlaps)],
    PS2 = ps2_boundaries$PS2[subjectHits(overlaps)]
  )
  
  # Remove any ranges that have become invalid (start > end) due to the pmax/pmin operation
  valid_ranges <- start(combined_PS) <= end(combined_PS)
  combined_PS <- combined_PS[valid_ranges]
  combined_PS
}

calculateHapBlockSize <- function(LRS_phased_vcf_file_path, depth_cutoff = 20, GQ_cutoff = 30,
                                  PS_width_cutoff = 10000, boundary_cutoff = 2000){
  if (!file.exists(LRS_phased_vcf_file_path)) {
    stop("The LRS VCF file path does not exist: ", LRS_phased_vcf_file_path)
  }
  
  message("Importing LRS VCF file...")
  vcf <- readVcf(LRS_phased_vcf_file_path)
  vcf_granges <- SummarizedExperiment::rowRanges(vcf)
  vcf_metadata <- geno(vcf)
  
  proband_column <- 1
  if (length(proband_column) == 0) {
    stop("Proband column identifier not found: ", proband_column_identifier)
  }
  
  if (dim(vcf)[2] != 2) {
    stop("VCF does not contain exactly 2 samples.  Number of samples: ", length(samples(header(vcf))))
  }
  
  if (!"GT" %in% names(vcf_metadata)) {
    stop("The 'GT' (genotype) field is missing in the VCF metadata from LRS.")
  }
  if (!"DP" %in% names(vcf_metadata)) {
    stop("The 'DP' (depth) field is missing in the VCF metadata from LRS.")
  }
  if (!"GQ" %in% names(vcf_metadata)) {
    stop("The 'GQ' (genotype quality) field is missing in the VCF metadata from LRS.")
  }
  if (!"PS" %in% names(vcf_metadata)) {
    stop("The 'PS' (phasing set) field is missing in the VCF metadata from LRS.")
  }
  
  vcf_granges$phasing1 <- vcf_metadata$GT[, proband_column]
  vcf_granges$phasing2 <- vcf_metadata$GT[, -proband_column]
  
  vcf_granges$depth1 <- vcf_metadata$DP[, proband_column]
  vcf_granges$depth2 <- vcf_metadata$DP[, -proband_column]
  
  vcf_granges$GQ1 <- vcf_metadata$GQ[, proband_column]
  vcf_granges$GQ2 <- vcf_metadata$GQ[, -proband_column]
  
  #the following extracts the phasing set to which each phased variant is assigned to in the proband and in the parent
  #and ensures the phasing set is represented by a unique identifier in the columns PS1 (for proband) and PS2 (for parent)
  vcf_granges$PS1 <- vcf_metadata$PS[, proband_column]
  has_PS1 <- !is.na(vcf_granges$PS1)
  vcf_granges$PS1[has_PS1] <- paste0(seqnames(vcf_granges[has_PS1]), "_", vcf_granges$PS1[has_PS1])
  vcf_granges$PS2 <- vcf_metadata$PS[, -proband_column]
  has_PS2 <- !is.na(vcf_granges$PS2)
  vcf_granges$PS2[has_PS2] <- paste0(seqnames(vcf_granges[has_PS2]), "_", vcf_granges$PS2[has_PS2])
  
  message("Reconstructing haplotypes...")
  hap_granges <- getHaplotypes(vcf_granges)
  
  message("Calculating cumulative size of haplotype blocks...")
  hap_boundary_coordinate_ranges <- getHaplotypeBlockCoordinates(hap_granges)
  haplotype_boundary_coordinate_granges_filtered <- haplotype_boundary_coordinate_granges[
    which(width(haplotype_boundary_coordinate_granges) > PS_width_cutoff)]
  haplotype_boundary_coordinate_granges_filtered <- haplotype_boundary_coordinate_granges_filtered - boundary_cutoff
  
  candidate_variant_indices <- which(hap_granges$phasing1 %in% c("1|0", "0|1") & ranges_to_subset$phasing2 == "0/0")
  candidate_variant_ranges <- hap_granges[candidate_variant_indices]
  overlaps <- findOverlaps(candidate_variant_indices, hap_boundary_coordinate_ranges)
  phased_candidates <- length(unique(queryHits(overlaps)))

  c(sum(width(hap_boundary_coordinates)), sum(width(haplotype_boundary_coordinate_granges_filtered)), 
    length(phased_candidates), length(candidate_variant_indices))
}


###
###
library(VariantAnnotation)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: hap_block_sizes.R <run_directory>")
}
## set directory
current_dir <- args[1]
setwd(current_dir)

### use duoNovo output from duos to separate candidate variants based on classification
### PM duo first
duoNovo_input_filepath_pm <- list.files(pattern = "PM\\.vcf\\.gz$")
message("calculating haplotype block sizes from PM duo...")
hap_block_size_pm <- calculateHapBlockSize(duoNovo_input_filepath_pm)

### now PF duo
###
duoNovo_intput_filepath_pf <- list.files(pattern = "PF\\.vcf\\.gz$")
message("calculating haplotype block sizes from PF duo...")
hap_block_size_pf <- calculateHapBlockSize(duoNovo_input_filepath_pf)

sample_id <- sub("\\..*$", "", duoNovo_input_filepath_pm)
save(hap_block_size_pm, hap_block_size_pf,
     file = paste0(sample_id, "_pm_pf_hap_block_sizes.rda"))


