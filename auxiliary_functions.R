getQualFilteredVariantGranges <- function(vcf, depth_threshold, 
                                          GQ_threshold){ #vcf is read with the Variant Annotation R package
  v <- vcf
  v_granges <- rowRanges(v)
  
  proband_column <- grep("-0$", samples(header(v)))
  vcf_metadata <- geno(v)
  
  v_granges$phasing1 <- vcf_metadata[[1]][, proband_column]
  v_granges$phasing2 <- vcf_metadata[[1]][, -proband_column]
  
  v_granges$depth1 <- vcf_metadata[[3]][, proband_column]
  v_granges$depth2 <- vcf_metadata[[3]][, -proband_column]
  
  v_granges$GQ1 <- vcf_metadata[[5]][, proband_column]
  v_granges$GQ2 <- vcf_metadata[[5]][, -proband_column]
  
  v_granges$PS1 <- vcf_metadata[[7]][, proband_column]
  v_granges$PS2 <- vcf_metadata[[7]][, -proband_column]
  
  v_granges_filtered <- v_granges[which(v_granges$depth1 >= depth_threshold & v_granges$GQ1 >= GQ_threshold & 
                                          v_granges$depth2 >= depth_threshold & v_granges$GQ2 >= GQ_threshold)]
  v_granges_filtered
}

getQualFilteredVariantGrangesSR <- function(vcf, depth_threshold, 
                                            GQ_threshold){ #vcf is read with the Variant Annotation R package
  v <- vcf
  v_granges <- rowRanges(v)
  
  proband_column <- grep("-0$", samples(header(v)))
  vcf_metadata <- geno(v)
  
  v_granges$gt1 <- vcf_metadata$GT[, proband_column]
  v_granges$gt2 <- vcf_metadata$GT[, -proband_column]
  
  v_granges$depth1 <- vcf_metadata$DP[, proband_column]
  v_granges$depth2 <- vcf_metadata$DP[, -proband_column]
  
  v_granges$GQ1 <- vcf_metadata$GQ[, proband_column]
  v_granges$GQ2 <- vcf_metadata$GQ[, -proband_column]
  
  v_granges_filtered <- v_granges[which(v_granges$depth1 >= depth_threshold & v_granges$GQ1 >= GQ_threshold & 
                                          v_granges$depth2 >= depth_threshold & v_granges$GQ2 >= GQ_threshold)]
  v_granges_filtered
}

getHaplotypes <- function(phased_vcf_granges){
  granges <- phased_vcf_granges
  granges$hap11 <- NA
  granges$hap21 <- NA
  
  ###proband
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
  granges$hap21[grep("0/0",  granges$phasing1)] <- 0
  granges$hap21[grep("1/1", granges$phasing1)] <- 1
  
  ###parent
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
  granges$hap22[grep("0/0",  granges$phasing2)] <- 0
  granges$hap22[grep("1/1", granges$phasing2)] <- 1
  
  granges
}

getPhasingSetCoordinates <- function(haplotype_granges){
  hap_granges <- haplotype_granges
  
  granges_by_ps1 <- split(hap_granges, hap_granges$PS1)
  
  all_chrs <- sapply(granges_by_ps1, function(xx) unique(seqnames(xx)))
  indices_to_use <- which(lengths(all_chrs) == 1)
  start_coords <- sapply(granges_by_ps1[indices_to_use], function(xx) min(start(xx)))
  end_coords <- sapply(granges_by_ps1[indices_to_use], function(xx) max(end(xx)))
  all_chrs <- unlist(all_chrs[indices_to_use])
  
  ps1_boundaries <- GRanges(
    seqnames = all_chrs, ranges = IRanges(start = start_coords, end = end_coords),
    PS1 = names(granges_by_ps1)[indices_to_use]
  )
  
  granges_by_ps2 <- split(hap_granges, hap_granges$PS2)
  start_coords <- sapply(granges_by_ps2, function(xx) min(start(xx)))
  end_coords <- sapply(granges_by_ps2, function(xx) max(end(xx)))
  all_chrs <- sapply(granges_by_ps2, function(xx) unique(seqnames(xx)))
  
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