getHaplotypesForTrio <- function(phased_variant_granges_trio) {
  granges <- phased_variant_granges_trio
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
  
  ### Extracting haplotypes for the sibling
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
  
  ### Extracting haplotypes for the mother
  granges$hap13 <- NA
  granges$hap23 <- NA
  
  granges$hap13[grep("0\\|", granges$phasing3)] <- 0
  granges$hap13[grep("1\\|", granges$phasing3)] <- 1
  granges$hap13[grep("2\\|", granges$phasing3)] <- 2
  granges$hap13[grep("3\\|", granges$phasing3)] <- 3
  granges$hap13[grep("4\\|", granges$phasing3)] <- 4
  granges$hap13[grep("0/0", granges$phasing3)] <- 0
  granges$hap13[grep("1/1", granges$phasing3)] <- 1
  
  granges$hap23[grep("\\|0", granges$phasing3)] <- 0
  granges$hap23[grep("\\|1", granges$phasing3)] <- 1
  granges$hap23[grep("\\|2", granges$phasing3)] <- 2
  granges$hap23[grep("\\|3", granges$phasing3)] <- 3
  granges$hap23[grep("\\|4", granges$phasing3)] <- 4
  granges$hap23[grep("0/0", granges$phasing3)] <- 0
  granges$hap23[grep("1/1", granges$phasing3)] <- 1
  
  granges
}



getHaplotypeBlockCoordinatesForTrio <- function(haplotype_granges) {
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
  
  # Split haplotype ranges by PS1 and calculate coordinates
  granges_by_ps3 <- split(hap_granges, hap_granges$PS3)
  start_coords <- start(granges_by_ps3)
  start_coords <- sapply(start_coords, min)
  end_coords <- end(granges_by_ps3)
  end_coords <- sapply(end_coords, max)
  all_chrs <- seqnames(granges_by_ps3)  
  all_chrs <- unlist(unique(all_chrs)) #Assuming all sequences are from the same chromosome in each list element
  
  ps3_boundaries <- GRanges(
    seqnames = all_chrs, ranges = IRanges(start = start_coords, end = end_coords),
    PS3 = names(granges_by_ps3)
  )
  
  
  ## ------------ pairwise intersection: PS1 ∩ PS2 ------------------
  ov12      <- findOverlaps(ps1_boundaries, ps2_boundaries, type = "any")
  inter12   <- pintersect(ps1_boundaries[queryHits(ov12)],
                          ps2_boundaries[subjectHits(ov12)],
                          drop.nohit.ranges = TRUE)
  
  ## ------------ add the third sample: (PS1 ∩ PS2) ∩ PS3 ----------
  ov123     <- findOverlaps(inter12, ps3_boundaries, type = "any")
  inter123  <- pintersect(inter12[queryHits(ov123)],
                          ps3_boundaries[subjectHits(ov123)],
                          drop.nohit.ranges = TRUE)
  
  ## ------------ assemble one row per 3-way shared segment ---------
  combined_PS <- inter123                              
  mcols(combined_PS) <- DataFrame(
    PS1 = ps1_boundaries$PS1[
      queryHits(ov12)[queryHits(ov123)]           ],
    PS2 = ps2_boundaries$PS2[
      subjectHits(ov12)[queryHits(ov123)]         ],
    PS3 = ps3_boundaries$PS3[
      subjectHits(ov123)                          ]
  )
  
  valid_ranges <- start(combined_PS) <= end(combined_PS)
  combined_PS <- combined_PS[valid_ranges]
  combined_PS
}