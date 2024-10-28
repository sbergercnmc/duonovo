duoNovo <- function(haplotype_granges, phasing_set_granges, proband_phasing,
                    PS_width_cutoff, boundary_cutoff, distance_cutoff, 
                    short_read_vcf){
  hap_granges <- haplotype_granges
  combined_PS <- phasing_set_granges
  combined_PS <- combined_PS[which(width(combined_PS) > PS_width_cutoff)]
  
  potential_de_novo_indices_allele1 <- which(hap_granges$phasing1 == proband_phasing & hap_granges$phasing2 == "0/0")
  hap_granges1 <- hap_granges[potential_de_novo_indices_allele1]
  
  ###the following only retains candidates where the het call
  ###in proband and hom reference call in parent is confirmed with short reads
  v_granges_duo_SR <- short_read_vcf
  v_granges_duo_SR <- v_granges_duo_SR[which(v_granges_duo_SR$gt1 %in% c("1/0", "0/1") & 
                                               v_granges_duo_SR$gt2 == "0/0")]
  hap_granges1 <- hap_granges1[unique(queryHits(findOverlaps(hap_granges1, v_granges_duo_SR)))]
  ###
  
  overlaps <- findOverlaps(hap_granges1, combined_PS - boundary_cutoff)
  overlapping_indices <- split(queryHits(overlaps), subjectHits(overlaps))
  
  haplotypes <- vector("list", length(overlapping_indices))
  names(haplotypes) <- names(overlapping_indices)
  
  indices <- as.numeric(names(overlapping_indices))
  PS1_ids <- combined_PS$PS1[indices]
  PS2_ids <- combined_PS$PS2[indices]
  
  hap_granges_no_denovo <- hap_granges[-potential_de_novo_indices_allele1]
  for(i in 1:length(indices)) {
    #index <- indices[i]
    variant_granges <- hap_granges_no_denovo[
      unique(queryHits(findOverlaps(hap_granges_no_denovo, combined_PS[
        which(combined_PS$PS1 == PS1_ids[i] & combined_PS$PS2 == PS2_ids[i])])))]
    hap11 <- variant_granges$hap11
    hap12 <- variant_granges$hap12
    hap21 <- variant_granges$hap21
    hap22 <- variant_granges$hap22  
    #
    hap_mat <- matrix(c(hap11, hap12, hap21, hap22), nrow = length(hap11), ncol = 4, byrow = FALSE)
    colnames(hap_mat) <- c("hap11", "hap12", "hap21", "hap22")
    haplotypes[[i]] <- hap_mat
    haplotypes
  }
  
  hamming_distance_mat <- sapply(haplotypes, function(xx) 
    c(sum(xx[, "hap11"] != xx[, "hap12"], na.rm = TRUE), 
      sum(xx[, "hap11"] != xx[, "hap22"], na.rm = TRUE), 
      sum(xx[, "hap21"] != xx[, "hap12"], na.rm = TRUE), 
      sum(xx[, "hap21"] != xx[, "hap22"], na.rm = TRUE))
  )
  
  hamming_distance_mins <- colMins(hamming_distance_mat)
  number_zeros <- colSums(hamming_distance_mat == 0)
  
  inheritance_stringent_hap1 <- apply(hamming_distance_mat, 2,
                                      function(xx) (xx[1] == 0 | xx[2] == 0) & 
                                        xx[3] > distance_cutoff & xx[4] > distance_cutoff)
  
  
  inheritance_stringent_hap2 <- apply(hamming_distance_mat, 2,
                                      function(xx) (xx[3] == 0 | xx[4] == 0) & 
                                        xx[1] > distance_cutoff & xx[2] > distance_cutoff)
  
  hap1_inherited <- unlist(overlapping_indices[which(inheritance_stringent_hap1 == TRUE)])
  hap2_inherited <- unlist(overlapping_indices[which(inheritance_stringent_hap2 == TRUE)])
  uncertain <- unlist(overlapping_indices[which(hamming_distance_mins > 0)])
  
  if(proband_phasing == "1|0"){
    result_list <- GRangesList(de_novo = hap_granges1[hap1_inherited], not_de_novo = hap_granges1[hap2_inherited], 
                               uncertain = hap_granges1[uncertain])
  } else if (proband_phasing == "0|1"){
    result_list <- GRangesList(de_novo = hap_granges1[hap2_inherited], not_de_novo = hap_granges1[hap1_inherited], 
                               uncertain = hap_granges1[uncertain])
  }
  result_list
} 