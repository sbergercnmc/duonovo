# Load necessary libraries
setwd('/scratch/PMGRC_LRdata/trios/data')
library(VariantAnnotation)

# Check if an argument is provided (subdirectory path)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No directory specified", call. = FALSE)
}

subdirs <- c("PMGRC-107-107-0", "PMGRC-323-323-0", "PMGRC-423-423-0", "PMGRC-5-5-0",
             "PMGRC-191-191-0", "PMGRC-332-332-0", "PMGRC-43-43-0", "PMGRC-615-615-0",
             "PMGRC-291-291-0", "PMGRC-351-351-0", "PMGRC-445-445-0", "PMGRC-645-645-0",
             "PMGRC-296-296-0", "PMGRC-372-372-0", "PMGRC-482-482-0", "PMGRC-661-661-0",
             "PMGRC-320-320-0", "PMGRC-417-417-0", "PMGRC-522-522-0", "PMGRC-770-770-0")


# Convert args[1] to an integer for indexing
task_id <- as.integer(args[1])
subdirectory <- subdirs[task_id]

setwd(subdirectory)

# Read de novo vcf file (obtained from short read sequencing) in the specified subdirectory
vcf <- readVcf('trio.shortRead.denovos.vcf.gz', genome = "hg38")

genotypes <- geno(vcf)$GT

proband_column <- grep("-0$", colnames(genotypes))
father_column <- grep("-1$", colnames(genotypes)) 
mother_column <- grep("-2$", colnames(genotypes))

child_genotype <- genotypes[, proband_column]
father_genotype <- genotypes[, father_column]
mother_genotype <- genotypes[, mother_column]

depth <- geno(vcf)$DP
child_depth <- depth[, proband_column]
father_depth <- depth[, father_column]
mother_depth <- depth[, mother_column]

gq <- geno(vcf)$GQ
child_gq <- gq[, proband_column]
father_gq <- gq[, father_column]
mother_gq <- gq[, mother_column]

###arrange as metadata in the granges
trio_de_novo <- rowRanges(vcf)

trio_de_novo$child_genotype <- child_genotype
trio_de_novo$father_genotype <- father_genotype
trio_de_novo$mother_genotype <- mother_genotype

trio_de_novo$child_depth <- child_depth
trio_de_novo$father_depth <- father_depth
trio_de_novo$mother_depth <- mother_depth

trio_de_novo$child_gq <- child_gq
trio_de_novo$father_gq <- father_gq
trio_de_novo$mother_gq <- mother_gq

setwd('/home/lboukas/PMGRC_longread/duo/vcfs')
save(trio_de_novo, file = paste0("de_novo_candidates/ground_truth_de_novo_SR_", 
                                 basename(subdirectory), ".rda"))

