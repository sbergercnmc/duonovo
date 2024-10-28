setwd('/scratch/PMGRC_LRdata/trios/data')
library(VariantAnnotation)

#check if an argument is provided (subdirectory path)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No directory specified", call. = FALSE)
}

subdirs <- c("PMGRC-107-107-0", "PMGRC-323-323-0", "PMGRC-5-5-0",
             "PMGRC-191-191-0", "PMGRC-615-615-0",
             "PMGRC-291-291-0", "PMGRC-351-351-0", "PMGRC-645-645-0",
             "PMGRC-296-296-0", "PMGRC-372-372-0", "PMGRC-482-482-0", "PMGRC-661-661-0",
             "PMGRC-320-320-0", "PMGRC-417-417-0", "PMGRC-522-522-0", "PMGRC-770-770-0")


task_id <- as.integer(args[1])
subdirectory <- subdirs[task_id]

setwd(subdirectory)

vcf_in_dir <- c("duo_proband_father.longread.hiphase.vcf.gz", 
                "duo_proband_mother.longread.hiphase.vcf.gz", 
                "trio.longread.hiphase.vcf.gz")

vcfs <- lapply(vcf_in_dir, function(x) readVcf(x, genome = "hg38"))

save(vcfs, file = paste0("/home/lboukas/PMGRC_longread/duo/vcfs/", basename(subdirectory), "_vcfs.rda"))
