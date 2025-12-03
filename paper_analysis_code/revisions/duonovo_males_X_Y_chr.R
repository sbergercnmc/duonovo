## import spreadsheet with parental age info
library(readr)
trios <- read_csv(trio_info_directory)
trios <- trios[-which(trios$child_SAMPLEID %in% problematic_trios), ]
##

trios_male_proband <- trios[which(trios$child_SEX == "Male"), ]

sample_names <- trios_male_proband$child_SAMPLEID
pf_dn_grl_male <- GRangesList(
  setNames(
    replicate(length(sample_names), GRanges(), simplify = FALSE),
    sample_names
  ))

pm_dn_grl_male <- GRangesList(
  setNames(
    replicate(length(sample_names), GRanges(), simplify = FALSE),
    sample_names
  ))

for (i in 1:length(sample_names)){
  setwd(sample_names[i])
  all_data <- list.files()
  load(file = all_data[grep('de_novo_variant_granges', all_data)])
  pf_dn_grl_male[[i]] <- de_novos_pf
  pm_dn_grl_male[[i]] <- de_novos_pm
  
  setwd(data_directory)
}

x_denovos_pf <- sapply(pf_dn_grl_male, function(xx) length(which(seqnames(xx) %in% c("chrX"))))
x_denovos_pm <- sapply(pm_dn_grl_male, function(xx) length(which(seqnames(xx) %in% c("chrX"))))
y_denovos_pf <- sapply(pf_dn_grl_male, function(xx) length(which(seqnames(xx) %in% c("chrY"))))
y_denovos_pm <- sapply(pm_dn_grl_male, function(xx) length(which(seqnames(xx) %in% c("chrY"))))



