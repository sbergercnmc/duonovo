args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript false_positive_sources.R <data_directory> <figure_directory>")
}

data_directory   <- args[1]
figure_directory <- args[2]

# ensure data directory exists
if (!dir.exists(data_directory)) {
  stop(sprintf("Data directory not found: '%s'", data_directory))
}

setwd(data_directory)
all_dirs <- list.files()

ppv_mat_extended_pf <- matrix(NA, nrow = 6, length(all_dirs))
rownames(ppv_mat_extended_pf) <- c('giab_problematic_dn', 'giab_problematic_total_classified', 
                                'clustered_dn', 'clustered_total_classified', 
                                'exonic_intronic_dn', 'exonic_intronic_total_classified')
ppv_mat_extended_pm <- ppv_mat_extended_pf

### PF duos
for (i in 1:length(all_dirs)){
  setwd(all_dirs[i])
  all_data <- list.files()
  load(file = all_data[grep('de_novo_variant_granges', all_data)])
  
  ## with giab problematic
  de_novos_pf <- de_novos_pf_with_giab_problematic[!grepl("\\.", de_novos_pf_with_giab_problematic$parentValidation_gt)]
  keep <- which(de_novos_pf$parentValidation_GQ >= 30 & de_novos_pf$parentValidation_depth >= 20)
  de_novos_pf <- de_novos_pf[keep]
  
  #in_gnomad <- which(de_novos_pf$gnomad41_genome_AF > 0)
  #if(length(in_gnomad) > 0){
  #  de_novos_pf <- de_novos_pf[-in_gnomad]
  #}
  ppv_mat_extended_pf['giab_problematic_dn', i] <- table(de_novos_pf$parentValidation_gt)["0/0"]
  ppv_mat_extended_pf['giab_problematic_total_classified', i] <- sum(table(de_novos_pf$parentValidation_gt))
  
  ## with clustered
  de_novos_pf <- de_novos_pf_with_clustered[!grepl("\\.", de_novos_pf_with_clustered$parentValidation_gt)]
  keep <- which(de_novos_pf$parentValidation_GQ >= 30 & de_novos_pf$parentValidation_depth >= 20)
  de_novos_pf <- de_novos_pf[keep]
  
  #in_gnomad <- which(de_novos_pf$gnomad41_genome_AF > 0)
  #if(length(in_gnomad) > 0){
  #  de_novos_pf <- de_novos_pf[-in_gnomad]
  #}
  ppv_mat_extended_pf['clustered_dn', i] <- table(de_novos_pf$parentValidation_gt)["0/0"]
  ppv_mat_extended_pf['clustered_total_classified', i] <- sum(table(de_novos_pf$parentValidation_gt))
  
  ## exonic/intronic only
  de_novos_pf <- de_novos_pf_exonic_intronic_only[!grepl("\\.", de_novos_pf_exonic_intronic_only$parentValidation_gt)]
  keep <- which(de_novos_pf$parentValidation_GQ >= 30 & de_novos_pf$parentValidation_depth >= 20)
  de_novos_pf <- de_novos_pf[keep]
  
  in_gnomad <- which(de_novos_pf$gnomad41_genome_AF > 0)
  if(length(in_gnomad) > 0){
    de_novos_pf <- de_novos_pf[-in_gnomad]
  }
  ppv_mat_extended_pf['exonic_intronic_dn', i] <- table(de_novos_pf$parentValidation_gt)["0/0"]
  ppv_mat_extended_pf['exonic_intronic_total_classified', i] <- sum(table(de_novos_pf$parentValidation_gt))
  
  setwd(data_directory)
}

### PM duos
for (i in 1:length(all_dirs)){
  setwd(all_dirs[i])
  all_data <- list.files()
  load(file = all_data[grep('de_novo_variant_granges', all_data)])
  
  ## with giab problematic
  de_novos_pm <- de_novos_pm_with_giab_problematic[!grepl("\\.", de_novos_pm_with_giab_problematic$parentValidation_gt)]
  keep <- which(de_novos_pm$parentValidation_GQ >= 30 & de_novos_pm$parentValidation_depth >= 20)
  de_novos_pm <- de_novos_pm[keep]
  
  #in_gnomad <- which(de_novos_pm$gnomad41_genome_AF > 0)
  #if(length(in_gnomad) > 0){
  #  de_novos_pm <- de_novos_pm[-in_gnomad]
  #}
  ppv_mat_extended_pm['giab_problematic_dn', i] <- table(de_novos_pm$parentValidation_gt)["0/0"]
  ppv_mat_extended_pm['giab_problematic_total_classified', i] <- sum(table(de_novos_pm$parentValidation_gt))
  
  ## with clustered
  de_novos_pm <- de_novos_pm_with_clustered[!grepl("\\.", de_novos_pm_with_clustered$parentValidation_gt)]
  keep <- which(de_novos_pm$parentValidation_GQ >= 30 & de_novos_pm$parentValidation_depth >= 20)
  de_novos_pm <- de_novos_pm[keep]
  
  #in_gnomad <- which(de_novos_pm$gnomad41_genome_AF > 0)
  #if(length(in_gnomad) > 0){
  #  de_novos_pm <- de_novos_pm[-in_gnomad]
  #}
  ppv_mat_extended_pm['clustered_dn', i] <- table(de_novos_pm$parentValidation_gt)["0/0"]
  ppv_mat_extended_pm['clustered_total_classified', i] <- sum(table(de_novos_pm$parentValidation_gt))
  
  ## exonic/intronic only
  de_novos_pm <- de_novos_pm_exonic_intronic_only[!grepl("\\.", de_novos_pm_exonic_intronic_only$parentValidation_gt)]
  keep <- which(de_novos_pm$parentValidation_GQ >= 30 & de_novos_pm$parentValidation_depth >= 20)
  de_novos_pm <- de_novos_pm[keep]
  
  in_gnomad <- which(de_novos_pm$gnomad41_genome_AF > 0)
  if(length(in_gnomad) > 0){
    de_novos_pm <- de_novos_pm[-in_gnomad]
  }
  ppv_mat_extended_pm['exonic_intronic_dn', i] <- table(de_novos_pm$parentValidation_gt)["0/0"]
  ppv_mat_extended_pm['exonic_intronic_total_classified', i] <- sum(table(de_novos_pm$parentValidation_gt))
  
  setwd(data_directory)
}

ppv_all_pf_giab <- sum(ppv_mat_extended_pf['giab_problematic_dn', ], 
                       na.rm = TRUE)/sum(ppv_mat_extended_pf['giab_problematic_total_classified', ], na.rm = TRUE)
ppv_all_pf_clustered <- sum(ppv_mat_extended_pf['clustered_dn', ], 
                       na.rm = TRUE)/sum(ppv_mat_extended_pf['clustered_total_classified', ], na.rm = TRUE)
ppv_all_pf_exonic_intronic <- sum(ppv_mat_extended_pf['exonic_intronic_dn', ], 
                       na.rm = TRUE)/sum(ppv_mat_extended_pf['exonic_intronic_total_classified', ], na.rm = TRUE)

ppv_all_pm_giab <- sum(ppv_mat_extended_pm['giab_problematic_dn', ], 
                       na.rm = TRUE)/sum(ppv_mat_extended_pm['giab_problematic_total_classified', ], na.rm = TRUE)
ppv_all_pm_clustered <- sum(ppv_mat_extended_pm['clustered_dn', ], 
                            na.rm = TRUE)/sum(ppv_mat_extended_pm['clustered_total_classified', ], na.rm = TRUE)
ppv_all_pm_exonic_intronic <- sum(ppv_mat_extended_pm['exonic_intronic_dn', ], 
                                  na.rm = TRUE)/sum(ppv_mat_extended_pm['exonic_intronic_total_classified', ], na.rm = TRUE)


pdf(file = paste0(figure_directory, "/collective_ppv_exonic_intronic.pdf"), height = 3.2, width = 2, pointsize = 8)
df <- data.frame(
  Condition = factor(c("PF", "PM", "PF (ex.in)", "PM (ex.in)"),
                     levels = c("PF", "PM", "PF (ex.in)", "PM (ex.in)")),
  PPV       = c(ppv_all_pf_no_gnomad, ppv_all_pm_no_gnomad, 
                ppv_all_pf_exonic_intronic, ppv_all_pm_exonic_intronic)
)

ggplot(df, aes(x = Condition, y = PPV, fill = alpha("red", 0.75))) +
  geom_col() +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "",
    y     = "PPV",
    x     = ""
  ) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x    = element_text(angle = 45, hjust = 1))
dev.off()



