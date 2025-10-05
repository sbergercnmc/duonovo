args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript ppv_GQ40.R <data_directory> <figure_directory>")
}

data_directory   <- args[1]
figure_directory <- args[2]

# ensure data directory exists
if (!dir.exists(data_directory)) {
  stop(sprintf("Data directory not found: '%s'", data_directory))
}

setwd(data_directory)
all_dirs <- list.files()

ppv_pf_gq40_no_gnomad <- rep(NA, length(all_dirs))
for (i in 1:length(all_dirs)){
  setwd(all_dirs[i])
  all_data <- list.files()
  load(file = all_data[grep('de_novo_variant_granges', all_data)])
  de_novos_pf <- de_novos_pf[!grepl("\\.", de_novos_pf$parentValidation_gt)]
  keep <- which(de_novos_pf$GQ_proband >= 40 & de_novos_pf$GQ_parent >= 40 & de_novos_pf$parentValidation_GQ >= 40 & 
                  de_novos_pf$parentValidation_depth >= 20)
  dn_gq40 <- de_novos_pf[keep]
  in_gnomad <- which(dn_gq40$gnomad41_genome_AF > 0)
  if(length(in_gnomad) > 0){
    dn_gq40 <- dn_gq40[-in_gnomad]
  }
  ppv_pf_gq40_no_gnomad[i] <- table(dn_gq40$parentValidation_gt)["0/0"]/sum(table(dn_gq40$parentValidation_gt))
  
  setwd(data_directory)
}

ppv_pm_gq40_no_gnomad <- rep(NA, length(all_dirs))
for (i in 1:length(all_dirs)){
  setwd(all_dirs[i])
  all_data <- list.files()
  load(file = all_data[grep('de_novo_variant_granges', all_data)])
  de_novos_pm <- de_novos_pm[!grepl("\\.", de_novos_pm$parentValidation_gt)]
  keep <- which(de_novos_pm$GQ_proband >= 40 & de_novos_pm$GQ_parent >= 40 & de_novos_pm$parentValidation_GQ >= 40 & 
                  de_novos_pm$parentValidation_depth >= 20)
  dn_gq40 <- de_novos_pm[keep]
  in_gnomad <- which(dn_gq40$gnomad41_genome_AF > 0)
  if(length(in_gnomad) > 0){
    dn_gq40 <- dn_gq40[-in_gnomad]
  }
  ppv_pm_gq40_no_gnomad[i] <- table(dn_gq40$parentValidation_gt)["0/0"]/sum(table(dn_gq40$parentValidation_gt))
  
  setwd(data_directory)
}


number_called_pm_gq40 <- rep(NA, length(all_dirs))
for (i in 1:length(all_dirs)){
  setwd(all_dirs[i])
  all_data <- list.files()
  load(file = all_data[grep('de_novo_variant_granges', all_data)])
  de_novos_pm <- de_novos_pm[!grepl("\\.", de_novos_pm$parentValidation_gt)]
  keep <- which(de_novos_pm$GQ_proband >= 40 & de_novos_pm$GQ_parent >= 40 & de_novos_pm$parentValidation_GQ >= 40 & 
                  de_novos_pm$parentValidation_depth >= 20)
  dn_gq40 <- de_novos_pm[keep]
  in_gnomad <- which(dn_gq40$gnomad41_genome_AF > 0)
  if(length(in_gnomad) > 0){
    dn_gq40 <- dn_gq40[-in_gnomad]
  }
  number_called_pm_gq40[i] <- sum(table(dn_gq40$parentValidation_gt))
  
  setwd(data_directory)
}

number_called_pf_gq40 <- rep(NA, length(all_dirs))
for (i in 1:length(all_dirs)){
  setwd(all_dirs[i])
  all_data <- list.files()
  load(file = all_data[grep('de_novo_variant_granges', all_data)])
  de_novos_pf <- de_novos_pf[!grepl("\\.", de_novos_pf$parentValidation_gt)]
  keep <- which(de_novos_pf$GQ_proband >= 40 & de_novos_pf$GQ_parent >= 40 & de_novos_pf$parentValidation_GQ >= 40 & 
                  de_novos_pf$parentValidation_depth >= 20)
  dn_gq40 <- de_novos_pf[keep]
  in_gnomad <- which(dn_gq40$gnomad41_genome_AF > 0)
  if(length(in_gnomad) > 0){
    dn_gq40 <- dn_gq40[-in_gnomad]
  }
  number_called_pf_gq40[i] <- sum(table(dn_gq40$parentValidation_gt))
  
  setwd(data_directory)
}


### plots
pdf(file = paste0(figure_directory, "/ppv_no_gnomad_GQ40.pdf"), height = 3.5, width = 5.9, pointsize = 8)
par(mfrow = c(2, 1))
precision_pf <- ppv_pf_gq40_no_gnomad
plot(sort(precision_pf, decreasing = TRUE), ylim = c(0, 1), pch = 19, col = "dark orange", xlab = "Duo index", 
     ylab = "PPV (PF)", main = "Variants absent from gnomAD", cex = 0.5, bty = 'l', yaxt = 'n')

precision_pm <- ppv_pm_gq40_no_gnomad
precision_pm <- precision_pm[-which(is.na(precision_pm))]
plot(sort(precision_pm, decreasing = TRUE), ylim = c(0, 1), pch = 19, col = "dark orange", xlab = "Duo index", 
     ylab = "PPV (PM)", main = "Variants absent from gnomAD", cex = 0.5, bty = 'l', yaxt = 'n')
dev.off()

pdf(file = paste0(figure_directory, "/ppv_n_called_pf_GQ40.pdf"), 
    height = 2, width = 5.9, pointsize = 8)

number_called <- number_called_pf_gq40
plot(number_called[order(precision_pf, decreasing = TRUE)], pch = 19, type = 'h', lwd = 1,
     cex = 0.5, bty = 'l', col = "dark orange", ylim = c(0, max(number_called)), 
     yaxt = 'n', xaxt = 'n', ylab = "# de novos called", xlab = "father-proband duos")
axis(2, at = c(0, 15, 30))
dev.off()

pdf(file = paste0(figure_directory, "/ppv_n_called_pm_GQ40.pdf"), 
    height = 2, width = 5.9, pointsize = 8)

number_called <- number_called_pm_gq40[which(number_called_pm_gq40 > 0)]
plot(number_called[order(precision_pm, decreasing = TRUE)], pch = 19, type = 'h', lwd = 1,
     cex = 0.5, bty = 'l', col = "dark orange", ylim = c(0, max(number_called)), 
     yaxt = 'n', xaxt = 'n', ylab = "# de novos called", xlab = "father-proband duos")
axis(2, at = c(0, 5, 10))
dev.off()
###


###
### Collective PPV
ppv_pf_gq40 <- matrix(NA, nrow = 2, ncol = length(all_dirs))
rownames(ppv_pf_gq40) <- c('dn', 'all')
for (i in 1:length(all_dirs)){
  setwd(all_dirs[i])
  all_data <- list.files()
  load(file = all_data[grep('de_novo_variant_granges', all_data)])
  de_novos_pf <- de_novos_pf[!grepl("\\.", de_novos_pf$parentValidation_gt)]
  keep <- which(de_novos_pf$GQ_proband >= 40 & de_novos_pf$GQ_parent >= 40 & de_novos_pf$parentValidation_GQ >= 40 & 
                  de_novos_pf$parentValidation_depth >= 20)
  dn_gq40 <- de_novos_pf[keep]
  ppv_pf_gq40['dn', i] <- table(dn_gq40$parentValidation_gt)["0/0"]
  ppv_pf_gq40['all', i] <- sum(table(dn_gq40$parentValidation_gt))
  
  setwd(data_directory)
}

ppv_pm_gq40 <- matrix(NA, nrow = 2, ncol = length(all_dirs))
rownames(ppv_pm_gq40) <- c('dn', 'all')
for (i in 1:length(all_dirs)){
  setwd(all_dirs[i])
  all_data <- list.files()
  load(file = all_data[grep('de_novo_variant_granges', all_data)])
  de_novos_pm <- de_novos_pm[!grepl("\\.", de_novos_pm$parentValidation_gt)]
  keep <- which(de_novos_pm$GQ_proband >= 40 & de_novos_pm$GQ_parent >= 40 & de_novos_pm$parentValidation_GQ >= 40 & 
                  de_novos_pm$parentValidation_depth >= 20)
  dn_gq40 <- de_novos_pm[keep]
  ppv_pm_gq40['dn', i] <- table(dn_gq40$parentValidation_gt)["0/0"]
  ppv_pm_gq40['all', i] <- sum(table(dn_gq40$parentValidation_gt))
  
  setwd(data_directory)
}

ppv_all_pf <- sum(ppv_pf[1, ], na.rm = TRUE)/sum(ppv_pf[2, ], na.rm = TRUE)
ppv_all_pm <- sum(ppv_pm[1, ], na.rm = TRUE)/sum(ppv_pm[2, ], na.rm = TRUE)
ppv_all_pf_gq40 <- sum(ppv_pf_gq40[1, ], na.rm = TRUE)/sum(ppv_pf_gq40[2, ], na.rm = TRUE)
ppv_all_pm_gq40 <- sum(ppv_pm_gq40[1, ], na.rm = TRUE)/sum(ppv_pm_gq40[2, ], na.rm = TRUE)

pdf(file = paste0(figure_directory, "/collective_ppv_GQ40.pdf"), height = 3.2, width = 2, pointsize = 8)
df <- data.frame(
  Condition = factor(c("PF", "PF (GQ > 40)", "PM", "PM (GQ > 40)"),
                     levels = c("PF", "PF (GQ > 40)", "PM", "PM (GQ > 40)")),
  PPV       = c(ppv_all_pf, ppv_all_pf_gq40, ppv_all_pm, ppv_all_pm_gq40)
)

ggplot(df, aes(x = Condition, y = PPV, fill = alpha("red", 0.75))) +
  geom_col() +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(
    title = "",
    y     = "PPV",
    x     = ""
  ) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x    = element_text(angle = 45, hjust = 1))
dev.off()
