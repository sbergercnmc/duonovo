args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript ppv_npv.R <data_directory> <figure_directory> <trio_info_directory>")
}
data_directory   <- args[1]
figure_directory <- args[2]
trio_info_directory <- args[3]

# ensure data directory exists
if (!dir.exists(data_directory)) {
  stop(sprintf("Data directory not found: '%s'", data_directory))
}

setwd(data_directory)

## import spreadsheet with proband ancestry info
library(readr)
trios <- read_csv(trio_info_directory) 
euro_indices <- which(trios$child_ANCESTRY == "EUR")
eur_ids <- trios$child_SAMPLEID[euro_indices]
non_eur_ids <- trios$child_SAMPLEID[-euro_indices]
##

dirs_eur <- eur_ids
dirs_non_eur <- non_eur_ids

dirs_eur <- dirs_eur[-which(dirs_eur %in% problematic_trios)]
dirs_non_eur <- dirs_non_eur[-which(dirs_non_eur %in% problematic_trios)]

ppv_pf_eur <- matrix(NA, ncol = length(dirs_eur), nrow = 3)
ppv_pm_eur <- ppv_pf_eur
ppv_pf_non_eur <- matrix(NA, ncol = length(dirs_non_eur), nrow = 3)
ppv_pm_non_eur <- ppv_pf_non_eur
rownames(ppv_pf_eur) <- rownames(ppv_pm_eur) <- rownames(ppv_pf_non_eur) <- rownames(ppv_pm_non_eur) <- c('dn', 'assessed', 'total')

for (i in 1:length(dirs_eur)){
  setwd(dirs_eur[i])
  all_data <- list.files()
  load(file = all_data[grep('ppv', all_data)])
  
  ppv_pf_eur['dn', i] <- ppv_pf_alt[1]
  ppv_pf_eur['assessed', i] <- ppv_pf_alt[2]
  ppv_pf_eur['total', i] <- ppv_pf_alt[3]
  
  ppv_pm_eur['dn', i] <- ppv_pm_alt[1]
  ppv_pm_eur['assessed', i] <- ppv_pm_alt[2]
  ppv_pm_eur['total', i] <- ppv_pm_alt[3]
  
  setwd(data_directory)
}

for (i in 1:length(dirs_non_eur)){
  setwd(dirs_non_eur[i])
  all_data <- list.files()
  load(file = all_data[grep('ppv', all_data)])
  
  ppv_pf_non_eur['dn', i] <- ppv_pf_alt[1]
  ppv_pf_non_eur['assessed', i] <- ppv_pf_alt[2]
  ppv_pf_non_eur['total', i] <- ppv_pf_alt[3]
  
  ppv_pm_non_eur['dn', i] <- ppv_pm_alt[1]
  ppv_pm_non_eur['assessed', i] <- ppv_pm_alt[2]
  ppv_pm_non_eur['total', i] <- ppv_pm_alt[3]
  
  setwd(data_directory)
}


# collective PPV
ppv_all_pf_eur <- sum(ppv_pf_eur[1, ], na.rm = TRUE)/sum(ppv_pf_eur[2, ], na.rm = TRUE)
ppv_all_pm_eur <- sum(ppv_pm_eur[1, ], na.rm = TRUE)/sum(ppv_pm_eur[2, ], na.rm = TRUE)
ppv_all_pf_non_eur <- sum(ppv_pf_non_eur[1, ], na.rm = TRUE)/sum(ppv_pf_non_eur[2, ], na.rm = TRUE)
ppv_all_pm_non_eur <- sum(ppv_pm_non_eur[1, ], na.rm = TRUE)/sum(ppv_pm_non_eur[2, ], na.rm = TRUE)


pdf(file = paste0(figure_directory, "/collective_ppv_by_ancestry_revision.pdf"), height = 3.2, width = 2, pointsize = 8)
df <- data.frame(
  Condition = factor(c("PF EUR", "PM EUR", "PF non-EUR", "PM non-EUR"),
                     levels = c("PF EUR", "PM EUR", "PF non-EUR", "PM non-EUR")),
  PPV       = c(ppv_all_pf_eur, ppv_all_pm_eur, ppv_all_pf_non_eur, ppv_all_pm_non_eur)
)

ggplot(df, aes(x = Condition, y = PPV)) +
  geom_col(fill = "darkorange") +
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

