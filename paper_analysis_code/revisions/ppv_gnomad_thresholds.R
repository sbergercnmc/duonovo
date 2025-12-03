dirs_to_exclude <- which(all_dirs %in% problematic_trios)
dirs_to_use <- all_dirs[-dirs_to_exclude]
setwd(data_directory)

getPPVMatrix <- function(gnomad_threshold, GQ_threshold, duo_type = c("PF", "PM")){
  ppv_mat <- matrix(NA, nrow = 2, ncol = length(dirs_to_use))
  rownames(ppv_mat) <- c('dn', 'all')
  for (i in 1:length(dirs_to_use)){
    setwd(dirs_to_use[i])
    all_data <- list.files()
    load(file = all_data[grep('de_novo_variant_granges', all_data)])
    if (duo_type == "PF"){
      de_novos <- de_novos_pf
    } else if (duo_type == "PM"){
      de_novos <- de_novos_pm
    }
    de_novos <- de_novos[!grepl("\\.", de_novos$parentValidation_gt)]
    keep <- which(de_novos$GQ_proband >= GQ_threshold & de_novos$GQ_parent >= GQ_threshold & 
                    de_novos$parentValidation_GQ >= GQ_threshold & 
                    de_novos$parentValidation_depth >= 20)
    dn_gq <- de_novos[keep]
    in_gnomad <- which(dn_gq$gnomad41_genome_AF > gnomad_threshold)
    if(length(in_gnomad) > 0){
      dn_gq <- dn_gq[-in_gnomad]
    }
    ppv_mat['dn', i] <- length(dn_gq) - length(grep("1", dn_gq$parentValidation_gt))
    ppv_mat['all', i] <- length(dn_gq)  
    setwd(data_directory)
  }
  ppv_mat
}

ppv_mat_list_pf <- lapply(c(0, 0.001, 0.01, 1), function(xx) {
  ppv_mat <- getPPVMatrix(xx, GQ_threshold = 30, duo_type = "PF")
  ppv_mat
})

ppv_mat_list_pm <- lapply(c(0, 0.001, 0.01, 1), function(xx) {
  ppv_mat <- getPPVMatrix(xx, GQ_threshold = 30, duo_type = "PM")
  ppv_mat
})

ppv_mat_list_pf_GQ40 <- lapply(c(0, 0.001, 0.01, 1), function(xx) {
  ppv_mat <- getPPVMatrix(xx, GQ_threshold = 40, duo_type = "PF")
  ppv_mat
})

ppv_mat_list_pm_GQ40 <- lapply(c(0, 0.001, 0.01, 1), function(xx) {
  ppv_mat <- getPPVMatrix(xx, GQ_threshold = 40, duo_type = "PM")
  ppv_mat
})


ppv_all_pf_0 <- sum(ppv_mat_list_pf[[1]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pf[[1]][2, ], na.rm = TRUE)
ppv_all_pf_001 <- sum(ppv_mat_list_pf[[2]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pf[[2]][2, ], na.rm = TRUE)
ppv_all_pf_01 <- sum(ppv_mat_list_pf[[3]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pf[[3]][2, ], na.rm = TRUE)
ppv_all_pf <- sum(ppv_mat_list_pf[[4]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pf[[4]][2, ], na.rm = TRUE)

ppv_all_pm_0 <- sum(ppv_mat_list_pm[[1]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pm[[1]][2, ], na.rm = TRUE)
ppv_all_pm_001 <- sum(ppv_mat_list_pm[[2]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pm[[2]][2, ], na.rm = TRUE)
ppv_all_pm_01 <- sum(ppv_mat_list_pm[[3]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pm[[3]][2, ], na.rm = TRUE)
ppv_all_pm <- sum(ppv_mat_list_pm[[4]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pm[[4]][2, ], na.rm = TRUE)

ppv_all_pf_0_GQ40 <- sum(ppv_mat_list_pf_GQ40[[1]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pf_GQ40[[1]][2, ], na.rm = TRUE)
ppv_all_pf_001_GQ40 <- sum(ppv_mat_list_pf_GQ40[[2]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pf_GQ40[[2]][2, ], na.rm = TRUE)
ppv_all_pf_01_GQ40 <- sum(ppv_mat_list_pf_GQ40[[3]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pf_GQ40[[3]][2, ], na.rm = TRUE)
ppv_all_pf_GQ40 <- sum(ppv_mat_list_pf_GQ40[[4]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pf_GQ40[[4]][2, ], na.rm = TRUE)

ppv_all_pm_0_GQ40 <- sum(ppv_mat_list_pm_GQ40[[1]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pm_GQ40[[1]][2, ], na.rm = TRUE)
ppv_all_pm_001_GQ40 <- sum(ppv_mat_list_pm_GQ40[[2]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pm_GQ40[[2]][2, ], na.rm = TRUE)
ppv_all_pm_01_GQ40 <- sum(ppv_mat_list_pm_GQ40[[3]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pm_GQ40[[3]][2, ], na.rm = TRUE)
ppv_all_pm_GQ40 <- sum(ppv_mat_list_pm_GQ40[[4]][1, ], na.rm = TRUE)/sum(ppv_mat_list_pm_GQ40[[4]][2, ], na.rm = TRUE)


### Plots
###
df <- data.frame(
  Parent = factor(
    rep(c("PF", "PM"), each = 4),
    levels = c("PF", "PM")
  ),
  Cutoff = factor(
    rep(c("0", "0.001", "0.01", "None"), times = 2),
    levels = c("0", "0.001", "0.01", "None")
  ),
  PPV = c(
    ppv_all_pf_0,  ppv_all_pf_001,  ppv_all_pf_01,  ppv_all_pf,
    ppv_all_pm_0,  ppv_all_pm_001,  ppv_all_pm_01,  ppv_all_pm
  )
)

df$Condition <- paste(
  df$Parent,
  df$Cutoff
)

pdf(file = paste0(figure_directory, "/collective_ppv_gnomad_thresholds.pdf"), 
    height = 2.4, width = 3, pointsize = 8)

ggplot(df, aes(x = Condition, y = PPV, fill = Parent)) +
  geom_col() +
  scale_fill_manual(values = c(
    "PF" = alpha("red", 0.75),
    "PM" = "cornflowerblue"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "",
    y = "PPV",
    title = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )
dev.off()

### GQ > 40
###
df <- data.frame(
  Parent = factor(
    rep(c("PF", "PM"), each = 4),
    levels = c("PF", "PM")
  ),
  Cutoff = factor(
    rep(c("0", "0.001", "0.01", "None"), times = 2),
    levels = c("0", "0.001", "0.01", "None")
  ),
  PPV = c(
    ppv_all_pf_0_GQ40,  ppv_all_pf_001_GQ40,  ppv_all_pf_01_GQ40,  ppv_all_pf_GQ40,
    ppv_all_pm_0_GQ40,  ppv_all_pm_001_GQ40,  ppv_all_pm_01_GQ40,  ppv_all_pm_GQ40
  )
)

df$Condition <- paste(
  df$Parent,
  df$Cutoff
)

pdf(file = paste0(figure_directory, "/collective_ppv_gnomad_thresholds_GQ.pdf"), 
    height = 2.4, width = 3, pointsize = 8)

ggplot(df, aes(x = Condition, y = PPV, fill = Parent)) +
  geom_col() +
  scale_fill_manual(values = c(
    "PF" = alpha("red", 0.75),
    "PM" = "cornflowerblue"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "",
    y = "PPV",
    title = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )
dev.off()




