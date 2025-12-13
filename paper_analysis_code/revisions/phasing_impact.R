setwd(qc_data_directory)
all_dirs <- list.files()
all_dirs_no_problematic <- all_dirs[-which(all_dirs %in% problematic_trios)]

hap_blocks_pf <- matrix(NA, nrow = 3, ncol = length(all_dirs_no_problematic))
hap_blocks_pm <- hap_blocks_pf
rownames(hap_blocks_pf) <- rownames(hap_blocks_pm) <- c("total_size", "filtered_size", "fraction_phased_variants")

for (i in 1:length(all_dirs_no_problematic)){
  setwd(all_dirs_no_problematic[i])
  all_data <- list.files()
  load(file = all_data[grep('hap_block_sizes', all_data)])
  
  hap_blocks_pf[1:2, i] <- hap_block_size_pf[1:2]
  hap_blocks_pm[1:2, i] <- hap_block_size_pm[1:2]
  
  hap_blocks_pf[3, i] <- hap_block_size_pf[3]/hap_block_size_pf[4]
  hap_blocks_pm[3, i] <- hap_block_size_pm[3]/hap_block_size_pm[4]
  setwd(qc_data_directory)
}

###
### Plots
library(ggplot2)
library(gridExtra)
pdf(paste0(figure_directory, "/phasing_QC.pdf"), width = 3.5, height = 2.75)
size_rows <- c("total_size", "filtered_size")

df_size_pf <- data.frame(
  Metric = factor(rep(size_rows, times = ncol(hap_blocks_pf)), levels = size_rows),
  Parent = factor("PF", levels = c("PF", "PM")),
  Value  = as.vector(hap_blocks_pf[size_rows, , drop = FALSE])
)

df_size_pm <- data.frame(
  Metric = factor(rep(size_rows, times = ncol(hap_blocks_pm)), levels = size_rows),
  Parent = factor("PM", levels = c("PF", "PM")),
  Value  = as.vector(hap_blocks_pm[size_rows, , drop = FALSE])
)

df_size <- rbind(df_size_pf, df_size_pm)

p_size <- ggplot(df_size, aes(x = Metric, y = Value, fill = Parent)) +
  geom_boxplot(outlier.shape = NA,
               position = position_dodge(width = 0.8),
               width = 0.7,
               color = "gray30") +
  scale_fill_manual(values = c("PF" = "lightpink", "PM" = "salmon")) +
  scale_y_continuous() +
  theme_classic() +
  labs(x = "", y = "Size") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

## 2) Boxplots: fraction_phased_variants (PF vs PM side-by-side)
frac_row <- "fraction_phased_variants"

df_frac_pf <- data.frame(
  Metric = factor(rep(frac_row, times = ncol(hap_blocks_pf)), levels = frac_row),
  Parent = factor("PF", levels = c("PF", "PM")),
  Value  = as.vector(hap_blocks_pf[frac_row, , drop = FALSE])
)

df_frac_pm <- data.frame(
  Metric = factor(rep(frac_row, times = ncol(hap_blocks_pm)), levels = frac_row),
  Parent = factor("PM", levels = c("PF", "PM")),
  Value  = as.vector(hap_blocks_pm[frac_row, , drop = FALSE])
)

df_frac <- rbind(df_frac_pf, df_frac_pm)

p_frac <- ggplot(df_frac, aes(x = Metric, y = Value, fill = Parent)) +
  geom_boxplot(outlier.shape = NA,
               position = position_dodge(width = 0.8),
               width = 0.7,
               color = "gray30",
               fatten = 0.4) +
  scale_fill_manual(values = c("PF" = "lightpink", "PM" = "salmon")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  labs(x = "", y = "Fraction phased variants") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

grid.arrange(p_size, p_frac, nrow = 1)
dev.off()
