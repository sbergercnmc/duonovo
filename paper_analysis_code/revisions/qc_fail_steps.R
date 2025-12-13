
setwd(qc_data_directory)
all_dirs <- list.files()
all_dirs_no_problematic <- all_dirs[-which(all_dirs %in% problematic_trios)]

qc_fail_steps_pf <- matrix(NA, nrow = 8, ncol = length(all_dirs_no_problematic))
qc_fail_steps_pm <- qc_fail_steps_pf
rownames(qc_fail_steps_pf) <- rownames(qc_fail_steps_pm) <- c("hap block overlap", "small hap block", "boundary", 
                                                              "depth", "GQ", "depth and GQ", 
                                                              "problematic_QC_fail", "problematic_classified")

for (i in 1:length(all_dirs_no_problematic)){
  setwd(all_dirs_no_problematic[i])
  all_data <- list.files()
  load(file = all_data[grep('hap_block_impact', all_data)])
  
  qc_fail_steps_pf[1:8, i] <- prop.table(hap_block_impact_pf[1:8])
  qc_fail_steps_pm[1:8, i] <- prop.table(hap_block_impact_pm[1:8])
  setwd(qc_data_directory)
}

###
### Plots
library(ggplot2)

pdf(paste0(figure_directory, "/QC_fail_steps.pdf"), width = 3.5, height = 2.75)
steps <- rownames(qc_fail_steps_pf)

df_pf <- data.frame(
  Step   = factor(rep(steps, times = ncol(qc_fail_steps_pf)), levels = steps),
  Parent = factor("PF", levels = c("PF", "PM")),
  Value  = as.vector(qc_fail_steps_pf)
)

df_pm <- data.frame(
  Step   = factor(rep(steps, times = ncol(qc_fail_steps_pm)), levels = steps),
  Parent = factor("PM", levels = c("PF", "PM")),
  Value  = as.vector(qc_fail_steps_pm)
)

df_long <- rbind(df_pf, df_pm)

ggplot(df_long, aes(x = Step, y = Value, fill = Parent)) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "gray30",
    fatten = 0.4, linewidth = 0) +
  scale_fill_manual(values = c("PF" = "lightpink", "PM" = "salmon")) +
  coord_cartesian(ylim = c(0, 0.4)) +
  theme_classic() +
  labs(x = "", y = "Value") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()

####
####




