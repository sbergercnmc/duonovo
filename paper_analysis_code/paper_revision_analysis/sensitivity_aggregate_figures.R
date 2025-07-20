args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript sensitivity_aggregate_figures.R <data_directory> <figure_directory>")
}

data_directory   <- args[1]
figure_directory <- args[2]

# ensure figure directory exists
if (!dir.exists(data_directory)) {
  stop(sprintf("Data directory not found: '%s'", data_directory))
}

setwd(data_directory)
all_dirs <- list.files()

problematic_trios <- c('25-224964', '25-230606', '24-467599', '25-229554', 
                       '24-441582', '24-441816', '24-441599', '24-441864', 
                       'UCI-008', 'UCI-031', '25-224985', '25-224968', '25-224942')
all_dirs_no_problematic <- all_dirs[-problematic_trios]

sensitivity_mat_pf <- matrix(NA, ncol = length(all_dirs_no_problematic), nrow = 3)
rownames(sensitivity_mat_pf) <- c('dn', 'other_parent', 'uncertain')
sensitivity_mat_pm <- sensitivity_mat_pf

for (i in 1:length(all_dirs_no_problematic)){
  setwd(all_dirs_no_problematic[i])
  all_data <- list.files()
  load(file = all_data[grep('sensitivity', all_data)])
  
  sensitivity_mat_pf['dn', i] <- sens_pf['de_novo']
  sensitivity_mat_pf['other_parent', i] <- sens_pf['on_other_parent_haplotype']
  sensitivity_mat_pf['uncertain', i] <- sens_pf['uncertain']
  
  sensitivity_mat_pm['dn', i] <- sens_pm['de_novo']
  sensitivity_mat_pm['other_parent', i] <- sens_pm['on_other_parent_haplotype']
  sensitivity_mat_pm['uncertain', i] <- sens_pm['uncertain']
  
  setwd(data_directory)
}

aggregate_pf <- rowSums(sensitivity_mat_pf)
aggregate_pm <- rowSums(sensitivity_mat_pm)

candidate_percentages_f <- aggregate_pf[1:3]/sum(aggregate_pf)
candidate_percentages_m <- aggregate_pm[1:3]/sum(aggregate_pm)

pdf(file = paste0(figure_directory, "/pf_aggregate_sensitivity.pdf"), 
    height = 1.8, width = 2.4, pointsize = 8)
data <- data.frame(
  category = names(candidate_percentages_f),
  value = candidate_percentages_f
)
data$category <- fct_relevel(data$category, "uncertain", "other_parent", "dn")

ggplot(data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c(alpha("dark gray", 0.75), alpha("forest green", 0.84), "dark orange")) +
  labs(title = "",
       x = "",
       y = "Percentage",
       fill = "") +
  theme_classic() +
  theme(legend.position = "bottom")
dev.off()

pdf(file = paste0(figure_directory, "/pm_aggregate_sensitivity.pdf"), 
    height = 1.8, width = 2.4, pointsize = 8)
data <- data.frame(
  category = names(candidate_percentages_m),
  value = candidate_percentages_m
)
data$category <- fct_relevel(data$category, "uncertain", "other_parent", "dn")

ggplot(data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c(alpha("dark gray", 0.75), alpha("forest green", 0.84), "dark orange")) +
  labs(title = "",
       x = "",
       y = "Percentage",
       fill = "") +
  theme_classic() +
  theme(legend.position = "bottom")
dev.off()

