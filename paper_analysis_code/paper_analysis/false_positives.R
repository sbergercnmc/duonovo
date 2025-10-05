args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript false_positives.R <data_directory> <figure_directory>")
}

data_directory   <- args[1]
figure_directory <- args[2]

# ensure figure directory exists
if (!dir.exists(data_directory)) {
  stop(sprintf("Data directory not found: '%s'", data_directory))
}

setwd(data_directory)
all_dirs <- list.files()

### false pos
false_pos_mat_pf <- matrix(NA, ncol = length(all_dirs), nrow = 3)
false_pos_mat_pm <- false_pos_mat_pf
rownames(false_pos_mat_pf) <- rownames(false_pos_mat_pm) <- c('false_pos', "true_neg", "uncertain")

for (i in 1:length(all_dirs)){
  setwd(all_dirs[i])
  all_data <- list.files()
  load(file = all_data[grep('false_pos', all_data)])
  
  false_pos_mat_pf[1, i] <- false_pos_pf['dn']
  false_pos_mat_pm[1, i] <- false_pos_pm['dn']
  
  false_pos_mat_pf[2, i] <- false_pos_pf['ndn']
  false_pos_mat_pm[2, i] <- false_pos_pm['ndn']
  
  false_pos_mat_pf[3, i] <- false_pos_pf['uncertain']
  false_pos_mat_pm[3, i] <- false_pos_pm['uncertain']
  
  setwd('~/Desktop/duoNovo_results_R_objects/data')
}

aggregate_transmitted_f <- rowSums(false_pos_mat_pf)
percentages_transmitted_f <- aggregate_transmitted_f/sum(aggregate_transmitted_f)
per_sample_false_pos_rate_f <- apply(false_pos_mat_pf, 2, function(xx) xx[1]/sum(xx))

aggregate_transmitted_m <- rowSums(false_pos_mat_pm)
percentages_transmitted_m <- aggregate_transmitted_m/sum(aggregate_transmitted_m)
per_sample_false_pos_rate_m <- apply(false_pos_mat_pm, 2, function(xx) xx[1]/sum(xx))

pdf(file = paste0(figure_directory, "/pf_transmitted_aggregate.pdf"), height = 1.8, width = 2.4, pointsize = 8)
data <- data.frame(
  category = names(percentages_transmitted_f),
  value = percentages_transmitted_f
)
data$category <- fct_relevel(data$category, "uncertain", "dn", "ndn")

ggplot(data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c(alpha("dark gray", 0.75), "cornflowerblue", rgb(1,0,0,0.62))) +
  labs(title = "",
       x = "",
       y = "Percentage",
       fill = "") +
  theme_classic() +
  theme(legend.position = "bottom")
dev.off()

pdf(file = paste0(figure_directory, "/pm_transmitted_aggregate.pdf"), height = 1.8, width = 2.4, pointsize = 8)
data <- data.frame(
  category = names(percentages_transmitted_m),
  value = percentages_transmitted_m
)
data$category <- fct_relevel(data$category, "uncertain", "dn", "ndn")

ggplot(data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c(alpha("dark gray", 0.75), "cornflowerblue", rgb(1,0,0,0.62))) +
  labs(title = "",
       x = "",
       y = "Percentage",
       fill = "") +
  theme_classic() +
  theme(legend.position = "bottom")
dev.off()


### Individual duos
###
pdf(file = paste0(figure_directory, "/all_duos_false_positives.pdf"), 
    height = 2, width = 5.9, pointsize = 8)
par(mar = c(3, 4, 3, 0.1) + 0.1)
x_base <- seq(1, 350, length.out = 117)

offset <- 0.4
x_red  <- x_base - offset
x_blue <- x_base + offset

x_sep <- (x_base[-1] + x_base[-length(x_base)]) / 2

plot(x_red, 100 * per_sample_false_pos_rate_f,
     ylim = c(0, 0.005), pch = 19, col = rgb(1, 0, 0, 0.57),
     xlab = "trio index", cex = 0.5,
     ylab = "% false positive duo novo classifications",
     bty = 'l', xaxt = 'n', yaxt = 'n',
     xaxs = "i",
     xlim = range(x_base) + c(-2, 2))

points(x_blue, 100 * per_sample_false_pos_rate_m,
       pch = 19, col = "cornflowerblue", cex = 0.5)

axis(2, at = c(0, 0.005))
abline(v = x_sep, lty = "longdash", col = rgb(0, 0, 0, 0.35))
legend <- legend("topleft", legend = c("father-proband duos", "mother-proband duos"), 
                 pch = 20, bty = 'n', cex = 0.62, 
                 col = c(rgb(1,0,0,0.57), "cornflowerblue"))
dev.off()

