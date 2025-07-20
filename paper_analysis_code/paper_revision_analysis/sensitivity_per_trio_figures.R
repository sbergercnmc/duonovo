args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript sensitivity_per_trio_figures.R <data_directory> <figure_directory>")
}

data_directory   <- args[1]
figure_directory <- args[2]

# ensure figure directory exists
if (!dir.exists(data_directory)) {
  stop(sprintf("Data directory not found: '%s'", data_directory))
}

setwd(data_directory)
all_dirs <- list.files()

trio_dn_count <- sapply(1:length(all_dirs), function(xx) {
  setwd(all_dirs[xx])
  all_data <- list.files()
  load(file = all_data[grep("sensitivity", all_data)])
  setwd(data_directory)
  
  sens_both['total']
})

problematic_trios <- which(trio_dn_count >= 150)
all_dirs_no_problematic <- all_dirs[-problematic_trios]

### Sensitivity
sens <- rep(NA, length(all_dirs_no_problematic))
sens_no_filter <- sens

for (i in 1:length(all_dirs_no_problematic)){
  setwd(all_dirs_no_problematic[i])
  all_data <- list.files()
  load(file = all_data[grep('sensitivity', all_data)])
  
  sens[i] <- sens_both[1]/sens_both[4]
  sens_no_filter[i] <- sens_both_no_ad_filter[1]/sens_both_no_ad_filter[4]
  
  setwd(data_directory)
}

pdf(file = paste0(figure_directory, "/sensitivity.pdf"), height = 2.4, width = 6.5, pointsize = 8)
par(mfrow = c(1, 3))
hist(trio_dn_count, breaks = 50, xlim = c(0, 700), col = "red", border = "white", 
     xlab = "# de novos from trio", main = "")
plot(trio_dn_count, sens, xlab = "# de novos from full trio", ylab = "sensitivity (both duos combined)", 
     bty = 'l', yaxt = 'n', xaxt = 'n', ylim = c(0, 0.75))
axis(1, at = seq(100, 600, by = 100))
axis(2, at = c(0.25, 0.5, 0.75))

plot(sens_no_filter, sens, xlab = "sensitivity (w/out allele depth filter)", 
     ylab = "sensitivity (with allele depth filter)", 
     bty = 'l', yaxt = 'n', xaxt = 'n', ylim = c(0, 0.75), xlim = c(0, 0.75))
axis(1, at = c(0.25, 0.5, 0.75))
axis(2, at = c(0.25, 0.5, 0.75))
abline(0, 1, lwd = 1.75, col = 2)
dev.off()
