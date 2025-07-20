args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript ref_vs_alt_dn_cal_rate.R <data_directory> <figure_directory>")
}

data_directory   <- args[1]
figure_directory <- args[2]

# ensure figure directory exists
if (!dir.exists(data_directory)) {
  stop(sprintf("Data directory not found: '%s'", data_directory))
}

setwd(data_directory)
all_dirs <- list.files()

### compute REF vs ALT candidate allele de novo call rate
dn_call_rate <- matrix(NA, ncol = length(all_dirs), nrow = 4)
rownames(dn_call_rate) <- c('pf_alt', 'pf_ref', 'pm_alt', 'pm_ref')

for (i in 1:length(all_dirs)){
  setwd(all_dirs[i])
  all_data <- list.files()
  load(file = all_data[grep('call_rate', all_data)])
  
  dn_call_rate['pf_alt', i] <- dn_call_rate_alt_pf[3]
  dn_call_rate['pf_ref', i] <- dn_call_rate_ref_pf[3]
  dn_call_rate['pm_alt', i] <- dn_call_rate_alt_pm[3]
  dn_call_rate['pm_ref', i] <- dn_call_rate_ref_pm[3]
  
  setwd(data_directory)
}

pdf(file = paste0(figure_directory, "/duo_novo_performance_ref_vs_alt.pdf"), 
    height = 2.4, width = 7, pointsize = 8)
par(mfrow = c(1, 2))
plot(1:117, 100*dn_call_rate['pf_alt', ], 
     pch = 19, col = "salmon", xlab = "father-proband duos", cex = 0.8,
     ylab = "% classified as de novo", 
     bty = 'l', xaxt = 'n', yaxt = 'n', xlim = c(0.8, 117.2), ylim = c(0, 0.012))
points(1:117, 100*dn_call_rate['pf_ref', ], 
       pch = 19, col = "gray50", cex = 0.8)
axis(2, at = c(0, 0.01))
abline(v = seq(1.5, 116.5, by = 1), lty = "longdash", col = rgb(0,0,0,0.4))
legend <- legend("topleft", legend = c("ALT allele", "REF allele"), 
                 pch = 20, bty = 'n', cex = 0.62, 
                 col = c('salmon', "gray50"))

plot(1:117, 100*dn_call_rate['pm_alt', ], 
     pch = 19, col = "salmon", xlab = "father-proband duos", cex = 0.8,
     ylab = "% classified as de novo", 
     bty = 'l', xaxt = 'n', yaxt = 'n', xlim = c(0.8, 117.2), ylim = c(0, 0.0046))
points(1:117, 100*dn_call_rate['pm_ref', ], 
       pch = 19, col = "gray50", cex = 0.8)
axis(2, at = c(0, 0.004))
abline(v = seq(1.5, 116.5, by = 1), lty = "longdash", col = rgb(0,0,0,0.4))
dev.off()
