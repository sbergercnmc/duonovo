


# total count of de novo variants for all trios including problematic ones
trio_dn_count <- sapply(1:length(problematic_trios), function(xx) {
  setwd(problematic_trios[xx])
  all_data <- list.files()
  load(file = all_data[grep("sensitivity", all_data)])
  setwd(data_directory)
  
  sens_both['total']
})

### compute number of de novo classifications from duoNovo
n_denovo_pf <- rep(NA, length(problematic_trios))
n_denovo_pm <- n_denovo_pf

for (i in 1:length(problematic_trios)){
  setwd(problematic_trios[i])
  all_data <- list.files()
  load(file = all_data[grep('ppv', all_data)])
  
  n_denovo_pf[i] <- ppv_pf_alt[3]
  n_denovo_pm[i] <- ppv_pm_alt[3]

  setwd(data_directory)
}

pdf(file = paste0(figure_directory, "/duonovo_problematic_trios.pdf"), height = 2.2, width = 2.2, pointsize = 8)
plot(n_denovo_pf + n_denovo_pm, trio_dn_count, xlim = c(0, 650), ylim = c(0, 650), col = "darkorange", 
     cex = 0.5, bty = 'l', xaxt = 'n', yaxt = 'n', 
     xlab = "# de novo variants (duoNovo)", ylab = "# de novo variants (full trio genotyping)")
abline(0, 1, col = rgb(0,0,0,0.5), lty = "longdash")
axis(1, at = c(0, 300, 600))
axis(2, at = c(0, 300, 600))
dev.off()

