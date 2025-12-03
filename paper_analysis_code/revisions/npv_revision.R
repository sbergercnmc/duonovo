colnames(npv_pf_mat) <- colnames(npv_pm_mat) <- all_dirs
problematic_trios <- c('25-224964', '25-230606', '24-467599', '25-229554', 
                       '24-441582', '24-441816', '24-441599', '24-441864', 
                       'UCI-008', 'UCI-031', '25-224985', '25-224968', '25-224942')
indices_to_exclude <- which(colnames(npv_pf_mat) %in% problematic_trios)
npv_pf_mat <- npv_pf_mat[, -indices_to_exclude]
npv_pm_mat <- npv_pm_mat[, -indices_to_exclude]

pdf(file = paste0(figure_directory, "/father_proband_npv_scatterplot.pdf"), 
    height = 2.2, width = 2.2, pointsize = 8)
npv <- npv_pf_mat['duoNovo', ]
baseline_npv <- npv_pf_mat['naive', ]
plot(baseline_npv, npv, ylim = c(0.9995, 1), xlim = c(0.9995, 1), bty = 'l', xaxt = 'n', yaxt = 'n', cex = 0.75,
     xlab = "negative predictive value (naive)", ylab = "negative predictive value (duoNovo)", col = "dark orange")
abline(0, 1, lwd = 0.75, col = rgb(0,0,0,0.7))
axis(1, at = c(0.9995, 1))
axis(2, at = c(0.9995, 1))
dev.off()

pdf(file = paste0(figure_directory, "/mother_proband_npv_scatterplot.pdf"), 
    height = 2.2, width = 2.2, pointsize = 8)
npv <- npv_pm_mat['duoNovo', ]
baseline_npv <- npv_pm_mat['naive', ]
plot(baseline_npv, npv, ylim = c(0.9995, 1), xlim = c(0.9995, 1), bty = 'l', xaxt = 'n', yaxt = 'n', cex = 0.75,
     xlab = "negative predictive value (naive)", ylab = "negative predictive value (duoNovo)", col = "dark orange")
abline(0, 1, lwd = 0.75, col = rgb(0,0,0,0.7))
axis(1, at = c(0.9995, 1))
axis(2, at = c(0.9995, 1))
dev.off()


####
####
pdf(file = paste0(figure_directory, "/npv_revision.pdf"), 
    height = 5.5, width = 6, pointsize = 8)
par(mfrow = c(2, 1))
npv <- npv_pf_mat['duoNovo', ]
baseline_npv <- npv_pf_mat['naive', ]
plot(npv, pch = 19, 
     cex = 0.75, bty = 'l', col = "dark orange", ylim = c(0.999, 1), xlab = "father-proband duos", 
     ylab = "negative predictive value", 
     xaxt = 'n', yaxt = 'n', xlim = c(0.8, 104))
points(baseline_npv, cex = 0.6, bty = 'l', col = "cornflowerblue", pch = 19)
axis(2, at = c(0.999, 1))
abline(v = seq(1.5, 103.5, by = 1), col = rgb(0,0,0,0.4), lty = "longdash", lwd = 0.5)
legend <- legend("bottomright", legend = c("duoNovo", "naive baseline"), 
                 pch = 20, bty = 'n', cex = 0.62, 
                 col = c("orange", "cornflowerblue"))

npv <- npv_pm_mat['duoNovo', ]
baseline_npv <- npv_pm_mat['naive', ]
plot(npv, pch = 19, 
     cex = 0.75, bty = 'l', col = "dark orange", ylim = c(0.999, 1), xlab = "mother-proband duos", 
     ylab = "negative predictive value", 
     xaxt = 'n', yaxt = 'n', xlim = c(0.8, 104))
points(baseline_npv, cex = 0.6, bty = 'l', col = "cornflowerblue", pch = 19)
axis(2, at = c(0.999, 1))
abline(v = seq(1.5, 103.5, by = 1), col = rgb(0,0,0,0.4), lty = "longdash", lwd = 0.5)

dev.off()




