colnames(ppv_pf) <- colnames(ppv_pf_no_gnomad) <- colnames(ppv_pm) <- colnames(ppv_pm_no_gnomad) <- all_dirs
indices_to_exclude <- which(colnames(ppv_pf_no_gnomad) %in% problematic_trios)
ppv_pf_no_gnomad <- ppv_pf_no_gnomad[, -indices_to_exclude]
ppv_pm_no_gnomad <- ppv_pm_no_gnomad[, -indices_to_exclude]
ppv_pf <- ppv_pf[, -indices_to_exclude]
ppv_pm <- ppv_pm[, -indices_to_exclude]

### PPV plots 
pdf(file = paste0(figure_directory, "/ppv_no_gnomad_revision.pdf"), height = 3.5, width = 5.9, pointsize = 8)
par(mfrow = c(2, 1))
precision_pf <- ppv_pf_no_gnomad['dn', ]/ppv_pf_no_gnomad['assessed', ]
plot(sort(precision_pf, decreasing = TRUE), ylim = c(0, 1), pch = 19, col = "dark orange", xlab = "Duo index", 
     ylab = "PPV (PF)", main = "Variants absent from gnomAD", cex = 0.5, bty = 'l', yaxt = 'n')
axis(2, at = c(0, 1))

precision_pm <- ppv_pm_no_gnomad['dn', ]/ppv_pm_no_gnomad['assessed', ]
precision_pm <- precision_pm[-which(is.na(precision_pm))]
plot(sort(precision_pm, decreasing = TRUE), ylim = c(0, 1), pch = 19, col = "dark orange", xlab = "Duo index", 
     ylab = "PPV (PM)", main = "Variants absent from gnomAD", cex = 0.5, bty = 'l', yaxt = 'n')
axis(2, at = c(0, 1))
dev.off()

pdf(file = paste0(figure_directory, "/ppv_n_called_revision.pdf"), 
    height = 3.5, width = 5.9, pointsize = 8)
par(mfrow = c(2, 1))
number_called <- ppv_pf_no_gnomad['assessed', ]

plot(number_called[order(precision_pf, decreasing = TRUE)], pch = 19, type = 'h', lwd = 0.9,
     cex = 0.5, bty = 'l', col = "dark orange", ylim = c(0, max(number_called)), 
     yaxt = 'n', xaxt = 'n', ylab = "# de novos called", xlab = "father-proband duos")
axis(2, at = c(0, 15, 30))


number_called <- ppv_pm_no_gnomad['assessed', ]
number_called <- number_called[-which(is.na(number_called))]
plot(number_called[order(precision_pm, decreasing = TRUE)], pch = 19, type = 'h', lwd = 0.9,
     cex = 0.5, bty = 'l', col = "dark orange", ylim = c(0, max(number_called)), 
     yaxt = 'n', xaxt = 'n', ylab = "# de novos called", xlab = "mother-proband duos")
axis(2, at = c(0, 5, 10))
dev.off()

