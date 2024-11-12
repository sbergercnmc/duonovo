###
pdf(file = "duo_novo_performance_false_positive_sources.pdf", height = 4.4, width = 6.5, pointsize = 8)
par(mfrow = c(2, 3))
plot(1- duo_novo_performance_pos_f["precision", ], 1 - duo_novo_performance_pos_f2["precision", ], 
     xlab = "false discovery rate (all)", 
     ylab = "false discovery rate (within genes only)", col = "dark orange", bty = 'l', 
     xaxt = 'n', yaxt = 'n', pch = 19, xlim = c(0, 0.75), ylim = c(0, 0.75))
axis(1, at = c(0, 0.25, 0.5, 0.75))
axis(2, at = c(0, 0.25, 0.5, 0.75))
abline(0, 1, col = rgb(0,0,0,0.5))

plot(1- duo_novo_performance_pos_f["precision", ], 1 - duo_novo_performance_pos_f5["precision", ], 
     xlab = "false discovery rate (excluding mult PS)", 
     ylab = "false discovery rate (including mult PS)", col = "dark orange", bty = 'l', 
     xaxt = 'n', yaxt = 'n', pch = 19, xlim = c(0, 0.75), ylim = c(0, 0.75))
axis(1, at = c(0, 0.25, 0.5, 0.75))
axis(2, at = c(0, 0.25, 0.5, 0.75))
abline(0, 1, col = rgb(0,0,0,0.5))

plot(1- duo_novo_performance_pos_f["precision", ], 1 - duo_novo_performance_pos_f6["precision", ], 
     xlab = "false discovery rate (excluding boundary)", 
     ylab = "false discovery rate (including boundary)", col = "dark orange", bty = 'l', 
     xaxt = 'n', yaxt = 'n', pch = 19, xlim = c(0, 0.75), ylim = c(0, 0.75))
axis(1, at = c(0, 0.25, 0.5, 0.75))
axis(2, at = c(0, 0.25, 0.5, 0.75))
abline(0, 1, col = rgb(0,0,0,0.5))

###
plot(1- duo_novo_performance_pos_m["precision", ], 1 - duo_novo_performance_pos_m2["precision", ], 
     xlab = "false discovery rate (all)", 
     ylab = "false discovery rate (within genes only)", col = "dark orange", bty = 'l', 
     xaxt = 'n', yaxt = 'n', pch = 19, xlim = c(0, 0.75), ylim = c(0, 0.75))
axis(1, at = c(0, 0.25, 0.5, 0.75))
axis(2, at = c(0, 0.25, 0.5, 0.75))
abline(0, 1, col = rgb(0,0,0,0.5))

plot(1- duo_novo_performance_pos_m["precision", ], 1 - duo_novo_performance_pos_m5["precision", ], 
     xlab = "false discovery rate (excluding mult PS)", 
     ylab = "false discovery rate (including mult PS)", col = "dark orange", bty = 'l', 
     xaxt = 'n', yaxt = 'n', pch = 19, xlim = c(0, 0.8), ylim = c(0, 0.8))
axis(1, at = c(0, 0.25, 0.5, 0.75))
axis(2, at = c(0, 0.25, 0.5, 0.75))
abline(0, 1, col = rgb(0,0,0,0.5))

plot(1- duo_novo_performance_pos_m["precision", ], 1 - duo_novo_performance_pos_m6["precision", ], 
     xlab = "false discovery rate (excluding boundary)", 
     ylab = "false discovery rate (including boudnary)", col = "dark orange", bty = 'l', 
     xaxt = 'n', yaxt = 'n', pch = 19, xlim = c(0, 0.8), ylim = c(0, 0.8))
axis(1, at = c(0, 0.25, 0.5, 0.75))
axis(2, at = c(0, 0.25, 0.5, 0.75))
abline(0, 1, col = rgb(0,0,0,0.5))
dev.off()

######
plot(1- duo_novo_performance_pos_m["precision", ], 1 - duo_novo_performance_pos_m_hom["precision", ], 
     xlab = "false discovery rate (including het)", 
     ylab = "false discovery rate (hom only)", col = "dark orange", bty = 'l', 
     xaxt = 'n', yaxt = 'n', pch = 19, xlim = c(0, 1), ylim = c(0, 1))
axis(1, at = c(0, 0.25, 0.5, 0.75))
axis(2, at = c(0, 0.25, 0.5, 0.75))
abline(0, 1, col = rgb(0,0,0,0.5))

plot(1- duo_novo_performance_pos_f["precision", ], 1 - duo_novo_performance_pos_f_hom["precision", ], 
     xlab = "false discovery rate (including het)", 
     ylab = "false discovery rate (hom only)", col = "dark orange", bty = 'l', 
     xaxt = 'n', yaxt = 'n', pch = 19, xlim = c(0, 1), ylim = c(0, 1))
axis(1, at = c(0, 0.25, 0.5, 0.75))
axis(2, at = c(0, 0.25, 0.5, 0.75))
abline(0, 1, col = rgb(0,0,0,0.5))
######


