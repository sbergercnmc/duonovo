load(file = "duo_novo_performance_results.rda")

total_f <- duo_novo_performance_pos_f["total", ]
total_m <- duo_novo_performance_pos_m["total", ]


pdf(file = "duo_novo_ratio_father_mother_dnm.pdf", height = 2.2, width = 2.2, pointsize = 8)
#par(mar = c(4, 4, 1, 1) + 0.1)
plot(total_f/total_m,  pch = 19, 
     cex = 1, bty = 'l', col = "dark orange", xlim = c(0.8, 16.2), ylim = c(0, 5), xaxt = "n", 
     yaxt = "n", ylab = "ratio of de novo from father-proband to mother-proband", 
     xlab = "Trio index")
axis(2, at = c(0, 2, 4))
axis(1, at = c(1:16), cex.axis = 0.5)
dev.off()

abline(h = mean(total_f/total_m), lty = "longdash", col = rgb(0,0,0,0.7))


library(readr)
trio_info <- read_csv('C:\\Users\\lboukas\\trio_info.csv')
all.equal(trio_info$PMGRC_ID_proband, colnames(duo_novo_performance_pos_f))

pdf(file = "duo_novo_dnm_vs_age.pdf", height = 2.2, width = 2.2, pointsize = 8)
#par(mar = c(4, 4, 1, 1) + 0.1)
plot(trio_info$father_age, total_f, pch = 19, 
     cex = 1.05, bty = 'l', col = "dark orange",
     xlab = "parent age at birth", xaxt = 'n', yaxt = 'n',
     ylab = "# de novo classifications", xlim = c(0, 55), ylim = c(0, 55))
axis(1, at = c(0, 25, 50))
axis(2, at = c(0, 25, 50))


points(trio_info$mother_age, total_m, pch = 2, 
     cex = 1.05, bty = 'l', col = "deep pink")
legend("topleft", legend = c("father-proband", "mother-proband"), pch = c(19, 2), 
       bty = 'n', col = c("dark orange", "deep pink"))

#fit poisson regression models
fit_f <- glm(total_f ~ father_age, data = trio_info, family = poisson(link = "log"))
fit_m <- glm(total_m ~ mother_age, data = trio_info, family = poisson(link = "log"))

#add a smooth trend line
age_range <- seq(15, 55, by = 1)
pred_f <- predict(fit_f, newdata = data.frame(father_age = age_range), type = "response")
pred_m <- predict(fit_m, newdata = data.frame(mother_age = age_range), type = "response")

lines(age_range, pred_f, col = "dark orange", lwd = 2)
lines(age_range, pred_m, col = "deep pink", lwd = 2)
dev.off()




