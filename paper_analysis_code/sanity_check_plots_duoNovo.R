load(file = "C:\\Users\\lboukas\\duo_novo_performance_results.rda")

total_f <- ppv_pf_alt["total", ]
total_m <- ppv_pm_alt["total", ]


pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_ratio_father_mother_dnm.pdf", 
    height = 2.2, width = 2.2, pointsize = 8)
#par(mar = c(4, 4, 1, 1) + 0.1)
plot(total_f/total_m,  pch = 19, 
     cex = 1, bty = 'l', col = "dark orange", xlim = c(0.8, 40.2), ylim = c(0, max(total_f/total_m)),
     xaxt = "n", 
     yaxt = "n", ylab = "ratio of de novo from father-proband to mother-proband", 
     xlab = "All trios")
axis(2, at = c(0, 4, 8, 12))
abline(h = median(total_f/total_m), lty = "longdash", col = rgb(0,0,0,0.7))
#axis(1, at = c(1:40), cex.axis = 0.3)
dev.off()




library(readr)
trio_info_full <- read_csv('C:\\Users\\lboukas\\trio_info.csv')
trio_info_full <- trio_info_full[which(trio_info_full$`Study ID` %in% colnames(ppv_pf_alt)), ]
trio_info <- data.frame(ID1 = trio_info_full$`Family ID`, 
                        ID2 = trio_info_full$`Study ID`,
                        father_age = round(trio_info_full$`Father Age at Birth`), 
                        mother_age = round(trio_info_full$`Mother Age at Birth`), 
                        n_denovo_pf = sapply(trio_info_full$`Study ID`, function(xx) 
                          ppv_pf_alt['total', xx]), 
                        n_denovo_pm = sapply(trio_info_full$`Study ID`, function(xx) 
                          ppv_pm_alt['total', xx]))

pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_dnm_vs_age.pdf", 
    height = 2.2, width = 2.2, pointsize = 8)
par(mar = c(4, 4, 2, 2) + 0.1)
same_parent1 <- c('24-435405', '24-435607')
same_parent2 <- c('24-441639', '24-441781')
same_parent_indices <- which(trio_info$ID2 %in% c(same_parent1, same_parent2))
plot(trio_info$father_age[-same_parent_indices], trio_info$n_denovo_pf[-same_parent_indices], pch = 19, 
     cex = 1.05, bty = 'l', col = "dark orange",
     xlab = "parent age at birth", xaxt = 'n', yaxt = 'n',
     ylab = "# de novo classifications", xlim = c(0, 53), ylim = c(0, 53))

same_parent_indices1 <- which(trio_info$ID2 %in% c(same_parent1))
same_parent_indices2 <- which(trio_info$ID2 %in% c(same_parent2))
points(trio_info$father_age[same_parent_indices1], 
       trio_info$n_denovo_pf[same_parent_indices1], pch = 19, 
       cex = 1.05, bty = 'l', col = "red")
points(trio_info$father_age[same_parent_indices2], 
       trio_info$n_denovo_pf[same_parent_indices2], pch = 19, 
       cex = 1.05, bty = 'l', col = "blue")
axis(1, at = c(0, 25, 50))
axis(2, at = c(0, 25, 50))


points(trio_info$mother_age[-same_parent_indices], trio_info$n_denovo_pm[-same_parent_indices], pch = 2, 
     cex = 1.05, bty = 'l', col = "deep pink")
points(trio_info$mother_age[same_parent_indices1], 
       trio_info$n_denovo_pm[same_parent_indices1], pch = 2, 
       cex = 1.05, bty = 'l', col = "red")
points(trio_info$mother_age[same_parent_indices2], 
       trio_info$n_denovo_pm[same_parent_indices2], pch = 2, 
       cex = 1.05, bty = 'l', col = "blue")
legend("topleft", legend = c("father-proband", "mother-proband"), pch = c(19, 2), 
       bty = 'n', col = c("dark orange", "deep pink"))

fit_f <- glm(n_denovo_pf ~ father_age, data = trio_info, family = poisson(link = "log"))
fit_m <- glm(n_denovo_pm ~ mother_age, data = trio_info, family = poisson(link = "log"))

# Generate predicted values for a smooth trend line
age_range <- seq(15, 55, by = 1)
pred_f <- predict(fit_f, newdata = data.frame(father_age = age_range), type = "response")
pred_m <- predict(fit_m, newdata = data.frame(mother_age = age_range), type = "response")

lines(age_range, pred_f, col = "dark orange", lwd = 1)
lines(age_range, pred_m, col = "deep pink", lwd = 1)
dev.off()




