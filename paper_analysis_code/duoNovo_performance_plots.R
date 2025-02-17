pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_father_proband_ppv.pdf", 
    height = 2.2, width = 2.2, pointsize = 8)

precision <- ppv_pf_alt["dn", ]/ppv_pf_alt["assessed", ]
number_called <- ppv_pf_alt["assessed", ]

plot(jitter(number_called, 2), precision, pch = 19, 
     cex = 0.75, bty = 'l', col = "dark orange", xlim = c(0, max(number_called)), ylim = c(0, 1),
     yaxt = 'n', xaxt = 'n', ylab = "positive predictive value", xlab = "# de novos called")
axis(1, at = c(0, 20, 40))
axis(2, at = c(0, 0.45, 0.9))
dev.off()


pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_father_proband_ppv2.pdf", 
    height = 2.2, width = 2.2, pointsize = 8)

precision <- ppv_pf_alt["dn", ]/ppv_pf_alt["assessed", ]

plot(sort(precision, decreasing = TRUE), pch = 19, 
     cex = 0.5, bty = 'l', col = "dark orange", ylim = c(0, 1),
     yaxt = 'n', xaxt = 'n', ylab = "positive predictive value", xlab = "father-proband duos")
axis(2, at = c(0, 0.5, 1))
dev.off()

pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_father_proband_n_called.pdf", 
    height = 1.5, width = 2.2, pointsize = 8)

number_called <- ppv_pf_alt["assessed", ]

plot(number_called[order(precision, decreasing = TRUE)], pch = 19, type = 'h', lwd = 1,
     cex = 0.5, bty = 'l', col = "dark orange", ylim = c(0, max(number_called)), 
     yaxt = 'n', xaxt = 'n', ylab = "# de novos called", xlab = "father-proband duos")
axis(2, at = c(0, 20, 40))
dev.off()

############
pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_father_proband_npv.pdf", 
    height = 2.4, width = 4, pointsize = 8)
npv <- npv_pf_alt[1, ]
baseline_npv <- npv_pf_alt[2, ]
plot(npv, pch = 19, 
     cex = 0.75, bty = 'l', col = "dark orange", ylim = c(0.9995, 1), xlab = "father-proband duos", 
     ylab = "negative predictive value", 
     xaxt = 'n', yaxt = 'n', xlim = c(0.8, 40))
points(baseline_npv, cex = 0.75, bty = 'l', col = "cornflowerblue", pch = 19)
axis(2, at = c(0.9995, 1))
abline(v = seq(1.5, 39.5, by = 1), col = rgb(0,0,0,0.4), lty = "longdash")
legend <- legend("bottomright", legend = c("duoNovo", "naive baseline"), 
                 pch = 20, bty = 'n', cex = 0.62, 
                 col = c("orange", "cornflowerblue"))
dev.off()


###########
pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_mother_proband_ppv.pdf", 
    height = 2.2, width = 2.2, pointsize = 8)
precision <- ppv_mf_alt["dn", ]/ppv_pm_alt["assessed", ]
number_called <- ppv_pm_alt["assessed", ]

plot(jitter(number_called, 2), precision, pch = 19, 
     cex = 0.75, bty = 'l', col = "dark orange", xlim = c(0, max(number_called)), ylim = c(0, 1),
     yaxt = 'n', xaxt = 'n', ylab = "positive predictive value", xlab = "# de novos called")
axis(1, at = c(0, 8, 16))
axis(2, at = c(0, 0.45, 0.9))
dev.off()

pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_mother_proband_ppv2.pdf", 
    height = 2.2, width = 2.2, pointsize = 8)

precision <- ppv_mf_alt["dn", ]/ppv_pm_alt["assessed", ]

plot(sort(precision, decreasing = TRUE), pch = 19, 
     cex = 0.5, bty = 'l', col = "dark orange", ylim = c(0, 1),
     yaxt = 'n', xaxt = 'n', ylab = "positive predictive value", xlab = "father-proband duos")
axis(2, at = c(0, 0.5, 1))
dev.off()

pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_mother_proband_n_called.pdf", 
    height = 1.5, width = 2.2, pointsize = 8)

number_called <- ppv_pm_alt["assessed", ]

plot(number_called[order(precision, decreasing = TRUE)], pch = 19, type = 'h', lwd = 1,
     cex = 0.5, bty = 'l', col = "dark orange", ylim = c(0, max(number_called)), 
     yaxt = 'n', xaxt = 'n', ylab = "# de novos called", xlab = "father-proband duos")
axis(2, at = c(0, 5, 10))
dev.off()

##########
pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_mother_proband_npv.pdf", 
    height = 2.4, width = 4, pointsize = 8)
npv <- npv_pm_alt[1, ]
baseline_npv <- npv_pm_alt[2, ]
plot(npv, pch = 19, 
     cex = 0.75, bty = 'l', col = "dark orange", ylim = c(0.9995, 1), xlab = "mother-proband duos", 
     ylab = "negative predictive value", 
     xaxt = 'n', yaxt = 'n', xlim = c(0.8, 40))
points(baseline_npv, cex = 0.75, bty = 'l', col = "cornflowerblue", pch = 19)
axis(2, at = c(0.9995, 1))
abline(v = seq(1.5, 39.5, by = 1), col = rgb(0,0,0,0.4), lty = "longdash")
legend <- legend("top", legend = c("duoNovo", "naive baseline"), 
                 pch = 20, bty = 'n', cex = 0.62, 
                 col = c("orange", "cornflowerblue"))
dev.off()

##########
pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_false_positives.pdf", 
    height = 2.4, width = 6.5, pointsize = 8)
par(mar = c(3, 4, 3, 0.1) + 0.1)
plot(seq(1, 60, by = 1.5), 100*per_sample_false_pos_rate_f, 
     ylim = c(0, 0.001), pch = 19, col = rgb(1,0,0,0.57), xlab = "trio index", cex = 1.1,
     ylab = "% false positive duo novo classifications", 
     bty = 'l', xaxt = 'n', yaxt = 'n', xlim = c(0.8, 61))
points(seq(1.65, 61, by = 1.5), 100*per_sample_false_pos_rate_m, 
       pch = 19, col = "cornflowerblue", cex = 1.1)
#axis(1, at = seq(0.5, 79.5, by = 2),
#     cex.axis = 0.8, las = 1)
axis(2, at = c(0, 0.001))
abline(v = seq(2, 60, by = 1.5), lty = "longdash", col = rgb(0,0,0,0.4))
legend <- legend("topleft", legend = c("father-proband duos", "mother-proband duos"), 
                 pch = 20, bty = 'n', cex = 0.62, 
                 col = c(rgb(1,0,0,0.57), "cornflowerblue"))
dev.off()

##########
pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_ref_vs_alt.pdf", 
    height = 2.4, width = 7, pointsize = 8)
par(mfrow = c(1, 2))
plot(1:40, 100*dn_call_rate_alt_pf, 
     pch = 19, col = "salmon", xlab = "father-proband duos", cex = 0.8,
     ylab = "% classified as de novo", 
     bty = 'l', xaxt = 'n', yaxt = 'n', xlim = c(0.8, 40), ylim = c(0, 0.01))
points(1:40, 100*dn_call_rate_ref_pf, 
       pch = 19, col = "gray50", cex = 0.8)
axis(2, at = c(0, 0.01))
abline(v = seq(1.5, 39.5, by = 1), lty = "longdash", col = rgb(0,0,0,0.4))
legend <- legend("topleft", legend = c("ALT allele", "REF allele"), 
                 pch = 20, bty = 'n', cex = 0.62, 
                 col = c('salmon', "gray50"))

plot(1:40, 100*dn_call_rate_alt_pm, 
     ylim = c(0, 0.0035), pch = 19, col = "salmon", xlab = "mother-proband duos", cex = 0.8,
     ylab = "% classified as de novo", 
     bty = 'l', xaxt = 'n', yaxt = 'n', xlim = c(0.8, 40))
points(1:40, 100*dn_call_rate_ref_pm, 
       pch = 19, col = "gray50", cex = 0.8)
axis(2, at = c(0, 0.003))
abline(v = seq(1.5, 39.5, by = 1), lty = "longdash", col = rgb(0,0,0,0.4))
dev.off()



########
pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_vs_trio_de_novo.pdf", 
    height = 2.4, width = 4, pointsize = 8)
plot(1:40, 100*sensitivity, 
     ylim = c(0, 80), pch = 19, col = 'salmon', xlab = "trio index", cex = 0.8,
     ylab = "% of trio de novo variants", 
     bty = 'l', xaxt = 'n', yaxt = 'n', xlim = c(1, 41))
points(1:40, 100*false_ndn, 
       pch = 19, col = "gray50", cex = 0.8)
points(1:40, 100*other, 
        cex = 0.8, lwd = 1, col = "gray10")
#axis(1, at = seq(0.5, 79.5, by = 2),
#     cex.axis = 0.8, las = 1)
axis(2, at = c(0, 40, 80))
abline(v = seq(1.5, 39.5, by = 1), lty = "longdash", col = rgb(0,0,0,0.4))
legend <- legend("top", legend = c("de novo", "not de novo", "other"), 
                 pch = c(20, 20, 1), bty = 'n', cex = 0.62, 
                 col = c("salmon", "gray50", "gray10"))
dev.off()


########
########
######## Aggregate plots
library(ggplot2)
library(forcats)

pdf(file = "C:/Users/lboukas/duoNovo_figures/pf_transmitted_aggregate.pdf", height = 1.8, width = 2.4, pointsize = 8)
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

##
pdf(file = "C:/Users/lboukas/duoNovo_figures/pm_transmitted_aggregate.pdf", height = 1.8, width = 2.4, pointsize = 8)
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

pdf(file = "C:/Users/lboukas/duoNovo_figures/pf_candidates_aggregate.pdf", 
    height = 1.8, width = 2.4, pointsize = 8)
data <- data.frame(
  category = names(candidate_percentages_f),
  value = candidate_percentages_f
)
data$category <- fct_relevel(data$category, "uncertain_f", "ndn_f", "dn_f")

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

pdf(file = "C:/Users/lboukas/duoNovo_figures/pm_candidates_aggregate.pdf", 
    height = 1.8, width = 2.4, pointsize = 8)
data <- data.frame(
  category = names(candidate_percentages_m),
  value = candidate_percentages_m
)
data$category <- fct_relevel(data$category, "uncertain_m", "ndn_m", "dn_m")

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


######
###### Ancestry-specific plots
pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_father_proband_ppv2_ancestry.pdf", 
    height = 2.2, width = 4.4, pointsize = 8)
par(mfrow = c(1, 2))
precision_eur <- trio_info_full$ppv_pf[which(trio_info_full$`Proband 1KG Group` == "EUR")]

plot(sort(precision_eur, decreasing = TRUE), pch = 19, 
     cex = 0.5, bty = 'l', col = "deep pink", ylim = c(0, 1),
     yaxt = 'n', xaxt = 'n', ylab = "positive predictive value", xlab = "father-proband duos")
axis(2, at = c(0, 0.5, 1))

precision_non_eur <- trio_info_full$ppv_pf[-which(trio_info_full$`Proband 1KG Group` == "EUR")]

plot(sort(precision_non_eur, decreasing = TRUE), pch = 19, 
     cex = 0.5, bty = 'l', col = "forest green", ylim = c(0, 1),
     yaxt = 'n', xaxt = 'n', ylab = "positive predictive value", xlab = "father-proband duos")
axis(2, at = c(0, 0.5, 1))
dev.off()

pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_father_proband_n_called_ancestry.pdf", 
    height = 1.5, width = 4.4, pointsize = 8)
par(mfrow = c(1, 2))
y_lim <-  c(0, max(trio_info_full$n_dn_assessed_pf))
number_called <- trio_info_full$n_dn_assessed_pf[which(trio_info_full$`Proband 1KG Group` == "EUR")]

plot(number_called[order(precision_eur, decreasing = TRUE)], pch = 19, type = 'h', lwd = 1,
     cex = 0.5, bty = 'l', col = "deep pink", ylim = y_lim, 
     yaxt = 'n', xaxt = 'n', ylab = "# de novos called", xlab = "father-proband duos")
axis(2, at = c(0, 20, 40))

number_called <- trio_info_full$n_dn_assessed_pf[-which(trio_info_full$`Proband 1KG Group` == "EUR")]

plot(number_called[order(precision_non_eur, decreasing = TRUE)], pch = 19, type = 'h', lwd = 1,
     cex = 0.5, bty = 'l', col = "forest green", ylim = y_lim, 
     yaxt = 'n', xaxt = 'n', ylab = "# de novos called", xlab = "father-proband duos")
axis(2, at = c(0, 20, 40))
dev.off()


###
pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_mother_proband_ppv2_ancestry.pdf", 
    height = 2.2, width = 4.4, pointsize = 8)
par(mfrow = c(1, 2))
precision_eur <- trio_info_full$ppv_pm[which(trio_info_full$`Proband 1KG Group` == "EUR")]

plot(sort(precision_eur, decreasing = TRUE), pch = 19, 
     cex = 0.5, bty = 'l', col = "deep pink", ylim = c(0, 1),
     yaxt = 'n', xaxt = 'n', ylab = "positive predictive value", xlab = "mother-proband duos")
axis(2, at = c(0, 0.5, 1))

precision_non_eur <- trio_info_full$ppv_pm[-which(trio_info_full$`Proband 1KG Group` == "EUR")]

plot(sort(precision_non_eur, decreasing = TRUE), pch = 19, 
     cex = 0.5, bty = 'l', col = "forest green", ylim = c(0, 1),
     yaxt = 'n', xaxt = 'n', ylab = "positive predictive value", xlab = "mother-proband duos")
axis(2, at = c(0, 0.5, 1))
dev.off()

pdf(file = "C:/Users/lboukas/duoNovo_figures/duo_novo_performance_mother_proband_n_called_ancestry.pdf", 
    height = 1.5, width = 4.4, pointsize = 8)
par(mfrow = c(1, 2))
y_lim <-  c(0, max(trio_info_full$n_dn_assessed_pm))

number_called <- trio_info_full$n_dn_assessed_pm[which(trio_info_full$`Proband 1KG Group` == "EUR")]

plot(number_called[order(precision_eur, decreasing = TRUE)], pch = 19, type = 'h', lwd = 1,
     cex = 0.5, bty = 'l', col = "deep pink", ylim = y_lim, 
     yaxt = 'n', xaxt = 'n', ylab = "# de novos called", xlab = "mother-proband duos")
axis(2, at = c(0, 5, 10))

number_called <- trio_info_full$n_dn_assessed_pm[-which(trio_info_full$`Proband 1KG Group` == "EUR")]

plot(number_called[order(precision_non_eur, decreasing = TRUE)], pch = 19, type = 'h', lwd = 1,
     cex = 0.5, bty = 'l', col = "forest green", ylim = y_lim, 
     yaxt = 'n', xaxt = 'n', ylab = "# de novos called", xlab = "mother-proband duos")
axis(2, at = c(0, 5, 10))
dev.off()










