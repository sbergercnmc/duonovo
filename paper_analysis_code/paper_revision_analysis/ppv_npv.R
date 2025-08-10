args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript ppv_npv.R <data_directory> <figure_directory>")
}

data_directory   <- args[1]
figure_directory <- args[2]

# ensure data directory exists
if (!dir.exists(data_directory)) {
  stop(sprintf("Data directory not found: '%s'", data_directory))
}

setwd(data_directory)
all_dirs <- list.files()

### compute PPV 
ppv_pf <- matrix(NA, ncol = length(all_dirs), nrow = 4)
ppv_pf_no_gnomad <- ppv_pf
ppv_pm <- matrix(NA, ncol = length(all_dirs), nrow = 4)
ppv_pm_no_gnomad <- ppv_pm
rownames(ppv_pf) <- rownames(ppv_pf_no_gnomad) <- rownames(ppv_pm) <- rownames(ppv_pm_no_gnomad) <- c('dn', 'assessed', 'total', 'all_candidates')

for (i in 1:length(all_dirs)){
  setwd(all_dirs[i])
  all_data <- list.files()
  load(file = all_data[grep('ppv', all_data)])
  
  ppv_pf['dn', i] <- ppv_pf_alt[1]
  ppv_pf['assessed', i] <- ppv_pf_alt[2]
  ppv_pf['total', i] <- ppv_pf_alt[3]
  ppv_pf['all_candidates', i] <- ppv_pf_alt[4]
  
  ppv_pf_no_gnomad['dn', i] <- ppv_pf_alt_no_gnomad[1]
  ppv_pf_no_gnomad['assessed', i] <- ppv_pf_alt_no_gnomad[2]
  ppv_pf_no_gnomad['total', i] <- ppv_pf_alt_no_gnomad[3]
  ppv_pf_no_gnomad['all_candidates', i] <- ppv_pf_alt_no_gnomad[4]
  
  ppv_pm['dn', i] <- ppv_pm_alt[1]
  ppv_pm['assessed', i] <- ppv_pm_alt[2]
  ppv_pm['total', i] <- ppv_pm_alt[3]
  ppv_pm['all_candidates', i] <- ppv_pm_alt[4]
  
  ppv_pm_no_gnomad['dn', i] <- ppv_pm_alt_no_gnomad[1]
  ppv_pm_no_gnomad['assessed', i] <- ppv_pm_alt_no_gnomad[2]
  ppv_pm_no_gnomad['total', i] <- ppv_pm_alt_no_gnomad[3]
  ppv_pm_no_gnomad['all_candidates', i] <- ppv_pm_alt_no_gnomad[4]
  
  setwd(data_directory)
}

ppv_all_pf_no_gnomad <- sum(ppv_pf_no_gnomad[1, ], na.rm = TRUE)/sum(ppv_pf_no_gnomad[2, ], na.rm = TRUE)
ppv_all_pm_no_gnomad <- sum(ppv_pm_no_gnomad[1, ], na.rm = TRUE)/sum(ppv_pm_no_gnomad[2, ], na.rm = TRUE)
ppv_all_pf <- sum(ppv_pf[1, ], na.rm = TRUE)/sum(ppv_pf[2, ], na.rm = TRUE)
ppv_all_pm <- sum(ppv_pm[1, ], na.rm = TRUE)/sum(ppv_pm[2, ], na.rm = TRUE)


### compute NPV
npv_pf_mat <- matrix(NA, ncol = length(all_dirs), nrow = 2)
npv_pm_mat <- matrix(NA, ncol = length(all_dirs), nrow = 2)
rownames(npv_pf_mat) <- rownames(npv_pm_mat) <- c('duoNovo', 'naive')
for (i in 1:length(all_dirs)){
  setwd(all_dirs[i])
  all_data <- list.files()
  load(file = all_data[grep('npv', all_data)])
  
  npv_pf_mat['duoNovo', i] <- npv_pf[1]
  npv_pf_mat['naive', i] <- npv_pf[2]
  
  npv_pm_mat['duoNovo', i] <- npv_pm[1]
  npv_pm_mat['naive', i] <- npv_pm[2]
  
  setwd(data_directory)
}


### PPV plots 
pdf(file = paste0(figure_directory, "/ppv_no_gnomad.pdf"), height = 3.5, width = 5.9, pointsize = 8)
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

pdf(file = paste0(figure_directory, "/ppv_n_called.pdf"), 
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

# collective PPV
pdf(file = paste0(figure_directory, "/collective_ppv.pdf"), height = 2.8, width = 2.2, pointsize = 8)
df <- data.frame(
  Condition = factor(c("PF no gnomAD", "PM no gnomAD", "PF", "PM"),
                     levels = c("PF no gnomAD", "PM no gnomAD", "PF", "PM")),
  PPV       = c(ppv_all_pf_no_gnomad, ppv_all_pm_no_gnomad, ppv_all_pf, ppv_all_pm)
)

ggplot(df, aes(x = Condition, y = PPV, fill = Condition)) +
  geom_col() +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  labs(
    title = "",
    y     = "PPV",
    x     = ""
  ) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x    = element_text(angle = 45, hjust = 1))
dev.off()

### NPV plots
pdf(file = paste0(figure_directory, "/npv.pdf"), 
    height = 5.5, width = 6, pointsize = 8)
par(mfrow = c(2, 1))
npv <- npv_pf_mat['duoNovo', ]
baseline_npv <- npv_pf_mat['naive', ]
plot(npv, pch = 19, 
     cex = 0.75, bty = 'l', col = "dark orange", ylim = c(0.998, 1), xlab = "father-proband duos", 
     ylab = "negative predictive value", 
     xaxt = 'n', yaxt = 'n', xlim = c(0.8, 117))
points(baseline_npv, cex = 0.6, bty = 'l', col = "cornflowerblue", pch = 19)
axis(2, at = c(0.998, 0.999, 1))
abline(v = seq(1.5, 116.5, by = 1), col = rgb(0,0,0,0.4), lty = "longdash", lwd = 0.5)
legend <- legend("bottomright", legend = c("duoNovo", "naive baseline"), 
                 pch = 20, bty = 'n', cex = 0.62, 
                 col = c("orange", "cornflowerblue"))

npv <- npv_pm_mat['duoNovo', ]
baseline_npv <- npv_pm_mat['naive', ]
plot(npv, pch = 19, 
     cex = 0.75, bty = 'l', col = "dark orange", ylim = c(0.998, 1), xlab = "mother-proband duos", 
     ylab = "negative predictive value", 
     xaxt = 'n', yaxt = 'n', xlim = c(0.8, 117))
points(baseline_npv, cex = 0.6, bty = 'l', col = "cornflowerblue", pch = 19)
axis(2, at = c(0.998, 0.999, 1))
abline(v = seq(1.5, 116.5, by = 1), col = rgb(0,0,0,0.4), lty = "longdash", lwd = 0.5)

dev.off()

pdf(file = paste0(figure_directory, "/father_proband_npv_hist.pdf"), 
    height = 2.2, width = 2.2, pointsize = 8)
npv <- npv_pf_mat['duoNovo', ]
baseline_npv <- npv_pf_mat['naive', ]
hist(npv - baseline_npv, breaks = 35, xlab = "duoNovo NPV - naive NPV", freq = FALSE, 
     col = "dark orange", main = "", xlim = c(-1e-5, 7e-4), xaxt = 'n')
axis(1, at = c(0, 7e-4))
dev.off()

pdf(file = paste0(figure_directory, "/mother_proband_npv_hist.pdf"), 
    height = 2.2, width = 2.2, pointsize = 8)
npv <- npv_pm_mat['duoNovo', ]
baseline_npv <- npv_pm_mat['naive', ]
hist(npv - baseline_npv, breaks = 35, xlab = "duoNovo NPV - naive NPV", freq = FALSE, 
     col = "dark orange", main = "", xlim = c(-1e-5, 7e-4), xaxt = 'n')
axis(1, at = c(0, 7e-4))
dev.off()





