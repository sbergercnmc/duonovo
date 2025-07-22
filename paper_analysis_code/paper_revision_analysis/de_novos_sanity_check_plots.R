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

## import spreadsheet with parental age info
library(readr)
trios <- read_csv('~/Downloads/trio_info_updated.csv')
##


### use the PPV matrix to get total count of de novos
ppv_pf <- matrix(NA, ncol = length(all_dirs), nrow = 3)
ppv_pm <- matrix(NA, ncol = length(all_dirs), nrow = 3)
rownames(ppv_pf) <- rownames(ppv_pm) <- c('dn', 'assessed', 'total')
colnames(ppv_pf) <- colnames(ppv_pm) <- all_dirs

for (i in 1:length(all_dirs)){
  setwd(all_dirs[i])
  all_data <- list.files()
  load(file = all_data[grep('ppv', all_data)])
  
  ppv_pf['dn', i] <- ppv_pf_alt[1]
  ppv_pf['assessed', i] <- ppv_pf_alt[2]
  ppv_pf['total', i] <- ppv_pf_alt[3]
  
  ppv_pm['dn', i] <- ppv_pm_alt[1]
  ppv_pm['assessed', i] <- ppv_pm_alt[2]
  ppv_pm['total', i] <- ppv_pm_alt[3]
  
  setwd(data_directory)
}

trios$n_denovo_pf <- ppv_pf["total", trios$child_SAMPLEID]
trios$n_denovo_pm <- ppv_pm["total", trios$child_SAMPLEID]

### total de novo counts PF vs PM
total_f <- ppv_pf["total", ]
total_m <- ppv_pm["total", ]

pdf(paste0(figure_directory, "/duo_novo_ratio_father_mother_dnm.pdf"), 
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

### de novos vs parent age
pdf(paste0(figure_directory, "/duo_novo_dnm_vs_age.pdf"), 
    height = 2.2, width = 2.2, pointsize = 8)
plot(trios$father_age_at_birth, trios$n_denovo_pf, pch = 19, 
     cex = 1.05, bty = 'l', col = "dark orange",
     xlab = "parent age at birth", xaxt = 'n', yaxt = 'n',
     ylab = "# de novo classifications", xlim = c(0, 52), ylim = c(0, 52))
points(trios$mother_age_at_birth, trios$n_denovo_pm, pch = 2, 
       cex = 1.05, bty = 'l', col = "deep pink")
legend("topleft", legend = c("father-proband", "mother-proband"), pch = c(19, 2), 
       bty = 'n', col = c("dark orange", "deep pink"))

fit_f <- glm(n_denovo_pf ~ father_age_at_birth, data = trios, family = poisson(link = "log"))
fit_m <- glm(n_denovo_pm ~ mother_age_at_birth, data = trios, family = poisson(link = "log"))

# Generate predicted values for a smooth trend line
age_range <- seq(15, 55, by = 1)
pred_f <- predict(fit_f, newdata = data.frame(father_age_at_birth = age_range), type = "response")
pred_m <- predict(fit_m, newdata = data.frame(mother_age_at_birth = age_range), type = "response")

lines(age_range, pred_f, col = "dark orange", lwd = 1)
lines(age_range, pred_m, col = "deep pink", lwd = 1)
dev.off()


