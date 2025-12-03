args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript ppv_npv.R <data_directory> <figure_directory> <trio_info_directory>")
}
data_directory   <- args[1]
figure_directory <- args[2]
trio_info_directory <- args[3]

# ensure data directory exists
if (!dir.exists(data_directory)) {
  stop(sprintf("Data directory not found: '%s'", data_directory))
}

setwd(data_directory)
all_dirs <- list.files()

## import spreadsheet with parental age info
library(readr)
trios <- read_csv(trio_info_directory)
trios <- trios[-which(trios$child_SAMPLEID %in% problematic_trios), ]
##

### the following lines use the ppv_pf and ppv_pm matrices from the ppv_revision script
### total de novo counts PF vs PM
trios$n_denovo_pf <- ppv_pf["total", trios$child_SAMPLEID]
trios$n_denovo_pm <- ppv_pm["total", trios$child_SAMPLEID]

total_f <- ppv_pf["total", ]
total_m <- ppv_pm["total", ]

pdf(paste0(figure_directory, "/duo_novo_ratio_father_mother_dnm_revision.pdf"), 
    height = 2.2, width = 4.4, pointsize = 8)
#par(mar = c(4, 4, 1, 1) + 0.1)
plot(total_f/total_m,  pch = 19, 
     cex = 1, bty = 'l', col = "dark orange", xlim = c(0.8, 104.2), ylim = c(0, max(total_f/total_m)),
     xaxt = "n", 
     yaxt = "n", ylab = "ratio of de novo from father-proband to mother-proband", 
     xlab = "All trios")
axis(2, at = c(1, 6, 11))
abline(h = median(total_f/total_m), lty = "longdash", col = rgb(0,0,0,0.7))
#axis(1, at = c(1:40), cex.axis = 0.3)
dev.off()

### de novos vs parent age
pdf(paste0(figure_directory, "/duo_novo_dnm_vs_age_revision.pdf"), 
    height = 3, width = 2.2, pointsize = 8)
plot(trios$father_age_at_birth, trios$n_denovo_pf, pch = 19, 
     cex = 0.5, bty = 'l', col = "dark orange",
     xlab = "parent age at birth", 
     ylab = "# de novo classifications", xlim = c(0, 52), ylim = c(0, 52))
points(trios$mother_age_at_birth, trios$n_denovo_pm, pch = 2, 
       cex = 0.5, bty = 'l', col = "deep pink")
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


