# Regression of GlcCer vs. PD 
### Pull ATP10D, GBA1, and LRRK2 allele status from bfiles ahead of time

## Set up
library(data.table)
library(dplyr)
library(pROC)

setwd("~/Desktop/glccer_temp/")

## Data prep
### Load gene data
ATP10D <- fread("ATP10D_for_cov.raw")
colnames(ATP10D) <- c("FID","rs13106975","rs6827604","rs11729411")

GBA_LRRK2 <- fread("GBA_LRRK2.raw")

### Load GlcCer data
gc <- fread("Plasma_GlcCer_C24.txt")
colnames(gc) <- c("FID","IID","GlcCer")

### Load covariate data
cov <- fread("pdhc_covar_2.txt")

### merge all into one dataframe
merge1 <- merge(cov, ATP10D, by="FID", all.x = TRUE)
merge2 <- merge(merge1, GBA_LRRK2, by="FID", all.x = TRUE)
merge3 <- merge(merge2, gc, by="FID", all.x = TRUE)
data <- merge3[,c("FID","Sex","Age","Status","AJ","wbc","GBA_N370S","LRRK2_G2019S","rs13106975","rs6827604","rs11729411","GlcCer","PC1","PC2","PC3","PC4","PC5")]

### Remove SWEDD individuals, GBA N370S carriers, and LRRK2 G2019S carriers
data <- data[data$Status != 3,]
data <- data[data$GBA_N370S == 0,]
data <- data[data$LRRK2_G2019S == 0,]

### Recode disease status 
data$Status[data$Status == 1] <- 0
data$Status[data$Status == 2] <- 1

### Remove NAs
data <- data[!is.na(data$GlcCer),]

## Regressions and AUC analysis
### Some prep for AUC
### Rank the samples based on GlcCer levels and sort by rank
### Replace N with total sample count for the analysis
data$rank <- (rank(data$GlcCer)/N)
data <- data[order(data$rank), ]

### Basic analysis
plot(x=data$GlcCer, y = data$Status)
title(main = "Disease status x Total GlcCer - no cov")

fit_basic <-  glm(Status ~ GlcCer, data = data, family = binomial)
summary(fit_basic)

lines(data$GlcCer, fit_basic$fitted.values)
par(pty = "s")
roc(data$Status, fit_basic$fitted.values, plot=TRUE, legacy.axes=TRUE, xlab = "False Positive Rate", 
    ylab = "True Positive Rate", col = "royalblue2", lwd = 4, print.auc = TRUE)

roc.info.basic <- roc(data$Status, fit_basic$fitted.values, plot=TRUE, legacy.axes=TRUE)
roc.df.basic <- data.frame(tpp=roc.info.basic$sensitivities*100, fpp=(1-roc.info.basic$specificities)*100, thresholds=roc.info.basic$thresholds)
head(roc.df.basic) # at Inf-, meaning all samples are classified as PD, the tpp is 100 bc all PD cases were correctly identified \
# but the fpp is also 100 bc all the controls were incorrectly identified
tail(roc.df.basic)
# we want to find the threshold that maximizes tpp and minimizes fpp 

### Covariate analysis
plot(x=data$GlcCer, y = data$Status)
title(main = "Disease status x Total GlcCer - w covariates")

fit_covs <- glm(Status ~ GlcCer + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5, data = data, family = "binomial")
summary(fit_covs)

lines(data$GlcCer, fit_covs$fitted.values)
par(pty = "s")
roc(data$Status, fit_covs$fitted.values, plot=TRUE, legacy.axes=TRUE, xlab = "False Positive Rate", 
    ylab = "True Positive Rate", col = "maroon3", lwd = 4, print.auc = TRUE)

roc.info.covs <- roc(data$Status, fit_basic$fitted.values, plot=TRUE, legacy.axes=TRUE)
roc.df.covs <- data.frame(tpp=roc.info.covs$sensitivities*100, fpp=(1-roc.info.covs$specificities)*100, thresholds=roc.info.covs$thresholds)
head(roc.df.covs)
tail(roc.df.covs)

### ATP10D analysis
### substitute in top SNP for the given isoform 
plot(x=data$GlcCer, y = data$Status)
title(main = "Disease status x Total GlcCer - w ATP10D adj")
fit_ATP10D <- glm(Status ~ GlcCer + Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + rs6827604, data = data, family = "binomial")
summary(fit_ATP10D)

lines(data$GlcCer, fit_ATP10D$fitted.values)
par(pty = "s")
roc(data$Status, fit_ATP10D$fitted.values, plot=TRUE, legacy.axes=TRUE, xlab = "False Positive Rate", 
    ylab = "True Positive Rate", col = "darkgoldenrod2", lwd = 4, print.auc = TRUE)

roc.info.ATP10D <- roc(data$Status, fit_basic$fitted.values, plot=TRUE, legacy.axes=TRUE)
roc.df.ATP10D <- data.frame(tpp=roc.info.ATP10D$sensitivities*100, fpp=(1-roc.info.ATP10D$specificities)*100, thresholds=roc.info.ATP10D$thresholds)
head(roc.df.ATP10D)
tail(roc.df.ATP10D)

### Combine AUC plots 
tiff("AUC_plot_Total-GlcCer_noGBA-LRRK2.png", width = 6, height = 6, units = "in", res = 300)

roc(data$Status, fit_basic$fitted.values, plot=TRUE, legacy.axes=TRUE, xlab = "False Positive Rate", 
    ylab = "True Positive Rate", col = "royalblue2", lwd = 4, print.auc = TRUE, print.auc.x = 0.25)
plot.roc(data$Status, fit_covs$fitted.values, col = "maroon3", lwd = 4, print.auc = TRUE, add = TRUE, print.auc.y = 0.4, print.auc.x = 0.25)
plot.roc(data$Status, fit_ATP10D$fitted.values, col = "darkgoldenrod2", lwd = 4, print.auc = TRUE, add = TRUE, print.auc.y = 0.3, print.auc.x = 0.25)
legend("bottomright", legend=c("no_cov","age_sex_5PCs","ATP10D"),col=c("royalblue2","maroon3","darkgoldenrod2"), lwd = 4, text.width = 0.28)

dev.off()






