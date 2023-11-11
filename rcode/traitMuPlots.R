# November 10, 1023
# Aim of this code is to create mu plots for the five trait models
rm(list=ls())
options(stringsAsFactors = FALSE)

load("output/heightDummyIntGrandZ.Rdata")
htModelFit <- rstan::extract(mdlHt)
sumHt <- summary(mdlHt)$summary
col4table <- c("mean","sd","2.5%","25%","50%","75%","97.5%","Rhat")
col4fig <- c("mean","sd","25%","50%","75%","Rhat")

mu_params <- c( 
  "mu_grand",
  "b_tranE",
  "b_tranlat",
 # "muPhenoSp",
  "muForceSp", 
  "muChillSp",
  "muPhotoSp",
  "betaTraitxForce",
  "betaTraitxChill",
  "betaTraitxPhoto")

meanHt <- sumHt[mu_params, col4table]

rownames(meanHt) = c( 
  #"Root trait intercept", "Lambda",
  "Grand mean",
  "Transect",
  "Transect by latitude",
#  "Budburst slope",
  "Forcing",
  "Photoperiod",
  "Chilling",
  "Trait-forcing effect",
  "Trait-chilling effect",
  "Trait-photoperiod effect"
)

load("output/lmaDummyIntGrandZ.Rdata")
lmaModelFit <- rstan::extract(mdlLMA)
sumLMA <- summary(mdlLMA)$summary

meanLMA <- sumLMA[mu_params, col4table]

rownames(meanLMA) = c( 
  #"Root trait intercept", "Lambda",
  "Grand mean",
  "Transect",
  "Transect by latitude",
 # "Budburst slope",
  "Forcing",
  "Photoperiod",
  "Chilling",
  "Trait-forcing effect",
  "Trait-chilling effect",
  "Trait-photoperiod effect"
)

load("output/dbhDummyIntGrandZ.Rdata")
dbhModelFit <- rstan::extract(mdlDBH)
sumDBH <- summary(mdlDBH)$summary

meanDBH <- sumDBH[mu_params, col4table]

rownames(meanDBH) = c( 
  #"Root trait intercept", "Lambda",
  "Grand mean",
  "Transect",
  "Transect by latitude",
 # "Budburst slope",
  "Forcing",
  "Photoperiod",
  "Chilling",
  "Trait-forcing effect",
  "Trait-chilling effect",
  "Trait-photoperiod effect"
)

load("output/ssdDummyIntGrandZ.Rdata")
ssdModelFit <- rstan::extract(mdlSSD)
sumSSD <- summary(mdlSSD)$summary

meanSSD <- sumSSD[mu_params, col4table]

rownames(meanSSD) = c( 
  #"Root trait intercept", "Lambda",
  "Grand mean",
  "Transect",
  "Transect by latitude",
 # "Budburst slope",
  "Forcing",
  "Photoperiod",
  "Chilling",
  "Trait-forcing effect",
  "Trait-chilling effect",
  "Trait-photoperiod effect"
)

load("output/cnDummyIntGrandZ.Rdata")
cnModelFit <- rstan::extract(mdlCN)
sumCN <- summary(mdlCN)$summary

meanCN <- sumCN[mu_params, col4table]

rownames(meanCN) = c( 
  #"Root trait intercept", "Lambda",
  "Grand mean",
  "Transect",
  "Transect by latitude",
  #"Budburst slope",
  "Forcing",
  "Photoperiod",
  "Chilling",
  "Trait-forcing effect",
  "Trait-chilling effect",
  "Trait-photoperiod effect"
)

pdf(file.path( "figures/triatMuPlot.pdf"), width = 7, height = 5)
par(mfrow = c(1,5), mar = c(5, 8, 2, 1))
# Upper panel: bud burst
plot(seq(-10,  10,length.out = nrow(meanHt)), 1:nrow(meanHt),
  type = "n", xlab = "Estimated change in budburst day",
  ylab = "", yaxt = "n")

axis(2, at = nrow(meanHt):1, labels = rownames(meanHt), las = 1, cex.axis = 0.8)
points(meanHt[, 'mean'],
  nrow(meanHt):1,
  pch = 16,
  col = "cyan4",
  cex = 1.5)
arrows(meanHt[, "97.5%"], nrow(meanHt):1, meanHt[, "2.5%"], nrow(meanHt):1,
  len = 0, col = "black")
abline(v = 0, lty = 3)

# LMA
par(mar = c(5, 1, 2, 1))

plot(seq(-10,  10,length.out = nrow(meanLMA)), 1:nrow(meanLMA),
  type = "n", xlab = "Estimated change in budburst day",
  ylab = "", yaxt = "n")
abline(v = 0, lty = 3)

points(meanLMA[, 'mean'],
  nrow(meanLMA):1,
  pch = 16,
  col = "maroon",
  cex = 1.5)
arrows(meanLMA[, "97.5%"], nrow(meanLMA):1, meanLMA[, "2.5%"], nrow(meanLMA):1,
  len = 0, col = "black")

# DBH
plot(seq(-10,  10,length.out = nrow(meanDBH)), 1:nrow(meanDBH),
  type = "n", xlab = "Estimated change in budburst day",
  ylab = "", yaxt = "n")
abline(v = 0, lty = 3)

points(meanDBH[, 'mean'],
  nrow(meanDBH):1,
  pch = 16,
  col = "goldenrod",
  cex = 1.5)
arrows(meanDBH[, "97.5%"], nrow(meanDBH):1, meanDBH[, "2.5%"], nrow(meanDBH):1,
  len = 0, col = "black")

# SSD

plot(seq(-10,  10,length.out = nrow(meanSSD)), 1:nrow(meanSSD),
  type = "n", xlab = "Estimated change in budburst day",
  ylab = "", yaxt = "n")
abline(v = 0, lty = 3)

points(meanSSD[, 'mean'],
  nrow(meanSSD):1,
  pch = 16,
  col = "purple4",
  cex = 1.5)
arrows(meanSSD[, "97.5%"], nrow(meanSSD):1, meanSSD[, "2.5%"], nrow(meanSSD):1,
  len = 0, col = "black")

# CN
plot(seq(-10,  10,length.out = nrow(meanCN)), 1:nrow(meanCN),
  type = "n", xlab = "Estimated change in budburst day",
  ylab = "", yaxt = "n")
abline(v = 0, lty = 3)

points(meanCN[, 'mean'],
  nrow(meanCN):1,
  pch = 16,
  col = "darkolivegreen",
  cex = 1.5)
arrows(meanCN[, "97.5%"], nrow(meanCN):1, meanCN[, "2.5%"], nrow(meanCN):1,
  len = 0, col = "black")
dev.off()
