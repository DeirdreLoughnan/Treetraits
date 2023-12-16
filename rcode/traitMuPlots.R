# November 10, 1023
# Aim of this code is to create mu plots for the five trait models
rm(list=ls())
options(stringsAsFactors = FALSE)

require(rstan)
library(forcats)
library(ggdist)
library(reshape2)
require(cowplot)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/Treetraits") 
}  else{
  setwd("/home/deirdre/Treetraits") # for midge
}

col4table <- c("mean","sd","2.5%","25%","50%","75%","97.5%","Rhat")
col4fig <- c("mean","sd","25%","50%","75%","Rhat")

mu_trait <- c( 
 # "mu_grand",
  "b_tranE",
  "b_tranlat")

mu_pheno <- c(
 # "muPhenoSp",
  "muForceSp", 
  "muChillSp",
  "muPhotoSp")

mu_bTC <- c(
  "betaTraitxForce",
  "betaTraitxChill",
  "betaTraitxPhoto")

quantile590 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.05, 0.9))
  return(returnQuanilte)
}


### mu plot of trait part

load("output/heightDummyIntGrandZ25.Rdata")
postHt<- rstan::extract(mdlHt)
sumHt <- summary(mdlHt)$summary

bTranHt <- apply(data.frame(postHt$b_tranE), 2, quantile590 )
bTranLatHt <- apply(data.frame(postHt$b_tranlat), 2, quantile590 )

traitHt <- sumHt[mu_trait, col4table]
traitHt <- cbind(traitHt, bTran_quan, bTranLat_quan)

load("output/dbhDummyIntGrandZ25.Rdata")
postDBH<- rstan::extract(mdlDBH)
sumDBH <- summary(mdlDBH)$summary

bTranDBH <- apply(data.frame(postDBH$b_tranE), 2, quantile590 )
bTranLatDBH <- apply(data.frame(postDBH$b_tranlat), 2, quantile590 )

traitDBH <- sumDBH[mu_trait, col4table]
traitDBH <- cbind(traitDBH, bTran_quan, bTranLat_quan)

rownames(traitDBH) = c( 
  "Transect",
  "Transect x latitude"
)

load("output/lmaDummyIntGrandZ25.Rdata")
postLMA<- rstan::extract(mdlLMA)
sumLMA <- summary(mdlLMA)$summary

bTranLMA <- apply(data.frame(postLMA$b_tranE), 2, quantile590 )
bTranLatLMA <- apply(data.frame(postLMA$b_tranlat), 2, quantile590 )

traitLMA <- sumLMA[mu_trait, col4table]
traitLMA <- cbind(traitLMA, bTran_quan, bTranLat_quan)

rownames(traitLMA) = c( 
  "Transect",
  "Transect x latitude"
)

load("output/cnDummyIntGrandZ25.Rdata")
postCN<- rstan::extract(mdl)
sumCN <- summary(mdl)$summary

quantile590 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.05, 0.9))
  return(returnQuanilte)
}

bTranCN <- apply(data.frame(postCN$b_tranE), 2, quantile590 )
bTranLatCN <- apply(data.frame(postCN$b_tranlat), 2, quantile590 )

traitCN <- sumCN[mu_trait, col4table]
traitCN <- cbind(traitCN, bTran_quan, bTranLat_quan)

rownames(traitCN) = c( 
  "Transect",
  "Transect x latitude"
)

load("output/ssdDummyIntGrandZ25.Rdata")
postSSD<- rstan::extract(mdlSSD)
sumSSD <- summary(mdlSSD)$summary

bTranSSD <- apply(data.frame(postSSD$b_tranE), 2, quantile590 )
bTranLatSSD <- apply(data.frame(postSSD$b_tranlat), 2, quantile590 )

traitSSD <- sumSSD[mu_trait, col4table]
traitSSD <- cbind(traitSSD, bTran_quan, bTranLat_quan)

rownames(traitSSD) = c( 
  "Transect",
  "Transect x latitude"
)

pdf(file.path( "../figures/traitMuPlot.pdf"), width = 15, height = 3)
par(mfrow = c(1,5), mar = c(5, 8, 2, 0.05))
# Upper panel: bud burst
plot(seq(-1,1,length.out = nrow(traitHt)), 1:nrow(traitHt),
  type = "n", xlab = "",
  ylab = "", yaxt = "n", cex.lab = 2)

axis(2, at = nrow(meanHt):1, labels = rownames(meanHt), las = 1, cex.axis = 1)
points(meanHt[, 'mean'],
  nrow(meanHt):1,
  pch = 16,
  col = "darkslategray4",
  cex = 2)
arrows(meanHt[, "97.5%"], nrow(meanHt):1, meanHt[, "2.5%"], nrow(meanHt):1,
  len = 0, col = "black")
abline(v = 0, lty = 3)
text(-10, 8.5, label = "a)", cex = 1.25)
# LMA
par(mar = c(5, 1, 2, 0.25))

plot(seq(-10,  10,length.out = nrow(meanLMA)), 1:nrow(meanLMA),
  type = "n", xlab = "",
  ylab = "", yaxt = "n", cex.lab = 2)
abline(v = 0, lty = 3)

points(meanLMA[, 'mean'],
  nrow(meanLMA):1,
  pch = 16,
  col = "darkolivegreen",
  cex = 2)
arrows(meanLMA[, "97.5%"], nrow(meanLMA):1, meanLMA[, "2.5%"], nrow(meanLMA):1,
  len = 0, col = "black")
text(-10, 8.5, label = "b)", cex = 1.25)

# DBH
plot(seq(-10,  10,length.out = nrow(meanDBH)), 1:nrow(meanDBH),
  type = "n", xlab = "Estimated change in budburst day",
  ylab = "", yaxt = "n", cex.lab = 1.5)
abline(v = 0, lty = 3)

points(meanDBH[, 'mean'],
  nrow(meanDBH):1,
  pch = 16,
  col = "goldenrod",
  cex = 2)
arrows(meanDBH[, "97.5%"], nrow(meanDBH):1, meanDBH[, "2.5%"], nrow(meanDBH):1,
  len = 0, col = "black")
text(-10, 8.5, label = "c)", cex = 1.25)

# SSD

plot(seq(-10,  10,length.out = nrow(meanSSD)), 1:nrow(meanSSD),
  type = "n", xlab = "",
  ylab = "", yaxt = "n", cex.lab = 2)
abline(v = 0, lty = 3)

points(meanSSD[, 'mean'],
  nrow(meanSSD):1,
  pch = 16,
  col = "maroon",
  cex = 2)
arrows(meanSSD[, "97.5%"], nrow(meanSSD):1, meanSSD[, "2.5%"], nrow(meanSSD):1,
  len = 0, col = "black")
text(-10, 8.5, label = "d)", cex = 1.25)

# CN
plot(seq(-10,  10,length.out = nrow(meanCN)), 1:nrow(meanCN),
  type = "n", xlab = "",
  ylab = "", yaxt = "n")
abline(v = 0, lty = 3, cex.lab = 2)

points(meanCN[, 'mean'],
  nrow(meanCN):1,
  pch = 16,
  col = "purple4",
  cex = 2)
arrows(meanCN[, "97.5%"], nrow(meanCN):1, meanCN[, "2.5%"], nrow(meanCN):1,
  len = 0, col = "black")
text(-10, 8.5, label = "e)", cex = 1.25)
dev.off()
### mu plot of cues

phenoHt <- sumHt[mu_pheno, col4table]
phenoHt <- cbind(phenoHt, bTran_quan, bTranLat_quan)

phenoDBH <- sumDBH[mu_pheno, col4table]
phenoDBH <- cbind(phenoDBH, bTran_quan, bTranLat_quan)

phenoSSD <- sumSSD[mu_pheno, col4table]
phenoSSD <- cbind(phenoSSD, bTran_quan, bTranLat_quan)

phenoCN <- sumCN[mu_pheno, col4table]
phenoCN <- cbind(phenoCN, bTran_quan, bTranLat_quan)

phenoLMA <- sumLMA[mu_pheno, col4table]
phenoLMA <- cbind(phenoLMA, bTran_quan, bTranLat_quan)

### mu plot of betaTraitxCue

meanHt <- sumHt[mu_p, col4table]

rownames(meanHt) = c( 
  #"Root trait intercept", "Lambda",
  "Grand mean",
  "Transect",
  "Transect x latitude",
#  "Budburst slope",
  "Forcing",
  "Chilling",
  "Photoperiod",
  "Trait x forcing",
  "Trait x chill",
  "Trait x photo"
)

load("output/mdl2023/z-scored/lmaDummyIntGrandZ.Rdata")
lmaModelFit <- rstan::extract(mdlLMA)
sumLMA <- summary(mdlLMA)$summary

meanLMA <- sumLMA[mu_params, col4table]

rownames(meanLMA) = c( 
  #"Root trait intercept", "Lambda",
  "Grand mean",
  "Transect",
  "Transect by latitude",
  #  "Budburst slope",
  "Forcing",
  "Chilling",
  "Photoperiod",
  "Trait-forcing effect",
  "Trait-chilling effect",
  "Trait-photoperiod effect"
)

load("output/mdl2023/z-scored/dbhDummyIntGrandZ.Rdata")
dbhModelFit <- rstan::extract(mdlDBH)
sumDBH <- summary(mdlDBH)$summary

meanDBH <- sumDBH[mu_params, col4table]

rownames(meanDBH) = c( 
  #"Root trait intercept", "Lambda",
  "Grand mean",
  "Transect",
  "Transect by latitude",
  #  "Budburst slope",
  "Forcing",
  "Chilling",
  "Photoperiod",
  "Trait-forcing effect",
  "Trait-chilling effect",
  "Trait-photoperiod effect"
)

load("output/mdl2023/z-scored/ssdDummyIntGrandZ.Rdata")
ssdModelFit <- rstan::extract(mdlSSD)
sumSSD <- summary(mdlSSD)$summary

meanSSD <- sumSSD[mu_params, col4table]

rownames(meanSSD) = c( 
  #"Root trait intercept", "Lambda",
  "Grand mean",
  "Transect",
  "Transect by latitude",
  #  "Budburst slope",
  "Forcing",
  "Chilling",
  "Photoperiod",
  "Trait-forcing effect",
  "Trait-chilling effect",
  "Trait-photoperiod effect"
)

load("output/mdl2023/z-scored/cnDummyIntGrandZ.Rdata")
cnModelFit <- rstan::extract(mdlCN)
sumCN <- summary(mdlCN)$summary

meanCN <- sumCN[mu_params, col4table]

rownames(meanCN) = c( 
  #"Root trait intercept", "Lambda",
  "Grand mean",
  "Transect",
  "Transect by latitude",
  #  "Budburst slope",
  "Forcing",
  "Chilling",
  "Photoperiod",
  "Trait-forcing effect",
  "Trait-chilling effect",
  "Trait-photoperiod effect"
)

pdf(file.path( "../figures/traitMuPlot.pdf"), width = 15, height = 3)
par(mfrow = c(1,5), mar = c(5, 8, 2, 0.05))
# Upper panel: bud burst
plot(seq(-10,  10,length.out = nrow(meanHt)), 1:nrow(meanHt),
     type = "n", xlab = "",
     ylab = "", yaxt = "n", cex.lab = 2)

axis(2, at = nrow(meanHt):1, labels = rownames(meanHt), las = 1, cex.axis = 1)
points(meanHt[, 'mean'],
       nrow(meanHt):1,
       pch = 16,
       col = "darkslategray4",
       cex = 2)
arrows(meanHt[, "97.5%"], nrow(meanHt):1, meanHt[, "2.5%"], nrow(meanHt):1,
       len = 0, col = "black")
abline(v = 0, lty = 3)
text(-10, 8.5, label = "a)", cex = 1.25)
# LMA
par(mar = c(5, 1, 2, 0.25))

plot(seq(-10,  10,length.out = nrow(meanLMA)), 1:nrow(meanLMA),
     type = "n", xlab = "",
     ylab = "", yaxt = "n", cex.lab = 2)
abline(v = 0, lty = 3)

points(meanLMA[, 'mean'],
       nrow(meanLMA):1,
       pch = 16,
       col = "darkolivegreen",
       cex = 2)
arrows(meanLMA[, "97.5%"], nrow(meanLMA):1, meanLMA[, "2.5%"], nrow(meanLMA):1,
       len = 0, col = "black")
text(-10, 8.5, label = "b)", cex = 1.25)

# DBH
plot(seq(-10,  10,length.out = nrow(meanDBH)), 1:nrow(meanDBH),
     type = "n", xlab = "Estimated change in budburst day",
     ylab = "", yaxt = "n", cex.lab = 1.5)
abline(v = 0, lty = 3)

points(meanDBH[, 'mean'],
       nrow(meanDBH):1,
       pch = 16,
       col = "goldenrod",
       cex = 2)
arrows(meanDBH[, "97.5%"], nrow(meanDBH):1, meanDBH[, "2.5%"], nrow(meanDBH):1,
       len = 0, col = "black")
text(-10, 8.5, label = "c)", cex = 1.25)

# SSD

plot(seq(-10,  10,length.out = nrow(meanSSD)), 1:nrow(meanSSD),
     type = "n", xlab = "",
     ylab = "", yaxt = "n", cex.lab = 2)
abline(v = 0, lty = 3)

points(meanSSD[, 'mean'],
       nrow(meanSSD):1,
       pch = 16,
       col = "maroon",
       cex = 2)
arrows(meanSSD[, "97.5%"], nrow(meanSSD):1, meanSSD[, "2.5%"], nrow(meanSSD):1,
       len = 0, col = "black")
text(-10, 8.5, label = "d)", cex = 1.25)

# CN
plot(seq(-10,  10,length.out = nrow(meanCN)), 1:nrow(meanCN),
     type = "n", xlab = "",
     ylab = "", yaxt = "n")
abline(v = 0, lty = 3, cex.lab = 2)

points(meanCN[, 'mean'],
       nrow(meanCN):1,
       pch = 16,
       col = "purple4",
       cex = 2)
arrows(meanCN[, "97.5%"], nrow(meanCN):1, meanCN[, "2.5%"], nrow(meanCN):1,
       len = 0, col = "black")
text(-10, 8.5, label = "e)", cex = 1.25)
dev.off()
