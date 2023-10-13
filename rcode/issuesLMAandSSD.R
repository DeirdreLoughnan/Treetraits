# Started Friday Oct 13 ~~ oooh sppoooky ~~

# aim of this code is to try and figure out why the lma and ssd models are not including the trait data
# interesting that both these traits are on very similar scales --- ie very small values (0.01-0.06)


## Here is my prposed plan of action
# 1. Posterior predictive checks
# 2. Run a quick lmer model of just the trait part ---> is it an issue with Stan?
# 3. Run it quickly on the trait only part ---> is it an issue with the priors?

## 1. Posterior predictive checks ---- starting with lma model
rm(list=ls())
options(stringsAsFactors = FALSE)

library(tidyr)
library(plyr)
library(dplyr)
library(reshape2)
library(rstan)
library(bayesplot)# nice posterior check plots 
library(shinystan)
library(stringr)
library(truncnorm)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  else{
  setwd("/home/deirdre/pheno_bc") # for midge
}

load("output/lmaDummyInt.Rdata")
post<- rstan::extract(mdl)

postLMA <- data.frame(post)

sumer <- data.frame(summary(mdl)$summary)


stanBMuSp <- sumer[1:47,]
pdf("figures/stanBmuSpesti.pdf")
hist(stanBMuSp$mean)
dev.off()
## pairs plot:
# png("figures/simulatedPairs.png")
pdf("figures/lmaPairsTrt.pdf")
pairs(mdl, pars = c("b_tranlat", "b_tranE", "muSp", "sigma_sp", "sigma_traity", "lp__"))
dev.off()

pdf("figures/lmaPairsPheno.pdf")
pairs(mdl, pars = c("muForceSp", "muChillSp", "muPhotoSp", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
dev.off()

## Compare histograms:

par(mfrow=c(1,1))
##Compare results to simulated values
# Trait parameters
hist(rnorm(1000, 0, 0.01), ylim = c(0,800))
hist(postLMA$b_tranlat, main = "b_tranLat", add = T, col = "pink")
abline(v = 0, col="red", lwd=3, lty=2)

hist(rnorm(1000, 0, 0.05), ylim = c(0,800), xlim = c(-1,1))
hist(postLMA$b_tranE, main = "b_tranE", add = T, col = "pink")
abline(v = 0, col="red", lwd=3, lty=2)

hist(rnorm(1000, 0, 0.05), ylim = c(0,800), xlim = c(-1,1))
hist(postLMA$muSp, main = "muSp", add = T, col = "pink")
abline(v = 0.05, col="red", lwd=3, lty=2)

hist(rnorm(1000, 0, 0.5), ylim = c(0,800), xlim = c(-1,1))
hist(postLMA$sigma_sp, main = "sigma_sp", add = T, col = "pink")
abline(v = 0, col="red", lwd=3, lty=2)

hist(rnorm(1000, 0, 0.5), ylim = c(0,800), xlim = c(-1,1))
hist(postLMA$sigma_traity, main = "sigma_triaty", add = T, col = "pink")
abline(v = 0, col="red", lwd=3, lty=2)

# phenology parameters
hist(postLMA$muPhenoSp, main = "muPhenoSp", add = T, col = "pink")
abline(v = 40, col="red", lwd=3, lty=2)

hist(postLMA$muForceSp, main = "muForceSp", add = T, col = "pink")
abline(v = -15, col="red", lwd=3, lty=2)

hist(postLMA$muChillSp, main = "muChillSp", add = T, col = "pink")
abline(v = -15, col="red", lwd=3, lty=2)

hist(postLMA$muPhotoSp, main = "muPhotoSp", add = T, col = "pink")
abline(v = -15, col="red", lwd=3, lty=2)

hist(postLMA$betaTraitxForce, main = "betaTraitxForce", add = T, col = "pink")
abline(v = 0, col="red", lwd=3, lty=2)
#
hist(postLMA$betaTraitxChill, main = "betaTraitxChill", add = T, col = "pink")
abline(v = 0, col="red", lwd=3, lty=2)
#
hist(postLMA$betaTraitxPhoto, main = "betaTraitxPhoto"), add = T, col = "pink")
abline(v = 0, col="red", lwd=3, lty=2)

hist(postLMA$sigmaPhenoSp, main = "sigmapheno_y", add = T, col = "pink")
abline(v = 5, col="red", lwd=3, lty=2)

hist(postLMA$sigmapheno_y, main = "sigmapheno_y", add = T, col = "pink")
abline(v = 10, col="red", lwd=3, lty=2)

hist(postLMA$sigmaChillSp, main = "sigmaChillSp", add = T, col = "pink")
abline(v = 5, col="red", lwd=3, lty=2)

hist(postLMA$sigmaForceSp, main = "sigmaForceSp", add = T, col = "pink")
abline(v = 5, col="red", lwd=3, lty=2)

hist(postLMA$sigmaPhotoSp, main = "sigmaPhotoSp", add = T, col = "pink")
abline(v = 5, col="red", lwd=3, lty=2)

# png("figures/simulatedPairs.png")
#pairs(mdl.LMA, pars = c("muForceSp", "muChillSp", "muPhotoSp", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
# dev.off()

#############################################################################################

## Simulate data using prior values:


