# Started May 20, 2022 by Deirdre

# Prior predictive checks should be done after you have thought about what plots you make with your data and what remotely viable priors should be. 
#1. Simulate the data ~ 10000 
#2. Look at plots of your simulated data and see if the range of values makes any sense; ie none are in the millions, all on a scale of what might be viable

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

Nrep <- 15 # rep per trait
Ntransect <- 2
Nspp <- 40 # number of species 

# First making a data frame for the test trait data
Ntrt <- Nspp * Ntransect * Nrep # total number of traits observations
Ntrt

n_spec <-Nspp # same number of species as teh traits section fo the model 
#Number of repeat observations per species
nRep <- 15
#Overall number of pbservations (rows)
Nph <- n_spec * nRep # for phenology model

#Prior Predictive Check (Run 1000 times and plot results)
#------------------------------------------------------------

#Number fo prior check itterations 
nRepPrior <- 300

# ppc for traits portion
priorCheckTrait <- data.frame(matrix(NA, Ntrt*nRepPrior, 3))
names(priorCheckTrait) <- c("simRep", "rep", "species")
priorCheckTrait$simRep <- rep(1:nRepPrior, each = Ntrt)
priorCheckTrait$rep <- rep(1:Ntrt, times = nRepPrior)
priorCheckTrait$species <- rep(1:Nspp, each = nRepPrior)
priorCheckTrait$transect <- rep(1:Ntransect, each = nRepPrior)

# #Make this the name of the full vector of ht per species values - alphaTraitSp 
# priorCheckTrait$alphaTraitSp <-  rep(rep(trt.dat$mu_grand_sp, times = nRepPrior))
# height <- trtPheno[complete.cases(trtPheno$ht),]
# 
# specieslist <- sort(unique(trtPheno$species))
# sitelist <- sort(unique(trtPheno$transect))

specieslist <- seq(1,Nspp, 1)
sitelist <- seq(1,Ntransect, 1)

ht.data <- list(n_spec = length(specieslist),
                n_site = length(sitelist),
                prior_mu_grand_mu = 20,
                prior_mu_grand_sigma = 10,
                prior_sigma_sp_mu = 4,
                prior_sigma_sp_sigma = 5,
                prior_sigma_site_mu = 2,
                prior_sigma_site_sigma = 5,
                prior_sigma_traity_mu = 3,
                prior_sigma_traity_sigma = 5,
                ## Phenology
                Nph = Nph,
                prior_muForceSp_mu = -15,
                prior_muForceSp_sigma = 10, #wider
                prior_muChillSp_mu = -15,
                prior_muChillSp_sigma = 10,#wider
                prior_muPhotoSp_mu = -15,
                prior_muPhotoSp_sigma = 10,#wider
                prior_muPhenoSp_mu = 40,
                prior_muPhenoSp_sigma = 10,#wider
                prior_sigmaForceSp_mu = 5,
                prior_sigmaForceSp_sigma = 5,
                prior_sigmaChillSp_mu = 5,#wider
                prior_sigmaChillSp_sigma = 5, #wider
                prior_sigmaPhotoSp_mu = 5,
                prior_sigmaPhotoSp_sigma = 5,
                prior_sigmaPhenoSp_mu = 5, #wider
                prior_sigmaPhenoSp_sigma = 5, #wider
                prior_betaTraitxForce_mu = 0,
                prior_betaTraitxForce_sigma = 1,
                prior_betaTraitxChill_mu = 0,
                prior_betaTraitxChill_sigma = 1,
                prior_betaTraitxPhoto_mu = 0,
                prior_betaTraitxPhoto_sigma = 1,
                prior_sigmaphenoy_mu = 10,
                prior_sigmaphenoy_sigma = 5 #wider
)



for (ir in 1:nRepPrior){
  # Parameter Values
  #ir <- 1
  
  muGrand <- rtruncnorm(1, a = 0, mean = ht.data$prior_mu_grand_mu, sd = ht.data$prior_mu_grand_sigma)
  sigmaSp <- rtruncnorm(1, a = 0, mean = ht.data$prior_sigma_sp_mu, sd = ht.data$prior_sigma_sp_sigma)
  sigmatransect <- rtruncnorm(1, a = 0, mean = ht.data$prior_sigma_site_mu, sd = ht.data$prior_sigma_site_sigma)
  
  alphaTraitSp <- rnorm(Nspp, 0, sigmaSp)
  priorCheckTrait$alphaTraitSp[priorCheckTrait$simRep == ir] <- rep(alphaTraitSp, each = nRep)
  
  muSp <- rnorm(Nspp, 0, sigmaSp)
  priorCheckTrait$muSp[priorCheckTrait$simRep == ir] <- rep(muSp, each = nRep)
  
  mutransect <- rnorm(Ntransect, 0, sigmatransect)
  priorCheckTrait$mutransect[priorCheckTrait$simRep == ir] <- rep(mutransect, each = nRep)
  
  #general varience
  priorCheckTrait$sigmaTrait_y[priorCheckTrait$simRep == ir] <- rnorm(ht.data$prior_sigma_traity_mu, ht.data$prior_sigma_traity_sigma)
  priorCheckTrait$e[priorCheckTrait$simRep == ir] <- rnorm(Ntrt, 0, priorCheckTrait$sigmaTrait_y)
  
  priorCheckTrait$yTraiti <- muGrand + priorCheckTrait$muSp + priorCheckTrait$mutransect + priorCheckTrait$e
}# end simulating new priors, from here vectorize code

#Final values
priorCheckTrait$muGrandSp <- muGrand + priorCheckTrait$muSp

priorCheckTraityTraiti <- priorCheckTrait[complete.cases(priorCheckTrait$yTraiti),]

png("figures/density_Trait_Prior_joint_ht.png")
plot(density(priorCheckTraityTraiti$yTraiti))
dev.off()

png("figures/GrandSp_PlotPrior_joint_ht.png")
plot(priorCheckTrait$yTraiti ~ priorCheckTrait$muGrandSp, xlab = "muGrandSp", ylab = "Trait")
dev.off()

png("figures/MuSp_PlotPrior_joint_ht.png")
plot(priorCheckTrait$yTraiti ~ priorCheckTrait$muSp, xlab = "MuSp", ylab = "Trait")
dev.off()

png("figures/Mutransect_PlotPrior_joint_ht.png")
plot(priorCheckTrait$yTraiti ~ priorCheckTrait$mutransect, xlab = "Mutransect", ylab = "Trait")
dev.off()
#####################################################################################

#Make a data frame for input simulation data
priorCheckPheno <- data.frame(matrix(NA, Nph*nRepPrior, 3))
names(priorCheckPheno) <- c("simRep","rep","species")

priorCheckPheno$simRep <- rep(1:nRepPrior, each = Nph)

priorCheckPheno$rep <- rep(c(1:Nph), times = nRepPrior)
priorCheckPheno$species <- rep(rep(c(1:n_spec), each = nRep), times = nRepPrior)

#Simulate SLA data per species
muGrandSp <- muGrand + muSp
#Make this the name of the full vector of sla per species values - alphaTraitSp 
priorCheckPheno$alphaTraitSp <-  rep(rep(muGrandSp, times = nRepPrior)) # use the mean mu_grand_sp, grand mean + transect, not transect

#Simulate cues (z scored)
priorCheckPheno$forcei <- rnorm(Nph, 1, 1)
priorCheckPheno$photoi <- rnorm(Nph, 1, 1) # less photoperiod 
priorCheckPheno$chilli <- rnorm(Nph, 1, 1) #more chilling

for (ir in 1:nRepPrior){
  # Parameter Values
  # ir <- 1
  
  #Species means
  sigmaPhenoSp <- rtruncnorm(1, a = 0, mean = ht.data$prior_sigmaPhenoSp_mu, sd = ht.data$prior_sigmaPhenoSp_sigma)
  muPhenoSp <- rnorm(1, ht.data$prior_muPhenoSp_mu, ht.data$prior_muPhenoSp_sigma)
  alphaPhenoSp <- rnorm(n_spec, muPhenoSp, sigmaPhenoSp)
  priorCheckPheno$alphaPhenoSp[priorCheckPheno$simRep == ir] <- rep(alphaPhenoSp, each = nRep)
  
  #Cue effects
  priorCheckPheno$betaTraitxForce[priorCheckPheno$simRep == ir] <- rnorm(1,ht.data$prior_betaTraitxForce_mu,ht.data$prior_betaTraitxForce_sigma)
  priorCheckPheno$betaTraitxPhoto[priorCheckPheno$simRep == ir] <- rnorm(1,ht.data$prior_betaTraitxPhoto_mu,ht.data$prior_betaTraitxPhoto_sigma)
  priorCheckPheno$betaTraitxChill[priorCheckPheno$simRep == ir] <- rnorm(1,ht.data$prior_betaTraitxChill_mu,ht.data$prior_betaTraitxChill_sigma)
  
  #Species level slopes sans trait data
  muForceSp <- rnorm(1,ht.data$prior_muForceSp_mu,  ht.data$prior_muForceSp_sigma)
  sigmaForceSp <- rtruncnorm(1, a = 0, mean = ht.data$prior_sigmaForceSp_mu,sd = ht.data$prior_sigmaForceSp_sigma)
  alphaForceSp <- rnorm(n_spec, muForceSp, sigmaForceSp)
  priorCheckPheno$alphaForceSp[priorCheckPheno$simRep == ir] <- rep(alphaForceSp, each = nRep)
  
  muPhotoSp <- rnorm(1, ht.data$prior_muPhotoSp_mu, ht.data$prior_muPhotoSp_sigma)
  sigmaPhotoSp <- rtruncnorm(1, a = 0, mean = ht.data$prior_sigmaPhotoSp_mu, sd = ht.data$prior_sigmaPhotoSp_sigma )
  alphaPhotoSp <- rnorm(n_spec, muPhotoSp, sigmaPhotoSp)
  priorCheckPheno$alphaPhotoSp[priorCheckPheno$simRep == ir] <- rep(alphaPhotoSp, each = nRep)
  
  muChillSp <-  rnorm(1,ht.data$prior_sigmaChillSp_mu,ht.data$prior_sigmaChillSp_sigma)
  sigmaChillSp <- rtruncnorm(1, a = 0, mean = ht.data$prior_sigmaChillSp_mu,sd = ht.data$prior_sigmaChillSp_sigma)
  alphaChillSp <- rnorm(n_spec, muChillSp, sigmaChillSp)
  priorCheckPheno$alphaChillSp[priorCheckPheno$simRep == ir] <- rep(alphaChillSp, each = nRep)
  
  
  #general varience
  priorCheckPheno$sigmapheno_y[priorCheckPheno$simRep == ir] <- rtruncnorm(ht.data$prior_sigma_traity_mu,  a = 0, ht.data$prior_sigma_traity_sigma)
  priorCheckPheno$e[priorCheckPheno$simRep == ir] <- rnorm(Nph, 0, priorCheckPheno$sigmapheno_y)
  
}# end simulating new priors, from here vectorize code
#slopes for each cue, combining trait and non-trait aspect of the slope.


priorCheckPheno$betaForceSp <- priorCheckPheno$alphaForceSp + priorCheckPheno$betaTraitxForce *  priorCheckPheno$alphaTraitSp

priorCheckPheno$betaPhotoSp <-  priorCheckPheno$alphaPhotoSp + priorCheckPheno$betaTraitxPhoto * priorCheckPheno$alphaTraitSp

priorCheckPheno$betaChillSp <-  priorCheckPheno$alphaChillSp + priorCheckPheno$betaTraitxChill * priorCheckPheno$alphaTraitSp

#Run full model to get mean simulated y values
priorCheckPheno$yMu <-  priorCheckPheno$alphaPhenoSp +  priorCheckPheno$betaForceSp * priorCheckPheno$forcei +  priorCheckPheno$betaPhotoSp* priorCheckPheno$photoi + priorCheckPheno$betaChillSp * priorCheckPheno$chilli

#Final values
priorCheckPheno$yPhenoi <- priorCheckPheno$yMu + priorCheckPheno$e

head(priorCheckPheno)

plot(priorCheckPheno$betaForceSp ~ priorCheckPheno$alphaTraitSp )
priorCheckPheno_posF <- priorCheckPheno[priorCheckPheno$betaForceSp > 0,]
plot(priorCheckPheno_posF$betaForceSp ~ priorCheckPheno_posF$alphaTraitSp )

png("figures/densityYPrior_joint_ht.png")
plot(density(priorCheckPheno$yPhenoi))
dev.off()

png("figures/photoPlotPrior_joint_ht.png")
plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$photoi, xlab = "Photoperiod", ylab = "Phenological Date")
dev.off()

png("figures/forcingPlotPrior_joint_ht.png")
plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$forcei, xlab = "Forcing", ylab = "Phenological Date")
dev.off()

png("figures/chillingPlotPrior_joint_ht.png")
plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$chilli, xlab = "Chillina", ylab = "Phenological Date")
dev.off()

