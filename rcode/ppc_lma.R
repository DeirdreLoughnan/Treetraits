#Started May 6, 2022 by Deirdre

# Aim of this code is to perform several ppc for the trait pheno model for NAm tree traits
# This model includes transect as a dummy variable but otherwise is the same as the traitors model

rm(list=ls())
options(stringsAsFactors = FALSE)

library(tidyr)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(rstan)
library(bayesplot)# nice posterior check plots 
library(shinystan)
library(stringr)
library(truncnorm)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/Treetraits")
}  else{
  setwd("/home/deirdre/Treetraits") # for midge
}

trtPheno <- read.csv("input/trtPhenoDummy.csv")
pheno.t <- read.csv("input/phenoDataWChill.csv")
# #########################################################

load("output/mdl2023/z-scored/lmaDummyIntGrandZ.Rdata")
postLMA<- data.frame(rstan::extract(mdlLMA))

write.csv(postLMA, "output/postLMA.csv", row.names = F)

cueTrt <- postLMA[, colnames(postLMA) %in% c("mu_grand","b_tranE", "b_tranlat","sigma_sp", "sigma_traity")]

mcmc_intervals(cueTrt) 

cueEffects <- postLMA[, colnames(postLMA) %in% c("muPhenoSp", "muForceSp", "muChillSp", "muPhotoSp", "sigmapheno_y")]

mcmc_intervals(cueEffects) 
  # + theme_classic() + 
  # labs(title = "main intercept, cue slopes and general error")

###############

postLMA_muSp <- postLMA[,colnames(postLMA) %in% grep( "muSp", colnames(postLMA), value = TRUE)]
colnames(postLMA_muSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postLMA_muSp) + 
  geom_vline(xintercept = mean(postLMA_muSp), linetype="dotted", color = "grey")  +
  theme_classic() 
  
postLMA_alpaForceSp <- postLMA[,colnames(postLMA) %in% grep( "alphaForceSp", colnames(postLMA), value = TRUE)]
colnames(postLMA_alpaForceSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postLMA_alpaForceSp) + 
  geom_vline(xintercept = mean(postLMA$muForceSp), linetype="dotted", color = "grey")  +
  theme_classic() + 
  labs(subtitle = paste0("Mean muForceSp was ", round(mean(postLMA$muForceSp),3)),
       title = "muForceSp - species forcing slopes no trait")

postLMA_betaForceSp <- postLMA[,colnames(postLMA) %in% grep( "betaForceSp", colnames(postLMA), value = TRUE)]
colnames(postLMA_betaForceSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postLMA_betaForceSp) + 
  theme_classic() + 
  labs(title = "betaForceSp - Species forcing slopes with trait value")

postLMA_alphaChillSp <- postLMA[,colnames(postLMA) %in% grep( "alphaChillSp", colnames(postLMA), value = TRUE)]
colnames(postLMA_alphaChillSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postLMA_alphaChillSp) + 
  geom_vline(xintercept = mean(postLMA$muChillSp), linetype="dotted", color = "grey")  +
  theme_classic() + 
  labs(subtitle = paste0("Mean muChillSp was ", round(mean(postLMA$muChillSp),3)),
       title = "alphaChillSp - Species chill slopes no trait")

postLMA_betaChillSp <- postLMA[,colnames(postLMA) %in% grep( "betaChillSp", colnames(postLMA), value = TRUE)]
colnames(postLMA_betaChillSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postLMA_betaChillSp) + 
 # theme_classic() + 
  labs(title = "betaChillSp - Species chilling slopes with trait value")

postLMA_alphaPhotoSp <- postLMA[,colnames(postLMA) %in% grep( "alphaPhotoSp", colnames(postLMA), value = TRUE)]
colnames(postLMA_alphaPhotoSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postLMA_alphaPhotoSp) + 
  geom_vline(xintercept = mean(postLMA$muPhotoSp), linetype="dotted", color = "grey")  +
  theme_classic() + 
  labs(subtitle = paste0("Mean muPhotoSp was ", round(mean(postLMA$muPhotoSp),3)),
       title = "muPhotoSp - Species photo period slopes no trait")

postLMA_betaPhotoSp <- postLMA[,colnames(postLMA) %in% grep( "betaPhotoSp", colnames(postLMA), value = TRUE)]
colnames(postLMA_betaPhotoSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postLMA_betaPhotoSp) + 
  theme_classic() + 
  labs(title = "betaPhotoSp - Species photoperiod slopes with trait value")

postLMA_betaTraitx <- postLMA[,colnames(postLMA) %in% grep( "betaTraitx", colnames(postLMA), value = TRUE)]

mcmc_intervals(postLMA_betaTraitx) + 
  theme_classic() + 
  labs(title = "effect's of traits on cue slopes")


#######################################################

#png("figures/simPosteriorHist.png")
# priors means

mu.grand <- 0.5 # the grand mean of the lma model
sigma.species <- 1 # we want to keep the variaiton across spp. high
sigma.transect <- 5
b.tranE <- 0
b.tranlat <- 0
sigma_traity <- 5

sigmaPhenoSp <- 5
muPhenoSp <- 40
betaTraitxForce <- 0
betaTraitxPhoto <- 0
betaTraitxChill <- 0
muForceSp <- -15
sigmaForceSp <- 5
muPhotoSp <- -15
sigmaPhotoSp <- 5
muChillSp <- -15
sigmaChillSp <- 5
sigmapheno_y <- 10


postLMA <- rstan::extract(mdlLMA)


##Compare results to simulated values
pdf("figures/ppcFig/histCompLMATrait.pdf")
par(mfrow=c(2,3))
hist(postLMA$mu_grand, main = paste("muPhenoSp is " , signif(mu.grand,3), sep = ""), xlim = c(0,0.15), col=rgb(0,0,1,1/4))
hist(rnorm(1000, 0.05,0.02), col=rgb(1,0,1,1/4), add = T)
abline(v = mu.grand, col="red", lwd=3, lty=2)

hist(postLMA$b_tranE, main = paste("muForceSp is " , signif(muForceSp,3), sep = ""),  xlim = c(-.5,.5),col=rgb(0,0,1,1/4))
hist(rnorm(1000, 0,.1), col=rgb(1,0,1,1/4), add = T)
abline(v = b.tranE, col="red", lwd=3, lty=2)

hist(postLMA$b_tranlat, main = paste("muChillSp is " , signif(muChillSp,3), sep = ""), col=rgb(0,0,1,1/4),  xlim = c(-.5,.5))
hist(rnorm(1000, 0,0.1), col=rgb(1,0,1,1/4), add = T)
abline(v = b.tranlat, col="red", lwd=3, lty=2)

hist(postLMA$sigma_sp, main = paste("muPhotoSp is " , signif(muPhotoSp,3), sep = ""), col=rgb(0,0,1,1/4), xlim =c(-0.10,0.1))
hist(rnorm(1000, 0.01,0.05), col=rgb(1,0,1,1/4), add = T)
abline(v = sigma.species, col="red", lwd=3, lty=2)

hist(postLMA$sigma_traity, main = paste("sigmapheno_y is " , signif(sigmapheno_y,3), sep = ""), col=rgb(0,0,1,1/4),  xlim = c(-0.5,0.5))
hist(rnorm(1000, 0.05,0.1), col=rgb(1,0,1,1/4), add = T)
abline(v = sigmapheno_y, col="red", lwd=3, lty=2)

dev.off()

pdf("figures/ppcFig/histCompLMAPheno.pdf")
par(mfrow=c(1,1))
hist(postLMA$muPhenoSp, main = "muPhenoSp" , col=rgb(0,0,1,1/4), xlim = c(5,70))
hist(rnorm(1000, 40,10), col=rgb(1,0,1,1/4), add = T)
abline(v = muPhenoSp, col="red", lwd=3, lty=2)

hist(postLMA$muForceSp, main = paste("muForceSp is " , signif(muForceSp,3), sep = ""),  xlim = c(-75,75),col=rgb(0,0,1,1/4))
hist(rnorm(1000, -15,10), col=rgb(1,0,1,1/4), add = T)
abline(v = muForceSp, col="red", lwd=3, lty=2)

hist(postLMA$muChillSp, main = paste("muChillSp is " , signif(muChillSp,3), sep = ""),  xlim = c(-75,75), col=rgb(0,0,1,1/4))
hist(rnorm(1000, -15,10), col=rgb(1,0,1,1/4), add = T)
abline(v = muChillSp, col="red", lwd=3, lty=2)

hist(postLMA$muPhotoSp, main = paste("muPhotoSp is " , signif(muPhotoSp,3), sep = ""), col=rgb(0,0,1,1/4),  xlim = c(-75,75))
hist(rnorm(1000, -15,10), col=rgb(1,0,1,1/4), add = T)
abline(v = muPhotoSp, col="red", lwd=3, lty=2)

hist(postLMA$sigmapheno_y, main = paste("sigmapheno_y is " , signif(sigmapheno_y,3), sep = ""), col=rgb(0,0,1,1/4),  xlim = c(-10,30))
hist(rnorm(1000, 10,1), col=rgb(1,0,1,1/4), add = T)
abline(v = sigmapheno_y, col="red", lwd=3, lty=2)

hist(postLMA$betaTraitxForce, main = "betaTraitxForce", col=rgb(0,0,1,1/4),  xlim = c(-5,5))
hist(rnorm(1000, 0,1), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)
#
hist(postLMA$betaTraitxChill, main ="betaTraitxChill", col=rgb(0,0,1,1/4),  xlim = c(-5,5))
hist(rnorm(1000, 0,1), col=rgb(1,0,1,1/4), add = T)
abline(v =0, col="red", lwd=3, lty=2)
#
hist(postLMA$betaTraitxPhoto, main = paste("betaTraitxPhoto is " , signif(betaTraitxPhoto,3), sep = ""), col=rgb(0,0,1,1/4),  xlim = c(-5,5))
hist(rnorm(1000, 0,2), col=rgb(1,0,1,1/4), add = T)
abline(v = betaTraitxPhoto, col="red", lwd=3, lty=2)

hist(postLMA$sigmaChillSp, main = paste("sigmaChillSp is " , signif(sigmaChillSp,3), sep = ""), col=rgb(0,0,1,1/4), xlim = c(-10,20))
hist(rnorm(1000, 5,5), col=rgb(1,0,1,1/4), add = T)
abline(v = sigmaChillSp, col="red", lwd=3, lty=2)

hist(postLMA$sigmaForceSp, main = paste("sigmaForceSp is " , signif(sigmaForceSp,3), sep = ""), col=rgb(0,0,1,1/4), xlim = c(-15,15))
hist(rnorm(1000, 5,5), col=rgb(1,0,1,1/4), add = T)
abline(v = sigmaForceSp, col="red", lwd=3, lty=2)

hist(postLMA$sigmaPhotoSp, main = paste("sigmaPhotoSp is " , signif(sigmaPhotoSp,3), sep = ""), col=rgb(0,0,1,1/4), xlim = c(-15,15))
hist(rnorm(1000, 5,5), col=rgb(1,0,1,1/4), add = T)
abline(v = sigmaPhotoSp, col="red", lwd=3, lty=2)

dev.off()
#pairs(mdl.LMA, pars = c("muForceSp", "muChillSp", "muPhotoSp", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
# dev.off()

  #Prior Predictive Check (Run 1000 times and plot results)
  #------------------------------------------------------------
  
  #Number fo prior check itterations 


  nRepPrior <- 100

  Nrep <- 15 # rep per trait
  nRep <- 15
  Npop <- 8
  Ntran <- 2
  Nspp <- 40 # number of species 

  Ntrt <- Nspp * Npop * Nrep# total number of traits observations
  Ntrt

  # ppc for traits portion
  priorCheckTrait <- data.frame(matrix(NA, Ntrt*nRepPrior, 3))
  names(priorCheckTrait) <- c("simRep", "rep", "species")
  priorCheckTrait$simRep <- rep(1:nRepPrior, each = Ntrt)
  priorCheckTrait$rep <- rep(1:Ntrt, times = nRepPrior)
  priorCheckTrait$species <- rep(1:Nspp, each = nRep)
  priorCheckTrait$pop <- rep(1:Npop, each = Nspp*Nrep)
  priorCheckTrait$tran <- rep(1:Ntran, each = 4*Nrep*Nspp)
  
  priorCheckTrait$dumE <- priorCheckTrait$tran 
  priorCheckTrait$dumE[priorCheckTrait$dumE == "1"] <- 0
  priorCheckTrait$dumE[priorCheckTrait$dumE == "2"] <- 1
  priorCheckTrait$dumE <- as.numeric(priorCheckTrait$dumE)

  
  # Add latitudes
  priorCheckTrait$lat <- priorCheckTrait$pop
  priorCheckTrait$lat[priorCheckTrait$lat == "1"] <- 54.7824
  priorCheckTrait$lat[priorCheckTrait$lat == "2"] <- 52.1417
  priorCheckTrait$lat[priorCheckTrait$lat == "3"] <- 50.8837
  priorCheckTrait$lat[priorCheckTrait$lat == "4"] <- 49.0646
  priorCheckTrait$lat[priorCheckTrait$lat == "5"] <- 45.9310
  priorCheckTrait$lat[priorCheckTrait$lat == "6"] <- 44.92466697
  priorCheckTrait$lat[priorCheckTrait$lat == "7"] <- 43.99837498
  priorCheckTrait$lat[priorCheckTrait$lat == "8"] <- 42.5315
  priorCheckTrait$lat <- as.numeric(priorCheckTrait$lat)

  priorCheckTrait$lat.z <- (priorCheckTrait$lat-mean(priorCheckTrait$lat,na.rm=TRUE))/(sd(priorCheckTrait$lat,na.rm=TRUE))

  #Make this the name of the full vector of sla per species values - alphaTraitSp 
  #priorCheckTrait$alphaTraitSp <-  rep(rep(trt.dat$mu_grand_sp, times = nRepPrior))
  leafMass <- trtPheno[complete.cases(trtPheno$lma),]
  
  specieslist <- sort(unique(trtPheno$species))
  sitelist <- sort(unique(trtPheno$transect))

  mu.grand <- 0.5 # the grand mean of the lma model
  sigma.species <- 1 # we want to keep the variation across spp. high
  sigma.transect <- 5
  b.tranE <- 0
  b.tranlat <- 0
  sigma_traity <- 5
  sigmapheno_y <-10
  
lma.data <- list(yTraiti = leafMass$lma,
                   N = nrow(leafMass),
                   n_spec = length(specieslist),
                   trait_species = as.numeric(as.factor(leafMass$species)),
                   n_tran = length(unique(leafMass$transect)),
                   tranE = leafMass$transect,
                   lati <- leafMass$latitude,

                   prior_muGrand_mu = -1.5,
                   prior_muGrand_sigma = 5, 
                   prior_sigma_sp_mu = 0,
                   prior_sigma_sp_sigma = 5,
                   prior_b_tranlat_mu = 0,
                   prior_b_tranlat_sigma = 1,
                   prior_b_tranE_mu = 0,
                   prior_b_tranE_sigma = 5,
                   prior_sigma_traity_mu = 0,
                   prior_sigma_traity_sigma = 5,
                   ## Phenology
                   Nph = nrow(pheno.t),
                   phenology_species = as.numeric(as.factor(pheno.t$species)),
                   yPhenoi = pheno.t$bb,
                   forcei = pheno.t$force.z2,
                   chilli = pheno.t$chillport.z2,
                   photoi = pheno.t$photo.z2,
                   prior_muForceSp_mu = 0,
                   prior_muForceSp_sigma = 10, 
                   prior_muChillSp_mu = 0,
                   prior_muChillSp_sigma = 10,
                   prior_muPhotoSp_mu = 0,
                   prior_muPhotoSp_sigma = 10,
                   prior_muPhenoSp_mu = 40,
                   prior_muPhenoSp_sigma = 20,
                   prior_sigmaForceSp_mu = 5,
                   prior_sigmaForceSp_sigma = 5,
                   prior_sigmaChillSp_mu = 5,
                   prior_sigmaChillSp_sigma = 5, 
                   prior_sigmaPhotoSp_mu = 5,
                   prior_sigmaPhotoSp_sigma = 5,
                   prior_sigmaPhenoSp_mu = 5, 
                   prior_sigmaPhenoSp_sigma = 5, 
                   prior_betaTraitxForce_mu = 0,
                   prior_betaTraitxForce_sigma = 3,
                   prior_betaTraitxChill_mu = 0,
                   prior_betaTraitxChill_sigma = 3,
                   prior_betaTraitxPhoto_mu = 0,
                   prior_betaTraitxPhoto_sigma = 3,
                   prior_sigmaphenoy_mu = 10,
                   prior_sigmaphenoy_sigma = 5 
) 
  
  for (ir in 1:nRepPrior){
    # Parameter Values

    
    #muGrand <- rtruncnorm(1, a = 0, mean = lma.data$prior_muGrand_mu, sd = lma.data$prior_muGrand_sigma)
    sigmaSp <- rtruncnorm(1, a = 0, mean = lma.data$prior_sigma_sp_mu, sd = lma.data$prior_sigma_sp_sigma)
    b_tranlat <- rtruncnorm(1, a = 0, mean = lma.data$prior_b_tranlat_mu, sd = lma.data$prior_b_tranlat_sigma)
    b.tranE <- rtruncnorm(1, a = 0, mean = lma.data$prior_b_tranE_mu, sd = lma.data$prior_b_tranE_sigma)
  
    muSp <- rnorm(Nspp, -3, 0.5)
    priorCheckTrait$muSp[priorCheckTrait$simRep == ir] <- rep(muSp, each = nRep)
    
    #general varience
    priorCheckTrait$sigma_traity[priorCheckTrait$simRep == ir] <- rtruncnorm(lma.data$prior_sigma_traity_mu,  a = 0, lma.data$prior_sigma_traity_sigma)
    priorCheckTrait$e[priorCheckTrait$simRep == ir] <- rnorm(Ntrt, 0, lma.data$prior_sigma_traity_mu)
    
    priorCheckTrait$yTraiti <-  priorCheckTrait$muSp + b_tranlat *(priorCheckTrait$lat.z * priorCheckTrait$dumE) + b.tranE * priorCheckTrait$dumE + priorCheckTrait$e
  }# end simulating new priors, from here vectorize code
  
  
  priorCheckTrait$Traiti <- priorCheckTrait[complete.cases(priorCheckTrait$yTraiti),]
  
  #png("figures/density_Trait_Prior_joint_lma.png")

  plot(density(priorCheckTrait$yTraiti))
  #dev.off()
  
  #png("figures/GrandSp_PlotPrior_joint_lma.png")
  plot(priorCheckTrait$yTraiti ~ priorCheckTrait$muSp, xlab = "ytrait", ylab = "Trait")
  #dev.off()
  
  #png("figures/MuSp_PlotPrior_joint_lma.png")
  plot(priorCheckTrait$yTraiti ~ priorCheckTrait$muSp, xlab = "MuSp", ylab = "Trait")
  #dev.off()
  
  #png("figures/Mutransect_PlotPrior_joint_lma.png")
  plot(priorCheckTrait$yTraiti ~ priorCheckTrait$dumE, xlab = "Mutransect", ylab = "Trait")
  #dev.off()
  #####################################################################################
  # Chilling - included as chill portions: 8 values between 45 and 80
  # Forcing included as continuous and accounting for thermoperiodicity

  
  n_spec <- 40 # same number of species as teh traits section fo the model 
  #Number of repeat observations per species
  nRep <- 5
  nChill <- 2 # sm high low, mp high low
  nPhoto <- 2 # high low
  nForce <- 2 # high amd low
  nPop <- 2
  nSites <- 4
  #Overall number of pbservations (rows)
  Nph <- n_spec * nRep*nChill * nForce * nPhoto * nPop*nSites # for phenology model
  
  
  priorCheckPheno <- data.frame(matrix(NA, Nph*nRepPrior, 2))
  names(priorCheckPheno) <- c("simRep","rep")
  priorCheckPheno$simRep <- rep(1:nRepPrior, each = Nph)
  priorCheckPheno$rep <- c(1:nRep)
  priorCheckPheno$species <- rep(c(1:n_spec), each = nRep, times = nRepPrior)
  priorCheckPheno$sites <- rep(1:nSites, each = n_spec * nRep*nChill * nForce * nPhoto * nPop, times = nRepPrior)
  priorCheckPheno$tran <- rep(1:nPop, each = n_spec * nRep*nChill * nForce * nPhoto * nSites, times = nRepPrior)
  
  #Simulate SLA data per species
  muGrandSp <-muSp
  #Make this the name of the full vector of sla per species values - alphaTraitSp 
  priorCheckPheno$alphaTraitSp <-  rep(rep(muGrandSp, each = nRep*nChill * nForce * nPhoto * nPop), times = nRepPrior)

  priorCheckPheno$alphaTraitSp <-  rep(rep(muGrandSp, each = nRep*nChill * nForce * nPhoto * nPop*nSites), times = nRepPrior)
  
  
  forceTrt <- c("LF", "HF")
  priorCheckPheno$force <- rep(forceTrt, each =  nRep, times = n_spec * nChill *  nPhoto * nPop *nRepPrior)
  priorCheckPheno$forcei <- priorCheckPheno$force
  priorCheckPheno$forcei[which(priorCheckPheno$force == "HF" & priorCheckPheno$tran == "1")] <- "15"
  priorCheckPheno$forcei[which(priorCheckPheno$force == "HF" & priorCheckPheno$tran == "2")] <- "13.33"
  priorCheckPheno$forcei[which(priorCheckPheno$force == "LF" & priorCheckPheno$tran == "1")] <- "10"
  priorCheckPheno$forcei[which(priorCheckPheno$force == "LF" & priorCheckPheno$tran == "2")] <- "8.33"
  priorCheckPheno$forcei <- as.numeric(priorCheckPheno$forcei)
  
  photoTrt <- c("LP", "HP")
  priorCheckPheno$photo <- rep(photoTrt, each =  nRep, times = n_spec * nChill *  nForce * nPop *nRepPrior)
  priorCheckPheno$photoi <- priorCheckPheno$photo
  priorCheckPheno$photoi[which(priorCheckPheno$photo == "HP" )] <- "12"
  priorCheckPheno$photoi[which(priorCheckPheno$photo == "LP")] <- "8"
  priorCheckPheno$photoi <- as.numeric(priorCheckPheno$photoi)
  
  chillTrt <- c("LC", "HC")
  chilling <- c( 55.09418, 75.32653, 54.95115, 74.67285, 56.62430, 82.06390,44.63226, 79.17552)
  priorCheckPheno$chill <- rep(chillTrt, each =  nRep, times = n_spec * nPhoto *  nForce  *nRepPrior)
  priorCheckPheno$chilli <- as.numeric(rep(chilling, each =  n_spec * nRep*nChill * nForce * nPhoto, times = nRepPrior))
  
  priorCheckPheno$force.z <- (priorCheckPheno$forcei-mean(priorCheckPheno$forcei,na.rm=TRUE))/(sd(priorCheckPheno$forcei,na.rm=TRUE))
  priorCheckPheno$photo.z <- (priorCheckPheno$photoi-mean(priorCheckPheno$photoi,na.rm=TRUE))/(sd(priorCheckPheno$photoi,na.rm=TRUE))
  priorCheckPheno$chill.z <- (priorCheckPheno$chilli-mean(priorCheckPheno$chilli,na.rm=TRUE))/(sd(priorCheckPheno$chilli,na.rm=TRUE))
  
  
  for (ir in 1:nRepPrior){
    # Parameter Values
    # ir <- 1
    
    sigmaPhenoSp <- rtruncnorm(1, a = 0, mean = lma.data $prior_sigmaPhenoSp_mu, sd = lma.data $prior_sigmaPhenoSp_sigma)
    muPhenoSp <- rnorm(1, lma.data $prior_muPhenoSp_mu, lma.data$prior_muPhenoSp_sigma)
    alphaPhenoSp <- rnorm(n_spec, muPhenoSp, sigmaPhenoSp)
    priorCheckPheno$alphaPhenoSp[priorCheckPheno$sim == ir] <- rep(alphaPhenoSp, each = nRep)
    
    #Cue effects
    priorCheckPheno$betaTraitxForce[priorCheckPheno$simRep == ir] <- rnorm(1,lma.data$prior_betaTraitxForce_mu,lma.data$prior_betaTraitxForce_sigma)
    priorCheckPheno$betaTraitxPhoto[priorCheckPheno$simRep == ir] <- rnorm(1,lma.data$prior_betaTraitxPhoto_mu,lma.data$prior_betaTraitxPhoto_sigma)
    priorCheckPheno$betaTraitxChill[priorCheckPheno$simRep == ir] <- rnorm(1,lma.data$prior_betaTraitxChill_mu,lma.data$prior_betaTraitxChill_sigma)
    
    #Species level slopes sans trait data
    muForceSp <- rnorm(1,lma.data$prior_muForceSp_mu,  lma.data$prior_muForceSp_sigma)
    sigmaForceSp <- rtruncnorm(1, a = 0,mean = lma.data$prior_sigmaForceSp_mu,sd = lma.data$prior_sigmaForceSp_sigma)
    alphaForceSp <- rnorm(n_spec, muForceSp, sigmaForceSp)
    priorCheckPheno$alphaForceSp[priorCheckPheno$simRep == ir] <- rep(alphaForceSp, each = nRep)
    
    muPhotoSp <- rnorm(1, lma.data$prior_muPhotoSp_mu, lma.data$prior_muPhotoSp_sigma)
    sigmaPhotoSp <- rtruncnorm(1, a = 0,mean = lma.data$prior_sigmaPhotoSp_mu, sd = lma.data$prior_sigmaPhotoSp_sigma )
    alphaPhotoSp <- rnorm(n_spec, muPhotoSp, sigmaPhotoSp)
    priorCheckPheno$alphaPhotoSp[priorCheckPheno$simRep == ir] <- rep(alphaPhotoSp, each = nRep)
    
    muChillSp <-  rnorm(1,lma.data$prior_muChillSp_mu,lma.data$prior_sigmaChillSp_sigma)
    sigmaChillSp <- rtruncnorm(1, a = 0,mean = lma.data$prior_sigmaChillSp_mu,sd = lma.data$prior_sigmaChillSp_sigma)
    alphaChillSp <- rnorm(n_spec, muChillSp, sigmaChillSp)
    priorCheckPheno$alphaChillSp[priorCheckPheno$simRep == ir] <- rep(alphaChillSp, each = nRep)
    
    
    #general varience
    priorCheckPheno$sigmapheno_y[priorCheckPheno$simRep == ir] <- rtruncnorm(lma.data$prior_sigmaphenoy_mu,  a = 0, lma.data$prior_sigmaphenoy_sigma)
    priorCheckPheno$e[priorCheckPheno$simRep == ir] <- rnorm(Nph, 0, sigmapheno_y)
    
  }# end simulating new priors, from here vectorize code
  #slopes for each cue, combining trait and non-trait aspect of the slope.
  
  
  priorCheckPheno$betaForceSp <- priorCheckPheno$alphaForceSp + priorCheckPheno$betaTraitxForce *  priorCheckPheno$alphaTraitSp
  
  priorCheckPheno$betaPhotoSp <-  priorCheckPheno$alphaPhotoSp + priorCheckPheno$betaTraitxPhoto * priorCheckPheno$alphaTraitSp
  
  priorCheckPheno$betaChillSp <-  priorCheckPheno$alphaChillSp + priorCheckPheno$betaTraitxChill * priorCheckPheno$alphaTraitSp
  
  #Run full model to get mean simulated y values
  priorCheckPheno$yMu <-  priorCheckPheno$alphaPhenoSp +  priorCheckPheno$betaForceSp * priorCheckPheno$force.z +  priorCheckPheno$betaPhotoSp* priorCheckPheno$photo.z + priorCheckPheno$betaChillSp * priorCheckPheno$chill.z
  
  #Final values
  priorCheckPheno$yPhenoi <- priorCheckPheno$yMu + priorCheckPheno$e
  
  plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$photoi, xlab = "Photoperiod", ylab = "Phenological Date")
  
  #png("figures/pheno_LMA_betaForceSp.png")
  plot(priorCheckPheno$betaForceSp ~ priorCheckPheno$alphaTraitSp )
  #dev.off()
  
  priorCheckPheno_posF <- priorCheckPheno[priorCheckPheno$betaForceSp > 0,]
  plot(priorCheckPheno_posF$betaForceSp ~ priorCheckPheno_posF$alphaTraitSp )
  
  #png("figures/densityYPrior_joint_LMA.png")
  plot(density(priorCheckPheno$yPhenoi))
