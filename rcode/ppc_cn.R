#Started May 6, 2022 by Deirdre

# Aim of this code is to perform several ppc for the trait pheno model for NAm tree traits
# This model includes transect as a dummy variable but otherwise is the same as the traitors model

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

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  else{
  setwd("/home/deirdre/pheno_bc") # for midge
}

# Start by getting the pheno data
dl <- read.csv("input/dl_allbb.csv")

temp <- str_split_fixed(dl$trt, "_", 3); head(temp)
dl$chill<- temp[,1]
dl$photo <- temp[,2]
dl$force <- temp[,3]

dl.chill <- read.csv("input/chilling_values_Hope_Smithers.csv")

dl.wchill <- merge(dl, dl.chill, by = c("population","chill"))
dl.wchill$lab3 <- dl.wchill$lab2
dl.wchill$lab2 <- paste(dl.wchill$species, dl.wchill$population, dl.wchill$rep, sep = "_")

df <- read.csv("input/df_dxb_prepped_data.csv")
df.chill <- read.csv("input/chilling_values_eastern.csv")
df.wchill <- merge(df, df.chill, by =c("population","chill"))
df.wchill <- df.wchill[, c("population", "chill","force","photo","lab2", "bb","species", "treatment","Chill_portions","Utah_Model")]

# mergeing the my data with DF
pheno <- rbind.fill(dl.wchill, df.wchill)
pheno$force.n <- pheno$force
pheno$force.n[pheno$force.n == "HF"] <- "1"
pheno$force.n[pheno$force.n == "LF"] <- "0"
pheno$force.n <- as.numeric(pheno$force.n)

pheno$photo.n <- pheno$photo
pheno$photo.n[pheno$photo.n == "HP"] <- "1"
pheno$photo.n[pheno$photo.n == "LP"] <- "0"
pheno$photo.n <- as.numeric(pheno$photo.n)

# pheno$site.n <- pheno$population
# pheno$site.n[pheno$site.n == "sm"] <- "1"
# pheno$site.n[pheno$site.n == "mp"] <- "2"
# pheno$site.n[pheno$site.n == "HF"] <- "3"
# pheno$site.n[pheno$site.n == "SH"] <- "4"
# pheno$site.n <- as.numeric(pheno$site.n)

pheno$transect <- pheno$population
pheno$transect[pheno$transect == "sm"] <- "0"
pheno$transect[pheno$transect == "mp"] <- "0"
pheno$transect[pheno$transect == "kl"] <- "0"
pheno$transect[pheno$transect == "af"] <- "0"
pheno$transect[pheno$transect == "GR"] <- "1"
pheno$transect[pheno$transect == "HF"] <- "1"
pheno$transect[pheno$transect == "SH"] <- "1"
pheno$transect[pheno$transect == "WM"] <- "1"

# head(pheno)
# #add dummy/ site level effects:
# pheno <- pheno %>%
#   mutate ( site2 = if_else(site.n == 2, 1, 0),
#            site3 = if_else(site.n == 3, 1, 0),
#            site4 = if_else(site.n == 4, 1, 0))

# standardize the 0/1 and standardize sites? 
pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)

# pheno$site2.z2 <- (pheno$site2-mean(pheno$site2,na.rm=TRUE))/(sd(pheno$site2,na.rm=TRUE)*2)
# pheno$site3.z2 <- (pheno$site3-mean(pheno$site3,na.rm=TRUE))/(sd(pheno$site3,na.rm=TRUE)*2)
# pheno$site4.z2 <- (pheno$site4-mean(pheno$site4,na.rm=TRUE))/(sd(pheno$site4,na.rm=TRUE)*2)

#going to split it into analysis of terminal bb and lateral bb
# Starting with the terminal buds:
#pheno.term <- pheno[,c("tbb", "chill.n", "force.n", "photo.n", "site.n", "species", "lab2")]
pheno.term <- pheno[,c("bb", "force.z2", "photo.z2", "population", "species", "lab2","Utah_Model","Chill_portions","chillport.z2","transect")] #"site2.z2", "site3.z2","site4.z2")]
pheno.t <- pheno.term[complete.cases(pheno.term$bb), ] # none,great!
length(unique(pheno.t$species))
#pheno.t <- pheno.term[complete.cases(pheno.term$bb), ] # 1780 rows data 
pheno.t$species <- tolower(pheno.t$species)
pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 47 bc two species occur in both transects

# Now get the trait data and subset to only include spp we have pheno data for:
setwd("..//Treetraits")
trtData <- read.csv("data/allTrt.csv", stringsAsFactors = FALSE)
head(trtData)

phenoSp <- sort(unique(pheno.t$species))

trtPheno <- trtData[trtData$species %in% phenoSp, ]
length(unique(trtPheno$species))

trtPheno$transect <- trtPheno$site
trtPheno$transect[trtPheno$transect == "sm"] <- "0"
trtPheno$transect[trtPheno$transect == "mp"] <- "0"
trtPheno$transect[trtPheno$transect == "kl"] <- "0"
trtPheno$transect[trtPheno$transect == "af"] <- "0"
trtPheno$transect[trtPheno$transect == "GR"] <- "1"
trtPheno$transect[trtPheno$transect == "HF"] <- "1"
trtPheno$transect[trtPheno$transect == "SH"] <- "1"
trtPheno$transect[trtPheno$transect == "WM"] <- "1"

#########################################################

fit <- readRDS("output/cn_stanfit_np.RDS")


post<- rstan::extract(fit)

postCN <- data.frame(post)

cueEffects <- postCN[, colnames(postCN) %in% c("muPhenoSp", "muForceSp", "muChillSp", "muPhotoSp", "sigmapheno_y")]

mcmc_intervals(cueEffects) + 
  theme_classic() + 
  labs(title = "main intercept, cue slopes and general error")

###############
postCN_alpaForceSp <- postCN[,colnames(postCN) %in% grep( "alphaForceSp", colnames(postCN), value = TRUE)]
colnames(postCN_alpaForceSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postCN_alpaForceSp) + 
  geom_vline(xintercept = mean(postCN$muForceSp), linetype="dotted", color = "grey")  +
  theme_classic() + 
  labs(subtitle = paste0("Mean muForceSp was ", round(mean(postCN$muForceSp),3)),
       title = "muForceSp - species forcing slopes no trait")

postCN_betaForceSp <- postCN[,colnames(postCN) %in% grep( "betaForceSp", colnames(postCN), value = TRUE)]
colnames(postCN_betaForceSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postCN_betaForceSp) + 
  theme_classic() + 
  labs(title = "betaForceSp - Species forcing slopes with trait value")

#Different species slopes for chilling, without the effect of trait
postCN_alphaChillSp <- postCN[,colnames(postCN) %in% grep( "alphaChillSp", colnames(postCN), value = TRUE)]
colnames(postCN_alphaChillSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postCN_alphaChillSp) + 
  geom_vline(xintercept = mean(postCN$muChillSp), linetype="dotted", color = "grey")  +
  theme_classic() + 
  labs(subtitle = paste0("Mean muChillSp was ", round(mean(postCN$muChillSp),3)),
       title = "alphaChillSp - Species chill slopes no trait")

#Different species slopes for forcing, with the effect of trait
postCN_betaChillSp <- postCN[,colnames(postCN) %in% grep( "betaChillSp", colnames(postCN), value = TRUE)]
colnames(postCN_betaChillSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postCN_betaChillSp) + 
  theme_classic() + 
  labs(title = "betaChillSp - Species chilling slopes with trait value")

#Different species slopes for photoperiod, without the effect of trait
postCN_alphaPhotoSp <- postCN[,colnames(postCN) %in% grep( "alphaPhotoSp", colnames(postCN), value = TRUE)]
colnames(postCN_alphaPhotoSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postCN_alphaPhotoSp) + 
  geom_vline(xintercept = mean(postCN$muPhotoSp), linetype="dotted", color = "grey")  +
  theme_classic() + 
  labs(subtitle = paste0("Mean muPhotoSp was ", round(mean(postCN$muPhotoSp),3)),
       title = "muPhotoSp - Species photo period slopes no trait")

#Different species slopes for forcing, with the effect of trait
postCN_betaPhotoSp <- postCN[,colnames(postCN) %in% grep( "betaPhotoSp", colnames(postCN), value = TRUE)]
colnames(postCN_betaPhotoSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postCN_betaPhotoSp) + 
  theme_classic() + 
  labs(title = "betaPhotoSp - Species photoperiod slopes with trait value")


#Different species slopes for forcing only the effect of trait
postCN_betaTraitx <- postCN[,colnames(postCN) %in% grep( "betaTraitx", colnames(postCN), value = TRUE)]

mcmc_intervals(postCN_betaTraitx) + 
  theme_classic() + 
  labs(title = "effect's of traits on cue slopes")

# require(bayesplot)
# y <- pheno$bb 
# yrep <-  postCN[,colnames(postCN) %in% grep( "y_hat", colnames(postCN), value = TRUE)]
# yrepM <- colMeans(yrep)
# 
# ppc_dens_overlay(y, yrepM[1:50,])

#################################################################################
Nrep <- 15 # rep per trait
Ntransect <- 2
Nspp <- 40 # number of species 

# First making a data frame for the test trait data
Ntrt <- Nspp * Ntransect * Nrep # total number of traits observations
Ntrt

mu.grand <- 30 # the grand mean of the height model
sigma.species <- 10 # we want to keep the variaiton across spp. high
sigma.transect <- 4
sigmaTrait_y <- 4

#make a dataframe for height
trt.dat <- data.frame(matrix(NA, Ntrt, 1))
names(trt.dat) <- c("rep")
trt.dat$rep <- c(1:Nrep)
trt.dat$transect <- rep(c(1:Ntransect), each = Nspp)
trt.dat$species <- rep(1:Nspp, Ntransect)

# now generating the species trait data, here it is for height
#the alphaTraitSp in Faiths original code:
alphaTraitSp <- rnorm(Nspp, 0, sigma.species)
trt.dat$alphaTraitSp <- rep(alphaTraitSp, Ntransect) #adding ht data for ea. sp

#now generating the effects of transect
mutransect <- rnorm(Ntransect, 0, sigma.transect) #intercept for each transect
trt.dat$mutransect <- rep(mutransect, each = Nspp) # generate data for ea transect

# general variance
trt.var <- 5 #sigmaTrait_y in the stan code
trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)

# generate yhat - heights -  for this first trt model
#trt.dat$yTraiti <- mu.grand + trt.dat$alphaTraitSp + trt.dat$mutransect + trt.dat$trt.er

for (i in 1:Ntrt){
  trt.dat$mu_grand_sp[i] <-  trt.dat$alphaTraitSp[i] +  mu.grand
}

for (i in 1:Ntrt){
  trt.dat$yTraiti[i] <-  trt.dat$alphaTraitSp[i] + trt.dat$mutransect[i] +  mu.grand
}

#Add grand mean to trait values 
alphaTraitspFull <-  alphaTraitSp + mu.grand
trt.dat$alphaTraitspFull <-trt.dat$alphaTraitSp + mu.grand

hist( trt.dat$yTraiti)
hist(alphaTraitspFull)


#### Pheno data generation ##############################
#All Negative betatraitX Values 
#--------------------------------------
n_spec <-Nspp # same number of species as teh traits section fo the model 
#Number of repeat observations per species
nRep <- 15
#Overall number of pbservations (rows)
Nph <- n_spec * nRep # for phenology model

#Make a data frame for input phenology simulation data
pheno.dat <- data.frame(matrix(NA, Nph, 2))
names(pheno.dat) <- c("rep","species")
pheno.dat$rep <- c(1:Nph)
pheno.dat$species <- rep(c(1:n_spec), each = nRep)

#Simulate mean SLA offset data per species (not mean value)
#meanSLA <- rnorm(n_spec, 0, 5)
#Make this the name of the full vector of sla per species values - alphaTraitSp 
pheno.dat$alphaTraitSp <- rep(alphaTraitSp + mu.grand, each = nRep)

#Simulate cues (z scored)
pheno.dat$forcei <- rnorm(Nph, 1, 1)

pheno.dat$photoi <- rnorm(Nph, 1, 0.5) # less photoperiod 
pheno.dat$chilli <- rnorm(Nph, 1, 1) #more chilling

# Parameter Values
#Species means
sigmaPhenoSp <- 10
muPhenoSp <- 40
alphaPhenoSp <- rnorm(n_spec, muPhenoSp, sigmaPhenoSp)
pheno.dat$alphaPhenoSp <- rep(alphaPhenoSp, each = nRep)

#Cue effects
betaTraitxForce <- -0.5 
betaTraitxPhoto <- -0.2
betaTraitxChill <- -0.4

#Species level slopes sans trait data
muForceSp <- -0.4
sigmaForceSp <- 4
alphaForceSp <- rnorm(n_spec, muForceSp, sigmaForceSp)
pheno.dat$alphaForceSp <- rep(alphaForceSp, each = nRep)

muPhotoSp <- -0.05
sigmaPhotoSp <- 4
alphaPhotoSp <- rnorm(n_spec, muPhotoSp, sigmaPhotoSp)
pheno.dat$alphaPhotoSp <- rep(alphaPhotoSp, each = nRep)

muChillSp <- -0.6
sigmaChillSp <- 4
alphaChillSp <- rnorm(n_spec, muChillSp, sigmaChillSp)
pheno.dat$alphaChillSp <- rep(alphaChillSp, each = nRep)

#general varience
sigmapheno_y <- 5
pheno.dat$e <- rnorm(Nph, 0, sigmapheno_y)

#slopes for each cue, combining trait and non-trait aspect of the slope.
#Make columns to put data 
pheno.dat$betaForceSp <- NA
pheno.dat$betaPhotoSp <- NA
pheno.dat$betaChillSp <- NA

for (iSp in 1:n_spec){
  
  #iSp <- 1
  #Select species data of interest 
  
  #Forcing
  betaForceSp <- alphaForceSp[iSp] + betaTraitxForce * alphaTraitspFull[iSp]
  pheno.dat$betaForceSp[pheno.dat$species == iSp] <- betaForceSp
  
  #chilling
  betaChillSp <- alphaChillSp[iSp] + betaTraitxChill * alphaTraitspFull[iSp]
  pheno.dat$betaChillSp[pheno.dat$species == iSp] <- betaChillSp
  
  #photoperiod
  betaPhotoSp <- alphaPhotoSp[iSp] + betaTraitxPhoto * alphaTraitspFull[iSp]
  pheno.dat$betaPhotoSp[pheno.dat$species == iSp] <- betaPhotoSp
  
  
}

#Run full model to get mean simulated y values
for (i in 1:Nph){
  pheno.dat$yMu[i] <-  pheno.dat$alphaPhenoSp[i] +  pheno.dat$betaForceSp[i] * pheno.dat$forcei[i] +  pheno.dat$betaPhotoSp[i] * pheno.dat$photoi[i] + pheno.dat$betaChillSp[i] * pheno.dat$chilli[i]
}

#Final values
pheno.dat$yPhenoi <- pheno.dat$yMu + pheno.dat$e

#What does the data look like?
plot(density(pheno.dat$yPhenoi))

#Plot trait values againt forcing slopes

head(pheno.dat)

plot(pheno.dat$betaForceSp ~ pheno.dat$alphaTraitSp)
plot(pheno.dat$alphaForceSp ~ pheno.dat$alphaTraitSp)

#######################################################

# png("figures/simPosteriorHist.png")
# par(mfrow=c(3,4))
#Compare results to simulated values
# hist(postCN$muPhenoSp, main = paste("muPhenoSp is " , signif(muPhenoSp,3), sep = ""), xlim = c(0,100))
# abline(v = muPhenoSp, col="red", lwd=3, lty=2)
# 
# hist(postCN$muForceSp, main = paste("muForceSp is " , signif(muForceSp,3), sep = ""))
# abline(v = muForceSp, col="red", lwd=3, lty=2)
# 
# hist(postCN$muChillSp, main = paste("muChillSp is " , signif(muChillSp,3), sep = ""))
# abline(v = muChillSp, col="red", lwd=3, lty=2)
# 
# hist(postCN$muPhotoSp, main = paste("muPhotoSp is " , signif(muPhotoSp,3), sep = ""))
# abline(v = muPhotoSp, col="red", lwd=3, lty=2)
# 
# hist(postCN$sigmapheno_y, main = paste("sigmapheno_y is " , signif(sigmapheno_y,3), sep = ""))
# abline(v = sigmapheno_y, col="red", lwd=3, lty=2)
# 
# plot(density(postCN$betaTraitxForce), main = paste("betaTraitxForce is " , signif(betaTraitxForcePos,3), sep = ""))
# abline(v = betaTraitxForcePos, col="red", lwd=3, lty=2)
# # 
# hist(postCN$betaTraitxChill, main = paste("betaTraitxChill is " , signif(betaTraitxChill,3), sep = ""))
# abline(v = betaTraitxChill, col="red", lwd=3, lty=2)
# # 
# hist(postCN$betaTraitxPhoto, main = paste("betaTraitxPhoto is " , signif(betaTraitxPhoto,3), sep = ""))
# abline(v = betaTraitxPhoto, col="red", lwd=3, lty=2)
# 
# hist(postCN$sigmaChillSp, main = paste("sigmaChillSp is " , signif(sigmaChillSp,3), sep = ""))
# abline(v = sigmaChillSp, col="red", lwd=3, lty=2)
# 
# hist(postCN$sigmaForceSp, main = paste("sigmaForceSp is " , signif(sigmaForceSp,3), sep = ""))
# abline(v = sigmaForceSp, col="red", lwd=3, lty=2)
# 
# hist(postCN$sigmaPhotoSp, main = paste("sigmaPhotoSp is " , signif(sigmaPhotoSp,3), sep = ""))
# abline(v = sigmaPhotoSp, col="red", lwd=3, lty=2)

# png("figures/simulatedPairs.png")
pairs(mdl.cn, pars = c("muForceSp", "muChillSp", "muPhotoSp", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__")) 
# dev.off()

  #Prior Predictive Check (Run 1000 times and plot results)
  #------------------------------------------------------------
  
  #Number fo prior check itterations 
  nRepPrior <- 300
  
  # ppc for traits portion
  priorCheckTrait <- data.frame(matrix(NA, Ntrt*nRepPrior, 3))
  names(priorCheckTrait) <- c("simRep", "rep", "species")
  priorCheckTrait$simRep <- rep(1:nRepPrior, each = Ntrt)
  priorCheckTrait$rep <- rep(1:Ntrt, times = nRepPrior)
  priorCheckTrait$species <- rep(1:Nspp, each = nRep)
  priorCheckTrait$transect <- rep(1:Ntransect, each = nRep)
  
  #traitSLA <- rnorm(Ntrt, 20, 5)
  
  #Make this the name of the full vector of sla per species values - alphaTraitSp 
  #priorCheckTrait$alphaTraitSp <-  rep(rep(trt.dat$mu_grand_sp, times = nRepPrior))
  carbNit <- trtPheno[complete.cases(trtPheno$C.N),]
  
  specieslist <- sort(unique(trtPheno$species))
  sitelist <- sort(unique(trtPheno$transect))
  leafMass <- trtPheno[complete.cases(trtPheno$lma),]
  
  cn.data <- list(yTraiti = carbNit$C.N,
                  N = nrow(carbNit),
                  n_spec = length(specieslist),
                  trait_species = as.numeric(as.factor(carbNit$species)),
                  n_site = length(sitelist),
                  site = as.numeric(as.factor(carbNit$transect)),
                  prior_mu_grand_mu = 20,
                  prior_mu_grand_sigma = 5, #widened
                  prior_sigma_sp_mu = 10,
                  prior_sigma_sp_sigma = 5,
                  prior_sigma_site_mu = 5,
                  prior_sigma_site_sigma = 2,
                  prior_sigma_traity_mu = 5,
                  prior_sigma_traity_sigma = 2,
                  ## Phenology
                  Nph = nrow(pheno.t),
                  phenology_species = as.numeric(as.factor(pheno.t$species)),
                  yPhenoi = pheno.t$bb,
                  forcei = pheno.t$force.z2,
                  chilli = pheno.t$chillport.z2,
                  photoi = pheno.t$photo.z2,
                  prior_muForceSp_mu = -15,
                  prior_muForceSp_sigma = 15, #10 #wider
                  prior_muChillSp_mu = -15,
                  prior_muChillSp_sigma = 15, #10 #wider
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
    
    muGrand <- rnorm(1,  mean = cn.data$prior_mu_grand_mu, sd = cn.data$prior_mu_grand_sigma)
    sigmaSp <- rnorm(1,  mean = cn.data$prior_sigma_sp_mu, sd = cn.data$prior_sigma_sp_sigma)
    sigmatransect <- rnorm(1, mean = cn.data$prior_sigma_site_mu, sd = cn.data$prior_sigma_site_sigma)
    
    alphaTraitSp <- rnorm(Nspp, 0, sigma.species)
    priorCheckTrait$alphaTraitSp[priorCheckTrait$simRep == ir] <- rep(alphaTraitSp, each = nRep)
    
    muSp <- rnorm(Nspp, 0, sigma.species)
    priorCheckTrait$muSp[priorCheckTrait$simRep == ir] <- rep(muSp, each = nRep)
    
    mutransect <- rnorm(Ntransect, 0, sigmatransect)
    priorCheckTrait$mutransect[priorCheckTrait$simRep == ir] <- rep(mutransect, each = nRep)
    
    #general varience
    priorCheckTrait$sigmaTrait_y[priorCheckTrait$simRep == ir] <- rnorm(cn.data$prior_sigma_traity_mu, cn.data$prior_sigma_traity_sigma)
    priorCheckTrait$e[priorCheckTrait$simRep == ir] <- rnorm(Ntrt, 0, sigmaTrait_y)
    
    priorCheckTrait$yTraiti <- muGrand + priorCheckTrait$muSp + priorCheckTrait$mutransect + priorCheckTrait$e
  }# end simulating new priors, from here vectorize code
  
  #Final values
  priorCheckTrait$muGrandSp <- muGrand + priorCheckTrait$muSp
  
  
  png("figures/density_Trait_Prior_joint_cn.png")
  plot(density(priorCheckTrait$yTraiti))
  dev.off()
  
  png("figures/GrandSp_PlotPrior_joint_cn.png")
  plot(priorCheckTrait$yTraiti ~ priorCheckTrait$muGrandSp, xlab = "muGrandSp", ylab = "Trait")
  dev.off()
  
  png("figures/MuSp_PlotPrior_joint_cn.png")
  plot(priorCheckTrait$yTraiti ~ priorCheckTrait$muSp, xlab = "MuSp", ylab = "Trait")
  dev.off()
  
  png("figures/Mutransect_PlotPrior_joint_cn.png")
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
    #ir <- 1
    
    #Species means
    sigmaPhenoSp <- rnorm(1, a = 0, mean = cn.data$prior_sigmaPhenoSp_mu, sd = cn.data$prior_sigmaPhenoSp_sigma)
    muPhenoSp <- rnorm(1, cn.data$prior_muPhenoSp_mu, cn.data$prior_muPhenoSp_sigma)
    alphaPhenoSp <- rnorm(n_spec, muPhenoSp, sigmaPhenoSp)
    priorCheckPheno$alphaPhenoSp[priorCheckPheno$simRep == ir] <- rep(alphaPhenoSp, each = nRep)
    
    #Cue effects
    priorCheckPheno$betaTraitxForce[priorCheckPheno$simRep == ir] <- rnorm(1,cn.data$prior_betaTraitxForce_mu,cn.data$prior_betaTraitxForce_sigma)
    priorCheckPheno$betaTraitxPhoto[priorCheckPheno$simRep == ir] <- rnorm(1,cn.data$prior_betaTraitxPhoto_mu,cn.data$prior_betaTraitxPhoto_sigma)
    priorCheckPheno$betaTraitxChill[priorCheckPheno$simRep == ir] <- rnorm(1,cn.data$prior_betaTraitxChill_mu,cn.data$prior_betaTraitxChill_sigma)
    
    #Species level slopes sans trait data
    muForceSp <- rnorm(1,cn.data$prior_muForceSp_mu,  cn.data$prior_muForceSp_sigma)
    sigmaForceSp <- rnorm(1, a = 0,mean = cn.data$prior_sigmaForceSp_mu,sd = cn.data$prior_sigmaForceSp_sigma)
    alphaForceSp <- rnorm(n_spec, muForceSp, sigmaForceSp)
    priorCheckPheno$alphaForceSp[priorCheckPheno$simRep == ir] <- rep(alphaForceSp, each = nRep)
    
    muPhotoSp <- rnorm(1, cn.data$prior_muPhotoSp_mu, cn.data$prior_muPhotoSp_sigma)
    sigmaPhotoSp <- rnorm(1, a = 0,mean = cn.data$prior_sigmaPhotoSp_mu, sd = cn.data$prior_sigmaPhotoSp_sigma )
    alphaPhotoSp <- rnorm(n_spec, muPhotoSp, sigmaPhotoSp)
    priorCheckPheno$alphaPhotoSp[priorCheckPheno$simRep == ir] <- rep(alphaPhotoSp, each = nRep)
    
    muChillSp <-  rnorm(1,cn.data$prior_sigmaChillSp_mu,cn.data$prior_sigmaChillSp_sigma)
    sigmaChillSp <- rnorm(1, a = 0,mean = cn.data$prior_sigmaChillSp_mu,sd = cn.data$prior_sigmaChillSp_sigma)
    alphaChillSp <- rnorm(n_spec, muChillSp, sigmaChillSp)
    priorCheckPheno$alphaChillSp[priorCheckPheno$simRep == ir] <- rep(alphaChillSp, each = nRep)
    
    
    #general varience
    priorCheckPheno$sigmapheno_y[priorCheckPheno$simRep == ir] <- rnorm(cn.data$prior_sigmaphenoy_mu,  a = 0, cn.data$prior_sigmaphenoy_sigma)
    priorCheckPheno$e[priorCheckPheno$simRep == ir] <- rnorm(Nph, 0, sigmapheno_y)
    
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
  
  png("figures/densityYPrior_joint.png")
  plot(density(priorCheckPheno$yPhenoi))
  dev.off()
  
  png("figures/photoPlotPrior_joint.png")
  plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$photoi, xlab = "Photoperiod", ylab = "Phenological Date")
  dev.off()
  
  png("figures/forcingPlotPrior_joint.png")
  plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$forcei, xlab = "Forcing", ylab = "Phenological Date")
  dev.off()
  
  png("figures/chillingPlotPrior_joint.png")
  plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$chilli, xlab = "Chillina", ylab = "Phenological Date")
  dev.off()
  
}

if(BayesSweave == TRUE){
  #For the BayesClass sweave documents 
  setwd("/home/faith/Documents/github/bayes2020/Projects/Faith/traitorsModel")
}
