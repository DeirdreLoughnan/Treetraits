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
library(truncnorm)

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

pheno$transect <- pheno$population
pheno$transect[pheno$transect == "sm"] <- "0"
pheno$transect[pheno$transect == "mp"] <- "0"
pheno$transect[pheno$transect == "kl"] <- "0"
pheno$transect[pheno$transect == "af"] <- "0"
pheno$transect[pheno$transect == "GR"] <- "1"
pheno$transect[pheno$transect == "HF"] <- "1"
pheno$transect[pheno$transect == "SH"] <- "1"
pheno$transect[pheno$transect == "WM"] <- "1"

# standardize the 0/1 and standardize sites? 
pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)

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

trtPheno$latitude <- trtPheno$site
trtPheno$latitude[trtPheno$latitude == "sm"] <- 54.7824
trtPheno$latitude[trtPheno$latitude == "mp"] <- 49.0646
trtPheno$latitude[trtPheno$latitude == "kl"] <- 52.1417
trtPheno$latitude[trtPheno$latitude == "af"] <- 50.8837
trtPheno$latitude[trtPheno$latitude == "GR"] <- 43.99837498
trtPheno$latitude[trtPheno$latitude == "HF"] <- 42.5315
trtPheno$latitude[trtPheno$latitude == "SH"] <- 45.9310
trtPheno$latitude[trtPheno$latitude == "WM"] <- 44.92466697
#########################################################

#fit <- readRDS("output/lma_stanfit.RDS")
load("output/lmaDummyInt.Rdata")
post<- rstan::extract(mdl)

postLMA <- data.frame(post)

cueEffects <- postLMA[, colnames(postLMA) %in% c("muPhenoSp", "muForceSp", "muChillSp", "muPhotoSp", "sigmapheno_y")]

mcmc_intervals(cueEffects) + 
  theme_classic() + 
  labs(title = "main intercept, cue slopes and general error")

###############
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

#Different species slopes for chilling, without the effect of trait
postLMA_alphaChillSp <- postLMA[,colnames(postLMA) %in% grep( "alphaChillSp", colnames(postLMA), value = TRUE)]
colnames(postLMA_alphaChillSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postLMA_alphaChillSp) + 
  geom_vline(xintercept = mean(postLMA$muChillSp), linetype="dotted", color = "grey")  +
  theme_classic() + 
  labs(subtitle = paste0("Mean muChillSp was ", round(mean(postLMA$muChillSp),3)),
       title = "alphaChillSp - Species chill slopes no trait")

#Different species slopes for forcing, with the effect of trait
postLMA_betaChillSp <- postLMA[,colnames(postLMA) %in% grep( "betaChillSp", colnames(postLMA), value = TRUE)]
colnames(postLMA_betaChillSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postLMA_betaChillSp) + 
  theme_classic() + 
  labs(title = "betaChillSp - Species chilling slopes with trait value")

#Different species slopes for photoperiod, without the effect of trait
postLMA_alphaPhotoSp <- postLMA[,colnames(postLMA) %in% grep( "alphaPhotoSp", colnames(postLMA), value = TRUE)]
colnames(postLMA_alphaPhotoSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postLMA_alphaPhotoSp) + 
  geom_vline(xintercept = mean(postLMA$muPhotoSp), linetype="dotted", color = "grey")  +
  theme_classic() + 
  labs(subtitle = paste0("Mean muPhotoSp was ", round(mean(postLMA$muPhotoSp),3)),
       title = "muPhotoSp - Species photo period slopes no trait")

#Different species slopes for forcing, with the effect of trait
postLMA_betaPhotoSp <- postLMA[,colnames(postLMA) %in% grep( "betaPhotoSp", colnames(postLMA), value = TRUE)]
colnames(postLMA_betaPhotoSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postLMA_betaPhotoSp) + 
  theme_classic() + 
  labs(title = "betaPhotoSp - Species photoperiod slopes with trait value")

#Different species slopes for forcing only the effect of trait
postLMA_betaTraitx <- postLMA[,colnames(postLMA) %in% grep( "betaTraitx", colnames(postLMA), value = TRUE)]

mcmc_intervals(postLMA_betaTraitx) + 
  theme_classic() + 
  labs(title = "effect's of traits on cue slopes")


#################################################################################
Nrep <- 10# rep per trait
Npop <- 8
Ntran <- 2
Nspp <- 10# number of species with traits (making this 20 just for speed for now)

# First making a data frame for the test trait data
Ntrt <- Nspp * Npop * Nrep# total number of traits observations
Ntrt

#make a dataframe for height
trt.dat <- data.frame(matrix(NA, Ntrt, 1))
names(trt.dat) <- c("rep")
trt.dat$rep <- c(1:Nrep)
trt.dat$species <- rep(1:Nspp, each = Nrep)
trt.dat$pop <- rep(1:Npop, each = Nspp*Nrep)
trt.dat$tran <- rep(1:Ntran, each = 4*Nrep*Nspp)


lati <- rnorm(8, 50, 5)
trt.dat$lat <- rep(lati, each = Nrep*Nspp)
trt.dat$lat <- as.numeric(trt.dat$lat)

mu.tranE <- 4
trt.dat$dumE <- as.numeric(ifelse(trt.dat$tran == "1","0","1"))
trt.dat$mutranE <- mu.tranE*trt.dat$dumE

mu.tranlat  <- 2
# sigma.tranlat = 1
#alpha.tranlat <- rnorm(Npop,mu.tranlat, sigma.tranlat)
#trt.dat$alpha.tranlat <- rep(alpha.tranlat, each = Nrep*Nspp)
trt.dat$alpha.tranlat <- mu.tranlat*(trt.dat$dumE*trt.dat$lat)

# trt.dat$species <- rep(1:Nspp, Nstudy)

# now generating the species trait data, here it is for height
#mu.grand <- 10
sigma.species <- 0.05 # we want to keep the variaiton across spp. high

#the alphaTraitSp in Faiths original code:
mu.trtsp <- rnorm(Nspp, 0.01, sigma.species)
trt.dat$mu.trtsp <- rep(mu.trtsp, each = Nrep) #adding ht data for ea. sp

# general variance
trt.var <- 1 #sigma_traity in the stan code
trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)

# generate yhat - heights -  for this first trt model
for (i in 1:Ntrt){
  trt.dat$yTraiti[i] <- 
    trt.dat$mu.trtsp[i] + trt.dat$trt.er[i] + trt.dat$mutranE[i]+ mu.tranlat * (trt.dat$dumE[i]*trt.dat$lat[i])  
}


# for (i in 1:Ntrt){
#   trt.dat$yTraiti[i] <- 
#     trt.dat$mu.trtsp[i] + trt.dat$trt.er[i] + trt.dat$mutranE[i]+ trt.dat$alpha.tranlat[i] * (trt.dat$dumE[i]*trt.dat$lat[i])  
# }

#####################################################################
##### Phenology test data ###########################

Nchill <- 2 # sm high low, mp high low
Nphoto <- 2 # high low
Nforce <- 2 # high amd low

Nrep <- 10
Nspp <- 10
Npop <- 2
Nph <- Nchill*Nphoto*Nforce*Nrep*Npop*Nspp #*Nforce - bc N force and N chill are completely confounded

pheno.dat <- data.frame(matrix(NA, Nph, 2))
names(pheno.dat) <- c("rep","species")
pheno.dat$rep <- c(1:Nrep)
pheno.dat$species <- rep(c(1:Nspp), each = Nrep)
pheno.dat$pop <- rep(1:Npop, each = Nspp*Nrep)

forceTrt <- c("LF", "HF")
forcei <- rnorm(Nforce, 1, 2)
pheno.dat$force <- rep(forceTrt, each = Nchill*Nphoto*Npop*Nrep*Nspp)
pheno.dat$forcei <- rep(forcei, each = Nchill*Nphoto*Npop*Nrep*Nspp)

chillTrt <- c("LC", "HC")
chilli <- rnorm(Nchill, 1, 5)
pheno.dat$chill <- rep(rep(chillTrt, each = Nphoto*Npop*Nrep*Nspp), times = Nforce)
pheno.dat$chilli <- rep(rep(chilli, each = Nphoto*Npop*Nrep*Nspp), times = Nforce)

photoTrt <- c("LP", "HP")
photoi <- rnorm(Nphoto, 1, 2)
pheno.dat$photo <- rep(rep(photoTrt, each = Npop*Nrep*Nspp), times = Nforce*Nchill)
pheno.dat$photoi <- rep(rep(photoi, each = Npop*Nrep*Nspp), times = Nforce*Nchill)

# sanity check that all treatments are the same
pheno.dat$label <- paste(pheno.dat$chill, pheno.dat$force, pheno.dat$photo, sep ="_")
pheno.dat$label <- paste(pheno.dat$chilli, pheno.dat$forcei, pheno.dat$photoi, sep ="_")
temp <- subset(pheno.dat, label == "HC_HF_HP"); dim(temp)

mu.force = -10 
sigma.force = 1
alpha.force.sp <- rnorm(Nspp, mu.force, sigma.force)
pheno.dat$alphaForceSp <- rep(alpha.force.sp, each = Nrep)

mu.photo = -15
sigma.photo = 1
alpha.photo.sp <- rnorm(Nspp, mu.photo, sigma.photo)
pheno.dat$alphaPhotoSp <- rep(alpha.photo.sp, each = Nrep)

mu.chill = -14
sigma.chill = 1
alpha.chill.sp <- rnorm(Nspp, mu.chill, sigma.chill)
pheno.dat$alphaChillSp <- rep(alpha.chill.sp, each = Nrep)

mu.pheno.sp = 80
sigma.pheno.sp = 30
alphaPhenoSp <- rnorm(Nspp, mu.pheno.sp, sigma.pheno.sp)
pheno.dat$alphaPhenoSp <- rep(alphaPhenoSp, each = Nrep)

sigma.pheno.y = 3
pheno.dat$ePhen <- rnorm(Nph, 0, sigma.pheno.y)

betaTraitxForce <- 0.3 
betaTraitxPhoto <- -0.2
betaTraitxChill <- -0.4

pheno.datTrait <- merge(pheno.dat, unique(trt.dat[,c("species","mu.trtsp")]), by = "species")
#head(pheno.datTrait,50)

for (i in 1:Nph){
  pheno.datTrait$betaForceSp[i] <-  pheno.datTrait$alphaForceSp[i] + (betaTraitxForce *  pheno.datTrait$mu.trtsp[i])
  
  pheno.datTrait$betaPhotoSp[i]<- pheno.datTrait$alphaPhotoSp[i] + (betaTraitxPhoto*  pheno.datTrait$mu.trtsp[i])
  
  pheno.datTrait$betaChillSp[i] <-pheno.datTrait$alphaChillSp[i] + (betaTraitxChill* pheno.datTrait$mu.trtsp[i])
}

for (i in 1:Nph){
  pheno.datTrait$yMu[i] <-  pheno.datTrait$alphaPhenoSp[i] +  pheno.datTrait$betaForceSp[i] * pheno.datTrait$forcei[i] +  pheno.datTrait$betaPhotoSp[i] * pheno.datTrait$photoi[i] + pheno.datTrait$betaChillSp[i] * pheno.datTrait$chilli[i]
}

pheno.datTrait$yPhenoi <- pheno.datTrait$yMu + pheno.datTrait$ePhen
dim(pheno.dat)


hist( trt.dat$yTraiti)
hist(alphaTraitspFull)


#What does the data look like?
plot(density(pheno.datTrait$yPhenoi))

#Plot trait values againt forcing slopes

head(pheno.dat)

plot(pheno.datTrait$betaForceSp ~ pheno.datTrait$mu.trtsp)
plot(pheno.dat$alphaForceSp ~ pheno.dat$mu.trtsp)

plot(pheno.dat$betaChillSp ~ pheno.dat$alphaTraitSp)
plot(pheno.dat$alphaChillSp ~ pheno.dat$alphaTraitSp)

plot(pheno.dat$betaPhotoSp ~ pheno.dat$alphaTraitSp)
plot(pheno.dat$alphaPhotoSp ~ pheno.dat$alphaTraitSp)

#######################################################

#png("figures/simPosteriorHist.png")
par(mfrow=c(1,1))
##Compare results to simulated values
hist(postLMA$muPhenoSp, main = paste("muPhenoSp is " , signif(muPhenoSp,3), sep = ""), xlim = c(0,100))
abline(v = muPhenoSp, col="red", lwd=3, lty=2)

hist(postLMA$muForceSp, main = paste("muForceSp is " , signif(muForceSp,3), sep = ""))
abline(v = muForceSp, col="red", lwd=3, lty=2)

hist(postLMA$muChillSp, main = paste("muChillSp is " , signif(muChillSp,3), sep = ""))
abline(v = muChillSp, col="red", lwd=3, lty=2)

hist(postLMA$muPhotoSp, main = paste("muPhotoSp is " , signif(muPhotoSp,3), sep = ""))
abline(v = muPhotoSp, col="red", lwd=3, lty=2)

hist(postLMA$sigmapheno_y, main = paste("sigmapheno_y is " , signif(sigmapheno_y,3), sep = ""))
abline(v = sigmapheno_y, col="red", lwd=3, lty=2)

hist(postLMA$betaTraitxForce, main = paste("betaTraitxForce is " , signif(betaTraitxForce,3), sep = ""))
abline(v = betaTraitxForce, col="red", lwd=3, lty=2)
#
hist(postLMA$betaTraitxChill, main = paste("betaTraitxChill is " , signif(betaTraitxChill,3), sep = ""))
abline(v = betaTraitxChill, col="red", lwd=3, lty=2)
#
hist(postLMA$betaTraitxPhoto, main = paste("betaTraitxPhoto is " , signif(betaTraitxPhoto,3), sep = ""))
abline(v = betaTraitxPhoto, col="red", lwd=3, lty=2)

hist(postLMA$sigmaChillSp, main = paste("sigmaChillSp is " , signif(sigmaChillSp,3), sep = ""))
abline(v = sigmaChillSp, col="red", lwd=3, lty=2)

hist(postLMA$sigmaForceSp, main = paste("sigmaForceSp is " , signif(sigmaForceSp,3), sep = ""))
abline(v = sigmaForceSp, col="red", lwd=3, lty=2)

hist(postLMA$sigmaPhotoSp, main = paste("sigmaPhotoSp is " , signif(sigmaPhotoSp,3), sep = ""))
abline(v = sigmaPhotoSp, col="red", lwd=3, lty=2)

# png("figures/simulatedPairs.png")
#pairs(mdl.LMA, pars = c("muForceSp", "muChillSp", "muPhotoSp", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
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
  priorCheckTrait$species <- rep(1:Nspp, each = Nrep)
  priorCheckTrait$transect <- rep(1:Ntran, each = Nrep)
  
  #traitSLA <- rnorm(Ntrt, 20, 5)
  
  #Make this the name of the full vector of sla per species values - alphaTraitSp 
  #priorCheckTrait$alphaTraitSp <-  rep(rep(trt.dat$mu_grand_sp, times = nRepPrior))
  leafMass <- trtPheno[complete.cases(trtPheno$lma),]
  
  specieslist <- sort(unique(trtPheno$species))
  sitelist <- sort(unique(trtPheno$transect))

lma.data <- list(yTraiti = leafMass$lma,
                 N =nrow(leafMass),
                 n_spec = length(specieslist),
                 trait_species = as.numeric(as.factor(leafMass$species)),
                 lati = leafMass$lat,
                 tranE = as.numeric(leafMass$dumE),
                 Nph = nrow(pheno.t),
                 phenology_species = as.numeric(as.factor(pheno.t$species)),
                 yPhenoi = pheno.t$bb,
                 forcei = pheno.t$force.z2,
                 chilli = pheno.t$chillport.z2,
                 photoi = pheno.t$photo.z2,
                 prior_mu_sp_mu = 0.5,
                 prior_mu_sp_sigma = 5, #widened
                 prior_bTranE_mu = 0,
                 prior_bTranE_sigma = 5,
                 prior_sigma_sp_mu = 1,
                 prior_sigma_sp_sigma = 5,
                 prior_sigma_traity_mu = 5,
                 prior_sigma_traity_sigma = 2,
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
    
   # muGrand <- rtruncnorm(1, a = 0, mean = lma.data$prior_mu_grand_mu, sd = lma.data$prior_mu_grand_sigma)
    sigmaSp <- rtruncnorm(1, a = 0, mean = lma.data$prior_sigma_sp_mu, sd = lma.data$prior_sigma_sp_sigma)
    #sigmaSite <- rtruncnorm(1, a = 0, mean = lma.data$prior_sigma_site_mu, sd = lma.data$prior_sigma_site_sigma)
    alphaTraitSp <- rnorm(Nspp, 0, sigma.species)
    priorCheckTrait$alphaTraitSp[priorCheckTrait$simRep == ir] <- rep(alphaTraitSp, each = Nrep)
    
    muSp <- rnorm(Nspp, 0, sigma.species)
    priorCheckTrait$muSp[priorCheckTrait$simRep == ir] <- rep(muSp, each = Nrep)
    
    mutranlat <- rnorm(Ntran, lma.data$prior_bTranE_mu, lma.data$prior_bTranE_sigma)
    priorCheckTrait$mutranlat[priorCheckTrait$simRep == ir] <- rep(mutranlat, each = nRep)
    
    #general varience
    priorCheckTrait$sigmaTrait_y[priorCheckTrait$simRep == ir] <- rtruncnorm(lma.data$prior_sigma_traity_mu,  a = 0, lma.data$prior_sigma_traity_sigma)
    priorCheckTrait$e[priorCheckTrait$simRep == ir] <- rnorm(Ntrt, 0, sigmaTrait_y)
    
    priorCheckTrait$yTraiti <- priorCheckTrait$muSp + priorCheckTrait$mutranlat + priorCheckTrait$e
  }# end simulating new priors, from here vectorize code
  
  #Final values
  
  priorCheckTraityTraiti <- priorCheckTrait[complete.cases(priorCheckTrait$yTraiti),]
  
  png("figures/density_Trait_Prior_joint_lma.png")
  plot(density(priorCheckTraityTraiti$yTraiti))
  dev.off()
  
  png("figures/GrandSp_PlotPrior_joint_lma.png")
  plot(priorCheckTrait$yTraiti ~ priorCheckTrait$muGrandSp, xlab = "muGrandSp", ylab = "Trait")
  dev.off()
  
  png("figures/MuSp_PlotPrior_joint_lma.png")
  plot(priorCheckTrait$yTraiti ~ priorCheckTrait$muSp, xlab = "MuSp", ylab = "Trait")
  dev.off()
  
  png("figures/Mutransect_PlotPrior_joint_lma.png")
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
    #Species means
    sigmaPhenoSp <- rtruncnorm(1, a = 0, mean = lma.data$prior_sigmaPhenoSp_mu, sd = lma.data$prior_sigmaPhenoSp_sigma)
    muPhenoSp <- rnorm(1, lma.data$prior_muPhenoSp_mu, lma.data$prior_muPhenoSp_sigma)
    alphaPhenoSp <- rnorm(n_spec, muPhenoSp, sigmaPhenoSp)
    priorCheckPheno$alphaPhenoSp[priorCheckPheno$simRep == ir] <- rep(alphaPhenoSp, each = nRep)
    
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
    
    muChillSp <-  rnorm(1,lma.data$prior_sigmaChillSp_mu,lma.data$prior_sigmaChillSp_sigma)
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
  priorCheckPheno$yMu <-  priorCheckPheno$alphaPhenoSp +  priorCheckPheno$betaForceSp * priorCheckPheno$forcei +  priorCheckPheno$betaPhotoSp* priorCheckPheno$photoi + priorCheckPheno$betaChillSp * priorCheckPheno$chilli
  
  #Final values
  priorCheckPheno$yPhenoi <- priorCheckPheno$yMu + priorCheckPheno$e
  
  head(priorCheckPheno)
  
  plot(priorCheckPheno$betaForceSp ~ priorCheckPheno$alphaTraitSp )
  priorCheckPheno_posF <- priorCheckPheno[priorCheckPheno$betaForceSp > 0,]
  plot(priorCheckPheno_posF$betaForceSp ~ priorCheckPheno_posF$alphaTraitSp )
  
  png("figures/densityYPrior_joint_lma.png")
  plot(density(priorCheckPheno$yPhenoi))
  dev.off()
  
  png("figures/photoPlotPrior_joint_lma.png")
  plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$photoi, xlab = "Photoperiod", ylab = "Phenological Date")
  dev.off()
  
  png("figures/forcingPlotPrior_joint_lma.png")
  plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$forcei, xlab = "Forcing", ylab = "Phenological Date")
  dev.off()
  
  png("figures/chillingPlotPrior_joint_lma.png")
  plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$chilli, xlab = "Chillina", ylab = "Phenological Date")
  dev.off()
  


if(BayesSweave == TRUE){
  #For the BayesClass sweave documents 
  setwd("/home/faith/Documents/github/bayes2020/Projects/Faith/traitorsModel")
}
