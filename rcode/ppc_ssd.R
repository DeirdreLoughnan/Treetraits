#Started May 6, 2022 by Deirdre

# Aim of this code is to perform several ppc for the trait pheno model for NAm tree traits
# This model includes transect as a dummy variable but otherwise is the same as the traitors model

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

pheno$transect <- as.numeric(pheno$transect)

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

trtPheno$latitude <- as.numeric(trtPheno$latitude)
trtPheno$transect <- as.numeric(trtPheno$transect)
#########################################################

load("output/ssdDummyIntGrand.Rdata")

post<- rstan::extract(mdlSSD)

postSSD <- data.frame(post)

cueTrt <- postSSD[, colnames(postSSD) %in% c("mu_grand","b_tranE", "b_tranlat","sigma_sp", "sigma_traity")]

mcmc_intervals(cueTrt) 

cueEffects <- postSSD[, colnames(postSSD) %in% c("muPhenoSp", "muForceSp", "muChillSp", "muPhotoSp", "sigmapheno_y")]

mcmc_intervals(cueEffects) 
# + theme_classic() + 
# labs(title = "main intercept, cue slopes and general error")

###############

postSSD_muSp <- postSSD[,colnames(postSSD) %in% grep( "muSp", colnames(postSSD), value = TRUE)]
colnames(postSSD_muSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postSSD_muSp) + 
  geom_vline(xintercept = mean(postSSD_muSp), linetype="dotted", color = "grey")  +
  theme_classic() 

postSSD_alpaForceSp <- postSSD[,colnames(postSSD) %in% grep( "alphaForceSp", colnames(postSSD), value = TRUE)]
colnames(postSSD_alpaForceSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postSSD_alpaForceSp) + 
  geom_vline(xintercept = mean(postSSD$muForceSp), linetype="dotted", color = "grey")  +
  theme_classic() + 
  labs(subtitle = paste0("Mean muForceSp was ", round(mean(postSSD$muForceSp),3)),
    title = "muForceSp - species forcing slopes no trait")

postSSD_betaForceSp <- postSSD[,colnames(postSSD) %in% grep( "betaForceSp", colnames(postSSD), value = TRUE)]
colnames(postSSD_betaForceSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postSSD_betaForceSp) + 
  theme_classic() + 
  labs(title = "betaForceSp - Species forcing slopes with trait value")

#Different species slopes for chilling, without the effect of trait
postSSD_alphaChillSp <- postSSD[,colnames(postSSD) %in% grep( "alphaChillSp", colnames(postSSD), value = TRUE)]
colnames(postSSD_alphaChillSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postSSD_alphaChillSp) + 
  geom_vline(xintercept = mean(postSSD$muChillSp), linetype="dotted", color = "grey")  +
  theme_classic() + 
  labs(subtitle = paste0("Mean muChillSp was ", round(mean(postSSD$muChillSp),3)),
    title = "alphaChillSp - Species chill slopes no trait")

#Different species slopes for forcing, with the effect of trait
postSSD_betaChillSp <- postSSD[,colnames(postSSD) %in% grep( "betaChillSp", colnames(postSSD), value = TRUE)]
colnames(postSSD_betaChillSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postSSD_betaChillSp) + 
  # theme_classic() + 
  labs(title = "betaChillSp - Species chilling slopes with trait value")

#Different species slopes for photoperiod, without the effect of trait
postSSD_alphaPhotoSp <- postSSD[,colnames(postSSD) %in% grep( "alphaPhotoSp", colnames(postSSD), value = TRUE)]
colnames(postSSD_alphaPhotoSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postSSD_alphaPhotoSp) + 
  geom_vline(xintercept = mean(postSSD$muPhotoSp), linetype="dotted", color = "grey")  +
  theme_classic() + 
  labs(subtitle = paste0("Mean muPhotoSp was ", round(mean(postSSD$muPhotoSp),3)),
    title = "muPhotoSp - Species photo period slopes no trait")

#Different species slopes for forcing, with the effect of trait
postSSD_betaPhotoSp <- postSSD[,colnames(postSSD) %in% grep( "betaPhotoSp", colnames(postSSD), value = TRUE)]
colnames(postSSD_betaPhotoSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postSSD_betaPhotoSp) + 
  theme_classic() + 
  labs(title = "betaPhotoSp - Species photoperiod slopes with trait value")

#Different species slopes for forcing only the effect of trait
postSSD_betaTraitx <- postSSD[,colnames(postSSD) %in% grep( "betaTraitx", colnames(postSSD), value = TRUE)]

mcmc_intervals(postSSD_betaTraitx) + 
  theme_classic() + 
  labs(title = "effect's of traits on cue slopes")


#################################################################################

## Generating test data:
# Model has transect as a dummy variable (0 and 1) and latitude as a continuous variable
# Nrep <- 15 # rep per trait
# Npop <- 8
# Ntran <- 2
# Nspp <- 40 # number of species 
# 
# # First making a data frame for the test trait data
# Ntrt <- Nspp * Npop * Nrep# total number of traits observations
# Ntrt
# 
# # prior values
# mu.grand <- 0.5 # the grand mean of the SSD model
# sigma.species <- 1 # we want to keep the variation across spp. high
# sigma.transect <- 5
# b.tranE <- 0
# b.tranlat <- 0
# sigma_traity <- 5
# 
# #make a dataframe for SSD
# trt.dat <- data.frame(matrix(NA, Ntrt, 1))
# names(trt.dat) <- c("rep")
# trt.dat$rep <- c(1:Nrep)
# trt.dat$species <- rep(1:Nspp, each = Nrep)
# trt.dat$pop <- rep(1:Npop, each = Nspp*Nrep)
# trt.dat$tran <- rep(1:Ntran, each = 4*Nrep*Nspp)
# 
# 
# lati <- rnorm(8, 50, 5)
# trt.dat$lat <- rep(lati, each = Nrep*Nspp)
# trt.dat$lat <- as.numeric(trt.dat$lat)
# 
# mu.tranE <- 4
# trt.dat$dumE <- as.numeric(ifelse(trt.dat$tran == "1","0","1"))
# trt.dat$mutranE <- mu.tranE*trt.dat$dumE
# 
# mu.tranlat  <- 2
# trt.dat$alpha.tranlat <- mu.tranlat*(trt.dat$dumE*trt.dat$lat)
# 
# sigma.species <- 5 # we want to keep the variation across spp. high
# #the alphaTraitSp in Faiths original code:
# mu.trtsp <- rnorm(Nspp, 0, sigma.species)
# trt.dat$mu.trtsp <- rep(mu.trtsp, each = Nrep) #adding ht data for ea. sp
# 
# # general variance
# trt.var <- 5 #sigmaTrait_y in the stan code
# trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)
# 
# # generate yhat - SSDs -  for this first trt model
# #trt.dat$yTraiti <- mu.grand + trt.dat$alphaTraitSp + trt.dat$mutransect + trt.dat$trt.er
# mu.grand.sp <- mu.trtsp + mu.grand
# 
# for (i in 1:Ntrt){
#   trt.dat$mu.grand.sp[i] <-  trt.dat$mu.trtsp[i] +  mu.grand
# }
# 
# for (i in 1:Ntrt){
#   trt.dat$yTraiti[i] <- trt.dat$mu.grand.sp[i]  + trt.dat$mutranE[i]+ mu.tranlat * (trt.dat$dumE[i]*trt.dat$lat[i]) + trt.dat$trt.er[i]
# }
# 
# hist( trt.dat$yTraiti) # bimodal - does split across the two transects
# temp <- subset(trt.dat, dumE == "1")
# hist(temp$yTraiti) 
# 
# hist(trt.dat$mu.grand.sp)
# 
# 
# #### Pheno data generation ##############################
# #All Negative betatraitX Values 
# #--------------------------------------
# Nspp <- 40 # same number of species as teh traits section fo the model 
# #Number of repeat observations per species
# nRep <- 15
# nChill <- 2 # sm high low, mp high low
# nPhoto <- 2 # high low
# nForce <- 2 # high amd low
# nPop <- 2
# #Overall number of pbservations (rows)
# Nph <- n_spec * nRep*nChill * nForce * nPhoto * nPop # for phenology model
# 
# #Make a data frame for input phenology simulation data
# pheno.dat <- data.frame(matrix(NA, Nph, 2))
# names(pheno.dat) <- c("rep","species")
# pheno.dat$rep <- c(1:Nph)
# pheno.dat$species <- rep(c(1:n_spec), each = nRep)
# 
# #Simulate mean SLA offset data per species (not mean value)
# #meanSLA <- rnorm(n_spec, 0, 5)
# #Make this the name of the full vector of sla per species values - alphaTraitSp 
# pheno.dat$alphaTraitSp <- rep(mu.trtsp + mu.grand, each = nRep)
# 
# forceTrt <- c("LF", "HF")
# forcei <- rnorm(Nforce, 1, 2)
# pheno.dat$force <- rep(forceTrt, each = Nchill*Nphoto*Npop*Nrep*Nspp)
# pheno.dat$forcei <- rep(forcei, each = Nchill*Nphoto*Npop*Nrep*Nspp)
# 
# chillTrt <- c("LC", "HC")
# chilli <- rnorm(Nchill, 1, 5)
# pheno.dat$chill <- rep(rep(chillTrt, each = Nphoto*Npop*Nrep*Nspp), times = Nforce)
# pheno.dat$chilli <- rep(rep(chilli, each = Nphoto*Npop*Nrep*Nspp), times = Nforce)
# 
# photoTrt <- c("LP", "HP")
# photoi <- rnorm(Nphoto, 1, 2)
# pheno.dat$photo <- rep(rep(photoTrt, each = Npop*Nrep*Nspp), times = Nforce*Nchill)
# pheno.dat$photoi <- rep(rep(photoi, each = Npop*Nrep*Nspp), times = Nforce*Nchill)
# # #Simulate cues (z scored)
# # pheno.dat$forcei <- rnorm(Nph, 1, 1)
# # 
# # pheno.dat$photoi <- rnorm(Nph, 1, 0.5) # less photoperiod 
# # pheno.dat$chilli <- rnorm(Nph, 1, 1) #more chilling
# 
# # Parameter Values
# #Species means
# sigmaPhenoSp <- 5
# muPhenoSp <- 40
# alphaPhenoSp <- rnorm(n_spec, muPhenoSp, sigmaPhenoSp)
# pheno.dat$alphaPhenoSp <- rep(alphaPhenoSp, each = nRep)
# 
# #Cue effects
# betaTraitxForce <- 0
# betaTraitxPhoto <- 0
# betaTraitxChill <- 0
# 
# #Species level slopes sans trait data
# muForceSp <- -15
# sigmaForceSp <- 5
# alphaForceSp <- rnorm(n_spec, muForceSp, sigmaForceSp)
# pheno.dat$alphaForceSp <- rep(alphaForceSp, each = nRep)
# 
# muPhotoSp <- -15
# sigmaPhotoSp <- 5
# alphaPhotoSp <- rnorm(n_spec, muPhotoSp, sigmaPhotoSp)
# pheno.dat$alphaPhotoSp <- rep(alphaPhotoSp, each = nRep)
# 
# muChillSp <- -15
# sigmaChillSp <- 5
# alphaChillSp <- rnorm(n_spec, muChillSp, sigmaChillSp)
# pheno.dat$alphaChillSp <- rep(alphaChillSp, each = nRep)
# 
# #general varience
# sigmapheno_y <- 10
# pheno.dat$e <- rnorm(Nph, 0, sigmapheno_y)
# 
# #slopes for each cue, combining trait and non-trait aspect of the slope.
# #Make columns to put data 
# pheno.dat$betaForceSp <- NA
# pheno.dat$betaPhotoSp <- NA
# pheno.dat$betaChillSp <- NA
# 
# for (iSp in 1:n_spec){
#   
#   #iSp <- 1
#   #Select species data of interest 
#   
#   #Forcing
#   betaForceSp <- alphaForceSp[iSp] + betaTraitxForce * mu.grand.sp[iSp]
#   pheno.dat$betaForceSp[pheno.dat$species == iSp] <- betaForceSp
#   
#   #chilling
#   betaChillSp <- alphaChillSp[iSp] + betaTraitxChill * mu.grand.sp[iSp]
#   pheno.dat$betaChillSp[pheno.dat$species == iSp] <- betaChillSp
#   
#   #photoperiod
#   betaPhotoSp <- alphaPhotoSp[iSp] + betaTraitxPhoto * mu.grand.sp[iSp]
#   pheno.dat$betaPhotoSp[pheno.dat$species == iSp] <- betaPhotoSp
#   
#   
# }
# 
# #Run full model to get mean simulated y values
# for (i in 1:Nph){
#   pheno.dat$yMu[i] <-  pheno.dat$alphaPhenoSp[i] +  pheno.dat$betaForceSp[i] * pheno.dat$forcei[i] +  pheno.dat$betaPhotoSp[i] * pheno.dat$photoi[i] + pheno.dat$betaChillSp[i] * pheno.dat$chilli[i]
# }
# 
# #Final values
# pheno.dat$yPhenoi <- pheno.dat$yMu + pheno.dat$e
# 
# #What does the data look like?
# plot(density(pheno.dat$yPhenoi))
# 
# #Plot trait values againt forcing slopes
# 
# head(pheno.dat)
# 
# plot(pheno.dat$betaForceSp ~ pheno.dat$alphaTraitSp)
# plot(pheno.dat$alphaForceSp ~ pheno.dat$alphaTraitSp)
# 
# plot(pheno.dat$betaChillSp ~ pheno.dat$alphaTraitSp)
# plot(pheno.dat$alphaChillSp ~ pheno.dat$alphaTraitSp)
# 
# plot(pheno.dat$betaPhotoSp ~ pheno.dat$alphaTraitSp)
# plot(pheno.dat$alphaPhotoSp ~ pheno.dat$alphaTraitSp)

#######################################################

#png("figures/simPosteriorHist.png")
# priors
mu.grand <- 1 # the grand mean of the SSD model
sigma.species <- 5 # we want to keep the variaiton across spp. high

b.tranE <- 0
b.tranlat <- 0
sigma_traity <- 2

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

postSSD <- rstan::extract(mdlSSD)
##Compare results to simulated values
par(mfrow=c(1,1))
##Compare results to simulated values

hist(postSSD$mu_grand, main = paste("muPhenoSp is " , signif(mu.grand,3), sep = ""), xlim = c(-1,1), col=rgb(0,0,1,1/4))
hist(rnorm(1000, 1,5), col=rgb(1,0,1,1/4), add = T)
abline(v = mu.grand, col="red", lwd=3, lty=2)
# try 0.5,0.5
hist(postSSD$b_tranE, main = paste("muForceSp is " , signif(muForceSp,3), sep = ""),  xlim = c(-5,5),col=rgb(0,0,1,1/4))
hist(rnorm(1000, 0,5), col=rgb(1,0,1,1/4), add = T)
abline(v = b.tranE, col="red", lwd=3, lty=2)

hist(postSSD$b_tranlat, main = paste("muChillSp is " , signif(muChillSp,3), sep = ""), col=rgb(0,0,1,1/4))
hist(rnorm(1000, 0,10), col=rgb(1,0,1,1/4), add = T)
abline(v = b.tranlat, col="red", lwd=3, lty=2)

hist(postSSD$sigma_sp, main = paste("muPhotoSp is " , signif(muPhotoSp,3), sep = ""), col=rgb(0,0,1,1/4), xlim =c(-10,10))
hist(rnorm(1000, 5,5), col=rgb(1,0,1,1/4), add = T)
abline(v = sigma.species, col="red", lwd=3, lty=2)
# try 0.5,0.5

hist(postSSD$sigma_traity, main = paste("sigmapheno_y is " , signif(sigmapheno_y,3), sep = ""), col=rgb(0,0,1,1/4),  xlim = c(-5,5))
hist(rnorm(1000, 0.5,0.5), col=rgb(1,0,1,1/4), add = T)
abline(v = sigmapheno_y, col="red", lwd=3, lty=2)
# try 0.5,0.5
hist(postSSD$muPhenoSp, main = paste("muPhenoSp is " , signif(muPhenoSp,3), sep = ""), xlim = c(0,100), col=rgb(0,0,1,1/4))
hist(rnorm(1000, 40,10), col=rgb(1,0,1,1/4), add = T)
abline(v = muPhenoSp, col="red", lwd=3, lty=2)

hist(postSSD$muForceSp, main = paste("muForceSp is " , signif(muForceSp,3), sep = ""),  xlim = c(-75,75),col=rgb(0,0,1,1/4))
hist(rnorm(1000, -15,10), col=rgb(1,0,1,1/4), add = T)
abline(v = muForceSp, col="red", lwd=3, lty=2)

hist(postSSD$muChillSp, main = paste("muChillSp is " , signif(muChillSp,3), sep = ""),  xlim = c(-75,75), col=rgb(0,0,1,1/4))
hist(rnorm(1000, -15,10), col=rgb(1,0,1,1/4), add = T)
abline(v = muChillSp, col="red", lwd=3, lty=2)

hist(postSSD$muPhotoSp, main = paste("muPhotoSp is " , signif(muPhotoSp,3), sep = ""), col=rgb(0,0,1,1/4),  xlim = c(-75,75))
hist(rnorm(1000, -15,10), col=rgb(1,0,1,1/4), add = T)
abline(v = muPhotoSp, col="red", lwd=3, lty=2)

hist(postSSD$sigmapheno_y, main = paste("sigmapheno_y is " , signif(sigmapheno_y,3), sep = ""), col=rgb(0,0,1,1/4),  xlim = c(-10,30))
hist(rnorm(1000, 10,5), col=rgb(1,0,1,1/4), add = T)
abline(v = sigmapheno_y, col="red", lwd=3, lty=2)

hist(postSSD$betaTraitxForce, main = paste("betaTraitxForce is " , signif(betaTraitxForce,3), sep = ""), col=rgb(0,0,1,1/4),  xlim = c(-5,5))
hist(rnorm(1000, 0,1), col=rgb(1,0,1,1/4), add = T)
abline(v = betaTraitxForce, col="red", lwd=3, lty=2)
#
hist(postSSD$betaTraitxChill, main = paste("betaTraitxChill is " , signif(betaTraitxChill,3), sep = ""), col=rgb(0,0,1,1/4),  xlim = c(-5,5))
hist(rnorm(1000, 0,1), col=rgb(1,0,1,1/4), add = T)
abline(v = betaTraitxChill, col="red", lwd=3, lty=2)
#
hist(postSSD$betaTraitxPhoto, main = paste("betaTraitxPhoto is " , signif(betaTraitxPhoto,3), sep = ""), col=rgb(0,0,1,1/4),  xlim = c(-5,5))
hist(rnorm(1000, 0,1), col=rgb(1,0,1,1/4), add = T)
abline(v = betaTraitxPhoto, col="red", lwd=3, lty=2)


hist(postSSD$sigmaChillSp, main = paste("sigmaChillSp is " , signif(sigmaChillSp,3), sep = ""), col=rgb(0,0,1,1/4), xlim = c(-10,20))
hist(rnorm(1000, 5,10), col=rgb(1,0,1,1/4), add = T)
abline(v = sigmaChillSp, col="red", lwd=3, lty=2)
# try 5,10
hist(postSSD$sigmaForceSp, main = paste("sigmaForceSp is " , signif(sigmaForceSp,3), sep = ""), col=rgb(0,0,1,1/4),  xlim = c(-10,10))
hist(rnorm(1000, 5,5), col=rgb(1,0,1,1/4), add = T)
abline(v = sigmaForceSp, col="red", lwd=3, lty=2)

hist(postSSD$sigmaPhotoSp, main = paste("sigmaPhotoSp is " , signif(sigmaPhotoSp,3), sep = ""), col=rgb(0,0,1,1/4),  xlim = c(-15,15))
hist(rnorm(1000, 5,5), col=rgb(1,0,1,1/4), add = T)
abline(v = sigmaPhotoSp, col="red", lwd=3, lty=2)

#pairs(mdl.SSD, pars = c("muForceSp", "muChillSp", "muPhotoSp", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
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

priorCheckTrait$lat.z <- (priorCheckTrait$lat-mean(priorCheckTrait$lat,na.rm=TRUE))/(sd(priorCheckTrait$lat,na.rm=TRUE)*2)

#Make this the name of the full vector of sla per species values - alphaTraitSp 
#priorCheckTrait$alphaTraitSp <-  rep(rep(trt.dat$mu_grand_sp, times = nRepPrior))
stemDen <- trtPheno[complete.cases(trtPheno$ssd),]

specieslist <- sort(unique(trtPheno$species))
sitelist <- sort(unique(trtPheno$transect))

mu.grand <- 1 # the grand mean of the SSD model
sigma.species <- 5 # we want to keep the variation across spp. high

b.tranE <- 0
b.tranlat <- 0
sigma_traity <- 2

SSD.data <- list(yTraiti = stemDen$SSD,
  N = nrow(stemDen),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(stemDen$species)),
  n_tran = length(unique(stemDen$transect)),
  tranE = stemDen$transect,
  lati <- stemDen$latitude,
  prior_mu_grand_mu = 1,
  prior_mu_grand_sigma = 5, 
  prior_sigma_sp_mu = 5,
  prior_sigma_sp_sigma = 5,
  prior_b_tranlat_mu = 0,
  prior_b_tranlat_sigma = 10,
  prior_b_tranE_mu = 0,
  prior_b_tranE_sigma = 5,
  prior_sigma_traity_mu = 2,
  prior_sigma_traity_sigma = 5,
  ## Phenology
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2,
  prior_muForceSp_mu = -15,
  prior_muForceSp_sigma = 10, 
  prior_muChillSp_mu = -15,
  prior_muChillSp_sigma = 10,
  prior_muPhotoSp_mu = -15,
  prior_muPhotoSp_sigma = 10,
  prior_muPhenoSp_mu = 40,
  prior_muPhenoSp_sigma = 10,
  prior_sigmaForceSp_mu = 5,
  prior_sigmaForceSp_sigma = 5,
  prior_sigmaChillSp_mu = 5,
  prior_sigmaChillSp_sigma = 5, 
  prior_sigmaPhotoSp_mu = 5,
  prior_sigmaPhotoSp_sigma = 5,
  prior_sigmaPhenoSp_mu = 5, 
  prior_sigmaPhenoSp_sigma = 5, 
  prior_betaTraitxForce_mu = 0,
  prior_betaTraitxForce_sigma = 1,
  prior_betaTraitxChill_mu = 0,
  prior_betaTraitxChill_sigma = 1,
  prior_betaTraitxPhoto_mu = 0,
  prior_betaTraitxPhoto_sigma = 1,
  prior_sigmaphenoy_mu = 5,
  prior_sigmaphenoy_sigma = 5 
) 

for (ir in 1:nRepPrior){
  # Parameter Values
  
  
  muGrand <- rtruncnorm(1, a = 0, mean = SSD.data$prior_mu_grand_mu, sd = SSD.data$prior_mu_grand_sigma)
  sigmaSp <- rtruncnorm(1, a = 0, mean = SSD.data$prior_sigma_sp_mu, sd = SSD.data$prior_sigma_sp_sigma)
  b_tranlat <- rtruncnorm(1, a = 0, mean = SSD.data$prior_b_tranlat_mu, sd = SSD.data$prior_b_tranlat_sigma)
  b.tranE <- rtruncnorm(1, a = 0, mean = SSD.data$prior_b_tranE_mu, sd = SSD.data$prior_b_tranE_sigma)
  
  muSp <- rnorm(Nspp, 0, SSD.data$prior_sigma_sp_mu)
  priorCheckTrait$muSp[priorCheckTrait$simRep == ir] <- rep(muSp, each = nRep)
  
  #general varience
  priorCheckTrait$sigma_traity[priorCheckTrait$simRep == ir] <- rtruncnorm(SSD.data$prior_sigma_traity_mu,  a = 0, SSD.data$prior_sigma_traity_sigma)
  priorCheckTrait$e[priorCheckTrait$simRep == ir] <- rnorm(Ntrt, 0, SSD.data$prior_sigma_traity_mu)
  
  priorCheckTrait$yTraiti <- muGrand + priorCheckTrait$muSp + b_tranlat *(priorCheckTrait$lat.z * priorCheckTrait$tran) + b.tranE * priorCheckTrait$tran + priorCheckTrait$e
}# end simulating new priors, from here vectorize code

#Final values
priorCheckTrait$muGrandSp <- muGrand + priorCheckTrait$muSp

priorCheckTrait$Traiti <- priorCheckTrait[complete.cases(priorCheckTrait$yTraiti),]

#png("figures/density_Trait_Prior_joint_SSD.png")
plot(density(priorCheckTrait$yTraiti))
#dev.off()

#png("figures/GrandSp_PlotPrior_joint_SSD.png")
plot(priorCheckTrait$yTraiti ~ priorCheckTrait$muGrandSp, xlab = "muGrandSp", ylab = "Trait")
#dev.off()

#png("figures/MuSp_PlotPrior_joint_SSD.png")
plot(priorCheckTrait$yTraiti ~ priorCheckTrait$muSp, xlab = "MuSp", ylab = "Trait")
#dev.off()

#png("figures/Mutransect_PlotPrior_joint_SSD.png")
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

# priorCheckPheno$simRep <- rep(1:nRepPrior, each = Nph)
# priorCheckPheno$rep <- rep(c(1:nRep), times = n_spec * nChill * nForce * nPhoto * nPop*nSites*nRepPrior)
# priorCheckPheno$species <- rep(rep(c(1:n_spec), each = nRep*nChill * nForce * nPhoto * nPop*nSites), times = nRepPrior)
# priorCheckPheno$pop <- rep(1:nPop, each = nRep*nChill * nForce * nPhoto*nSites,  times = n_spec*nRepPrior)
# priorCheckPheno$sites <- rep(1:nsites, each = nRep*nChill * nForce * nPhoto * nPop,  times = n_spec*nRepPrior)
# 
#Simulate SLA data per species
muGrandSp <- muGrand + muSp
#Make this the name of the full vector of sla per species values - alphaTraitSp 
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

for (ir in 1:nRepPrior){
  # Parameter Values
  # ir <- 1
  
  sigmaPhenoSp <- rtruncnorm(1, a = 0, mean = SSD.data$prior_sigmaPhenoSp_mu, sd = SSD.data$prior_sigmaPhenoSp_sigma)
  muPhenoSp <- rnorm(1, SSD.data$prior_muPhenoSp_mu, SSD.data$prior_muPhenoSp_sigma)
  alphaPhenoSp <- rnorm(n_spec, muPhenoSp, sigmaPhenoSp)
  priorCheckPheno$alphaPhenoSp[priorCheckPheno$sim == ir] <- rep(alphaPhenoSp, each = nRep)
  
  #Cue effects
  priorCheckPheno$betaTraitxForce[priorCheckPheno$simRep == ir] <- rnorm(1,SSD.data$prior_betaTraitxForce_mu,SSD.data$prior_betaTraitxForce_sigma)
  priorCheckPheno$betaTraitxPhoto[priorCheckPheno$simRep == ir] <- rnorm(1,SSD.data$prior_betaTraitxPhoto_mu,SSD.data$prior_betaTraitxPhoto_sigma)
  priorCheckPheno$betaTraitxChill[priorCheckPheno$simRep == ir] <- rnorm(1,SSD.data$prior_betaTraitxChill_mu,SSD.data$prior_betaTraitxChill_sigma)
  
  #Species level slopes sans trait data
  muForceSp <- rnorm(1,SSD.data$prior_muForceSp_mu,  SSD.data$prior_muForceSp_sigma)
  sigmaForceSp <- rtruncnorm(1, a = 0,mean = SSD.data$prior_sigmaForceSp_mu,sd = SSD.data$prior_sigmaForceSp_sigma)
  alphaForceSp <- rnorm(n_spec, muForceSp, sigmaForceSp)
  priorCheckPheno$alphaForceSp[priorCheckPheno$simRep == ir] <- rep(alphaForceSp, each = nRep)
  
  muPhotoSp <- rnorm(1, SSD.data$prior_muPhotoSp_mu, SSD.data$prior_muPhotoSp_sigma)
  sigmaPhotoSp <- rtruncnorm(1, a = 0,mean = SSD.data$prior_sigmaPhotoSp_mu, sd = SSD.data$prior_sigmaPhotoSp_sigma )
  alphaPhotoSp <- rnorm(n_spec, muPhotoSp, sigmaPhotoSp)
  priorCheckPheno$alphaPhotoSp[priorCheckPheno$simRep == ir] <- rep(alphaPhotoSp, each = nRep)
  
  muChillSp <-  rnorm(1,SSD.data$prior_sigmaChillSp_mu,SSD.data$prior_sigmaChillSp_sigma)
  sigmaChillSp <- rtruncnorm(1, a = 0,mean = SSD.data$prior_sigmaChillSp_mu,sd = SSD.data$prior_sigmaChillSp_sigma)
  alphaChillSp <- rnorm(n_spec, muChillSp, sigmaChillSp)
  priorCheckPheno$alphaChillSp[priorCheckPheno$simRep == ir] <- rep(alphaChillSp, each = nRep)
  
  
  #general varience
  priorCheckPheno$sigmapheno_y[priorCheckPheno$simRep == ir] <- rtruncnorm(SSD.data$prior_sigmaphenoy_mu,  a = 0, SSD.data$prior_sigmaphenoy_sigma)
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

#png("figures/densityYPrior_joint_SSD.png")
plot(density(priorCheckPheno$yPhenoi))
#dev.off()

# png("figures/photoPlotPrior_joint_SSD.png")
plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$photoi, xlab = "Photoperiod", ylab = "Phenological Date")
#dev.off()

#png("figures/forcingPlotPrior_joint_SSD.png")
plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$forcei, xlab = "Forcing", ylab = "Phenological Date")
# dev.off()

# png("figures/chillingPlotPrior_joint_SSD.png")
plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$chilli, xlab = "Chilling", ylab = "Phenological Date")
# dev.off()



# rm(list=ls())
# options(stringsAsFactors = FALSE)
# 
# library(tidyr)
# library(plyr)
# library(dplyr)
# library(reshape2)
# library(rstan)
# library(bayesplot)# nice posterior check plots 
# library(shinystan)
# library(stringr)
# 
# if(length(grep("deirdreloughnan", getwd()) > 0)) { 
#   setwd("~/Documents/github/pheno_bc") 
# }  else{
#   setwd("/home/deirdre/pheno_bc") # for midge
# }
# 
# # Start by getting the pheno data
# dl <- read.csv("input/dl_allbb.csv")
# 
# temp <- str_split_fixed(dl$trt, "_", 3); head(temp)
# dl$chill<- temp[,1]
# dl$photo <- temp[,2]
# dl$force <- temp[,3]
# 
# dl.chill <- read.csv("input/chilling_values_Hope_Smithers.csv")
# 
# dl.wchill <- merge(dl, dl.chill, by = c("population","chill"))
# dl.wchill$lab3 <- dl.wchill$lab2
# dl.wchill$lab2 <- paste(dl.wchill$species, dl.wchill$population, dl.wchill$rep, sep = "_")
# 
# df <- read.csv("input/df_dxb_prepped_data.csv")
# df.chill <- read.csv("input/chilling_values_eastern.csv")
# df.wchill <- merge(df, df.chill, by =c("population","chill"))
# df.wchill <- df.wchill[, c("population", "chill","force","photo","lab2", "bb","species", "treatment","Chill_portions","Utah_Model")]
# 
# # mergeing the my data with DF
# pheno <- rbind.fill(dl.wchill, df.wchill)
# pheno$force.n <- pheno$force
# pheno$force.n[pheno$force.n == "HF"] <- "1"
# pheno$force.n[pheno$force.n == "LF"] <- "0"
# pheno$force.n <- as.numeric(pheno$force.n)
# 
# pheno$photo.n <- pheno$photo
# pheno$photo.n[pheno$photo.n == "HP"] <- "1"
# pheno$photo.n[pheno$photo.n == "LP"] <- "0"
# pheno$photo.n <- as.numeric(pheno$photo.n)
# 
# # pheno$site.n <- pheno$population
# # pheno$site.n[pheno$site.n == "sm"] <- "1"
# # pheno$site.n[pheno$site.n == "mp"] <- "2"
# # pheno$site.n[pheno$site.n == "HF"] <- "3"
# # pheno$site.n[pheno$site.n == "SH"] <- "4"
# # pheno$site.n <- as.numeric(pheno$site.n)
# 
# pheno$transect <- pheno$population
# pheno$transect[pheno$transect == "sm"] <- "0"
# pheno$transect[pheno$transect == "mp"] <- "0"
# pheno$transect[pheno$transect == "kl"] <- "0"
# pheno$transect[pheno$transect == "af"] <- "0"
# pheno$transect[pheno$transect == "GR"] <- "1"
# pheno$transect[pheno$transect == "HF"] <- "1"
# pheno$transect[pheno$transect == "SH"] <- "1"
# pheno$transect[pheno$transect == "WM"] <- "1"
# 
# # head(pheno)
# # #add dummy/ site level effects:
# # pheno <- pheno %>%
# #   mutate ( site2 = if_else(site.n == 2, 1, 0),
# #            site3 = if_else(site.n == 3, 1, 0),
# #            site4 = if_else(site.n == 4, 1, 0))
# 
# # standardize the 0/1 and standardize sites? 
# pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
# pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
# pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)
# 
# # pheno$site2.z2 <- (pheno$site2-mean(pheno$site2,na.rm=TRUE))/(sd(pheno$site2,na.rm=TRUE)*2)
# # pheno$site3.z2 <- (pheno$site3-mean(pheno$site3,na.rm=TRUE))/(sd(pheno$site3,na.rm=TRUE)*2)
# # pheno$site4.z2 <- (pheno$site4-mean(pheno$site4,na.rm=TRUE))/(sd(pheno$site4,na.rm=TRUE)*2)
# 
# #going to split it into analysis of terminal bb and lateral bb
# # Starting with the terminal buds:
# #pheno.term <- pheno[,c("tbb", "chill.n", "force.n", "photo.n", "site.n", "species", "lab2")]
# pheno.term <- pheno[,c("bb", "force.z2", "photo.z2", "population", "species", "lab2","Utah_Model","Chill_portions","chillport.z2","transect")] #"site2.z2", "site3.z2","site4.z2")]
# pheno.t <- pheno.term[complete.cases(pheno.term$bb), ] # none,great!
# length(unique(pheno.t$species))
# #pheno.t <- pheno.term[complete.cases(pheno.term$bb), ] # 1780 rows data 
# pheno.t$species <- tolower(pheno.t$species)
# pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
# sort(unique(pheno.t$species.fact)) # 47 bc two species occur in both transects
# 
# # Now get the trait data and subset to only include spp we have pheno data for:
# setwd("..//Treetraits")
# trtData <- read.csv("data/allTrt.csv", stringsAsFactors = FALSE)
# head(trtData)
# 
# phenoSp <- sort(unique(pheno.t$species))
# 
# trtPheno <- trtData[trtData$species %in% phenoSp, ]
# length(unique(trtPheno$species))
# 
# trtPheno$transect <- trtPheno$site
# trtPheno$transect[trtPheno$transect == "sm"] <- "0"
# trtPheno$transect[trtPheno$transect == "mp"] <- "0"
# trtPheno$transect[trtPheno$transect == "kl"] <- "0"
# trtPheno$transect[trtPheno$transect == "af"] <- "0"
# trtPheno$transect[trtPheno$transect == "GR"] <- "1"
# trtPheno$transect[trtPheno$transect == "HF"] <- "1"
# trtPheno$transect[trtPheno$transect == "SH"] <- "1"
# trtPheno$transect[trtPheno$transect == "WM"] <- "1"
# 
# #########################################################
# 
# fit <- readRDS("output/ssd_stanfit.RDS")
# 
# 
# post<- rstan::extract(fit)
# 
# postssd <- data.frame(post)
# 
# cueEffects <- postssd[, colnames(postssd) %in% c("muPhenoSp", "muForceSp", "muChillSp", "muPhotoSp", "sigmapheno_y")]
# 
# mcmc_intervals(cueEffects) + 
#   theme_classic() + 
#   labs(title = "main intercept, cue slopes and general error")
# 
# ###############
# postssd_alpaForceSp <- postssd[,colnames(postssd) %in% grep( "alphaForceSp", colnames(postssd), value = TRUE)]
# colnames(postssd_alpaForceSp) <- levels(as.factor(trtPheno$species))
# 
# mcmc_intervals(postssd_alpaForceSp) + 
#   geom_vline(xintercept = mean(postssd$muForceSp), linetype="dotted", color = "grey")  +
#   theme_classic() + 
#   labs(subtitle = paste0("Mean muForceSp was ", round(mean(postssd$muForceSp),3)),
#        title = "muForceSp - species forcing slopes no trait")
# 
# postssd_betaForceSp <- postssd[,colnames(postssd) %in% grep( "betaForceSp", colnames(postssd), value = TRUE)]
# colnames(postssd_betaForceSp) <- levels(as.factor(trtPheno$species))
# 
# mcmc_intervals(postssd_betaForceSp) + 
#   theme_classic() + 
#   labs(title = "betaForceSp - Species forcing slopes with trait value")
# 
# #Different species slopes for chilling, without the effect of trait
# postssd_alphaChillSp <- postssd[,colnames(postssd) %in% grep( "alphaChillSp", colnames(postssd), value = TRUE)]
# colnames(postssd_alphaChillSp) <- levels(as.factor(trtPheno$species))
# 
# mcmc_intervals(postssd_alphaChillSp) + 
#   geom_vline(xintercept = mean(postssd$muChillSp), linetype="dotted", color = "grey")  +
#   theme_classic() + 
#   labs(subtitle = paste0("Mean muChillSp was ", round(mean(postssd$muChillSp),3)),
#        title = "alphaChillSp - Species chill slopes no trait")
# 
# #Different species slopes for forcing, with the effect of trait
# postssd_betaChillSp <- postssd[,colnames(postssd) %in% grep( "betaChillSp", colnames(postssd), value = TRUE)]
# colnames(postssd_betaChillSp) <- levels(as.factor(trtPheno$species))
# 
# mcmc_intervals(postssd_betaChillSp) + 
#   theme_classic() + 
#   labs(title = "betaChillSp - Species chilling slopes with trait value")
# 
# #Different species slopes for photoperiod, without the effect of trait
# postssd_alphaPhotoSp <- postssd[,colnames(postssd) %in% grep( "alphaPhotoSp", colnames(postssd), value = TRUE)]
# colnames(postssd_alphaPhotoSp) <- levels(as.factor(trtPheno$species))
# 
# mcmc_intervals(postssd_alphaPhotoSp) + 
#   geom_vline(xintercept = mean(postssd$muPhotoSp), linetype="dotted", color = "grey")  +
#   theme_classic() + 
#   labs(subtitle = paste0("Mean muPhotoSp was ", round(mean(postssd$muPhotoSp),3)),
#        title = "muPhotoSp - Species photo period slopes no trait")
# 
# #Different species slopes for forcing, with the effect of trait
# postssd_betaPhotoSp <- postssd[,colnames(postssd) %in% grep( "betaPhotoSp", colnames(postssd), value = TRUE)]
# colnames(postssd_betaPhotoSp) <- levels(as.factor(trtPheno$species))
# 
# mcmc_intervals(postssd_betaPhotoSp) + 
#   theme_classic() + 
#   labs(title = "betaPhotoSp - Species photoperiod slopes with trait value")
# 
# 
# #Different species slopes for forcing only the effect of trait
# postssd_betaTraitx <- postssd[,colnames(postssd) %in% grep( "betaTraitx", colnames(postssd), value = TRUE)]
# 
# mcmc_intervals(postssd_betaTraitx) + 
#   theme_classic() + 
#   labs(title = "effect's of traits on cue slopes")
# 
# # require(bayesplot)
# # y <- pheno$bb 
# # yrep <-  postssd[,colnames(postssd) %in% grep( "y_hat", colnames(postssd), value = TRUE)]
# # yrepM <- colMeans(yrep)
# # 
# # ppc_dens_overlay(y, yrepM[1:50,])
# 
# #################################################################################
# Nrep <- 15 # rep per trait
# Ntransect <- 2
# Nspp <- 40 # number of species 
# 
# # First making a data frame for the test trait data
# Ntrt <- Nspp * Ntransect * Nrep # total number of traits observations
# Ntrt
# 
# mu.grand <- 1 # the grand mean of the height model
# sigma.species <- 5 # we want to keep the variaiton across spp. high
# sigma.transect <- 2
# sigmaTrait_y <- 2
# 
# #make a dataframe for height
# trt.dat <- data.frame(matrix(NA, Ntrt, 1))
# names(trt.dat) <- c("rep")
# trt.dat$rep <- c(1:Nrep)
# trt.dat$transect <- rep(c(1:Ntransect), each = Nspp)
# trt.dat$species <- rep(1:Nspp, Ntransect)
# 
# # now generating the species trait data, here it is for height
# #the alphaTraitSp in Faiths original code:
# alphaTraitSp <- rnorm(Nspp, 0, sigma.species)
# trt.dat$alphaTraitSp <- rep(alphaTraitSp, Ntransect) #adding ht data for ea. sp
# 
# #now generating the effects of transect
# mutransect <- rnorm(Ntransect, 0, sigma.transect) #intercept for each transect
# trt.dat$mutransect <- rep(mutransect, each = Nspp) # generate data for ea transect
# 
# # general variance
# trt.var <- 2 #sigmaTrait_y in the stan code
# trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)
# 
# # generate yhat - heights -  for this first trt model
# #trt.dat$yTraiti <- mu.grand + trt.dat$alphaTraitSp + trt.dat$mutransect + trt.dat$trt.er
# 
# for (i in 1:Ntrt){
#   trt.dat$mu_grand_sp[i] <-  trt.dat$alphaTraitSp[i] +  mu.grand
# }
# 
# for (i in 1:Ntrt){
#   trt.dat$yTraiti[i] <-  trt.dat$alphaTraitSp[i] + trt.dat$mutransect[i] +  mu.grand
# }
# 
# #Add grand mean to trait values 
# alphaTraitspFull <-  alphaTraitSp + mu.grand
# trt.dat$alphaTraitspFull <-trt.dat$alphaTraitSp + mu.grand
# 
# hist( trt.dat$yTraiti)
# hist(alphaTraitspFull)
# 
# 
# #### Pheno data generation ##############################
# #All Negative betatraitX Values 
# #--------------------------------------
# n_spec <-Nspp # same number of species as teh traits section fo the model 
# #Number of repeat observations per species
# nRep <- 15
# #Overall number of pbservations (rows)
# Nph <- n_spec * nRep # for phenology model
# 
# #Make a data frame for input phenology simulation data
# pheno.dat <- data.frame(matrix(NA, Nph, 2))
# names(pheno.dat) <- c("rep","species")
# pheno.dat$rep <- c(1:Nph)
# pheno.dat$species <- rep(c(1:n_spec), each = nRep)
# 
# #Simulate mean SLA offset data per species (not mean value)
# #meanSLA <- rnorm(n_spec, 0, 5)
# #Make this the name of the full vector of sla per species values - alphaTraitSp 
# pheno.dat$alphaTraitSp <- rep(alphaTraitSp + mu.grand, each = nRep)
# 
# #Simulate cues (z scored)
# pheno.dat$forcei <- rnorm(Nph, 1, 1)
# 
# pheno.dat$photoi <- rnorm(Nph, 1, 0.5) # less photoperiod 
# pheno.dat$chilli <- rnorm(Nph, 1, 1) #more chilling
# 
# # Parameter Values
# #Species means
# sigmaPhenoSp <- 10
# muPhenoSp <- 40
# alphaPhenoSp <- rnorm(n_spec, muPhenoSp, sigmaPhenoSp)
# pheno.dat$alphaPhenoSp <- rep(alphaPhenoSp, each = nRep)
# 
# #Cue effects
# betaTraitxForce <- -0.5 
# betaTraitxPhoto <- -0.2
# betaTraitxChill <- -0.4
# 
# #Species level slopes sans trait data
# muForceSp <- -0.4
# sigmaForceSp <- 4
# alphaForceSp <- rnorm(n_spec, muForceSp, sigmaForceSp)
# pheno.dat$alphaForceSp <- rep(alphaForceSp, each = nRep)
# 
# muPhotoSp <- -0.05
# sigmaPhotoSp <- 4
# alphaPhotoSp <- rnorm(n_spec, muPhotoSp, sigmaPhotoSp)
# pheno.dat$alphaPhotoSp <- rep(alphaPhotoSp, each = nRep)
# 
# muChillSp <- -0.6
# sigmaChillSp <- 4
# alphaChillSp <- rnorm(n_spec, muChillSp, sigmaChillSp)
# pheno.dat$alphaChillSp <- rep(alphaChillSp, each = nRep)
# 
# #general varience
# sigmapheno_y <- 5
# pheno.dat$e <- rnorm(Nph, 0, sigmapheno_y)
# 
# #slopes for each cue, combining trait and non-trait aspect of the slope.
# #Make columns to put data 
# pheno.dat$betaForceSp <- NA
# pheno.dat$betaPhotoSp <- NA
# pheno.dat$betaChillSp <- NA
# 
# for (iSp in 1:n_spec){
#   
#   #iSp <- 1
#   #Select species data of interest 
#   
#   #Forcing
#   betaForceSp <- alphaForceSp[iSp] + betaTraitxForce * alphaTraitspFull[iSp]
#   pheno.dat$betaForceSp[pheno.dat$species == iSp] <- betaForceSp
#   
#   #chilling
#   betaChillSp <- alphaChillSp[iSp] + betaTraitxChill * alphaTraitspFull[iSp]
#   pheno.dat$betaChillSp[pheno.dat$species == iSp] <- betaChillSp
#   
#   #photoperiod
#   betaPhotoSp <- alphaPhotoSp[iSp] + betaTraitxPhoto * alphaTraitspFull[iSp]
#   pheno.dat$betaPhotoSp[pheno.dat$species == iSp] <- betaPhotoSp
#   
#   
# }
# 
# #Run full model to get mean simulated y values
# for (i in 1:Nph){
#   pheno.dat$yMu[i] <-  pheno.dat$alphaPhenoSp[i] +  pheno.dat$betaForceSp[i] * pheno.dat$forcei[i] +  pheno.dat$betaPhotoSp[i] * pheno.dat$photoi[i] + pheno.dat$betaChillSp[i] * pheno.dat$chilli[i]
# }
# 
# #Final values
# pheno.dat$yPhenoi <- pheno.dat$yMu + pheno.dat$e
# 
# #What does the data look like?
# plot(density(pheno.dat$yPhenoi))
# 
# #Plot trait values againt forcing slopes
# 
# head(pheno.dat)
# 
# plot(pheno.dat$betaForceSp ~ pheno.dat$alphaTraitSp)
# plot(pheno.dat$alphaForceSp ~ pheno.dat$alphaTraitSp)
# 
# plot(pheno.dat$betaChillSp ~ pheno.dat$alphaTraitSp)
# plot(pheno.dat$alphaChillSp ~ pheno.dat$alphaTraitSp)
# 
# plot(pheno.dat$betaPhotoSp ~ pheno.dat$alphaTraitSp)
# plot(pheno.dat$alphaPhotoSp ~ pheno.dat$alphaTraitSp)
# 
# #######################################################
# 
# # png("figures/simPosteriorHist.png")
# # par(mfrow=c(3,4))
# #Compare results to simulated values
# hist(postssd$muPhenoSp, main = paste("muPhenoSp is " , signif(muPhenoSp,3), sep = ""), xlim = c(0,100))
# abline(v = muPhenoSp, col="red", lwd=3, lty=2)
# 
# hist(postssd$muForceSp, main = paste("muForceSp is " , signif(muForceSp,3), sep = ""))
# abline(v = muForceSp, col="red", lwd=3, lty=2)
# 
# hist(postssd$muChillSp, main = paste("muChillSp is " , signif(muChillSp,3), sep = ""))
# abline(v = muChillSp, col="red", lwd=3, lty=2)
# 
# hist(postssd$muPhotoSp, main = paste("muPhotoSp is " , signif(muPhotoSp,3), sep = ""))
# abline(v = muPhotoSp, col="red", lwd=3, lty=2)
# 
# hist(postssd$sigmapheno_y, main = paste("sigmapheno_y is " , signif(sigmapheno_y,3), sep = ""))
# abline(v = sigmapheno_y, col="red", lwd=3, lty=2)
# 
# hist(postssd$betaTraitxForce, main = paste("betaTraitxForce is " , signif(betaTraitxForce,3), sep = ""))
# abline(v = betaTraitxForce, col="red", lwd=3, lty=2)
# #
# hist(postssd$betaTraitxChill, main = paste("betaTraitxChill is " , signif(betaTraitxChill,3), sep = ""))
# abline(v = betaTraitxChill, col="red", lwd=3, lty=2)
# #
# hist(postssd$betaTraitxPhoto, main = paste("betaTraitxPhoto is " , signif(betaTraitxPhoto,3), sep = ""))
# abline(v = betaTraitxPhoto, col="red", lwd=3, lty=2)
# 
# hist(postssd$sigmaChillSp, main = paste("sigmaChillSp is " , signif(sigmaChillSp,3), sep = ""))
# abline(v = sigmaChillSp, col="red", lwd=3, lty=2)
# 
# hist(postssd$sigmaForceSp, main = paste("sigmaForceSp is " , signif(sigmaForceSp,3), sep = ""))
# abline(v = sigmaForceSp, col="red", lwd=3, lty=2)
# 
# hist(postssd$sigmaPhotoSp, main = paste("sigmaPhotoSp is " , signif(sigmaPhotoSp,3), sep = ""))
# abline(v = sigmaPhotoSp, col="red", lwd=3, lty=2)
# 
# # png("figures/simulatedPairs.png")
# #pairs(mdl.ssd, pars = c("muForceSp", "muChillSp", "muPhotoSp", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__")) 
# # dev.off()
# 
#   #Prior Predictive Check (Run 1000 times and plot results)
#   #------------------------------------------------------------
#   
#   #Number fo prior check itterations 
#   nRepPrior <- 300
#   
#   # ppc for traits portion
#   priorCheckTrait <- data.frame(matrix(NA, Ntrt*nRepPrior, 3))
#   names(priorCheckTrait) <- c("simRep", "rep", "species")
#   priorCheckTrait$simRep <- rep(1:nRepPrior, each = Ntrt)
#   priorCheckTrait$rep <- rep(1:Ntrt, times = nRepPrior)
#   priorCheckTrait$species <- rep(1:Nspp, each = nRep)
#   priorCheckTrait$transect <- rep(1:Ntransect, each = nRep)
#   
#   #traitSLA <- rnorm(Ntrt, 20, 5)
#   
#   #Make this the name of the full vector of sla per species values - alphaTraitSp 
#   #priorCheckTrait$alphaTraitSp <-  rep(rep(trt.dat$mu_grand_sp, times = nRepPrior))
#   stemDen <- trtPheno[complete.cases(trtPheno$ssd),]
#   
#   specieslist <- sort(unique(trtPheno$species))
#   sitelist <- sort(unique(trtPheno$transect))
#   leafMass <- trtPheno[complete.cases(trtPheno$ssd),]
#   
#   ssd.data <- list(yTraiti = stemDen$ssd,
#                    N = nrow(stemDen),
#                    n_spec = length(specieslist),
#                    trait_species = as.numeric(as.factor(stemDen$species)),
#                    n_site = length(sitelist),
#                    site = as.numeric(as.factor(stemDen$transect)),
#                    prior_mu_grand_mu = 1,
#                    prior_mu_grand_sigma = 5, #widened
#                    prior_sigma_sp_mu = 4, #10
#                    prior_sigma_sp_sigma = 5,
#                    prior_sigma_site_mu = 2, #5
#                    prior_sigma_site_sigma = 5, #2
#                    prior_sigma_traity_mu = 2, #5
#                    prior_sigma_traity_sigma = 5,
#                    ## Phenology
#                    Nph = nrow(pheno.t),
#                    phenology_species = as.numeric(as.factor(pheno.t$species)),
#                    yPhenoi = pheno.t$bb,
#                    forcei = pheno.t$force.z2,
#                    chilli = pheno.t$chillport.z2,
#                    photoi = pheno.t$photo.z2,
#                    prior_muForceSp_mu = -15,
#                    prior_muForceSp_sigma = 10, #wider
#                    prior_muChillSp_mu = -15,
#                    prior_muChillSp_sigma = 10,#wider
#                    prior_muPhotoSp_mu = -15,
#                    prior_muPhotoSp_sigma = 10,#wider
#                    prior_muPhenoSp_mu = 40,
#                    prior_muPhenoSp_sigma = 10,#wider
#                    prior_sigmaForceSp_mu = 5,
#                    prior_sigmaForceSp_sigma = 5,
#                    prior_sigmaChillSp_mu = 5,#wider
#                    prior_sigmaChillSp_sigma = 5, #wider
#                    prior_sigmaPhotoSp_mu = 5,
#                    prior_sigmaPhotoSp_sigma = 5,
#                    prior_sigmaPhenoSp_mu = 5, #wider
#                    prior_sigmaPhenoSp_sigma = 5, #wider
#                    prior_betaTraitxForce_mu = 0,
#                    prior_betaTraitxForce_sigma = 1,
#                    prior_betaTraitxChill_mu = 0,
#                    prior_betaTraitxChill_sigma = 1,
#                    prior_betaTraitxPhoto_mu = 0,
#                    prior_betaTraitxPhoto_sigma = 1,
#                    prior_sigmaphenoy_mu = 10,
#                    prior_sigmaphenoy_sigma = 5 #wider
#   )
#   
#   for (ir in 1:nRepPrior){
#     # Parameter Values
#     #ir <- 1
#     
#     muGrand <- rtruncnorm(1, a = 0, mean = ssd.data$prior_mu_grand_mu, sd = ssd.data$prior_mu_grand_sigma)
#     sigmaSp <- rtruncnorm(1, a = 0, mean = ssd.data$prior_sigma_sp_mu, sd = ssd.data$prior_sigma_sp_sigma)
#     sigmatransect <- rtruncnorm(1, a = 0, mean = ssd.data$prior_sigma_site_mu, sd = ssd.data$prior_sigma_site_sigma)
#     
#     alphaTraitSp <- rnorm(Nspp, 0, sigma.species)
#     priorCheckTrait$alphaTraitSp[priorCheckTrait$simRep == ir] <- rep(alphaTraitSp, each = nRep)
#     
#     muSp <- rnorm(Nspp, 0, sigma.species)
#     priorCheckTrait$muSp[priorCheckTrait$simRep == ir] <- rep(muSp, each = nRep)
#     
#     mutransect <- rnorm(Ntransect, 0, sigmatransect)
#     priorCheckTrait$mutransect[priorCheckTrait$simRep == ir] <- rep(mutransect, each = nRep)
#     
#     #general varience
#     priorCheckTrait$sigmaTrait_y[priorCheckTrait$simRep == ir] <- rnorm(ssd.data$prior_sigma_traity_mu, ssd.data$prior_sigma_traity_sigma)
#     priorCheckTrait$e[priorCheckTrait$simRep == ir] <- rnorm(Ntrt, 0, sigmaTrait_y)
#     
#     priorCheckTrait$yTraiti <- muGrand + priorCheckTrait$muSp + priorCheckTrait$mutransect + priorCheckTrait$e
#   }# end simulating new priors, from here vectorize code
#   
#   #Final values
#   priorCheckTrait$muGrandSp <- muGrand + priorCheckTrait$muSp
#   
#   
#   png("figures/density_Trait_Prior_joint_ssd.png")
#   plot(density(priorCheckTrait$yTraiti))
#   dev.off()
#   
#   png("figures/GrandSp_PlotPrior_joint_ssd.png")
#   plot(priorCheckTrait$yTraiti ~ priorCheckTrait$muGrandSp, xlab = "muGrandSp", ylab = "Trait")
#   dev.off()
#   
#   png("figures/MuSp_PlotPrior_joint_ssd.png")
#   plot(priorCheckTrait$yTraiti ~ priorCheckTrait$muSp, xlab = "MuSp", ylab = "Trait")
#   dev.off()
#   
#   png("figures/Mutransect_PlotPrior_joint_ssd.png")
#   plot(priorCheckTrait$yTraiti ~ priorCheckTrait$mutransect, xlab = "Mutransect", ylab = "Trait")
#   dev.off()
#   #####################################################################################
#   
#   #Make a data frame for input simulation data
#   priorCheckPheno <- data.frame(matrix(NA, Nph*nRepPrior, 3))
#   names(priorCheckPheno) <- c("simRep","rep","species")
#   
#   priorCheckPheno$simRep <- rep(1:nRepPrior, each = Nph)
#   
#   priorCheckPheno$rep <- rep(c(1:Nph), times = nRepPrior)
#   priorCheckPheno$species <- rep(rep(c(1:n_spec), each = nRep), times = nRepPrior)
#   
#   
#   
#   #Simulate SLA data per species
#   muGrandSp <- muGrand + muSp
#   #Make this the name of the full vector of sla per species values - alphaTraitSp 
#   priorCheckPheno$alphaTraitSp <-  rep(rep(muGrandSp, times = nRepPrior)) # use the mean mu_grand_sp, grand mean + transect, not transect
#   
#   
#   #Simulate cues (z scored)
#   priorCheckPheno$forcei <- rnorm(Nph, 1, 1)
#   priorCheckPheno$photoi <- rnorm(Nph, 1, 1) # less photoperiod 
#   priorCheckPheno$chilli <- rnorm(Nph, 1, 1) #more chilling
#   
#   head(priorCheckPheno)
#   
#   for (ir in 1:nRepPrior){
#     # Parameter Values
#     #ir <- 1
#     
#     #Species means
#     sigmaPhenoSp <- rtruncnorm(1, a = 0, mean = ssd.data$prior_sigmaPhenoSp_mu, sd = ssd.data$prior_sigmaPhenoSp_sigma)
#     muPhenoSp <- rnorm(1, ssd.data$prior_muPhenoSp_mu, ssd.data$prior_muPhenoSp_sigma)
#     alphaPhenoSp <- rnorm(n_spec, muPhenoSp, sigmaPhenoSp)
#     priorCheckPheno$alphaPhenoSp[priorCheckPheno$simRep == ir] <- rep(alphaPhenoSp, each = nRep)
#     
#     #Cue effects
#     priorCheckPheno$betaTraitxForce[priorCheckPheno$simRep == ir] <- rnorm(1,ssd.data$prior_betaTraitxForce_mu,ssd.data$prior_betaTraitxForce_sigma)
#     priorCheckPheno$betaTraitxPhoto[priorCheckPheno$simRep == ir] <- rnorm(1,ssd.data$prior_betaTraitxPhoto_mu,ssd.data$prior_betaTraitxPhoto_sigma)
#     priorCheckPheno$betaTraitxChill[priorCheckPheno$simRep == ir] <- rnorm(1,ssd.data$prior_betaTraitxChill_mu,ssd.data$prior_betaTraitxChill_sigma)
#     
#     #Species level slopes sans trait data
#     muForceSp <- rnorm(1,ssd.data$prior_muForceSp_mu,  ssd.data$prior_muForceSp_sigma)
#     sigmaForceSp <- rtruncnorm(1, a = 0, mean = ssd.data$prior_sigmaForceSp_mu,sd = ssd.data$prior_sigmaForceSp_sigma)
#     alphaForceSp <- rnorm(n_spec, muForceSp, sigmaForceSp)
#     priorCheckPheno$alphaForceSp[priorCheckPheno$simRep == ir] <- rep(alphaForceSp, each = nRep)
#     
#     muPhotoSp <- rnorm(1, ssd.data$prior_muPhotoSp_mu, ssd.data$prior_muPhotoSp_sigma)
#     sigmaPhotoSp <- rtruncnorm(1, a = 0, mean = ssd.data$prior_sigmaPhotoSp_mu, sd = ssd.data$prior_sigmaPhotoSp_sigma )
#     alphaPhotoSp <- rnorm(n_spec, muPhotoSp, sigmaPhotoSp)
#     priorCheckPheno$alphaPhotoSp[priorCheckPheno$simRep == ir] <- rep(alphaPhotoSp, each = nRep)
#     
#     muChillSp <-  rnorm(1,ssd.data$prior_sigmaChillSp_mu,ssd.data$prior_sigmaChillSp_sigma)
#     sigmaChillSp <- rtruncnorm(1, a = 0, mean = ssd.data$prior_sigmaChillSp_mu,sd = ssd.data$prior_sigmaChillSp_sigma)
#     alphaChillSp <- rnorm(n_spec, muChillSp, sigmaChillSp)
#     priorCheckPheno$alphaChillSp[priorCheckPheno$simRep == ir] <- rep(alphaChillSp, each = nRep)
#     
#     
#     #general varience
#     priorCheckPheno$sigmapheno_y[priorCheckPheno$simRep == ir] <- rtruncnorm(ssd.data$prior_sigma_traity_mu,  a = 0, ssd.data$prior_sigma_traity_sigma)
#     priorCheckPheno$e[priorCheckPheno$simRep == ir] <- rnorm(Nph, 0, sigmapheno_y)
#     
#   }# end simulating new priors, from here vectorize code
#   #slopes for each cue, combining trait and non-trait aspect of the slope.
#   
#   
#   priorCheckPheno$betaForceSp <- priorCheckPheno$alphaForceSp + priorCheckPheno$betaTraitxForce *  priorCheckPheno$alphaTraitSp
#   
#   priorCheckPheno$betaPhotoSp <-  priorCheckPheno$alphaPhotoSp + priorCheckPheno$betaTraitxPhoto * priorCheckPheno$alphaTraitSp
#   
#   priorCheckPheno$betaChillSp <-  priorCheckPheno$alphaChillSp + priorCheckPheno$betaTraitxChill * priorCheckPheno$alphaTraitSp
#   
#   #Run full model to get mean simulated y values
#   priorCheckPheno$yMu <-  priorCheckPheno$alphaPhenoSp +  priorCheckPheno$betaForceSp * priorCheckPheno$forcei +  priorCheckPheno$betaPhotoSp* priorCheckPheno$photoi + priorCheckPheno$betaChillSp * priorCheckPheno$chilli
#   
#   #Final values
#   priorCheckPheno$yPhenoi <- priorCheckPheno$yMu + priorCheckPheno$e
#   
#   head(priorCheckPheno)
#   plot(priorCheckPheno$betaForceSp ~ priorCheckPheno$alphaTraitSp )
#   priorCheckPheno_posF <- priorCheckPheno[priorCheckPheno$betaForceSp > 0,]
#   plot(priorCheckPheno_posF$betaForceSp ~ priorCheckPheno_posF$alphaTraitSp )
#   
#   png("figures/densityYPrior_joint.png")
#   plot(density(priorCheckPheno$yPhenoi))
#   dev.off()
#   
#   png("figures/photoPlotPrior_joint.png")
#   plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$photoi, xlab = "Photoperiod", ylab = "Phenological Date")
#   dev.off()
#   
#   png("figures/forcingPlotPrior_joint.png")
#   plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$forcei, xlab = "Forcing", ylab = "Phenological Date")
#   dev.off()
#   
#   png("figures/chillingPlotPrior_joint.png")
#   plot(priorCheckPheno$yPhenoi ~ priorCheckPheno$chilli, xlab = "Chillina", ylab = "Phenological Date")
#   dev.off()
#   
# 
# 
# if(BayesSweave == TRUE){
#   #For the BayesClass sweave documents 
#   setwd("/home/faith/Documents/github/bayes2020/Projects/Faith/traitorsModel")
# }
