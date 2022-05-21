#Started May 6, 2022 by Deirdre

# Aim of this code is to perform several ppc for the trait pheno model for NAm tree traits
# This model includes transect as a dummy variable but otherwise is the same as the traitors model

rm(list=ls())
options(stringsAsFactors = FALSE)
## Set number of cores
options(mc.cores = 4)

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
  setwd("~/Documents/github/Treetraits") 
}  else{
  setwd("/home/deirdre/Treetraits") # for midge
}

## Posterior predictive check 1: does the model give us back values from simulated data:
# 1. Simulate some data:
# 2. Does the model give us good estimates of the data - pairs plots, mu plots, histograms with red bars for the given data values

Nrep <- 15 # rep per trait
Ntransect <- 2
Nspp <- 40 # number of species 

# First making a data frame for the test trait data
Ntrt <- Nspp * Ntransect * Nrep # total number of traits observations
Ntrt

mu.grand <- 20 # the grand mean of the height model
sigma.species <- 4 # we want to keep the variaiton across spp. high
sigma.transect <- 2
sigmaTrait_y <- 3

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
muPhenoSp <- 30
alphaPhenoSp <- rnorm(n_spec, muPhenoSp, sigmaPhenoSp)
pheno.dat$alphaPhenoSp <- rep(alphaPhenoSp, each = nRep)

#Cue effects
betaTraitxForce <- -0.5 
betaTraitxPhoto <- -0.2
betaTraitxChill <- -0.4

#Species level slopes sans trait data
muForceSp <- -15
sigmaForceSp <- 5
alphaForceSp <- rnorm(n_spec, muForceSp, sigmaForceSp)
pheno.dat$alphaForceSp <- rep(alphaForceSp, each = nRep)

muPhotoSp <- -15
sigmaPhotoSp <- 5
alphaPhotoSp <- rnorm(n_spec, muPhotoSp, sigmaPhotoSp)
pheno.dat$alphaPhotoSp <- rep(alphaPhotoSp, each = nRep)

muChillSp <- -15
sigmaChillSp <- 5
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
specieslist <- sort(unique(trt.dat$species))
sitelist <- sort(unique(trt.dat$transect))

ht.data <- list(yTraiti = trt.dat$yTraiti, 
                N = nrow(trt.dat),
                n_spec = length(specieslist),
                trait_species = as.numeric(as.factor(trt.dat$species)),
                n_site = length(sitelist),
                site = as.numeric(as.factor(trt.dat$transect)),
                prior_mu_grand_mu = 20,
                prior_mu_grand_sigma = 10,
                prior_sigma_sp_mu = 4,
                prior_sigma_sp_sigma = 5,
                prior_sigma_site_mu = 2,
                prior_sigma_site_sigma = 5,
                prior_sigma_traity_mu = 3,
                prior_sigma_traity_sigma = 5,
                ## Phenology
                Nph = nrow(pheno.dat),
                phenology_species = as.numeric(as.factor(pheno.dat$species)),
                yPhenoi = pheno.dat$yPhenoi,
                forcei = pheno.dat$forcei,
                chilli = pheno.dat$chilli,
                photoi = pheno.dat$photoi,
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

#ht.data$site
mdlSimHt <- stan("stan/jointMdl.stan",
               data = ht.data,
               iter = 6000,
               warmup = 3000,
               chains = 4)
# include = FALSE, pars = c("y_hat"))
save(mdlSimHt, file = "output/htSimPPC.Rda")


# plot(pheno.dat$betaForceSp ~ pheno.dat$alphaTraitSp)
# plot(pheno.dat$alphaForceSp ~ pheno.dat$alphaTraitSp)
# 
# plot(pheno.dat$betaChillSp ~ pheno.dat$alphaTraitSp)
# plot(pheno.dat$alphaChillSp ~ pheno.dat$alphaTraitSp)
# 
# plot(pheno.dat$betaPhotoSp ~ pheno.dat$alphaTraitSp)
# plot(pheno.dat$alphaPhotoSp ~ pheno.dat$alphaTraitSp)


#######################################################
postHt <- extract(mdlSimHt)
png("simPosteriorHist_ht.png")
par(mfrow=c(3,4))
#Compare results to simulated values
hist(postHt$muPhenoSp, main = paste("muPhenoSp is " , signif(muPhenoSp,3), sep = ""), xlim = c(0,100))
abline(v = muPhenoSp, col="red", lwd=3, lty=2)

hist(postHt$muForceSp, main = paste("muForceSp is " , signif(muForceSp,3), sep = ""))
abline(v = muForceSp, col="red", lwd=3, lty=2)

hist(postHt$muChillSp, main = paste("muChillSp is " , signif(muChillSp,3), sep = ""))
abline(v = muChillSp, col="red", lwd=3, lty=2)

hist(postHt$muPhotoSp, main = paste("muPhotoSp is " , signif(muPhotoSp,3), sep = ""))
abline(v = muPhotoSp, col="red", lwd=3, lty=2)

hist(postHt$sigmapheno_y, main = paste("sigmapheno_y is " , signif(sigmapheno_y,3), sep = ""))
abline(v = sigmapheno_y, col="red", lwd=3, lty=2)

hist(postHt$betaTraitxForce, main = paste("betaTraitxForce is " , signif(betaTraitxForce,3), sep = ""))
abline(v = betaTraitxForce, col="red", lwd=3, lty=2)
#
hist(postHt$betaTraitxChill, main = paste("betaTraitxChill is " , signif(betaTraitxChill,3), sep = ""))
abline(v = betaTraitxChill, col="red", lwd=3, lty=2)
#
hist(postHt$betaTraitxPhoto, main = paste("betaTraitxPhoto is " , signif(betaTraitxPhoto,3), sep = ""))
abline(v = betaTraitxPhoto, col="red", lwd=3, lty=2)

hist(postHt$sigmaChillSp, main = paste("sigmaChillSp is " , signif(sigmaChillSp,3), sep = ""), xlim = c(0,17))
abline(v = sigmaChillSp, col="red", lwd=3, lty=2)

hist(postHt$sigmaForceSp, main = paste("sigmaForceSp is " , signif(sigmaForceSp,3), sep = ""))
abline(v = sigmaForceSp, col="red", lwd=3, lty=2)

hist(postHt$sigmaPhotoSp, main = paste("sigmaPhotoSp is " , signif(sigmaPhotoSp,3), sep = ""), xlim = c(0,6))
abline(v = sigmaPhotoSp, col="red", lwd=3, lty=2)

dev.off()

png("figures/simulatedPairs.png")
pairs(mdlSimHt, pars = c("muForceSp", "muChillSp", "muPhotoSp", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto","sigmapheno_y", "lp__"))
dev.off()






####################################################################################

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

#########################################################

fit <- readRDS("output/height_stanfit.RDS")

png("figures/realPairs.png")
pairs(fit, pars = c("muForceSp", "muChillSp", "muPhotoSp", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto","sigmapheno_y", "lp__"))
dev.off()

post<- rstan::extract(fit)

postHt <- data.frame(post)

cueEffects <- postHt[, colnames(postHt) %in% c("muPhenoSp", "muForceSp", "muChillSp", "muPhotoSp", "sigmapheno_y")]

mcmc_intervals(cueEffects) + 
  theme_classic() + 
  labs(title = "main intercept, cue slopes and general error")

###############
postHt_alpaForceSp <- postHt[,colnames(postHt) %in% grep( "alphaForceSp", colnames(postHt), value = TRUE)]
colnames(postHt_alpaForceSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postHt_alpaForceSp) + 
  geom_vline(xintercept = mean(postHt$muForceSp), linetype="dotted", color = "grey")  +
  theme_classic() + 
  labs(subtitle = paste0("Mean muForceSp was ", round(mean(postHt$muForceSp),3)),
       title = "muForceSp - species forcing slopes no trait")

postHt_betaForceSp <- postHt[,colnames(postHt) %in% grep( "betaForceSp", colnames(postHt), value = TRUE)]
colnames(postHt_betaForceSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postHt_betaForceSp) + 
  theme_classic() + 
  labs(title = "betaForceSp - Species forcing slopes with trait value")

#Different species slopes for chilling, without the effect of trait
postHt_alphaChillSp <- postHt[,colnames(postHt) %in% grep( "alphaChillSp", colnames(postHt), value = TRUE)]
colnames(postHt_alphaChillSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postHt_alphaChillSp) + 
  geom_vline(xintercept = mean(postHt$muChillSp), linetype="dotted", color = "grey")  +
  theme_classic() + 
  labs(subtitle = paste0("Mean muChillSp was ", round(mean(postHt$muChillSp),3)),
       title = "alphaChillSp - Species chill slopes no trait")

#Different species slopes for forcing, with the effect of trait
postHt_betaChillSp <- postHt[,colnames(postHt) %in% grep( "betaChillSp", colnames(postHt), value = TRUE)]
colnames(postHt_betaChillSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postHt_betaChillSp) + 
  theme_classic() + 
  labs(title = "betaChillSp - Species chilling slopes with trait value")

#Different species slopes for photoperiod, without the effect of trait
postHt_alphaPhotoSp <- postHt[,colnames(postHt) %in% grep( "alphaPhotoSp", colnames(postHt), value = TRUE)]
colnames(postHt_alphaPhotoSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postHt_alphaPhotoSp) + 
  geom_vline(xintercept = mean(postHt$muPhotoSp), linetype="dotted", color = "grey")  +
  theme_classic() + 
  labs(subtitle = paste0("Mean muPhotoSp was ", round(mean(postHt$muPhotoSp),3)),
       title = "muPhotoSp - Species photo period slopes no trait")

#Different species slopes for forcing, with the effect of trait
postHt_betaPhotoSp <- postHt[,colnames(postHt) %in% grep( "betaPhotoSp", colnames(postHt), value = TRUE)]
colnames(postHt_betaPhotoSp) <- levels(as.factor(trtPheno$species))

mcmc_intervals(postHt_betaPhotoSp) + 
  theme_classic() + 
  labs(title = "betaPhotoSp - Species photoperiod slopes with trait value")


#Different species slopes for forcing only the effect of trait
postHt_betaTraitx <- postHt[,colnames(postHt) %in% grep( "betaTraitx", colnames(postHt), value = TRUE)]

mcmc_intervals(postHt_betaTraitx) + 
  theme_classic() + 
  labs(title = "effect's of traits on cue slopes")



##################################################################################
 # Plot the y vs the ypred
mdlOut <- load("output/ht_raw.Rda")

sumer <- summary(mdl.ht)$summary
post <- rstan::extract(mdl.ht)

require(bayesplot)
htData <- trtPheno[complete.cases(trtPheno$ht),]
y <- htData$ht
yrep <-  post$y_hat

ppc_dens_overlay(y, yrep[1:50,])

pheno.term$species <- tolower(pheno.term$species)

htTrtSplvl <- aggregate( htData["ht"],
                         htData[c("species")],
                         FUN = mean)

htPhenoSplvl <- aggregate( pheno.term["bb"],
                           pheno.term[c("species")],
                           FUN = mean)

muSp <- sumer[grep("muSp\\[", rownames(sumer))]
alphaForceSp <- sumer[grep("alphaForceSp\\[", rownames(sumer))]
alphaChillSp <- sumer[grep("alphaChillSp\\[", rownames(sumer))]
alphaPhotoSp <- sumer[grep("alphaPhotoSp\\[", rownames(sumer))]
alphaPhenoSp <- sumer[grep("alphaPhenoSp\\[", rownames(sumer))]

muGrandSp <- sumer[grep("mu_grand_sp\\[", rownames(sumer))]
betaForceSp <- sumer[grep("betaForceSp\\[", rownames(sumer))]
betaChillSp <- sumer[grep("betaChillSp\\[", rownames(sumer))]
betaPhotoSp <- sumer[grep("betaPhotoSp\\[", rownames(sumer))]

plot(muSp ~ muGrandSp)
plot(alphaForceSp ~ betaForceSp)
plot(alphaChillSp ~ betaChillSp)
plot(alphaPhotoSp ~ betaPhotoSp)

plot(alphaForceSp ~ alphaChillSp)
plot(alphaChillSp ~ alphaPhotoSp)
plot(alphaPhotoSp ~ alphaForceSp)

plot(betaForceSp ~ betaChillSp)
plot(betaChillSp ~ betaPhotoSp)
plot(betaPhotoSp ~ betaForceSp)

plot(alphaPhenoSp ~ alphaForceSp)
plot(alphaPhenoSp ~ alphaChillSp)
plot(alphaPhenoSp ~ alphaPhotoSp)

plot(htTrtSplvl$ht ~ muGrandSp)
plot(htTrtSplvl$ht ~ muSp)

plot(htSplvl$bb ~ betaForceSp)
plot(htSplvl$bb ~ betaChillSp)
plot(htSplvl$bb ~ betaPhotoSp)

plot(htSplvl$bb ~ alphaPhenoSp)
plot(htSplvl$bb ~ muGrandSp)

