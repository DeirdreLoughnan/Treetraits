# One more thesis chapter to go!

# Started Oct 10 2023 by D Loughnan
# aim of this code is to run the model with latitude as a continuous variable on the real trait data for my bc trait chapter

rm(list=ls())
options(stringsAsFactors = FALSE)

## Load libraries
library(rstan)
require(shinystan)
require(stringr)
library(plyr)
library(dplyr)


library(ggplot2)
library(reshape2)
# library(viridis)
#library(bayesplot)
#library(tidybayes)
library(gridExtra) 

options(mc.cores = 4)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/Treetraits") 
}  else{
  setwd("/home/deirdre/Treetraits") # for midge
}

trtPheno <- read.csv("input/trtPhenoZScore.csv", stringsAsFactors = FALSE)
pheno.t <- read.csv("input/bbPhenoZScore.csv", stringsAsFactors = FALSE)

specieslist <- sort(unique(trtPheno$species))
sitelist <- sort(unique(trtPheno$transect))

stem <- trtPheno[complete.cases(trtPheno$ssd),]

ssd.data <- list(yTraiti = stem$ssd, 
  N = nrow(stem),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(stem$species)),
  n_tran = length(unique(stem$transect)),
  lati = stem$latitude,
  tranE = as.numeric(stem$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

ssdZ.data <- list(yTraiti = stem$ssd, 
  N = nrow(stem),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(stem$species)),
  n_tran = length(unique(stem$transect)),
  lati = stem$latZ,
  tranE = as.numeric(stem$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

# mdl <- stan("stan/ssdDummyInt.stan",
#   data = ssd.data,
#   iter = 4000, warmup = 3000, chains=4,
#   include = FALSE, pars = c("y_hat")
# )

mdlSSD <- stan("stan/ssdDummyIntGrandNarrow.stan",
  data = ssd.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)
# Regarding the above model:
# A lot of the model below's priors for the traits still seem huge, what if we make the really tiny
# Also what if we make the prior for betaTraitxForce very small - does the prior also get small or does it stay the same ie is the prior defining that parameter? --- Yes that is very much the case

mdlSSD <- stan("stan/ssdDummyIntGrand.stan",
  data = ssd.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)
save(mdlSSD, file="output/ssdDummyIntGrand_np_contlat2.Rdata")


mdlSSD <- stan("stan/ssdDummyIntGrand.stan",
  data = ssdZ.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)
save(mdlSSD, file="output/ssdDummyIntGrand_zlat.Rdata")

# new priors for trait part of model:
postSSD<- rstan::extract(mdlSSD)

mu.grand <- 0.5 # the grand mean of the SSD model
sigma.species <- 0.5 # we want to keep the variaiton across spp. high

b.tranE <- 0
b.tranlat <- 0
sigma_traity <- 0.5

sigmaPhenoSp <- 5
muPhenoSp <- 40
betaTraitxForce <- 0
betaTraitxPhoto <- 0
betaTraitxChill <- 0
muForceSp <- -10
sigmaForceSp <- 5
muPhotoSp <- -5
sigmaPhotoSp <- 5
muChillSp <- -10
sigmaChillSp <- 5
sigmapheno_y <- 10

pdf("ssdGrandTrait.pdf")
par(mfrow = c(2,3))
hist(postSSD$mu_grand, main = "muGrand", xlim = c(-2,2), col=rgb(0,0,1,1/4))
hist(rnorm(1000, 0.5,0.5), col=rgb(1,0,1,1/4), add = T)
abline(v = mu.grand, col="red", lwd=3, lty=2)
# try 0.5,0.5
hist(postSSD$b_tranE, main = "b_tranE",  xlim = c(-5,5),col=rgb(0,0,1,1/4))
hist(rnorm(1000, 0,0.5), col=rgb(1,0,1,1/4), add = T)
abline(v = b.tranE, col="red", lwd=3, lty=2)

hist(postSSD$b_tranlat, main = "b_tranlat", col=rgb(0,0,1,1/4), xlim = c(-2,2))
hist(rnorm(1000, 0,5), col=rgb(1,0,1,1/4), add = T)
abline(v = b.tranlat, col="red", lwd=3, lty=2)

hist(postSSD$sigma_sp, main = "sigmaSp", col=rgb(0,0,1,1/4), xlim =c(-1,1))
hist(rnorm(1000, 0.5,0.5), col=rgb(1,0,1,1/4), add = T)
abline(v = sigma.species, col="red", lwd=3, lty=2)
# try 0.5,0.5

hist(postSSD$sigma_traity, main = "sigmaTraitY", col=rgb(0,0,1,1/4),  xlim = c(-2,2))
hist(rnorm(1000, 0.5,0.5), col=rgb(1,0,1,1/4), add = T)
abline(v = sigma_traity, col="red", lwd=3, lty=2)
# try 0.5,0.5
dev.off()

pdf("ssdGrandPheno.pdf")
par(mfrow = c(4,3))
hist(postSSD$muPhenoSp, main = "muPhenoSp", xlim = c(0,100), col=rgb(0,0,1,1/4))
hist(rnorm(1000, 40,10), col=rgb(1,0,1,1/4), add = T)
abline(v = muPhenoSp, col="red", lwd=3, lty=2)

hist(postSSD$muForceSp, main = "muForceSp", xlim = c(-75,75),col=rgb(0,0,1,1/4))
hist(rnorm(1000, -150,10), col=rgb(1,0,1,1/4), add = T)
abline(v = muForceSp, col="red", lwd=3, lty=2)

hist(postSSD$muChillSp, main = "muChillSp",  xlim = c(-75,75), col=rgb(0,0,1,1/4))
hist(rnorm(1000, -10,10), col=rgb(1,0,1,1/4), add = T)
abline(v = muChillSp, col="red", lwd=3, lty=2)

hist(postSSD$muPhotoSp, main = "muPhotoSp", col=rgb(0,0,1,1/4),  xlim = c(-75,75))
hist(rnorm(1000, -10,10), col=rgb(1,0,1,1/4), add = T)
abline(v = muPhotoSp, col="red", lwd=3, lty=2)

hist(postSSD$sigmapheno_y, main = "sigmaPhenoY", col=rgb(0,0,1,1/4),  xlim = c(-10,30))
hist(rnorm(1000, 10,5), col=rgb(1,0,1,1/4), add = T)
abline(v = sigmapheno_y, col="red", lwd=3, lty=2)

hist(postSSD$betaTraitxForce, main = "betaTraitxForce", col=rgb(0,0,1,1/4),  xlim = c(-5,5))
hist(rnorm(1000, 0,1), col=rgb(1,0,1,1/4), add = T)
abline(v = betaTraitxForce, col="red", lwd=3, lty=2)
#
hist(postSSD$betaTraitxChill, main = "betaTraitxChill", col=rgb(0,0,1,1/4),  xlim = c(-5,5))
hist(rnorm(1000, 0,1), col=rgb(1,0,1,1/4), add = T)
abline(v = betaTraitxChill, col="red", lwd=3, lty=2)
#
hist(postSSD$betaTraitxPhoto, main = "betaTraitxPhoto", col=rgb(0,0,1,1/4),  xlim = c(-5,5))
hist(rnorm(1000, 0,1), col=rgb(1,0,1,1/4), add = T)
abline(v = betaTraitxPhoto, col="red", lwd=3, lty=2)

hist(postSSD$sigmaChillSp, main = "sigmaChillSp", col=rgb(0,0,1,1/4), xlim = c(-10,20))
hist(rnorm(1000, 5,5), col=rgb(1,0,1,1/4), add = T)
abline(v = sigmaChillSp, col="red", lwd=3, lty=2)
# try 5,10
hist(postSSD$sigmaForceSp, main = "sigmaForceSp", col=rgb(0,0,1,1/4),  xlim = c(-10,10))
hist(rnorm(1000, 5,5), col=rgb(1,0,1,1/4), add = T)
abline(v = sigmaForceSp, col="red", lwd=3, lty=2)

hist(postSSD$sigmaPhotoSp, main = "sigmaPhotoSp", col=rgb(0,0,1,1/4),  xlim = c(-15,15))
hist(rnorm(1000, 5,5), col=rgb(1,0,1,1/4), add = T)
abline(v = sigmaPhotoSp, col="red", lwd=3, lty=2)
dev.off()
#pairs(mdl.SSD, pars = c("muForceSp", "muChillSp", "muPhotoSp", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
# dev.off()

sumer <- data.frame(summary(mdlSSD)$summary[c("mu_grand","b_tranE","b_tranlat", "muForceSp", "muChillSp", "muPhotoSp","muPhenoSp","betaTraitxForce", "betaTraitxChill","betaTraitxPhoto","sigma_traity" ,"sigma_sp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigmaPhenoSp","sigmapheno_y"),c("mean","2.5%","25%","50%", "75%","97.5%")])

sumer

hist(postSSD$betaTraitxForce)

# model summary output from z scored lat and new priors
# mean       X2.5.         X25.        X50.        X75.      X97.5.
# mu_grand         0.48410858  0.45131073  0.472986867  0.48390107  0.49509298  0.51801054
# b_tranE          0.03508660 -0.03822459  0.009840029  0.03562252  0.06039372  0.10825437
# b_tranlat       -0.03152138 -0.09057138 -0.051204501 -0.03118538 -0.01153042  0.02667818
# muForceSp       -4.89601335 -6.16400732 -5.322177988 -4.88880245 -4.45409780 -3.65876037
# muChillSp       -6.50497481 -8.42086997 -7.152231743 -6.50331783 -5.84827889 -4.56534914
# muPhotoSp       -1.68089901 -2.68195909 -2.060120881 -1.69217172 -1.30781107 -0.61869568
# muPhenoSp       27.05933771 24.10077560 26.118194783 27.07307300 28.02216107 29.90294945
# betaTraitxForce  0.24009186 -1.64313060 -0.446224572  0.23784542  0.88318001  2.24579671
# betaTraitxChill  0.11762155 -1.81669123 -0.569369238  0.12022128  0.81468523  2.00748446
# betaTraitxPhoto -0.08271429 -1.99491807 -0.777760798 -0.09063730  0.60187953  1.78358485
# sigma_traity     0.22126969  0.21294985  0.218246853  0.22120366  0.22430236  0.23016126
# sigma_sp         0.07426488  0.05528013  0.066741579  0.07348217  0.08106810  0.09722841
# sigmaForceSp     2.63364323  2.04659789  2.382982100  2.60546273  2.84539169  3.40564433
# sigmaChillSp     4.94666423  3.76014262  4.457071663  4.89436686  5.36788043  6.45476144
# sigmaPhotoSp     1.19926097  0.76788034  1.029104149  1.18464757  1.35387981  1.73249436
# sigmaPhenoSp     9.74656196  7.95425950  8.992181156  9.69896458 10.39372617 11.98836342
# sigmapheno_y     9.49989262  9.29795973  9.427498099  9.49835069  9.57084883  9.71004722

## without z scored lat but new priors:

# mean       X2.5.        X25.         X50.         X75.       X97.5.
# mu_grand         0.483576125  0.45068609  0.47198709  0.483603390  0.495134755  0.517079048
# b_tranE          0.405923995 -0.22780956  0.18387810  0.403827187  0.627007038  1.048199448
# b_tranlat       -0.007658589 -0.02217381 -0.01261965 -0.007637408 -0.002574076  0.006699118
# muForceSp       -4.903537510 -6.15643367 -5.34296224 -4.910011411 -4.477230788 -3.610107494
# muChillSp       -6.517674761 -8.38008272 -7.14927759 -6.520970704 -5.885056224 -4.620153264
# muPhotoSp       -1.702613053 -2.78691409 -2.04618219 -1.697200128 -1.343186694 -0.715734397
# muPhenoSp       27.023341073 24.16384303 26.07937445 27.002043548 27.996403512 29.901392327
# betaTraitxForce  0.241512711 -1.69340510 -0.40603783  0.242603210  0.928512498  2.144371272
# betaTraitxChill  0.104728344 -1.80397335 -0.55926910  0.082970409  0.757137982  2.049386850
# betaTraitxPhoto -0.057603843 -1.90991265 -0.70223002 -0.065487524  0.548241351  1.887014933
# sigma_traity     0.221283515  0.21311641  0.21844897  0.221218273  0.224127993  0.229515213
# sigma_sp         0.074280934  0.05635948  0.06682640  0.073461401  0.080716223  0.097173642
# sigmaForceSp     2.631100396  2.02538758  2.39125615  2.603735727  2.839575362  3.421499094
# sigmaChillSp     4.924515260  3.74607170  4.43283420  4.852105965  5.336432406  6.451621378
# sigmaPhotoSp     1.200705018  0.76414611  1.03553333  1.186115518  1.344534976  1.726014576
# sigmaPhenoSp     9.715331817  7.95841412  9.01184341  9.651556731 10.334582960 11.943436488
# sigmapheno_y     9.498085146  9.30290265  9.42702828  9.495614067  9.569685691  9.703852618