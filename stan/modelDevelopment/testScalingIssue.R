# Aim of this code is to test whether the scaling is actually an issue:
# If we make height values very very small do we run into the same issues? Is height actually also constrained by the prior?

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
# pheno$force.n[pheno$force.n == "HF"] <- "1"
# pheno$force.n[pheno$force.n == "LF"] <- "0"
# pheno$force.n <- as.numeric(pheno$force.n)

pheno$force.n[pheno$force.n == "HF" & pheno$population == "mp"] <- "15"
pheno$force.n[pheno$force.n == "HF" & pheno$population == "sm"] <- "15"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "mp"] <- "10"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "sm"] <- "10"

# (8*20+16*10)/24 #13.33
# (8*15+16*5)/24 #8.33

pheno$force.n[pheno$force.n == "HF" & pheno$population == "HF"] <- "13.33"
pheno$force.n[pheno$force.n == "HF" & pheno$population == "SH"] <- "13.33"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "HF"] <- "8.33"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "SH"] <- "8.33"
pheno$force.n <- as.numeric(pheno$force.n)
# pheno$photo.n <- pheno$photo
# pheno$photo.n[pheno$photo.n == "HP"] <- "1"
# pheno$photo.n[pheno$photo.n == "LP"] <- "0"
# pheno$photo.n <- as.numeric(pheno$photo.n)

pheno$photo.n <- pheno$photo
pheno$photo.n[pheno$photo.n == "HP"] <- "12"
pheno$photo.n[pheno$photo.n == "LP"] <- "8"
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

# head(pheno)
# #add dummy/ site level effects:
# pheno <- pheno %>%
#   mutate ( site2 = if_else(site.n == 2, 1, 0),
#            site3 = if_else(site.n == 3, 1, 0),
#            site4 = if_else(site.n == 4, 1, 0))

pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)

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
trtPheno$latZ <- (trtPheno$latitude-mean(trtPheno$latitude,na.rm=TRUE))/(sd(trtPheno$latitude,na.rm=TRUE)*2)

specieslist <- sort(unique(trtPheno$species))
sitelist <- sort(unique(trtPheno$transect))
######################################################################
# 1. Height 
height <- trtPheno[complete.cases(trtPheno$ht),]
height$ht.z2 <- (height$ht-mean(height$ht,na.rm=TRUE))/(sd(height$ht,na.rm=TRUE)*2)

ht.data <- list(yTraiti = height$ht, 
  N = nrow(height),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(height$species)),
  n_tran = length(unique(height$transect)),
  lati = height$latZ,
  tranE = as.numeric(height$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

htZ.data <- list(yTraiti = height$ht.z2, 
  N = nrow(height),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(height$species)),
  n_tran = length(unique(height$transect.z2)),
  lati = height$lat.z2,
  tranE = as.numeric(height$transect.z2),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

mdlHt <- stan("stan/heightDummyIntGrandWide.stan",
  data = htZ.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)

save(mdl, file="output/heightDummyInt.Rdata")
save(mdlHt, file="output/heightDummyIntGrandZ25.Rdata")
# ssm <- as.shinystan(mdl)
# launch_shinystan(ssm)

sumer <- data.frame(summary(mdlHt)$summary[c(
  "b_tranE","b_tranlat", "muForceSp", "muChillSp", "muPhotoSp","muPhenoSp","betaTraitxForce", "betaTraitxChill","betaTraitxPhoto","sigma_traity" ,"sigma_sp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigmaPhenoSp","sigmapheno_y"),c("mean","2.5%","25%","50%", "75%","97.5%")])
sumer

postHt <- rstan::extract(mdlHt)
par(mfrow = c(1,3))
hist(postHt$betaTraitxForce, main = "betaTraitxForce",  xlim = c(-100, 100))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postHt$betaTraitxChill, main = "betaTraitxChill",  xlim = c(-100, 100))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postHt$betaTraitxPhoto, main = "betaTraitxPhoto",  xlim = c(-100, 100))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

##### What happens if the height values were really really tiny:

ht.tiny <- list(yTraiti = (height$ht/100), 
  N = nrow(height),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(height$species)),
  n_tran = length(unique(height$transect)),
  lati = height$latZ,
  tranE = as.numeric(height$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

mdlHtTiny <- stan("stan/heightDummyIntGrandWide.stan",
  data = ht.tiny,
  iter = 3000, warmup = 2000, chains=4,
  include = FALSE, pars = c("y_hat")
)

sumerT <- data.frame(summary(mdlHtTiny)$summary[c(
  "b_tranE","b_tranlat", "muForceSp", "muChillSp", "muPhotoSp","muPhenoSp","betaTraitxForce", "betaTraitxChill","betaTraitxPhoto","sigma_traity" ,"sigma_sp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigmaPhenoSp","sigmapheno_y"),c("mean","2.5%","25%","50%", "75%","97.5%")])
sumerT

postHt <- rstan::extract(mdlHtTiny)

par(mfrow = c(1,3))
hist(postHt$betaTraitxForce, main = "betaTraitxForce",  xlim = c(-10, 10))
hist(rnorm(1000, 0,1), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postHt$betaTraitxChill, main = "betaTraitxChill",  xlim = c(-10, 10))
hist(rnorm(1000, 0,1), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postHt$betaTraitxPhoto, main = "betaTraitxPhoto",  xlim = c(-10, 10))
hist(rnorm(1000, 0,1), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)



## For context:
load("output/heightDummyIntGrandZ.Rdata")
