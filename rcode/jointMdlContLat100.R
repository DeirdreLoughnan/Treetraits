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
pheno$force.z <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE))
pheno$photo.z <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE))
pheno$chillport.z <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE))

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
trtPheno$transect <- as.numeric(trtPheno$transect)
trtPheno$transect.z2 <- (trtPheno$transect-mean(trtPheno$transect,na.rm=TRUE))/(sd(trtPheno$transect,na.rm=TRUE)*2) # Gelman --- 0/1 so use 2sd

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
trtPheno$lat.z <- (trtPheno$latitude-mean(trtPheno$latitude,na.rm=TRUE))/(sd(trtPheno$latitude,na.rm=TRUE))
trtPheno$lat.z2 <- (trtPheno$latitude-mean(trtPheno$latitude,na.rm=TRUE))/(sd(trtPheno$latitude,na.rm=TRUE)*2)
trtPheno$lat.100 <- (trtPheno$latitude/100)

specieslist <- sort(unique(trtPheno$species))
sitelist <- sort(unique(trtPheno$transect))

# write.csv(trtPheno ,"input/trtPhenoZScore.csv", row.names = F)
# write.csv(pheno.t,"input/bbPhenoZScore.csv", row.names = F)
height <- trtPheno[complete.cases(trtPheno$ht),]

ht.data <- list(yTraiti = height$ht, 
  N = nrow(height),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(height$species)),
  n_tran = length(unique(height$transect)),
  lati = height$latitude,
  tranE = as.numeric(height$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z,
  chilli = pheno.t$chillport.z,
  photoi = pheno.t$photo.z
)

mdlHt <- stan("stan/heightDummyIntGrand.stan",
  data = ht.data,
  iter = 6000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)

save(mdlHt, file="output/htContLatHundoLat.Rdata")

######################################################################
# 2. Leaf mass area
leafMA <- trtPheno[complete.cases(trtPheno$lma),]
leafMA$lma.z <- (leafMA$lma-mean(leafMA$lma,na.rm=TRUE))/(sd(leafMA$lma,na.rm=TRUE))*2

## put it on the same scale as height
leafMA$lma100 <- leafMA$lma*100
leafMA$lmalog <- log10(leafMA$lma)

lma <- list(yTraiti = leafMA$lma, 
  N = nrow(leafMA),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(leafMA$species)),
  n_tran = length(unique(leafMA$transect)),
  lati = leafMA$lat.z,
  tranE = as.numeric(leafMA$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z,
  chilli = pheno.t$chillport.z,
  photoi = pheno.t$photo.z
)

# mdlLMA <- stan("stan/lmaDummyIntGrand.stan",
#   data = lma,
#   iter = 4000, warmup = 3000, chains=4,
#   include = FALSE, pars = c("y_hat")
# )
# 
# save(mdlLMA, file="output/lmaContLatNone.Rdata")

lma.100 <- list(yTraiti = leafMA$lma100, 
  N = nrow(leafMA),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(leafMA$species)),
  n_tran = length(unique(leafMA$transect)),
  lati = leafMA$latitude,
  tranE = as.numeric(leafMA$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z,
  chilli = pheno.t$chillport.z,
  photoi = pheno.t$photo.z
)

lma.10 <- list(yTraiti = (leafMA$lma)*10, 
  N = nrow(leafMA),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(leafMA$species)),
  n_tran = length(unique(leafMA$transect)),
  lati = leafMA$lat.z,
  tranE = as.numeric(leafMA$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z,
  chilli = pheno.t$chillport.z,
  photoi = pheno.t$photo.z
)
## What about log 10? 
lma.log <- list(yTraiti = leafMA$lmalog, 
  N = nrow(leafMA),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(leafMA$species)),
  n_tran = length(unique(leafMA$transect.z2)),
  #lati = leafMA$latitude,
  lati = leafMA$lat.z2,
  tranE = as.numeric(leafMA$transect.z2),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

mdlLMA <- stan("stan/lmaDummyIntGrand.stan",
  data = lma.100,
  include = FALSE, pars = c("y_hat"),
  iter = 4000, chains= 4, warmup = 3000)
  #control = list(max_treedepth =11)) # may run better with 8000


save(mdlLMA, file="output/lmaContLatHundoLat.Rdata")

sumer100 <- data.frame(summary(mdlLMA)$summary[c("mu_grand",
  "b_tranE","b_tranlat", "muForceSp", "muChillSp", "muPhotoSp","muPhenoSp","betaTraitxForce", "betaTraitxChill","betaTraitxPhoto","sigma_traity" ,"sigma_sp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigmaPhenoSp","sigmapheno_y"),c("mean","2.5%","25%","50%", "75%","97.5%")])
sumer100

postLMA<- data.frame(rstan::extract(mdlLMA))

par(mfrow = c(1,1))
hist(postLMA$betaTraitxForce, main = "betaTraitxForce",  xlim = c(-50, 50))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postLMA$betaTraitxChill, main = "betaTraitxChill",  xlim = c(-200, 200))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postLMA$betaTraitxPhoto, main = "betaTraitxPhoto",  xlim = c(-200, 200))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)


sumer <- data.frame(summary(mdlLMAHundo)$summary[c("mu_grand","b_tranE","b_tranlat", "muForceSp", "muChillSp", "muPhotoSp","muPhenoSp","betaTraitxForce", "betaTraitxChill","betaTraitxPhoto","sigma_traity" ,"sigma_sp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigmaPhenoSp","sigmapheno_y"),c("mean","2.5%","25%","50%", "75%","97.5%")])
sumer
bMuSp <- summary(mdl)$summary[grep("b_muSp\\["),c("mean","2.5%","25%","50%", "75%","97.5%")]


####################################################
# 3. diameter at breast height
diam <- trtPheno[complete.cases(trtPheno$dbh),] 
diam$dbh.z2 <- (diam$dbh-mean(diam$dbh,na.rm=TRUE))/(sd(diam$dbh,na.rm=TRUE)*2)

dbh.data <- list(yTraiti = diam$dbh, 
  N = nrow(diam),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(diam$species)),
  n_tran = length(unique(diam$transect)),
  lati = diam$latitude,
  tranE = as.numeric(diam$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z,
  chilli = pheno.t$chillport.z,
  photoi = pheno.t$photo.z
)

mdlDBH <- stan("stan/diamDummyIntGrand.stan",
  data = dbh.data,
  iter = 6000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)

save(mdlDBH, file="output/dbhContLatHundoLat.Rdata")

sumer <- data.frame(summary(mdlDBH)$summary[c(
  "b_tranE","b_tranlat", "muForceSp", "muChillSp", "muPhotoSp","muPhenoSp","betaTraitxForce", "betaTraitxChill","betaTraitxPhoto","sigma_traity" ,"sigma_sp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigmaPhenoSp","sigmapheno_y"),c("mean","2.5%","25%","50%", "75%","97.5%")])
sumer

postDBH<- data.frame(rstan::extract(mdlDBH))

par(mfrow = c(1,3))
hist(postDBH$betaTraitxForce, main = "betaTraitxForce",  xlim = c(-200, 200))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postDBH$betaTraitxChill, main = "betaTraitxChill",  xlim = c(-200, 200))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postDBH$betaTraitxPhoto, main = "betaTraitxPhoto",  xlim = c(-200, 200))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

######################################################################
# 4. stem specific density
stem <- trtPheno[complete.cases(trtPheno$ssd),]

stem$ssd.z2 <- (stem$ssd-mean(stem$ssd,na.rm=TRUE))/(sd(stem$ssd,na.rm=TRUE)*2)
stem$ssd100 <- stem$ssd*100
stem$ssdlog <- log10(stem$ssd)

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
  forcei = pheno.t$force.z,
  chilli = pheno.t$chillport.z,
  photoi = pheno.t$photo.z
)


mdlSSD <- stan("stan/ssdDummyIntGrand.stan",
  data = ssd.data,
  iter = 6000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat"),
  control = list(max_treedepth =12)
)

save(mdlSSD, file="output/ssdContLatHundoLat.Rdata")

########## 
# ssd.data <- list(yTraiti = stem$ssdlog, 
#   N = nrow(stem),
#   n_spec = length(specieslist),
#   trait_species = as.numeric(as.factor(stem$species)),
#   n_tran = length(unique(stem$transect)),
#   lati = stem$lat.z,
#   tranE = as.numeric(stem$transect),
#   Nph = nrow(pheno.t),
#   phenology_species = as.numeric(as.factor(pheno.t$species)),
#   yPhenoi = pheno.t$bb,
#   forcei = pheno.t$force.z,
#   chilli = pheno.t$chillport.z,
#   photoi = pheno.t$photo.z
# )
# 
# 
# mdlSSD <- stan("stan/ssdDummyIntGrand.stan",
#   data = ssd.data,
#   iter = 4000, warmup = 3000, chains=4,
#   include = FALSE, pars = c("y_hat")
# )

# save(mdlSSD, file="output/ssdContLatLog.Rdata")

sumer <- data.frame(summary(mdlSSD)$summary[c(
  "b_tranE","b_tranlat", "muForceSp", "muChillSp", "muPhotoSp","muPhenoSp","betaTraitxForce", "betaTraitxChill","betaTraitxPhoto","sigma_traity" ,"sigma_sp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigmaPhenoSp","sigmapheno_y"),c("mean","2.5%","25%","50%", "75%","97.5%")])
sumer

postSSD<- data.frame(rstan::extract(mdlSSD))

par(mfrow = c(1,1))
hist(postSSD$betaTraitxForce, main = "betaTraitxForce",  xlim = c(-200, 200))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postSSD$betaTraitxChill, main = "betaTraitxChill",  xlim = c(-200, 200))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postSSD$betaTraitxPhoto, main = "betaTraitxPhoto",  xlim = c(-200, 200))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)
# looks the exact same as if you z-score the lma, so let's not since I am not sure it is a good idea

######################################################################
# 5. carbon to nitrogen ratio
carbNit <- trtPheno[complete.cases(trtPheno$C.N),]
carbNit$cn.z2 <- (carbNit$C.N-mean(carbNit$C.N,na.rm=TRUE))/(sd(carbNit$C.N,na.rm=TRUE)*2)
carbNit$cn.port <- (carbNit$C.N/100)
carbNit$cn.port <- (carbNit$C.N/100)

# cn.data <- list(yTraiti = carbNit$C.N, 
#   N = nrow(carbNit),
#   n_spec = length(specieslist),
#   trait_species = as.numeric(as.factor(carbNit$species)),
#   n_tran = length(unique(carbNit$transect)),
#   lati = carbNit$lat.z,
#   tranE = as.numeric(carbNit$transect),
#   Nph = nrow(pheno.t),
#   phenology_species = as.numeric(as.factor(pheno.t$species)),
#   yPhenoi = pheno.t$bb,
#   forcei = pheno.t$force.z,
#   chilli = pheno.t$chillport.z,
#   photoi = pheno.t$photo.z
# )

cn.data <- list(yTraiti = carbNit$cn.port, 
  N = nrow(carbNit),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(carbNit$species)),
  n_tran = length(unique(carbNit$transect)),
  lati = carbNit$lat.,
  tranE = as.numeric(carbNit$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z,
  chilli = pheno.t$chillport.z,
  photoi = pheno.t$photo.z
)

mdlCN <- stan("stan/lmaDummyIntGrand.stan",
  data = cn.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
  #,control = list(adapt_delta = 0.99, max_treedepth =12)
)

# mdlCN <- stan("stan/cnDummyIntGrand.stan",
#   data = cnZ.data,
#   iter = 4000, warmup = 3000, chains=4,
#   include = FALSE, pars = c("y_hat")
#   #,control = list(adapt_delta = 0.99, max_treedepth =12)
# )
sumer <- data.frame(summary(mdlCN)$summary[c(
  "b_tranE","b_tranlat", "muForceSp", "muChillSp", "muPhotoSp","muPhenoSp","betaTraitxForce", "betaTraitxChill","betaTraitxPhoto","sigma_traity" ,"sigma_sp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigmaPhenoSp","sigmapheno_y"),c("mean","2.5%","25%","50%", "75%","97.5%")])
sumer

postSSD<- data.frame(rstan::extract(mdlSSD))

par(mfrow = c(1,3))
hist(postSSD$betaTraitxForce, main = "betaTraitxForce",  xlim = c(-200, 200))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postSSD$betaTraitxChill, main = "betaTraitxChill",  xlim = c(-200, 200))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postSSD$betaTraitxPhoto, main = "betaTraitxPhoto",  xlim = c(-200, 200))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)
# looks the exact same as if you z-score the lma, so let's not since I am not sure it is a good idea


save(mdlCN, file="output/cnContLat.Rdata")

#####
######################################################################
# 5. leaf nitrogen content
nit <- trtPheno[complete.cases(trtPheno$per.N),]
#carbNit$n.z2 <- (carbNit$-mean(carbNit$C.N,na.rm=TRUE))/(sd(carbNit$C.N,na.rm=TRUE)*2)
# carbNit$n.port <- (carbNit$per.N/100)
# carbNit$nLog <- log10(carbNit$per.N)

nit.data <- list(yTraiti = nit$per.N, 
  N = nrow(nit),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(nit$species)),
  n_tran = length(unique(nit$transect)),
  lati = nit$latitude,
  tranE = as.numeric(nit$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z,
  chilli = pheno.t$chillport.z,
  photoi = pheno.t$photo.z
)

mdlPerN <- stan("stan/cnDummyIntGrand.stan",
  data = nit.data,
  iter = 6000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
  #,control = list(adapt_delta = 0.99, max_treedepth =12)
)

save(mdlPerN, file="output/lncContLatHundoLat.Rdata")

# load("output/nitDummyInt.Rdata")
#  ssm <- as.shinystan(mdl)
#  launch_shinystan(ssm)

sumer <- summary(mdlPerN)$summary

ModelFit <- rstan::extract(mdlPerN)

muSp <- data.frame(ModelFit$b_muSp)
muSpMean <- colMeans(muSp)

betaForceSp <- data.frame(ModelFit$betaForceSp)
betaForceSpMean <- colMeans(betaForceSp)

quantile2575 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.25, 0.75))
  return(returnQuanilte)
}

bf_quan <- apply(betaForceSp, 2, quantile2575) 
mu_quan <- apply(muSp, 2, quantile2575)

bfs <- rbind(betaForceSpMean, bf_quan)
bfs_t <- t(bfs)
bfs_df <- data.frame(bfs_t)
colnames(bfs_df)[colnames(bfs_df) == "X25."] <- "force25"
colnames(bfs_df)[colnames(bfs_df) == "X75."] <- "force75"
bfs_df$species <- specieslist

mg<- rbind(muSpMean, mu_quan)
mg_t <- t(mg)
mg_df <- data.frame(mg_t)
colnames(mg_df)[colnames(mg_df) == "X25."] <- "trait25"
colnames(mg_df)[colnames(mg_df) == "X75."] <- "trait75"
mg_df$species <- specieslist

muForceSp <- data.frame(ModelFit$muForceSp)
muForceSpMean <- colMeans(muForceSp)

betaTraitxForce <- data.frame(ModelFit$betaTraitxForce)
betaTraitxForceMean <- colMeans(betaTraitxForce)

# mg_df_east <- mg_df[mg_df$species %in% eastSp, ]
# mg_df_west <- mg_df[mg_df$species %in% westSp, ]
# 
# bfs_df_east <- bfs_df[bfs_df$species %in% eastSp, ]
# bfs_df_west <- bfs_df[bfs_df$species %in% westSp, ]

plot( x= mg_df$muSpMean, y = bfs_df$betaForceSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bfs_df$force25), max(bfs_df$force75)), ylab = "Species level forcing slope", xlab = "Trait value", cex.lab = 1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df[,"muSpMean"], # x mean
  bfs_df[,"force25"], # y 25
  mg_df[,"muSpMean"],
  bfs_df[,"force75"],
  length = 0, col= "#218380", lwd = 2
)

arrows(
  mg_df[,"trait25"], # x mean
  bfs_df[,"betaForceSpMean"], # y 25
  mg_df[,"trait75"], # x mean
  bfs_df[,"betaForceSpMean"],
  length = 0, col = "#218380", lwd = 2
)


mtext(side = 3, text = "C:N, Forcing", adj = 0, cex = 1.25)
for(j in 1:length(muForceSp[,1])){
  abline(a = muForceSp[j,], b = betaTraitxForceMean, col=alpha("#73d2de", 0.085))
}
abline(a=muForceSpMean, b=betaTraitxForceMean, col = "black")

