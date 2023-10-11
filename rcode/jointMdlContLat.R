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

specieslist <- sort(unique(trtPheno$species))
sitelist <- sort(unique(trtPheno$transect))
######################################################################
# 1. Height 
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
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

length(ht.data$tranE)
head(ht.data$tranE)


mdl <- stan("stan/heightDummyInt.stan",
  data = ht.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)

save(mdl, file="output/heightDummyInt.Rdata")

ssm <- as.shinystan(mdl)
launch_shinystan(ssm)

sumer <- summary(mdl)$summary

######################################################################
# 2. Leaf mass area
leafMA <- trtPheno[complete.cases(trtPheno$lma),]

lma.data <- list(yTraiti = leafMA$lma, 
  N = nrow(leafMA),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(leafMA$species)),
  n_tran = length(unique(leafMA$transect)),
  lati = leafMA$latitude,
  tranE = as.numeric(leafMA$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

mdl <- stan("stan/lmaDummyInt.stan",
  data = lma.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)

save(mdl, file="output/lmaDummyInt.Rdata")

######################################################################
# 3. diameter at breast height
diam <- trtPheno[complete.cases(trtPheno$dbh),]

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
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

mdl <- stan("stan/diamDummyInt.stan",
  data = dbh.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)

save(mdl, file="output/dbhDummyInt.Rdata")

######################################################################
# 4. stem specific density
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

mdl <- stan("stan/ssdDummyInt.stan",
  data = ssd.data,
  iter = 4000, warmup = 3000, chains=4,
  include = FALSE, pars = c("y_hat")
)

save(mdl, file="output/ssdDummyInt.Rdata")

######################################################################
# 5. carbon to nitrogen ratio
carbNit <- trtPheno[complete.cases(trtPheno$C.N),]

cn.data <- list(yTraiti = carbNit$C.N, 
  N = nrow(carbNit),
  n_spec = length(specieslist),
  trait_species = as.numeric(as.factor(carbNit$species)),
  n_tran = length(unique(carbNit$transect)),
  lati = carbNit$latitude,
  tranE = as.numeric(carbNit$transect),
  Nph = nrow(pheno.t),
  phenology_species = as.numeric(as.factor(pheno.t$species)),
  yPhenoi = pheno.t$bb,
  forcei = pheno.t$force.z2,
  chilli = pheno.t$chillport.z2,
  photoi = pheno.t$photo.z2
)

mdl <- stan("stan/cnDummyInt.stan",
  data = cn.data,
  iter = 6000, warmup = 5000, chains=4,
  include = FALSE, pars = c("y_hat")
)

save(mdl, file="output/cnDummyInt.Rdata")
