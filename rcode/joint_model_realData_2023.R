# started by DL Aug 30, 2023
# after reviewing the test data again, I am more confident in this joint model working. Let's see if it works on the real data!

# in this model:
# we are modeling the east and west individually
# site will be included as a dummy var in the trait model
# no site in pheno model
# pheno mdl similar to the phenology paper: cont var for f/c/p - 2z score BUT no interactions or phylogeny

rm(list=ls())
options(stringsAsFactors = FALSE)

## Load libraries
library(rstan)
require(stringr)
library(plyr)
library(ggplot2)
library(reshape2)
#library(viridis)
#library(forcats)
#library(ggdist)
#library(tidybayes)
library(gridExtra)
#require(shinystan)
library(dplyr)
# for arranging plots 
#library(patchwork) # another way of arranging plots 
#library(rethinking)

## Set number of cores
options(mc.cores = 4)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/Treetraits") 
}  else{
  setwd("/home/deirdre/Treetraits") # for midge
}

# Start by getting the pheno data
# dl <- read.csv("input/dl_allbb_mini.csv")
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
pheno <- read.csv("input/phenoDataWChill.csv")

pheno$force.n <- pheno$force
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

pheno$photo.n <- pheno$photo
pheno$photo.n[pheno$photo.n == "HP"] <- "12"
pheno$photo.n[pheno$photo.n == "LP"] <- "8"
pheno$photo.n <- as.numeric(pheno$photo.n)

pheno$site.n <- pheno$population
pheno$site.n[pheno$site.n == "sm"] <- "1"
pheno$site.n[pheno$site.n == "mp"] <- "2"
pheno$site.n[pheno$site.n == "HF"] <- "3"
pheno$site.n[pheno$site.n == "SH"] <- "4"
pheno$site.n <- as.numeric(pheno$site.n)

# west = 0, east = 1
pheno$transect <- pheno$population
pheno$transect[pheno$transect == "sm"] <- "0"
pheno$transect[pheno$transect == "mp"] <- "0"
pheno$transect[pheno$transect == "GR"] <- "1"
pheno$transect[pheno$transect == "HF"] <- "1"

# head(pheno)
# #add dummy/ site level effects:

# standardize the 0/1 and standardize sites? 
pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)

# pheno$site2.z2 <- (pheno$site2-mean(pheno$site2,na.rm=TRUE))/(sd(pheno$site2,na.rm=TRUE)*2)
# pheno$site3.z2 <- (pheno$site3-mean(pheno$site3,na.rm=TRUE))/(sd(pheno$site3,na.rm=TRUE)*2)
# pheno$site4.z2 <- (pheno$site4-mean(pheno$site4,na.rm=TRUE))/(sd(pheno$site4,na.rm=TRUE)*2)


pheno.term <- pheno[,c("bb", "force.z2", "photo.z2", "population", "species", "lab2","Utah_Model","Chill_portions","chillport.z2","transect", "site.n")] #"site2.z2", "site3.z2","site4.z2")]
pheno.t <- pheno.term[complete.cases(pheno.term$bb), ] # none,great!
length(unique(pheno.t$species))
#pheno.t <- pheno.term[complete.cases(pheno.term$bb), ] # 1780 rows data 
pheno.t$species <- tolower(pheno.t$species)
pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 47 bc two species occur in both transects

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Now get the trait data and subset to only include spp we have pheno data for:

trtData <- read.csv("input/allTrt.csv", stringsAsFactors = FALSE)


phenoSp <- sort(unique(pheno.t$species))

trtPheno <- trtData[trtData$species %in% phenoSp, ]
length(unique(trtPheno$species))


trtPheno$site.n <- trtPheno$site
trtPheno$site.n[trtPheno$site.n == "sm"] <- "1"
trtPheno$site.n[trtPheno$site.n == "mp"] <- "2"
trtPheno$site.n[trtPheno$site.n == "kl"] <- "3"
trtPheno$site.n[trtPheno$site.n == "af"] <- "4"
trtPheno$site.n[trtPheno$site.n == "GR"] <- "5"
trtPheno$site.n[trtPheno$site.n == "HF"] <- "6"
trtPheno$site.n[trtPheno$site.n == "SH"] <- "7"
trtPheno$site.n[trtPheno$site.n == "WM"] <- "8"

trtPheno$transect <- trtPheno$site
trtPheno$transect[trtPheno$transect == "sm"] <- "0"
trtPheno$transect[trtPheno$transect == "mp"] <- "0"
trtPheno$transect[trtPheno$transect == "kl"] <- "0"
trtPheno$transect[trtPheno$transect == "af"] <- "0"
trtPheno$transect[trtPheno$transect == "GR"] <- "1"
trtPheno$transect[trtPheno$transect == "HF"] <- "1"
trtPheno$transect[trtPheno$transect == "SH"] <- "1"
trtPheno$transect[trtPheno$transect == "WM"] <- "1"

trtPheno <- trtPheno %>%
  mutate ( site2 = if_else(site.n == 2, 1, 0),
           site3 = if_else(site.n == 3, 1, 0),
           site4 = if_else(site.n == 4, 1, 0),
           site5 = if_else(site.n == 5, 1, 0),
           site6 = if_else(site.n == 6, 1, 0),
           site7 = if_else(site.n == 7, 1, 0),
           site8 = if_else(site.n == 8, 1, 0))

#trtPheno <- read.csv("input/trtPhenoDummy.csv")
# remove the trait spp we don't have pheno data for:
phenoSp <- sort(unique(pheno.t$species))

trtPheno <- trtPheno[trtPheno$species %in% phenoSp, ]

specieslist <- sort(unique(trtPheno$species))
sitelist <- sort(unique(trtPheno$site.n))

# Now subset into the two transects and

# starting with height
# 1. Eastern species
# height <- trtPheno[complete.cases(trtPheno$ht),]
# 
# exPop <- c("sm", "mp")
# pheno.E <- pheno.t[!pheno.t$population %in% exPop, ]
# 
# excl <- c("kl", "af", "mp", "sm")
# htE <- height[!height$site %in% excl, ]
# 
# unique(pheno.E$population)
# unique(htE$site)
# 
# htE.data <- list(yTraiti = htE$ht,
#                  N = nrow(htE),
#                  n_spec = length(unique(htE$species)),
#                  trait_species = as.numeric(as.factor(htE$species)),
#                  n_pop = length(unique(htE$site.n)),
#                  pop2 = htE$site2,
#                  pop3 = htE$site3,
#                  pop4 = htE$site4,
#                  ## Phenology
#                  Nph = nrow(pheno.E),
#                  phenology_species = as.numeric(as.factor(pheno.E$species)),
#                  yPhenoi = pheno.E$bb,
#                  forcei = pheno.E$force.z2,
#                  chilli = pheno.E$chillport.z2,
#                  photoi = pheno.E$photo.z2
# )
# 
# mdl.htE <- stan("stan/heightPhenoDummy.stan",
#                 data = htE.data,
#                 iter = 4000,
#                 warmup = 3000,
#                 chains = 4, include = FALSE, pars = c("y_hat"))
#                 #control = list(adapt_delta =0.99))
# 
# # first run - no controls - no divergent transitions after warmup
# save(mdl.htE, file = "output/htEastSept2023.Rda")
# 
# load("output/htEastSept2023.Rda")
# sumerhtE <- summary(mdl.htE)$summary
# 
# pdf("pairsHtE.pdf")
# pairs(mdl.htE, pars = c("mu_grand", "muPhenoSp","muForceSp", "muChillSp", "muPhotoSp","b_pop2","b_pop3","b_pop4","sigma_sp","sigmaPhenoSp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigma_traity", "sigmapheno_y", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
# dev.off()
# 
# range(summary(mdl.htE)$summary[, "n_eff"]) # 354.734 9061.874
# range(summary(mdl.htE)$summary[, "Rhat"]) # 0.999, 1.0217
# 
# # 2. Western species
# 
# exPop <- c("SH", "HF")
# pheno.W <- pheno.t[!pheno.t$population %in% exPop, ]
# 
# excl <- c("HF","SH", "WM","GR")
# htW <- height[!height$site %in% excl, ]
# 
# unique(pheno.W$population)
# unique(htW$site)
# 
# htW.data <- list(yTraiti = htW$ht,
#                  N = nrow(htW),
#                  n_spec = length(unique(htW$species)),
#                  trait_species = as.numeric(as.factor(htW$species)),
#                  n_pop = length(unique(htW$site.n)),
#                  pop2 = htW$site2,
#                  pop3 = htW$site3,
#                  pop4 = htW$site4,
#                  ## Phenology
#                  Nph = nrow(pheno.W),
#                  phenology_species = as.numeric(as.factor(pheno.W$species)),
#                  yPhenoi = pheno.W$bb,
#                  forcei = pheno.W$force.z2,
#                  chilli = pheno.W$chillport.z2,
#                  photoi = pheno.W$photo.z2
# )
# 
# # mdl.htW <- stan("stan/heightPhenoDummy.stan",
# #                 data = htW.data,
# #                 iter = 4000,
# #                 warmup = 3000,
# #                 chains = 4, include = FALSE, pars = c("y_hat"))
# 
# mdl.htWncp <- stan("stan/heightPhenoDummyNcpPhoto.stan",
#                 data = htW.data,
#                 iter = 4000,
#                 warmup = 3000,
#                 chains = 4, include = FALSE, pars = c("y_hat"))
# 
# #control = list(adapt_delta =0.99))
# 
# # first run - no controls - ncp of photo fixed divergencies
# save(mdl.htWncp, file = "output/htWestSept2023NcpPhoto.Rda")
# 
# #load("output/htWestSept2023.Rda")
# sumerhtW <- summary(mdl.htW)$summary
# 
# pdf("figures/pairshtW_MuPara.pdf")
# pairs(mdl.htW, pars = c("mu_grand", "muPhenoSp","muForceSp", "muChillSp", "muPhotoSp","b_pop2","b_pop3","b_pop4", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
# dev.off()
# 
# pdf("figures/pairshtW_sigmaPara_ncp.pdf")
# pairs(mdl.htWncp, pars = c("sigma_sp","sigmaPhenoSp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigma_traity", "sigmapheno_y", "lp__"))
# dev.off()
# 
# range(summary(mdl.htW)$summary[, "n_eff"]) # 120.9332 6779.7469
# range(summary(mdl.htW)$summary[, "Rhat"]) # 0.9991601 1.0226981

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# 2. Leaf mass area:
# Eastern species
lma <- trtPheno[complete.cases(trtPheno$lma),]

exPop <- c("sm", "mp")
pheno.E <- pheno.t[!pheno.t$population %in% exPop, ]

excl <- c("kl", "af", "mp", "sm")
lmaE <- lma[!lma$site %in% excl, ]

unique(pheno.E$population)
unique(lmaE$site)

lmaE.data <- list(yTraiti = lmaE$lma,
                 N = nrow(lmaE),
                 n_spec = length(unique(lmaE$species)),
                 trait_species = as.numeric(as.factor(lmaE$species)),
                 n_pop = length(unique(lmaE$site.n)),
                 pop2 = lmaE$site2,
                 pop3 = lmaE$site3,
                 pop4 = lmaE$site4,
                 ## Phenology
                 Nph = nrow(pheno.E),
                 phenology_species = as.numeric(as.factor(pheno.E$species)),
                 yPhenoi = pheno.E$bb,
                 forcei = pheno.E$force.z2,
                 chilli = pheno.E$chillport.z2,
                 photoi = pheno.E$photo.z2
)

# mdl.lmaE <- stan("stan/lmaPhenoDummy.stan",
#                 data = lmaE.data,
#                 iter = 4000,
#                 warmup = 3000,
#                 chains = 4, include = FALSE, pars = c("y_hat"))
# #control = list(adapt_delta =0.99))
# 
# # first run - no controls - no divergent transitions after warmup
# save(mdl.lmaE, file = "output/lmaEastSept2023.Rda")
# 
# load("output/lmaEastSept2023.Rda")
# sumerlmaE <- summary(mdl.lmaE)$summary
# 
# pdf("pairslmaE.pdf")
# pairs(mdl.lmaE, pars = c("mu_grand", "muPhenoSp","muForceSp", "muChillSp", "muPhotoSp","b_pop2","b_pop3","b_pop4","sigma_sp","sigmaPhenoSp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigma_traity", "sigmapheno_y", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
# dev.off()
# 
# range(summary(mdl.lmaE)$summary[, "n_eff"]) # 354.734 9061.874
# range(summary(mdl.lmaE)$summary[, "Rhat"]) # 0.999, 1.0217

# 2. Western species

exPop <- c("SH", "HF")
pheno.W <- pheno.t[!pheno.t$population %in% exPop, ]

excl <- c("HF","SH", "WM","GR")
lmaW <- lma[!lma$site %in% excl, ]

unique(pheno.W$population)
unique(lmaW$site)

lmaW.data <- list(yTraiti = lmaW$lma,
                 N = nrow(lmaW),
                 n_spec = length(unique(lmaW$species)),
                 trait_species = as.numeric(as.factor(lmaW$species)),
                 n_pop = length(unique(lmaW$site.n)),
                 pop2 = lmaW$site2,
                 pop3 = lmaW$site3,
                 pop4 = lmaW$site4,
                 ## Phenology
                 Nph = nrow(pheno.W),
                 phenology_species = as.numeric(as.factor(pheno.W$species)),
                 yPhenoi = pheno.W$bb,
                 forcei = pheno.W$force.z2,
                 chilli = pheno.W$chillport.z2,
                 photoi = pheno.W$photo.z2
)

# mdl.lmaW <- stan("stan/lmaPhenoDummy.stan",
#                 data = lmaW.data,
#                 iter = 4000,
#                 warmup = 3000,
#                 chains = 4, include = FALSE, pars = c("y_hat"))
# 184 div trans, rhat = 1.06

mdl.lmaWncp <- stan("stan/lmaPhenoDummyNcpPhoto.stan",
                   data = lmaW.data,
                   iter = 4000,
                   warmup = 3000,
                   chains = 4, include = FALSE, pars = c("y_hat"))

#control = list(adapt_delta =0.99))

# first run - no controls - ncp of photo fixed divergencies
save(mdl.lmaWncp, file = "output/lmaWestSept2023NcpPhoto.Rda")

#load("output/lmaWestSept2023.Rda")
sumerlmaW <- summary(mdl.lmaW)$summary

# pdf("figures/pairslmaW_MuPara.pdf")
# pairs(mdl.lmaW, pars = c("mu_grand", "muPhenoSp","muForceSp", "muChillSp", "muPhotoSp","b_pop2","b_pop3","b_pop4", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
# dev.off()
# 
# pdf("figures/pairslmaW_sigmaPara_ncp.pdf")
# pairs(mdl.lmaWncp, pars = c("sigma_sp","sigmaPhenoSp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigma_traity", "sigmapheno_y", "lp__"))
# dev.off()

range(summary(mdl.lmaWncp)$summary[, "n_eff"]) # 732.0782 9192.9273
range(summary(mdl.lmaWncp)$summary[, "Rhat"]) #  0.9990993 1.0050612


#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# 3. Stem specific density:
# Eastern species
ssd <- trtPheno[complete.cases(trtPheno$ssd),]

exPop <- c("sm", "mp")
pheno.E <- pheno.t[!pheno.t$population %in% exPop, ]

excl <- c("kl", "af", "mp", "sm")
ssdE <- ssd[!ssd$site %in% excl, ]

unique(pheno.E$population)
unique(ssdE$site)

ssdE.data <- list(yTraiti = ssdE$ssd,
                  N = nrow(ssdE),
                  n_spec = length(unique(ssdE$species)),
                  trait_species = as.numeric(as.factor(ssdE$species)),
                  n_pop = length(unique(ssdE$site.n)),
                  pop2 = ssdE$site2,
                  pop3 = ssdE$site3,
                  pop4 = ssdE$site4,
                  ## Phenology
                  Nph = nrow(pheno.E),
                  phenology_species = as.numeric(as.factor(pheno.E$species)),
                  yPhenoi = pheno.E$bb,
                  forcei = pheno.E$force.z2,
                  chilli = pheno.E$chillport.z2,
                  photoi = pheno.E$photo.z2
)

# mdl.ssdE <- stan("stan/ssdPhenoDummy.stan",
#                  data = ssdE.data,
#                  iter = 4000,
#                  warmup = 3000,
#                  chains = 4, include = FALSE, pars = c("y_hat"))
# #control = list(adapt_delta =0.99))
# 
# # first run - no controls - no divergent transitions after warmup
# save(mdl.ssdE, file = "output/ssdEastSept2023.Rda")
# 
# load("output/ssdEastSept2023.Rda")
# sumerssdE <- summary(mdl.ssdE)$summary
# 
# pdf("pairsssdE.pdf")
# pairs(mdl.ssdE, pars = c("mu_grand", "muPhenoSp","muForceSp", "muChillSp", "muPhotoSp","b_pop2","b_pop3","b_pop4","sigma_sp","sigmaPhenoSp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigma_traity", "sigmapheno_y", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
# dev.off()
# 
# range(summary(mdl.ssdE)$summary[, "n_eff"]) # 354.734 9061.874
# range(summary(mdl.ssdE)$summary[, "Rhat"]) # 0.999, 1.0217

# 2. Western species

exPop <- c("SH", "HF")
pheno.W <- pheno.t[!pheno.t$population %in% exPop, ]

excl <- c("HF","SH", "WM","GR")
ssdW <- ssd[!ssd$site %in% excl, ]

unique(pheno.W$population)
unique(ssdW$site)

ssdW.data <- list(yTraiti = ssdW$ssd,
                  N = nrow(ssdW),
                  n_spec = length(unique(ssdW$species)),
                  trait_species = as.numeric(as.factor(ssdW$species)),
                  n_pop = length(unique(ssdW$site.n)),
                  pop2 = ssdW$site2,
                  pop3 = ssdW$site3,
                  pop4 = ssdW$site4,
                  ## Phenology
                  Nph = nrow(pheno.W),
                  phenology_species = as.numeric(as.factor(pheno.W$species)),
                  yPhenoi = pheno.W$bb,
                  forcei = pheno.W$force.z2,
                  chilli = pheno.W$chillport.z2,
                  photoi = pheno.W$photo.z2
)

# mdl.ssdW <- stan("stan/ssdPhenoDummy.stan",
#                  data = ssdW.data,
#                  iter = 4000,
#                  warmup = 3000,
#                  chains = 4, include = FALSE, pars = c("y_hat"))
# 155 div transitions, rhat = 1.07
mdl.ssdWncp <- stan("stan/ssdPhenoDummyNcpPhoto.stan",
                    data = ssdW.data,
                    iter = 4000,
                    warmup = 3000,
                    chains = 4, include = FALSE, pars = c("y_hat"))

#control = list(adapt_delta =0.99))

# first run - no controls - ncp of photo fixed divergencies
save(mdl.ssdWncp, file = "output/ssdWestSept2023NcpPhoto.Rda")

#load("output/ssdWestSept2023.Rda")
sumerssdW <- summary(mdl.ssdW)$summary

# pdf("figures/pairsssdW_MuPara.pdf")
# pairs(mdl.ssdW, pars = c("mu_grand", "muPhenoSp","muForceSp", "muChillSp", "muPhotoSp","b_pop2","b_pop3","b_pop4", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
# dev.off()
# 
# pdf("figures/pairsssdW_sigmaPara_ncp.pdf")
# pairs(mdl.ssdWncp, pars = c("sigma_sp","sigmaPhenoSp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigma_traity", "sigmapheno_y", "lp__"))
# dev.off()

range(summary(mdl.ssdWncp)$summary[, "n_eff"]) # 1071.936 10084.375
range(summary(mdl.ssdWncp)$summary[, "Rhat"]) #  0.9990951 1.0056746

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# 4. diameter at breast height
# Eastern species
dbh <- trtPheno[complete.cases(trtPheno$dbh),]

exPop <- c("sm", "mp")
pheno.E <- pheno.t[!pheno.t$population %in% exPop, ]

excl <- c("kl", "af", "mp", "sm")
dbhE <- dbh[!dbh$site %in% excl, ]

unique(pheno.E$population)
unique(dbhE$site)

dbhE.data <- list(yTraiti = dbhE$dbh,
                  N = nrow(dbhE),
                  n_spec = length(unique(dbhE$species)),
                  trait_species = as.numeric(as.factor(dbhE$species)),
                  n_pop = length(unique(dbhE$site.n)),
                  pop2 = dbhE$site2,
                  pop3 = dbhE$site3,
                  pop4 = dbhE$site4,
                  ## Phenology
                  Nph = nrow(pheno.E),
                  phenology_species = as.numeric(as.factor(pheno.E$species)),
                  yPhenoi = pheno.E$bb,
                  forcei = pheno.E$force.z2,
                  chilli = pheno.E$chillport.z2,
                  photoi = pheno.E$photo.z2
)

# mdl.dbhE <- stan("stan/dbhPhenoDummy.stan",
#                  data = dbhE.data,
#                  iter = 4000,
#                  warmup = 3000,
#                  chains = 4, include = FALSE, pars = c("y_hat"))
# #control = list(adapt_delta =0.99))
# 
# # first run - no controls - no divergent transitions after warmup
# save(mdl.dbhE, file = "output/dbhEastSept2023.Rda")
# 
# load("output/dbhEastSept2023.Rda")
# sumerdbhE <- summary(mdl.dbhE)$summary
# 
# pdf("pairsdbhE.pdf")
# pairs(mdl.dbhE, pars = c("mu_grand", "muPhenoSp","muForceSp", "muChillSp", "muPhotoSp","b_pop2","b_pop3","b_pop4","sigma_sp","sigmaPhenoSp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigma_traity", "sigmapheno_y", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
# dev.off()
# 
# range(summary(mdl.dbhE)$summary[, "n_eff"]) # 354.734 9061.874
# range(summary(mdl.dbhE)$summary[, "Rhat"]) # 0.999, 1.0217

# 2. Western species

exPop <- c("SH", "HF")
pheno.W <- pheno.t[!pheno.t$population %in% exPop, ]

excl <- c("HF","SH", "WM","GR")
dbhW <- dbh[!dbh$site %in% excl, ]

unique(pheno.W$population)
unique(dbhW$site)

dbhW.data <- list(yTraiti = dbhW$dbh,
                  N = nrow(dbhW),
                  n_spec = length(unique(dbhW$species)),
                  trait_species = as.numeric(as.factor(dbhW$species)),
                  n_pop = length(unique(dbhW$site.n)),
                  pop2 = dbhW$site2,
                  pop3 = dbhW$site3,
                  pop4 = dbhW$site4,
                  ## Phenology
                  Nph = nrow(pheno.W),
                  phenology_species = as.numeric(as.factor(pheno.W$species)),
                  yPhenoi = pheno.W$bb,
                  forcei = pheno.W$force.z2,
                  chilli = pheno.W$chillport.z2,
                  photoi = pheno.W$photo.z2
)

# mdl.dbhW <- stan("stan/dbhPhenoDummy.stan",
#                  data = dbhW.data,
#                  iter = 4000,
#                  warmup = 3000,
#                  chains = 4, include = FALSE, pars = c("y_hat"))

mdl.dbhWncp <- stan("stan/dbhPhenoDummyNcpPhoto.stan",
                    data = dbhW.data,
                    iter = 4000,
                    warmup = 3000,
                    chains = 4, include = FALSE, pars = c("y_hat"))

#control = list(adapt_delta =0.99))

# first run - no controls - ncp of photo fixed divergencies
save(mdl.dbhWncp, file = "output/dbhWestSept2023NcpPhoto.Rda")

#load("output/dbhWestSept2023.Rda")
sumerdbhW <- summary(mdl.dbhW)$summary

# pdf("figures/pairsdbhW_MuPara.pdf")
# pairs(mdl.dbhW, pars = c("mu_grand", "muPhenoSp","muForceSp", "muChillSp", "muPhotoSp","b_pop2","b_pop3","b_pop4", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
# dev.off()
# 
# pdf("figures/pairsdbhW_sigmaPara_ncp.pdf")
# pairs(mdl.dbhWncp, pars = c("sigma_sp","sigmaPhenoSp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigma_traity", "sigmapheno_y", "lp__"))
# dev.off()

range(summary(mdl.dbhWncp)$summary[, "n_eff"]) # 229.106 8988.451
range(summary(mdl.dbhWncp)$summary[, "Rhat"]) #0.9991407 1.0141438

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# 5. carbon:nitrogen:
# Eastern species
cn <- trtPheno[complete.cases(trtPheno$C.N),]

exPop <- c("sm", "mp")
pheno.E <- pheno.t[!pheno.t$population %in% exPop, ]

excl <- c("kl", "af", "mp", "sm")
cnE <- cn[!cn$site %in% excl, ]

unique(pheno.E$population)
unique(cnE$site)

cnE.data <- list(yTraiti = cnE$C.N,
                  N = nrow(cnE),
                  n_spec = length(unique(cnE$species)),
                  trait_species = as.numeric(as.factor(cnE$species)),
                  n_pop = length(unique(cnE$site.n)),
                  pop2 = cnE$site2,
                  pop3 = cnE$site3,
                  pop4 = cnE$site4,
                  ## Phenology
                  Nph = nrow(pheno.E),
                  phenology_species = as.numeric(as.factor(pheno.E$species)),
                  yPhenoi = pheno.E$bb,
                  forcei = pheno.E$force.z2,
                  chilli = pheno.E$chillport.z2,
                  photoi = pheno.E$photo.z2
)

# mdl.cnE <- stan("stan/cnPhenoDummy.stan",
#                  data = cnE.data,
#                  iter = 4000,
#                  warmup = 3000,
#                  chains = 4, include = FALSE, pars = c("y_hat"))
# #control = list(adapt_delta =0.99))
# 
# # first run - no controls - no divergent transitions after warmup
# save(mdl.cnE, file = "output/cnEastSept2023.Rda")
# 
# load("output/cnEastSept2023.Rda")
# sumercnE <- summary(mdl.cnE)$summary
# 
# pdf("pairscnE.pdf")
# pairs(mdl.cnE, pars = c("mu_grand", "muPhenoSp","muForceSp", "muChillSp", "muPhotoSp","b_pop2","b_pop3","b_pop4","sigma_sp","sigmaPhenoSp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigma_traity", "sigmapheno_y", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
# dev.off()

# range(summary(mdl.cnE)$summary[, "n_eff"]) # 354.734 9061.874
# range(summary(mdl.cnE)$summary[, "Rhat"]) # 0.999, 1.0217

# 2. Western species

exPop <- c("SH", "HF")
pheno.W <- pheno.t[!pheno.t$population %in% exPop, ]

excl <- c("HF","SH", "WM","GR")
cnW <- cn[!cn$site %in% excl, ]

unique(pheno.W$population)
unique(cnW$site)

cnW.data <- list(yTraiti = cnW$C.N,
                  N = nrow(cnW),
                  n_spec = length(unique(cnW$species)),
                  trait_species = as.numeric(as.factor(cnW$species)),
                  n_pop = length(unique(cnW$site.n)),
                  pop2 = cnW$site2,
                  pop3 = cnW$site3,
                  pop4 = cnW$site4,
                  ## Phenology
                  Nph = nrow(pheno.W),
                  phenology_species = as.numeric(as.factor(pheno.W$species)),
                  yPhenoi = pheno.W$bb,
                  forcei = pheno.W$force.z2,
                  chilli = pheno.W$chillport.z2,
                  photoi = pheno.W$photo.z2
)

# mdl.cnW <- stan("stan/cnPhenoDummy.stan",
#                  data = cnW.data,
#                  iter = 4000,
#                  warmup = 3000,
#                  chains = 4, include = FALSE, pars = c("y_hat"))

mdl.cnWncp <- stan("stan/cnPhenoDummyNcpPhoto.stan",
                    data = cnW.data,
                    iter = 4000,
                    warmup = 3000,
                    chains = 4, include = FALSE, pars = c("y_hat"))

#control = list(adapt_delta =0.99))

# first run - no controls - ncp of photo fixed divergencies
save(mdl.cnWncp, file = "output/cnWestSept2023NcpPhoto.Rda")

#load("output/cnWestSept2023.Rda")
sumercnW <- summary(mdl.cnW)$summary

# pdf("figures/pairscnW_MuPara.pdf")
# pairs(mdl.cnW, pars = c("mu_grand", "muPhenoSp","muForceSp", "muChillSp", "muPhotoSp","b_pop2","b_pop3","b_pop4", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
# dev.off()
# 
# pdf("figures/pairscnW_sigmaPara_ncp.pdf")
# pairs(mdl.cnWncp, pars = c("sigma_sp","sigmaPhenoSp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigma_traity", "sigmapheno_y", "lp__"))
# dev.off()

range(summary(mdl.cnW)$summary[, "n_eff"]) # 120.9332 6779.7469
range(summary(mdl.cnW)$summary[, "Rhat"]) # 0.9991601 1.0226981
