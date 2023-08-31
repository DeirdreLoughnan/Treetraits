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


# remove the trait spp we don't have pheno data for:
phenoSp <- sort(unique(pheno.t$species))

trtPheno <- trtPheno[trtPheno$species %in% phenoSp, ]

specieslist <- sort(unique(trtPheno$species))
sitelist <- sort(unique(trtPheno$site.n))

# Now subset into the two transects and

# starting with height
# 1. Eastern species
height <- trtPheno[complete.cases(trtPheno$ht),]

exPop <- c("HF", "SH")
pheno.W <- pheno.t[!pheno.t$population %in% exPop, ]

excl <- c("kl", "af", "mp", "sm")
htE <- height[height$site %in% excl, ]

htE.data <- list(yTraiti = htE$ht,
                 N = nrow(htE),
                 n_spec = length(unique(htE$species)),
                 trait_species = as.numeric(as.factor(htE$species)),
                 n_pop = length(unique(htE$site.n)),
                 pop2 = htE$site2,
                 pop3 = htE$site3,
                 pop4 = htE$site4,
                 ## Phenology
                 Nph = nrow(pheno.W),
                 phenology_species = as.numeric(as.factor(pheno.W$species)),
                 yPhenoi = pheno.W$bb,
                 forcei = pheno.W$force.z2,
                 chilli = pheno.W$chillport.z2,
                 photoi = pheno.W$photo.z2
)

mdl.htE <- stan("stan/heightPhenoDummy.stan",
                data = htE.data,
                iter = 4000,
                warmup = 3000,
                chains = 4, include = FALSE, pars = c("y_hat"))
                #control = list(adapt_delta =0.99))

# first run - no controls - 313 divergent transitions after warmup
save(mdl.htE, file = "output/htEastSept2023.Rda")
sumerhtE <- summary(mdl.htE)$summary

pdf("pairsHtE.pdf")
pairs(mdl.htE, pars = c("mu_grand", "muPhenoSp","muForceSp", "muChillSp", "muPhotoSp","b_pop2","b_pop3","b_pop4","sigma_sp","sigmaPhenoSp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigma_traity", "sigmapheno_y", "betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto", "lp__"))
dev.off()

range(summary(mdl.htE)$summary[, "n_eff"]) # 394.271, 13592.123
range(summary(mdl.htE)$summary[, "Rhat"])

