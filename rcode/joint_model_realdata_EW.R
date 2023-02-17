# Started April 22, 2022 by Deirdre

# The purpose of this code is to get the traitors joint model working for my trait phenology data! 

rm(list=ls())
options(stringsAsFactors = FALSE)

## Load libraries
library(rstan)
require(stringr)
library(plyr)
library(ggplot2)
library(reshape2)
library(viridis)
#library(bayesplot)
#library(tidybayes)
library(gridExtra)
#require(shinystan)
#library(dplyr)
 # for arranging plots 
#library(patchwork) # another way of arranging plots 
#library(rethinking)

## Set number of cores
options(mc.cores = 4)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  else{
  setwd("/home/deirdre/pheno_bc") # for midge
}

# Start by getting the pheno data
dl <- read.csv("input/dl_allbb_mini.csv")

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

pheno$site.n <- pheno$population
pheno$site.n[pheno$site.n == "sm"] <- "1"
pheno$site.n[pheno$site.n == "mp"] <- "2"
pheno$site.n[pheno$site.n == "HF"] <- "3"
pheno$site.n[pheno$site.n == "SH"] <- "4"
pheno$site.n <- as.numeric(pheno$site.n)

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
pheno.term <- pheno[,c("bb", "force.z2", "photo.z2", "population", "species", "lab2","Utah_Model","Chill_portions","chillport.z2","transect", "site.n")] #"site2.z2", "site3.z2","site4.z2")]
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


##########################################################

# Sites all, not just transects
specieslist <- sort(unique(trtPheno$species))
sitelist <- sort(unique(trtPheno$site.n))
height <- trtPheno[complete.cases(trtPheno$ht),]



exPop <- c("HF", "SH")
pheno.tW <- pheno.t[!pheno.t$population %in% exPop, ]

# head(pheno)
# #add dummy/ site level effects:
# 
# ht3.data <- list(yTraiti = ht3$ht,
#                  N = nrow(ht3),
#                  n_spec = length(unique(ht3$species)),
#                  trait_species = as.numeric(as.factor(ht3$species)),
#                  n_site = length(unique(ht3$site.n)),
#                  site = as.numeric(as.factor(ht3$site)),
#                  prior_mu_grand_mu = 20,
#                  prior_mu_grand_sigma = 10,
#                  prior_sigma_sp_mu = 4,
#                  prior_sigma_sp_sigma = 5,
#                  prior_sigma_site_mu = 2,
#                  prior_sigma_site_sigma = 10,
#                  prior_sigma_traity_mu = 3,
#                  prior_sigma_traity_sigma = 5,
#                  ## Phenology
#                  Nph = nrow(pheno.tW),
#                  phenology_species = as.numeric(as.factor(pheno.tW$species)),
#                  yPhenoi = pheno.tW$bb,
#                  forcei = pheno.tW$force.z2,
#                  chilli = pheno.tW$chillport.z2,
#                  photoi = pheno.tW$photo.z2,
#                  prior_muForceSp_mu = -15,
#                  prior_muForceSp_sigma = 10, #wider
#                  prior_muChillSp_mu = -15,
#                  prior_muChillSp_sigma = 10,#wider
#                  prior_muPhotoSp_mu = -15,
#                  prior_muPhotoSp_sigma = 10,#wider
#                  prior_muPhenoSp_mu = 40,
#                  prior_muPhenoSp_sigma = 10,#wider
#                  prior_sigmaForceSp_mu = 5,
#                  prior_sigmaForceSp_sigma = 5,
#                  prior_sigmaChillSp_mu = 5,#wider
#                  prior_sigmaChillSp_sigma = 5, #wider
#                  prior_sigmaPhotoSp_mu = 5,
#                  prior_sigmaPhotoSp_sigma = 5,
#                  prior_sigmaPhenoSp_mu = 5, #wider
#                  prior_sigmaPhenoSp_sigma = 5, #wider
#                  prior_betaTraitxForce_mu = 0,
#                  prior_betaTraitxForce_sigma = 1,
#                  prior_betaTraitxChill_mu = 0,
#                  prior_betaTraitxChill_sigma = 1,
#                  prior_betaTraitxPhoto_mu = 0,
#                  prior_betaTraitxPhoto_sigma = 1,
#                  prior_sigmaphenoy_mu = 10,
#                  prior_sigmaphenoy_sigma = 5 #wider
# )


#ht.data$site
# mdl.ncp <- stan("stan/jointMdl_Feb13_ncpPhoto.stan",
#                 data = ht3.data,
#                 iter = 6000,
#                 warmup = 4000,
#                 chains = 4, include = FALSE, pars = c("y_hat"),
#                 control = list(adapt_delta =0.99))
# 
# # pairs(mdl.ht3, pars = c("mu_grand", "mu_force", "mu_chill", "mu_photo","b_site2","b_site3","b_site4","sigma_a", "sigma_force", "sigma_chill", "sigma_force","sigma_y", "lp__")) 
# 
# save(mdl.ncp, file = "output/ht_westernNCP.Rda")
# sumer3 <- summary(mdl.ncp)$summary

 # 310 divergent transitions with original priors 
# Doubled sigma = 2187 div trans and bad RHat
## N effective?
# range(summary(mdl.ncp)$summary[, "n_eff"]) # 394.271, 13592.123
# range(summary(mdl.ncp)$summary[, "Rhat"])
# 
# ### Add species and study names to Stan object
# names(mdl.ht3)[grep(pattern = "^muSp", x = names(mdl.ht3))] <- paste(specieslist, sep = "")
# names(mdl.ht3)[grep(pattern = "^b_site", x = names(mdl.ht3))] <- paste(sitelist, sep = "")
# ##
# names(mdl.ht3)[grep(pattern = "^alphaForceSp", x = names(mdl.ht3))] <- paste(specieslist, sep = "")
# names(mdl.ht3)[grep(pattern = "^alphaChillSp", x = names(mdl.ht3))] <- paste(specieslist, sep = "")
# names(mdl.ht3)[grep(pattern = "^alphaPhotoSp", x = names(mdl.ht3))] <- paste(specieslist, sep = "")
# names(mdl.ht3)[grep(pattern = "^alphaPhenoSp", x = names(mdl.ht3))] <- paste(specieslist, sep = "")
# names(mdl.ht3)[grep(pattern = "^betaForceSp", x = names(mdl.ht3))] <- paste(specieslist, sep = "")
# names(mdl.ht3)[grep(pattern = "^betaChillSp", x = names(mdl.ht3))] <- paste(specieslist, sep = "")
# names(mdl.ht3)[grep(pattern = "^betaPhotoSp", x = names(mdl.ht3))] <- paste(specieslist, sep = "")
# names(mdl.ht3)[grep(pattern = "^betaPhenoSp", x = names(mdl.ht3))] <- paste(specieslist, sep = "")
# 
# saveRDS(object = mdl.ht3, file = "height_western_stanfit.RDS")

# 3. Just the eastern trait data
exPop <- c("sm", "mp")
pheno.tE <- pheno.t[!pheno.t$population %in% exPop, ]

# excl <- c("kl", "af", "mp", "sm")
# ht4 <- height[!height$site %in% excl, ]
# 
# ht4.data <- list(yTraiti = ht4$ht,
#                  N = nrow(ht4),
#                  n_spec = length(unique(ht4$species)),
#                  trait_species = as.numeric(as.factor(ht4$species)),
#                  n_site = length(unique(ht4$site.n)),
#                  site = as.numeric(as.factor(ht4$site.n)),
#                  prior_mu_grand_mu = 20,
#                  prior_mu_grand_sigma = 10,
#                  prior_sigma_sp_mu = 4,
#                  prior_sigma_sp_sigma = 5,
#                  prior_sigma_site_mu = 2,
#                  prior_sigma_site_sigma = 5,
#                  prior_sigma_traity_mu = 3,
#                  prior_sigma_traity_sigma = 5,
#                  ## Phenology
#                  Nph = nrow(pheno.tE),
#                  phenology_species = as.numeric(as.factor(pheno.tE$species)),
#                  yPhenoi = pheno.tE$bb,
#                  forcei = pheno.tE$force.z2,
#                  chilli = pheno.tE$chillport.z2,
#                  photoi = pheno.tE$photo.z2,
#                  prior_muForceSp_mu = -15,
#                  prior_muForceSp_sigma = 10, #wider
#                  prior_muChillSp_mu = -15,
#                  prior_muChillSp_sigma = 10,#wider
#                  prior_muPhotoSp_mu = -15,
#                  prior_muPhotoSp_sigma = 10,#wider
#                  prior_muPhenoSp_mu = 40,
#                  prior_muPhenoSp_sigma = 10,#wider
#                  prior_sigmaForceSp_mu = 5,
#                  prior_sigmaForceSp_sigma = 5,
#                  prior_sigmaChillSp_mu = 5,#wider
#                  prior_sigmaChillSp_sigma = 5, #wider
#                  prior_sigmaPhotoSp_mu = 5,
#                  prior_sigmaPhotoSp_sigma = 5,
#                  prior_sigmaPhenoSp_mu = 5, #wider
#                  prior_sigmaPhenoSp_sigma = 5, #wider
#                  prior_betaTraitxForce_mu = 0,
#                  prior_betaTraitxForce_sigma = 1,
#                  prior_betaTraitxChill_mu = 0,
#                  prior_betaTraitxChill_sigma = 1,
#                  prior_betaTraitxPhoto_mu = 0,
#                  prior_betaTraitxPhoto_sigma = 1,
#                  prior_sigmaphenoy_mu = 10,
#                  prior_sigmaphenoy_sigma = 5 #wider
# )
# 
# #ht.data$site
# mdl.ht4 <- stan("stan/jointMdl_Feb13_NoPhoto.stan",
#                 data = ht4.data,
#                 iter = 6000,
#                 warmup = 4000,
#                 chains = 4)
# # include = FALSE, pars = c("y_hat"))
# 
# save(mdl.ht4, file = "output/ht_eastern_FC.Rda")
# ## N effective?
# summary(mdl.ht4)$summary[, "n_eff"] # 394.271, 13592.123
# summary(mdl.ht4)$summary[, "Rhat"]
# ### Add species and study names to Stan object
# names(mdl.ht4)[grep(pattern = "^muSp", x = names(mdl.ht4))] <- paste(specieslist, sep = "")
# names(mdl.ht4)[grep(pattern = "^b_site", x = names(mdl.ht4))] <- paste(sitelist, sep = "")
# ##
# names(mdl.ht4)[grep(pattern = "^alphaForceSp", x = names(mdl.ht4))] <- paste(specieslist, sep = "")
# names(mdl.ht4)[grep(pattern = "^alphaChillSp", x = names(mdl.ht4))] <- paste(specieslist, sep = "")
# names(mdl.ht4)[grep(pattern = "^alphaPhotoSp", x = names(mdl.ht4))] <- paste(specieslist, sep = "")
# names(mdl.ht4)[grep(pattern = "^alphaPhenoSp", x = names(mdl.ht4))] <- paste(specieslist, sep = "")
# names(mdl.ht4)[grep(pattern = "^betaForceSp", x = names(mdl.ht4))] <- paste(specieslist, sep = "")
# names(mdl.ht4)[grep(pattern = "^betaChillSp", x = names(mdl.ht4))] <- paste(specieslist, sep = "")
# names(mdl.ht4)[grep(pattern = "^betaPhotoSp", x = names(mdl.ht4))] <- paste(specieslist, sep = "")
# names(mdl.ht4)[grep(pattern = "^betaPhenoSp", x = names(mdl.ht4))] <- paste(specieslist, sep = "")
# 
# saveRDS(object = mdl.ht4, file = "height_eastern_stanfit.RDS")

##########################################################
# LMA
# 8 sites instead of 2 transects
specieslist <- sort(unique(trtPheno$species))
sitelist <- sort(unique(trtPheno$transect))
leafMass <- trtPheno[complete.cases(trtPheno$lma),]


# 2. Just the western trait data
excl <- c("HF", "SH", "GR", "WM")
leafMass3 <- leafMass[!leafMass$site %in% excl, ]
unique(leafMass3$site)

lma3.data <- list(yTraiti = leafMass3$lma,
                  N = nrow(leafMass3),
                  n_spec = length(unique(leafMass3$species)),
                  trait_species = as.numeric(as.factor(leafMass3$species)),
                  n_site = length(unique(leafMass3$site.n)),
                  site = as.numeric(as.factor(leafMass3$site.n)),
                  prior_mu_grand_mu = 0.5,
                  prior_mu_grand_sigma = 5, #widened
                  prior_sigma_sp_mu = 10,
                  prior_sigma_sp_sigma = 5,
                  prior_sigma_site_mu = 5,
                  prior_sigma_site_sigma = 2,
                  prior_sigma_traity_mu = 5,
                  prior_sigma_traity_sigma = 2,
                  ## Phenology
                  Nph = nrow(pheno.tW),
                  phenology_species = as.numeric(as.factor(pheno.tW$species)),
                  yPhenoi = pheno.tW$bb,
                  forcei = pheno.tW$force.z2,
                  chilli = pheno.tW$chillport.z2,
                  photoi = pheno.tW$photo.z2,
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


mdl.lma3 <- stan("stan/jointMdl_Feb13_ncpPhoto.stan",
                 data = lma3.data,
                 iter = 6000,
                 warmup = 4000,
             chains = 4, include = FALSE, pars = c("y_hat"),
             control = list(max_treedepth = 15))
.# include = FALSE, pars = c("y_hat"),
# control = list(max_treedepth = 11))

save(mdl.lma3, file = "output/lma_western_ncp.Rda")
## N effective?
# range(summary(mdl.lma3)$summary[, "n_eff"] )# 394.271, 13592.123
# 
# ### Add species and study names to Stan object
# names(mdl.lma3)[grep(pattern = "^muSp", x = names(mdl.lma3))] <- paste(specieslist, sep = "")
# names(mdl.lma3)[grep(pattern = "^b_site", x = names(mdl.lma3))] <- paste(sitelist, sep = "")
# ##
# names(mdl.lma3)[grep(pattern = "^alphaForceSp", x = names(mdl.lma3))] <- paste(specieslist, sep = "")
# names(mdl.lma3)[grep(pattern = "^alphaChillSp", x = names(mdl.lma3))] <- paste(specieslist, sep = "")
# names(mdl.lma3)[grep(pattern = "^alphaPhotoSp", x = names(mdl.lma3))] <- paste(specieslist, sep = "")
# names(mdl.lma3)[grep(pattern = "^alphaPhenoSp", x = names(mdl.lma3))] <- paste(specieslist, sep = "")
# names(mdl.lma3)[grep(pattern = "^betaForceSp", x = names(mdl.lma3))] <- paste(specieslist, sep = "")
# names(mdl.lma3)[grep(pattern = "^betaChillSp", x = names(mdl.lma3))] <- paste(specieslist, sep = "")
# names(mdl.lma3)[grep(pattern = "^betaPhotoSp", x = names(mdl.lma3))] <- paste(specieslist, sep = "")
# names(mdl.lma3)[grep(pattern = "^betaPhenoSp", x = names(mdl.lma3))] <- paste(specieslist, sep = "")
# 
# saveRDS(object = mdl.lma3, file = "lma_western_stanfit.RDS")

# 2. Just the eastern trait data
excl <- c("sm", "mp","af","kl")
leafMass4 <- leafMass[!leafMass$site %in% excl, ]
unique(leafMass4$site)

# lma4.data <- list(yTraiti = leafMass4$lma,
#                   N = nrow(leafMass4),
#                   n_spec = length(unique(leafMass4$species)),
#                   trait_species = as.numeric(as.factor(leafMass4$species)),
#                   n_site = length(unique(leafMass4$site.n)),
#                   site = as.numeric(as.factor(leafMass4$site.n)),
#                   prior_mu_grand_mu = 0.5,
#                   prior_mu_grand_sigma = 5, #widened
#                   prior_sigma_sp_mu = 10,
#                   prior_sigma_sp_sigma = 5,
#                   prior_sigma_site_mu = 5,
#                   prior_sigma_site_sigma = 2,
#                   prior_sigma_traity_mu = 5,
#                   prior_sigma_traity_sigma = 2,
#                   ## Phenology
#                   Nph = nrow(pheno.tE),
#                   phenology_species = as.numeric(as.factor(pheno.tE$species)),
#                   yPhenoi = pheno.tE$bb,
#                   forcei = pheno.tE$force.z2,
#                   chilli = pheno.tE$chillport.z2,
#                   photoi = pheno.tE$photo.z2,
#                   prior_muForceSp_mu = -15,
#                   prior_muForceSp_sigma = 10, #wider
#                   prior_muChillSp_mu = -15,
#                   prior_muChillSp_sigma = 10,#wider
#                   prior_muPhotoSp_mu = -15,
#                   prior_muPhotoSp_sigma = 10,#wider
#                   prior_muPhenoSp_mu = 40,
#                   prior_muPhenoSp_sigma = 10,#wider
#                   prior_sigmaForceSp_mu = 5,
#                   prior_sigmaForceSp_sigma = 5,
#                   prior_sigmaChillSp_mu = 5,#wider
#                   prior_sigmaChillSp_sigma = 5, #wider
#                   prior_sigmaPhotoSp_mu = 5,
#                   prior_sigmaPhotoSp_sigma = 5,
#                   prior_sigmaPhenoSp_mu = 5, #wider
#                   prior_sigmaPhenoSp_sigma = 5, #wider
#                   prior_betaTraitxForce_mu = 0,
#                   prior_betaTraitxForce_sigma = 1,
#                   prior_betaTraitxChill_mu = 0,
#                   prior_betaTraitxChill_sigma = 1,
#                   prior_betaTraitxPhoto_mu = 0,
#                   prior_betaTraitxPhoto_sigma = 1,
#                   prior_sigmaphenoy_mu = 10,
#                   prior_sigmaphenoy_sigma = 5 #wider
# )
# 
# 
# mdl.lma4 <- stan("stan/jointMdl.stan",
#                  data = lma4.data,
#                  iter = 4000,
#                  warmup = 3000,
#                  chains = 4, include = FALSE, pars = c("y_hat"))
# # control = list(max_treedepth = 11))
# 
# save(mdl.lma4, file = "output/lma_eastern.Rda")
## N effective?
# range(summary(mdl.lma4)$summary[, "n_eff"] )# 394.271, 13592.123
# 
# ### Add species and study names to Stan object
# names(mdl.lma4)[grep(pattern = "^muSp", x = names(mdl.lma4))] <- paste(specieslist, sep = "")
# names(mdl.lma4)[grep(pattern = "^b_site", x = names(mdl.lma4))] <- paste(sitelist, sep = "")
# ##
# names(mdl.lma4)[grep(pattern = "^alphaForceSp", x = names(mdl.lma4))] <- paste(specieslist, sep = "")
# names(mdl.lma4)[grep(pattern = "^alphaChillSp", x = names(mdl.lma4))] <- paste(specieslist, sep = "")
# names(mdl.lma4)[grep(pattern = "^alphaPhotoSp", x = names(mdl.lma4))] <- paste(specieslist, sep = "")
# names(mdl.lma4)[grep(pattern = "^alphaPhenoSp", x = names(mdl.lma4))] <- paste(specieslist, sep = "")
# names(mdl.lma4)[grep(pattern = "^betaForceSp", x = names(mdl.lma4))] <- paste(specieslist, sep = "")
# names(mdl.lma4)[grep(pattern = "^betaChillSp", x = names(mdl.lma4))] <- paste(specieslist, sep = "")
# names(mdl.lma4)[grep(pattern = "^betaPhotoSp", x = names(mdl.lma4))] <- paste(specieslist, sep = "")
# names(mdl.lma4)[grep(pattern = "^betaPhenoSp", x = names(mdl.lma4))] <- paste(specieslist, sep = "")
# 
# saveRDS(object = mdl.lma4, file = "lma_eastern_stanfit.RDS")
# 

##########################################################
# ssd
stemDen <- trtPheno[complete.cases(trtPheno$ssd),]

# 2. western sites only
excl <- c("HF", "SH", "GR", "WM")
stemDen3 <- stemDen[!stemDen$site %in% excl, ]
unique(stemDen3$site)

ssd3.data <- list(yTraiti = stemDen3$ssd,
                  N = nrow(stemDen3),
                  n_spec = length(unique(stemDen3$species)),
                  trait_species = as.numeric(as.factor(stemDen3$species)),
                  n_site = length(unique(stemDen3$site.n)),
                  site = as.numeric(as.factor(stemDen3$site.n)),
                  prior_mu_grand_mu = 1,
                  prior_mu_grand_sigma = 5, #widened
                  prior_sigma_sp_mu = 4, #10
                  prior_sigma_sp_sigma = 5,
                  prior_sigma_site_mu = 2, #5
                  prior_sigma_site_sigma = 5, #2
                  prior_sigma_traity_mu = 2, #5
                  prior_sigma_traity_sigma = 5,
                  ## Phenology
                  Nph = nrow(pheno.tW),
                  phenology_species = as.numeric(as.factor(pheno.tW$species)),
                  yPhenoi = pheno.tW$bb,
                  forcei = pheno.tW$force.z2,
                  chilli = pheno.tW$chillport.z2,
                  photoi = pheno.tW$photo.z2,
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


mdl.ssd3 <- stan("stan/jointMdl_Feb13.stan",
                 data = ssd3.data,
                 iter = 6000,
                 warmup = 4000,
               chains = 4, include = FALSE, pars = c("y_hat"))
             #  control = list(max_treedepth = 15))

save(mdl.ssd3, file = "output/ssd_western_ncp.Rda")
## N effective?
range(summary(mdl.ssd3)$summary[, "n_eff"])
range(summary(mdl.ssd3)$summary[, "Rhat"])# 1013.934 8423.027
# 
# ### Add species and study names to Stan object
# names(mdl.ssd3)[grep(pattern = "^muSp", x = names(mdl.ssd3))] <- paste(specieslist, sep = "")
# names(mdl.ssd3)[grep(pattern = "^b_site", x = names(mdl.ssd3))] <- paste(sitelist, sep = "")
# ##
# names(mdl.ssd3)[grep(pattern = "^alphaForceSp", x = names(mdl.ssd3))] <- paste(specieslist, sep = "")
# names(mdl.ssd3)[grep(pattern = "^alphaChillSp", x = names(mdl.ssd3))] <- paste(specieslist, sep = "")
# names(mdl.ssd3)[grep(pattern = "^alphaPhotoSp", x = names(mdl.ssd3))] <- paste(specieslist, sep = "")
# names(mdl.ssd3)[grep(pattern = "^alphaPhenoSp", x = names(mdl.ssd3))] <- paste(specieslist, sep = "")
# names(mdl.ssd3)[grep(pattern = "^betaForceSp", x = names(mdl.ssd3))] <- paste(specieslist, sep = "")
# names(mdl.ssd3)[grep(pattern = "^betaChillSp", x = names(mdl.ssd3))] <- paste(specieslist, sep = "")
# names(mdl.ssd3)[grep(pattern = "^betaPhotoSp", x = names(mdl.ssd3))] <- paste(specieslist, sep = "")
# names(mdl.ssd3)[grep(pattern = "^betaPhenoSp", x = names(mdl.ssd3))] <- paste(specieslist, sep = "")
# 
# saveRDS(object = mdl.ssd3, file = "SSD_western_stanfit_np.RDS")

# 3. eastern sites
excl <- c("kl", "af", "mp", "sm")
stemDen4 <- stemDen[!stemDen$site %in% excl, ]
unique(stemDen4$site)

# ssd4.data <- list(yTraiti = stemDen4$ssd,
#                   N = nrow(stemDen4),
#                   n_spec = length(unique(stemDen4$species)),
#                   trait_species = as.numeric(as.factor(stemDen4$species)),
#                   n_site = length(stemDen4$site.n),
#                   site = as.numeric(as.factor(stemDen4$site.n)),
#                   prior_mu_grand_mu = 1,
#                   prior_mu_grand_sigma = 5, #widened
#                   prior_sigma_sp_mu = 4, #10
#                   prior_sigma_sp_sigma = 5,
#                   prior_sigma_site_mu = 2, #5
#                   prior_sigma_site_sigma = 5, #2
#                   prior_sigma_traity_mu = 2, #5
#                   prior_sigma_traity_sigma = 5,
#                   ## Phenology
#                   Nph = nrow(pheno.tE),
#                   phenology_species = as.numeric(as.factor(pheno.tE$species)),
#                   yPhenoi = pheno.tE$bb,
#                   forcei = pheno.tE$force.z2,
#                   chilli = pheno.tE$chillport.z2,
#                   photoi = pheno.tE$photo.z2,
#                   prior_muForceSp_mu = -15,
#                   prior_muForceSp_sigma = 10, #wider
#                   prior_muChillSp_mu = -15,
#                   prior_muChillSp_sigma = 10,#wider
#                   prior_muPhotoSp_mu = -15,
#                   prior_muPhotoSp_sigma = 10,#wider
#                   prior_muPhenoSp_mu = 40,
#                   prior_muPhenoSp_sigma = 10,#wider
#                   prior_sigmaForceSp_mu = 5,
#                   prior_sigmaForceSp_sigma = 5,
#                   prior_sigmaChillSp_mu = 5,#wider
#                   prior_sigmaChillSp_sigma = 5, #wider
#                   prior_sigmaPhotoSp_mu = 5,
#                   prior_sigmaPhotoSp_sigma = 5,
#                   prior_sigmaPhenoSp_mu = 5, #wider
#                   prior_sigmaPhenoSp_sigma = 5, #wider
#                   prior_betaTraitxForce_mu = 0,
#                   prior_betaTraitxForce_sigma = 1,
#                   prior_betaTraitxChill_mu = 0,
#                   prior_betaTraitxChill_sigma = 1,
#                   prior_betaTraitxPhoto_mu = 0,
#                   prior_betaTraitxPhoto_sigma = 1,
#                   prior_sigmaphenoy_mu = 10,
#                   prior_sigmaphenoy_sigma = 5 #wider
# )
# 
# 
# mdl.ssd4 <- stan("stan/jointMdl.stan",
#                  data = ssd4.data,
#                  iter = 4000,
#                  warmup = 3000,
#                  chains = 4,include = FALSE, pars = c("y_hat"))
# 
# save(mdl.ssd4, file = "output/ssd_eastern.Rda")
## N effective?
# summary(mdl.ssd4)$summary[, "n_eff"] # 1013.934 8443.047
# 
# ### Add species and study names to Stan object
# names(mdl.ssd4)[grep(pattern = "^muSp", x = names(mdl.ssd4))] <- paste(specieslist, sep = "")
# names(mdl.ssd4)[grep(pattern = "^b_site", x = names(mdl.ssd4))] <- paste(sitelist, sep = "")
# ##
# names(mdl.ssd4)[grep(pattern = "^alphaForceSp", x = names(mdl.ssd4))] <- paste(specieslist, sep = "")
# names(mdl.ssd4)[grep(pattern = "^alphaChillSp", x = names(mdl.ssd4))] <- paste(specieslist, sep = "")
# names(mdl.ssd4)[grep(pattern = "^alphaPhotoSp", x = names(mdl.ssd4))] <- paste(specieslist, sep = "")
# names(mdl.ssd4)[grep(pattern = "^alphaPhenoSp", x = names(mdl.ssd4))] <- paste(specieslist, sep = "")
# names(mdl.ssd4)[grep(pattern = "^betaForceSp", x = names(mdl.ssd4))] <- paste(specieslist, sep = "")
# names(mdl.ssd4)[grep(pattern = "^betaChillSp", x = names(mdl.ssd4))] <- paste(specieslist, sep = "")
# names(mdl.ssd4)[grep(pattern = "^betaPhotoSp", x = names(mdl.ssd4))] <- paste(specieslist, sep = "")
# names(mdl.ssd4)[grep(pattern = "^betaPhenoSp", x = names(mdl.ssd4))] <- paste(specieslist, sep = "")
# 
# saveRDS(object = mdl.ssd4, file = "SSD_eastern_stanfit_np.RDS")
# ##########################################################
# # dbh
# 
diam <- trtPheno[complete.cases(trtPheno$dbh),]


# 2. Just the western trait data
excl <- c("SH", "HF", "GR", "WM")
diam3 <- diam[!diam$site %in% excl, ]
unique(diam3$site)

# dbh3.data <- list(yTraiti = diam3$dbh,
#                   N = nrow(diam3),
#                   n_spec = length(unique(diam3$species)),
#                   trait_species = as.numeric(as.factor(diam3$species)),
#                   n_site = length(unique(diam3$site.n)),
#                   site = as.numeric(as.factor(diam3$site.n)),
#                   prior_mu_grand_mu = 15,
#                   prior_mu_grand_sigma = 5, #widened
#                   prior_sigma_sp_mu = 10,
#                   prior_sigma_sp_sigma = 5,
#                   prior_sigma_site_mu = 5,
#                   prior_sigma_site_sigma = 2,
#                   prior_sigma_traity_mu = 5,
#                   prior_sigma_traity_sigma = 2,
#                   ## Phenology
#                   Nph = nrow(pheno.tW),
#                   phenology_species = as.numeric(as.factor(pheno.tW$species)),
#                   yPhenoi = pheno.tW$bb,
#                   forcei = pheno.tW$force.z2,
#                   chilli = pheno.tW$chillport.z2,
#                   photoi = pheno.tW$photo.z2,
#                   prior_muForceSp_mu = -15,
#                   prior_muForceSp_sigma = 10, #wider
#                   prior_muChillSp_mu = -15,
#                   prior_muChillSp_sigma = 10,#wider
#                   prior_muPhotoSp_mu = -15,
#                   prior_muPhotoSp_sigma = 10,#wider
#                   prior_muPhenoSp_mu = 40,
#                   prior_muPhenoSp_sigma = 10,#wider
#                   prior_sigmaForceSp_mu = 5,
#                   prior_sigmaForceSp_sigma = 5,
#                   prior_sigmaChillSp_mu = 5,#wider
#                   prior_sigmaChillSp_sigma = 5, #wider
#                   prior_sigmaPhotoSp_mu = 5,
#                   prior_sigmaPhotoSp_sigma = 5,
#                   prior_sigmaPhenoSp_mu = 5, #wider
#                   prior_sigmaPhenoSp_sigma = 5, #wider
#                   prior_betaTraitxForce_mu = 0,
#                   prior_betaTraitxForce_sigma = 1,
#                   prior_betaTraitxChill_mu = 0,
#                   prior_betaTraitxChill_sigma = 1,
#                   prior_betaTraitxPhoto_mu = 0,
#                   prior_betaTraitxPhoto_sigma = 1,
#                   prior_sigmaphenoy_mu = 10,
#                   prior_sigmaphenoy_sigma = 5 #wider
# )
# 
# 
# mdl.dbh3 <- stan("stan/jointMdl_Feb13_ncpPhoto.stan",
#                  data = dbh3.data,
#                  iter = 6000,
#                  warmup = 3000,
#                  chains = 4, include = FALSE, pars = c("y_hat"),
#                  control = list(adapt_delta =0.99))
# # run for more iterations
# 
# save(mdl.dbh3, file = "output/dbh_western_ncp.Rda")
# ## N effective?
# # range(summary(mdl.dbh3)$summary[, "n_eff"]) # 394.371, 13593.133
# # 
# # ### Add species and study names to Stan object
# # names(mdl.dbh3)[grep(pattern = "^muSp", x = names(mdl.dbh3))] <- paste(specieslist, sep = "")
# # names(mdl.dbh3)[grep(pattern = "^b_site", x = names(mdl.dbh3))] <- paste(sitelist, sep = "")
# # ##
# # names(mdl.dbh3)[grep(pattern = "^alphaForceSp", x = names(mdl.dbh3))] <- paste(specieslist, sep = "")
# # names(mdl.dbh3)[grep(pattern = "^alphaChillSp", x = names(mdl.dbh3))] <- paste(specieslist, sep = "")
# # names(mdl.dbh3)[grep(pattern = "^alphaPhotoSp", x = names(mdl.dbh3))] <- paste(specieslist, sep = "")
# # names(mdl.dbh3)[grep(pattern = "^alphaPhenoSp", x = names(mdl.dbh3))] <- paste(specieslist, sep = "")
# # names(mdl.dbh3)[grep(pattern = "^betaForceSp", x = names(mdl.dbh3))] <- paste(specieslist, sep = "")
# # names(mdl.dbh3)[grep(pattern = "^betaChillSp", x = names(mdl.dbh3))] <- paste(specieslist, sep = "")
# # names(mdl.dbh3)[grep(pattern = "^betaPhotoSp", x = names(mdl.dbh3))] <- paste(specieslist, sep = "")
# # names(mdl.dbh3)[grep(pattern = "^betaPhenoSp", x = names(mdl.dbh))] <- paste(specieslist, sep = "")
# # 
# # saveRDS(object = mdl.dbh3, file = "dbh_western_stanfit.RDS")
# 
# # 3. Just the eastern trait data
# excl <- c("kl", "af", "sm", "mp")
# diam4 <- diam[!diam$site %in% excl, ]
# unique(diam4$site)
# 
# dbh4.data <- list(yTraiti = diam4$dbh,
#                   N = nrow(diam4),
#                   n_spec = length(unique(diam4$species)),
#                   trait_species = as.numeric(as.factor(diam4$species)),
#                   n_site = length(unique(diam4$site.n)),
#                   site = as.numeric(as.factor(diam4$site.n)),
#                   prior_mu_grand_mu = 15,
#                   prior_mu_grand_sigma = 5, #widened
#                   prior_sigma_sp_mu = 10,
#                   prior_sigma_sp_sigma = 5,
#                   prior_sigma_site_mu = 5,
#                   prior_sigma_site_sigma = 2,
#                   prior_sigma_traity_mu = 5,
#                   prior_sigma_traity_sigma = 2,
#                   ## Phenology
#                   Nph = nrow(pheno.tE),
#                   phenology_species = as.numeric(as.factor(pheno.tE$species)),
#                   yPhenoi = pheno.tE$bb,
#                   forcei = pheno.tE$force.z2,
#                   chilli = pheno.tE$chillport.z2,
#                   photoi = pheno.tE$photo.z2,
#                   prior_muForceSp_mu = -15,
#                   prior_muForceSp_sigma = 10, #wider
#                   prior_muChillSp_mu = -15,
#                   prior_muChillSp_sigma = 10,#wider
#                   prior_muPhotoSp_mu = -15,
#                   prior_muPhotoSp_sigma = 10,#wider
#                   prior_muPhenoSp_mu = 40,
#                   prior_muPhenoSp_sigma = 10,#wider
#                   prior_sigmaForceSp_mu = 5,
#                   prior_sigmaForceSp_sigma = 5,
#                   prior_sigmaChillSp_mu = 5,#wider
#                   prior_sigmaChillSp_sigma = 5, #wider
#                   prior_sigmaPhotoSp_mu = 5,
#                   prior_sigmaPhotoSp_sigma = 5,
#                   prior_sigmaPhenoSp_mu = 5, #wider
#                   prior_sigmaPhenoSp_sigma = 5, #wider
#                   prior_betaTraitxForce_mu = 0,
#                   prior_betaTraitxForce_sigma = 1,
#                   prior_betaTraitxChill_mu = 0,
#                   prior_betaTraitxChill_sigma = 1,
#                   prior_betaTraitxPhoto_mu = 0,
#                   prior_betaTraitxPhoto_sigma = 1,
#                   prior_sigmaphenoy_mu = 10,
#                   prior_sigmaphenoy_sigma = 5 #wider
# )
# 
# 
# mdl.dbh4 <- stan("stan/jointMdl.stan",
#                  data = dbh4.data,
#                  iter = 6000,
#                  warmup = 3000,
#                  chains = 4, include = FALSE, pars = c("y_hat"))
# # run for more iterations
# 
# save(mdl.dbh4, file = "output/dbh_eastern.Rda")
## N effective?
# range(summary(mdl.dbh4)$summary[, "n_eff"]) # 394.471, 13594.143
# 
# ### Add species and study names to Stan object
# names(mdl.dbh4)[grep(pattern = "^muSp", x = names(mdl.dbh4))] <- paste(specieslist, sep = "")
# names(mdl.dbh4)[grep(pattern = "^b_site", x = names(mdl.dbh4))] <- paste(sitelist, sep = "")
# ##
# names(mdl.dbh4)[grep(pattern = "^alphaForceSp", x = names(mdl.dbh4))] <- paste(specieslist, sep = "")
# names(mdl.dbh4)[grep(pattern = "^alphaChillSp", x = names(mdl.dbh4))] <- paste(specieslist, sep = "")
# names(mdl.dbh4)[grep(pattern = "^alphaPhotoSp", x = names(mdl.dbh4))] <- paste(specieslist, sep = "")
# names(mdl.dbh4)[grep(pattern = "^alphaPhenoSp", x = names(mdl.dbh4))] <- paste(specieslist, sep = "")
# names(mdl.dbh4)[grep(pattern = "^betaForceSp", x = names(mdl.dbh4))] <- paste(specieslist, sep = "")
# names(mdl.dbh4)[grep(pattern = "^betaChillSp", x = names(mdl.dbh4))] <- paste(specieslist, sep = "")
# names(mdl.dbh4)[grep(pattern = "^betaPhotoSp", x = names(mdl.dbh4))] <- paste(specieslist, sep = "")
# names(mdl.dbh4)[grep(pattern = "^betaPhenoSp", x = names(mdl.dbh))] <- paste(specieslist, sep = "")
# 
# saveRDS(object = mdl.dbh4, file = "dbh_eastern_stanfit.RDS")
# # ##########################################################
# # C:N

carbNit <- trtPheno[complete.cases(trtPheno$C.N),]


#2. western sites
excl <- c("HF", "SH", "GR", "WM")
carbNit3 <- carbNit[!carbNit$site %in% excl, ]
unique(carbNit3$site)

cn3.data <- list(yTraiti = carbNit3$C.N,
                 N = nrow(carbNit3),
                 n_spec = length(unique(carbNit3$species)),
                 trait_species = as.numeric(as.factor(carbNit3$species)),
                 n_site = length(unique(carbNit3$site.n)),
                 site = as.numeric(as.factor(carbNit3$site.n)),
                 prior_mu_grand_mu = 20,
                 prior_mu_grand_sigma = 5, #widened
                 prior_sigma_sp_mu = 10,
                 prior_sigma_sp_sigma = 5,
                 prior_sigma_site_mu = 5,
                 prior_sigma_site_sigma = 2,
                 prior_sigma_traity_mu = 5,
                 prior_sigma_traity_sigma = 2,
                 ## Phenology
                 Nph = nrow(pheno.tW),
                 phenology_species = as.numeric(as.factor(pheno.tW$species)),
                 yPhenoi = pheno.tW$bb,
                 forcei = pheno.tW$force.z2,
                 chilli = pheno.tW$chillport.z2,
                 photoi = pheno.tW$photo.z2,
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


mdl.cn3 <- stan("stan/jointMdl_Feb13_ncpPhoto.stan",
                data = cn3.data,
                iter = 6000,
                warmup = 4000,
                chains = 4, include = FALSE, pars = c("y_hat"))
              

# bluk ESS is low - run for more iterations

save(mdl.cn3, file = "output/cn_western_ncp.Rda")
## N effective?
# range(summary(mdl.cn3)$summary[, "n_eff"]) #
# 
# ### Add species and study names to Stan object
# names(mdl.cn3)[grep(pattern = "^muSp", x = names(mdl.cn3))] <- paste(specieslist, sep = "")
# names(mdl.cn3)[grep(pattern = "^b_site", x = names(mdl.cn3))] <- paste(sitelist, sep = "")
# ##
# names(mdl.cn3)[grep(pattern = "^alphaForceSp", x = names(mdl.cn3))] <- paste(specieslist, sep = "")
# names(mdl.cn3)[grep(pattern = "^alphaChillSp", x = names(mdl.cn3))] <- paste(specieslist, sep = "")
# names(mdl.cn3)[grep(pattern = "^alphaPhotoSp", x = names(mdl.cn3))] <- paste(specieslist, sep = "")
# names(mdl.cn3)[grep(pattern = "^alphaPhenoSp", x = names(mdl.cn3))] <- paste(specieslist, sep = "")
# names(mdl.cn3)[grep(pattern = "^betaForceSp", x = names(mdl.cn3))] <- paste(specieslist, sep = "")
# names(mdl.cn3)[grep(pattern = "^betaChillSp", x = names(mdl.cn3))] <- paste(specieslist, sep = "")
# names(mdl.cn3)[grep(pattern = "^betaPhotoSp", x = names(mdl.cn3))] <- paste(specieslist, sep = "")
# names(mdl.cn3)[grep(pattern = "^betaPhenoSp", x = names(mdl.cn3))] <- paste(specieslist, sep = "")
# 
# saveRDS(object = mdl.cn3, file = "cn_western_stanfit_np.RDS")

# 3. eastern sites only 
excl <- c("kl", "af", "sm", "mp")
carbNit4 <- carbNit[!carbNit$site %in% excl, ]
unique(carbNit4$site)

cn4.data <- list(yTraiti = carbNit4$C.N,
                 N = nrow(carbNit4),
                 n_spec = length(unique(carbNit4$species)),
                 trait_species = as.numeric(as.factor(carbNit4$species)),
                 n_site = length(unique(carbNit4$site.n)),
                 site = as.numeric(as.factor(carbNit4$site.n)),
                 prior_mu_grand_mu = 20,
                 prior_mu_grand_sigma = 5, #widened
                 prior_sigma_sp_mu = 10,
                 prior_sigma_sp_sigma = 5,
                 prior_sigma_site_mu = 5,
                 prior_sigma_site_sigma = 2,
                 prior_sigma_traity_mu = 5,
                 prior_sigma_traity_sigma = 2,
                 ## Phenology
                 Nph = nrow(pheno.tE),
                 phenology_species = as.numeric(as.factor(pheno.tE$species)),
                 yPhenoi = pheno.tE$bb,
                 forcei = pheno.tE$force.z2,
                 chilli = pheno.tE$chillport.z2,
                 photoi = pheno.tE$photo.z2,
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


mdl.cn4 <- stan("stan/jointMdl.stan",
                data = cn4.data,
                iter = 6000,
                warmup = 5000,
                chains = 4, include = FALSE, pars = c("y_hat"),
                control = list(adapt_delta =0.99))

# bluk ESS is low - run for more iterations

save(mdl.cn4, file = "output/cn_eastern.Rda")
## N effective?
sumcn4 <- summary(mdl.cn4)$summary
range(summary(mdl.cn4)$summary[, "n_eff"])
range(summary(mdl.cn4)$summary[, "Rhat"])#
# 
# ### Add species and study names to Stan object
# names(mdl.cn4)[grep(pattern = "^muSp", x = names(mdl.cn4))] <- paste(specieslist, sep = "")
# names(mdl.cn4)[grep(pattern = "^b_site", x = names(mdl.cn4))] <- paste(sitelist, sep = "")
# ##
# names(mdl.cn4)[grep(pattern = "^alphaForceSp", x = names(mdl.cn4))] <- paste(specieslist, sep = "")
# names(mdl.cn4)[grep(pattern = "^alphaChillSp", x = names(mdl.cn4))] <- paste(specieslist, sep = "")
# names(mdl.cn4)[grep(pattern = "^alphaPhotoSp", x = names(mdl.cn4))] <- paste(specieslist, sep = "")
# names(mdl.cn4)[grep(pattern = "^alphaPhenoSp", x = names(mdl.cn4))] <- paste(specieslist, sep = "")
# names(mdl.cn4)[grep(pattern = "^betaForceSp", x = names(mdl.cn4))] <- paste(specieslist, sep = "")
# names(mdl.cn4)[grep(pattern = "^betaChillSp", x = names(mdl.cn4))] <- paste(specieslist, sep = "")
# names(mdl.cn4)[grep(pattern = "^betaPhotoSp", x = names(mdl.cn4))] <- paste(specieslist, sep = "")
# names(mdl.cn4)[grep(pattern = "^betaPhenoSp", x = names(mdl.cn4))] <- paste(specieslist, sep = "")
# 
# saveRDS(object = mdl.cn4, file = "cn_eastern_stanfit_np.RDS")
# # plot ypred plots:
# y <- carbNit$C.N
# 
# postCN <- rstan::extract(mdl.cn)
# yrep <-  postCN$y_hat
# 
# pdf("figures/cn_ypred_trt.pdf")
# ppc_dens_overlay(y, yrep[1:100,])
# dev.off()

# sumht <- summary(mdl.8)$summary
# 
# col4table <- c("mean","sd","2.5%","50%","97.5%","Rhat")
# 
# mu_params <-   c("mu_grand",
#                  "muPhenoSp",
#                  "muForceSp", # muAlphaForce
#                  "muChillSp",
#                  "muPhotoSp",
#                  "1",
#                  "betaTraitxForce",
#                  "betaTraitxChill",
#                  "betaTraitxPhoto",
#                  "sigma_sp",
#                  "sigma_site",
#                  "sigma_traity",
#                  "sigmaPhenoSp",
#                  "sigmaForceSp",
#                  "sigmaChillSp",
#                  "sigmaPhotoSp", #sigma_alpha_photo
#                  "sigmapheno_y")
# esti <- sumht[mu_params, col4table]
# 
# #temp <- c(mugrandtrait, muStudy, muGrandSpname, betaForceSpname)
# rownames(esti) =c("mu_grand",
#                   "muPhenoSp",
#                   "muForceSp",
#                   "muChillSp",
#                   "muPhotoSp",
#                   "b_site",
#                   "betaTraitxForce",
#                   "betaTraitxChill",
#                   "betaTraitxPhoto",
#                   "sigma_sp",
#                   "sigma_site",
#                   "sigma_traity",
#                   "sigmaPhenoSp",
#                   "sigmaForceSp",
#                   "sigmaChillSp",
#                   "sigmaPhotoSp",
#                   "sigmapheno_y")
# 
# esti.table <- sumht[mu_params, col4table]
# row.names(esti.table) <- row.names(esti)
# esti.table
# 
# 
# # Why is the western model so bad?
# 
# 
# 
# 
# load("output/ht_western.Rda")
# post <- rstan::extract(mdl.ht3)
# 
# h1 <- hist(rnorm(1000, 20,10), col=rgb(0,0,1,1/4))
# hist(post$mu_grand, add = T,  col=rgb(1,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 40,10), col=rgb(0,0,1,1/4))
# hist(post$muPhenoSp, add = T,  col=rgb(1,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, -15,10), col=rgb(0,0,1,1/4))
# hist(post$muForceSp, add = T,  col=rgb(1,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, -15,10), col=rgb(0,0,1,1/4))
# hist(post$muChillSp, add = T,  col=rgb(1,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, -15,10), col=rgb(0,0,1,1/4))
# hist(post$muPhotoSp, add = T,  col=rgb(1,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,1), col=rgb(0,0,1,1/4))
# hist(post$betaTraitxForce, add = T,  col=rgb(1,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,1), col=rgb(0,0,1,1/4))
# hist(post$betaTraitxChill, add = T,  col=rgb(1,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,1), col=rgb(0,0,1,1/4))
# hist(post$betaTraitxPhoto, add = T,  col=rgb(1,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,5), col=rgb(0,0,1,1/4))
# hist(post$musite, add = T,  col=rgb(1,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 4,5), col=rgb(0,0,1,1/4))
# hist(post$sigma_sp,  add = T,  col=rgb(1,0,1,1/4))
# 
# # not great
# h1 <- hist(rnorm(1000, 2,10), col=rgb(0,0,1,1/4))
# hist(post$sigma_site,  add = T,  col=rgb(1,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 3,5), col=rgb(0,0,1,1/4))
# hist(post$sigma_traity,  add = T,  col=rgb(1,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 5,5), col=rgb(0,0,1,1/4))
# hist(post$sigmaForceSp,  add = T,  col=rgb(1,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 5,5), col=rgb(0,0,1,1/4))
# hist(post$sigmaChillSp,  add = T,  col=rgb(1,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 5,5), col=rgb(0,0,1,1/4))
# hist(post$sigmaPhotoSp,  add = T,  col=rgb(1,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 5,5), col=rgb(0,0,1,1/4))
# hist(post$sigmaPhenoSp,  add = T,  col=rgb(1,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 10,5), col=rgb(0,0,1,1/4))
# hist(post$sigmapheno_y,  add = T,  col=rgb(1,0,1,1/4))
# 
#
pdf("figures/ssdPairsWest2.pdf", width = 6, height = 6)
pairs(mdl.ssd3, pars = c("mu_grand",
                        "muPhenoSp",
                        "muForceSp", # muAlphaForce
                        "muChillSp",
                        "muPhotoSp",
                        "musite",
                        "betaTraitxForce",
                        "betaTraitxChill",
                        "betaTraitxPhoto","sigma_sp",
                        "sigma_site",
                        "sigma_traity",
                        "sigmaPhenoSp",
                        "sigmaForceSp",
                        "sigmaChillSp",
                        "sigmaPhotoSp", #sigma_alpha_photo
                        "sigmapheno_y", "lp__"))
dev.off()

pdf("figures/ssdPairsWest2.pdf", width = 6, height = 6)
pairs(mdl.ssd3, pars = c( "sigma_sp",
                        "sigma_site",
                        "sigma_traity",
                        # "sigmaPhenoSp",
                        # "sigmaForceSp",
                        # "sigmaChillSp",
                        "sigmaPhotoSp", #sigma_alpha_photo
                        "sigmapheno_y", "lp__"))
dev.off()
# par(mfrow = c(4,2))
# plot(post$mu_grand ~ post$muSp[,1])
# plot(post$muPhenoSp ~ post$b_site)
# plot(post$muForceSp ~ post$muPhotoSp)
# plot(post$muForceSp ~ post$muChillSp)
# plot(post$muPhotoSp ~ post$muChillSp)
# plot(post$betaTraitxChill ~ post$muChillSp)
# plot(post$betaTraitxForce ~ post$muForceSp)
# plot(post$betaTraitxPhoto ~ post$muPhotoSp)
# 
# par(mfrow = c(4,2))
# plot(post$sigma_sp ~ post$sigma_site)
# plot(post$sigma_sp ~ post$sigma_traity)
# plot(post$sigma_sp ~ post$sigmaPhenoSp)
# plot(post$sigmaPhenoSp ~ post$sigmaForceSp)
# plot(post$sigmaPhenoSp ~ post$sigmaChillSp)
# plot(post$sigmaPhenoSp ~ post$sigmaPhotoSp)
# plot(post$sigmaPhenoSp ~ post$sigmapheno_y)
# 
# 90+1530+1440+300+846
# 5617.34+1900-5996
# 
# 
