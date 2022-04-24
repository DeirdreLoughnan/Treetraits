# Started April 22, 2022 by Deirdre

# The purpose of this code is to get the traitors joint model working for my trait phenology data! 

rm(list=ls())
options(stringsAsFactors = FALSE)

## Load libraries
library(rstan)
require(shinystan)
require(stringr)
library(dplyr)
library(plyr)

library(ggplot2)
library(reshape2)
library(viridis)
library(bayesplot)
library(tidybayes)
library(gridExtra) # for arranging plots 
library(patchwork) # another way of arranging plots 
library(rethinking)

## Set number of cores
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

pheno$site.n <- pheno$population
pheno$site.n[pheno$site.n == "sm"] <- "1"
pheno$site.n[pheno$site.n == "mp"] <- "2"
pheno$site.n[pheno$site.n == "HF"] <- "3"
pheno$site.n[pheno$site.n == "SH"] <- "4"
pheno$site.n <- as.numeric(pheno$site.n)

head(pheno)
#add dummy/ site level effects:
pheno <- pheno %>%
  mutate ( site2 = if_else(site.n == 2, 1, 0),
           site3 = if_else(site.n == 3, 1, 0),
           site4 = if_else(site.n == 4, 1, 0))

# standardize the 0/1 and standardize sites? 
pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)

pheno$site2.z2 <- (pheno$site2-mean(pheno$site2,na.rm=TRUE))/(sd(pheno$site2,na.rm=TRUE)*2)
pheno$site3.z2 <- (pheno$site3-mean(pheno$site3,na.rm=TRUE))/(sd(pheno$site3,na.rm=TRUE)*2)
pheno$site4.z2 <- (pheno$site4-mean(pheno$site4,na.rm=TRUE))/(sd(pheno$site4,na.rm=TRUE)*2)

#going to split it into analysis of terminal bb and lateral bb
# Starting with the terminal buds:
#pheno.term <- pheno[,c("tbb", "chill.n", "force.n", "photo.n", "site.n", "species", "lab2")]
pheno.term <- pheno[,c("bb", "force.z2", "photo.z2", "population", "species", "lab2","Utah_Model","Chill_portions","chillport.z2", "site2.z2", "site3.z2","site4.z2")]
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

##########################################################
# Height 

specieslist <- sort(unique(trtPheno$species))
sitelist <- sort(unique(trtPheno$site))
height <- trtPheno[complete.cases(trtPheno$ht),]

ht.data <- list(yTraiti = height$ht,
                 N = nrow(height),
                 n_spec = length(specieslist),
                 trait_species = as.numeric(as.factor(height$species)),
                 n_site = length(sitelist),
                 site = as.numeric(as.factor(height$site)),
                prior_mu_grand_mu = 20,
                prior_mu_grand_sigma = 10,
                prior_sigma_sp_mu = 4,
                prior_sigma_sp_sigma = 5,
                prior_sigma_site_mu = 2,
                prior_sigma_site_sigma = 5,
                prior_sigma_traity_mu = 3,
                prior_sigma_traity_sigma = 5,
                 ## Phenology
                 Nph = nrow(pheno.t),
                 phenology_species = as.numeric(as.factor(pheno.t$species)),
                 yPhenoi = pheno.t$bb,
                 forcei = pheno.t$force.z2,
                 chilli = pheno.t$chillport.z2,
                 photoi = pheno.t$photo.z2,
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


mdl.ht <- stan("stan/jointMdl.stan",
                      data = ht.data,
                      iter = 2000,
                      warmup = 1000,
                      chains = 4,
                      include = FALSE, pars = c("y_hat"))

save(mdl.ht, file = "output/ht_raw.Rda")
## N effective?
summary(mdl.ht)$summary[, "n_eff"] # 394.271, 13592.123

### Add species and study names to Stan object
names(mdl.ht)[grep(pattern = "^muSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^musite", x = names(mdl.ht))] <- paste(sitelist, sep = "")
##
names(mdl.ht)[grep(pattern = "^alphaForceSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^alphaChillSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^alphaPhotoSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^alphaPhenoSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^betaForceSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^betaChillSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^betaPhotoSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^betaPhenoSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")

pdf(file = "SLA_estimates_37spp.pdf", onefile = TRUE)
plot(mdl.ht, pars = c("mu_grand", "muSp"))
plot(mdl.ht, pars = c("musite"))
plot(mdl.ht, pars = c("muPhenoSp", "alphaPhenoSp"))
plot(mdl.ht, pars = c("muForceSp", "alphaForceSp"))
plot(mdl.ht, pars = c("muChillSp", "alphaChillSp"))
plot(mdl.ht, pars = c("muPhotoSp", "alphaPhotoSp"))
plot(mdl.ht, pars = c("betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto"))
plot(mdl.ht, pars = c("betaTraitxForce","betaForceSp"))
plot(mdl.ht, pars = c("betaTraitxChill","betaChillSp"))
plot(mdl.ht, pars = c("betaTraitxPhoto","betaPhotoSp"))
plot(mdl.ht, pars = c("sigma_traity", "sigma_site", "sigma_sp", "sigmaPhenoSp", "sigmapheno_y"))
dev.off()

saveRDS(object = mdl.ht, file = "height_stanfit.RDS")

##########################################################
# LMA

specieslist <- sort(unique(trtPheno$species))
sitelist <- sort(unique(trtPheno$site))
leafMass <- trtPheno[complete.cases(trtPheno$lma),]

lma.data <- list(yTraiti = leafMass$lma,
                N = nrow(leafMass),
                n_spec = length(specieslist),
                trait_species = as.numeric(as.factor(leafMass$species)),
                n_site = length(sitelist),
                site = as.numeric(as.factor(leafMass$site)),
                prior_mu_grand_mu = 0.5,
                prior_mu_grand_sigma = 5, #widened
                prior_sigma_sp_mu = 10,
                prior_sigma_sp_sigma = 5,
                prior_sigma_site_mu = 5,
                prior_sigma_site_sigma = 2,
                prior_sigma_traity_mu = 5,
                prior_sigma_traity_sigma = 2,
                ## Phenology
                Nph = nrow(pheno.t),
                phenology_species = as.numeric(as.factor(pheno.t$species)),
                yPhenoi = pheno.t$bb,
                forcei = pheno.t$force.z2,
                chilli = pheno.t$chillport.z2,
                photoi = pheno.t$photo.z2,
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


mdl.lma <- stan("stan/jointMdl.stan",
               data = lma.data,
               iter = 2000,
               warmup = 1000,
               chains = 4,
               include = FALSE, pars = c("y_hat"))

save(mdl.lma, file = "output/lma_raw.Rda")
## N effective?
summary(mdl.lma)$summary[, "n_eff"] # 394.271, 13592.123

### Add species and study names to Stan object
names(mdl.lma)[grep(pattern = "^muSp", x = names(mdl.lma))] <- paste(specieslist, sep = "")
names(mdl.lma)[grep(pattern = "^musite", x = names(mdl.lma))] <- paste(sitelist, sep = "")
##
names(mdl.lma)[grep(pattern = "^alphaForceSp", x = names(mdl.lma))] <- paste(specieslist, sep = "")
names(mdl.lma)[grep(pattern = "^alphaChillSp", x = names(mdl.lma))] <- paste(specieslist, sep = "")
names(mdl.lma)[grep(pattern = "^alphaPhotoSp", x = names(mdl.lma))] <- paste(specieslist, sep = "")
names(mdl.lma)[grep(pattern = "^alphaPhenoSp", x = names(mdl.lma))] <- paste(specieslist, sep = "")
names(mdl.lma)[grep(pattern = "^betaForceSp", x = names(mdl.lma))] <- paste(specieslist, sep = "")
names(mdl.lma)[grep(pattern = "^betaChillSp", x = names(mdl.lma))] <- paste(specieslist, sep = "")
names(mdl.lma)[grep(pattern = "^betaPhotoSp", x = names(mdl.lma))] <- paste(specieslist, sep = "")
names(mdl.lma)[grep(pattern = "^betaPhenoSp", x = names(mdl.lma))] <- paste(specieslist, sep = "")

pdf(file = "LMA_estimates.pdf", onefile = TRUE)
plot(mdl.lma, pars = c("mu_grand", "muSp"))
plot(mdl.lma, pars = c("musite"))
plot(mdl.lma, pars = c("muPhenoSp", "alphaPhenoSp"))
plot(mdl.lma, pars = c("muForceSp", "alphaForceSp"))
plot(mdl.lma, pars = c("muChillSp", "alphaChillSp"))
plot(mdl.lma, pars = c("muPhotoSp", "alphaPhotoSp"))
plot(mdl.lma, pars = c("betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto"))
plot(mdl.lma, pars = c("betaTraitxForce","betaForceSp"))
plot(mdl.lma, pars = c("betaTraitxChill","betaChillSp"))
plot(mdl.lma, pars = c("betaTraitxPhoto","betaPhotoSp"))
plot(mdl.lma, pars = c("sigma_traity", "sigma_site", "sigma_sp", "sigmaPhenoSp", "sigmapheno_y"))
dev.off()

saveRDS(object = mdl.lma, file = "lma_stanfit.RDS")

##########################################################
# ssd
hist(trtPheno$ssd)
stemDen <- trtPheno[complete.cases(trtPheno$ssd),]

ssd.data <- list(yTraiti = stemDen$ssd,
                 N = nrow(stemDen),
                 n_spec = length(specieslist),
                 trait_species = as.numeric(as.factor(stemDen$species)),
                 n_site = length(sitelist),
                 site = as.numeric(as.factor(stemDen$site)),
                 prior_mu_grand_mu = 1,
                 prior_mu_grand_sigma = 5, #widened
                 prior_sigma_sp_mu = 10,
                 prior_sigma_sp_sigma = 5,
                 prior_sigma_site_mu = 5,
                 prior_sigma_site_sigma = 2,
                 prior_sigma_traity_mu = 5,
                 prior_sigma_traity_sigma = 2,
                 ## Phenology
                 Nph = nrow(pheno.t),
                 phenology_species = as.numeric(as.factor(pheno.t$species)),
                 yPhenoi = pheno.t$bb,
                 forcei = pheno.t$force.z2,
                 chilli = pheno.t$chillport.z2,
                 photoi = pheno.t$photo.z2,
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


mdl.ssd <- stan("stan/jointMdl.stan",
                data = ssd.data,
                iter = 2000,
                warmup = 1000,
                chains = 4,
                include = FALSE, pars = c("y_hat"))

save(mdl.ssd, file = "output/ssd_raw.Rda")
## N effective?
summary(mdl.ssd)$summary[, "n_eff"] # 394.271, 13592.123

### Add species and study names to Stan object
names(mdl.ssd)[grep(pattern = "^muSp", x = names(mdl.ssd))] <- paste(specieslist, sep = "")
names(mdl.ssd)[grep(pattern = "^musite", x = names(mdl.ssd))] <- paste(sitelist, sep = "")
##
names(mdl.ssd)[grep(pattern = "^alphaForceSp", x = names(mdl.ssd))] <- paste(specieslist, sep = "")
names(mdl.ssd)[grep(pattern = "^alphaChillSp", x = names(mdl.ssd))] <- paste(specieslist, sep = "")
names(mdl.ssd)[grep(pattern = "^alphaPhotoSp", x = names(mdl.ssd))] <- paste(specieslist, sep = "")
names(mdl.ssd)[grep(pattern = "^alphaPhenoSp", x = names(mdl.ssd))] <- paste(specieslist, sep = "")
names(mdl.ssd)[grep(pattern = "^betaForceSp", x = names(mdl.ssd))] <- paste(specieslist, sep = "")
names(mdl.ssd)[grep(pattern = "^betaChillSp", x = names(mdl.ssd))] <- paste(specieslist, sep = "")
names(mdl.ssd)[grep(pattern = "^betaPhotoSp", x = names(mdl.ssd))] <- paste(specieslist, sep = "")
names(mdl.ssd)[grep(pattern = "^betaPhenoSp", x = names(mdl.ssd))] <- paste(specieslist, sep = "")

pdf(file = "SSD_estimates.pdf", onefile = TRUE)
plot(mdl.ssd, pars = c("mu_grand", "muSp"))
plot(mdl.ssd, pars = c("musite"))
plot(mdl.ssd, pars = c("muPhenoSp", "alphaPhenoSp"))
plot(mdl.ssd, pars = c("muForceSp", "alphaForceSp"))
plot(mdl.ssd, pars = c("muChillSp", "alphaChillSp"))
plot(mdl.ssd, pars = c("muPhotoSp", "alphaPhotoSp"))
plot(mdl.ssd, pars = c("betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto"))
plot(mdl.ssd, pars = c("betaTraitxForce","betaForceSp"))
plot(mdl.ssd, pars = c("betaTraitxChill","betaChillSp"))
plot(mdl.ssd, pars = c("betaTraitxPhoto","betaPhotoSp"))
plot(mdl.ssd, pars = c("sigma_traity", "sigma_site", "sigma_sp", "sigmaPhenoSp", "sigmapheno_y"))
dev.off()

saveRDS(object = mdl.ssd, file = "SSD_stanfit.RDS")

##########################################################
# dbh

diam <- trtPheno[complete.cases(trtPheno$dbh),]

dbh.data <- list(yTraiti = diam$dbh,
                 N = nrow(diam),
                 n_spec = length(specieslist),
                 trait_species = as.numeric(as.factor(diam$species)),
                 n_site = length(sitelist),
                 site = as.numeric(as.factor(diam$site)),
                 prior_mu_grand_mu = 15,
                 prior_mu_grand_sigma = 5, #widened
                 prior_sigma_sp_mu = 10,
                 prior_sigma_sp_sigma = 5,
                 prior_sigma_site_mu = 5,
                 prior_sigma_site_sigma = 2,
                 prior_sigma_traity_mu = 5,
                 prior_sigma_traity_sigma = 2,
                 ## Phenology
                 Nph = nrow(pheno.t),
                 phenology_species = as.numeric(as.factor(pheno.t$species)),
                 yPhenoi = pheno.t$bb,
                 forcei = pheno.t$force.z2,
                 chilli = pheno.t$chillport.z2,
                 photoi = pheno.t$photo.z2,
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


mdl.dbh <- stan("stan/jointMdl.stan",
                data = dbh.data,
                iter = 2000,
                warmup = 1000,
                chains = 4,
                include = FALSE, pars = c("y_hat"))
# run for more iterations

save(mdl.diam, file = "output/dbh_raw.Rda")
## N effective?
summary(mdl.dbh)$summary[, "n_eff"] # 394.271, 13592.123

### Add species and study names to Stan object
names(mdl.dbh)[grep(pattern = "^muSp", x = names(mdl.dbh))] <- paste(specieslist, sep = "")
names(mdl.dbh)[grep(pattern = "^musite", x = names(mdl.dbh))] <- paste(sitelist, sep = "")
##
names(mdl.dbh)[grep(pattern = "^alphaForceSp", x = names(mdl.dbh))] <- paste(specieslist, sep = "")
names(mdl.dbh)[grep(pattern = "^alphaChillSp", x = names(mdl.dbh))] <- paste(specieslist, sep = "")
names(mdl.dbh)[grep(pattern = "^alphaPhotoSp", x = names(mdl.dbh))] <- paste(specieslist, sep = "")
names(mdl.dbh)[grep(pattern = "^alphaPhenoSp", x = names(mdl.dbh))] <- paste(specieslist, sep = "")
names(mdl.dbh)[grep(pattern = "^betaForceSp", x = names(mdl.dbh))] <- paste(specieslist, sep = "")
names(mdl.dbh)[grep(pattern = "^betaChillSp", x = names(mdl.dbh))] <- paste(specieslist, sep = "")
names(mdl.dbh)[grep(pattern = "^betaPhotoSp", x = names(mdl.dbh))] <- paste(specieslist, sep = "")
names(mdl.dbh)[grep(pattern = "^betaPhenoSp", x = names(mdl.dbh))] <- paste(specieslist, sep = "")

pdf(file = "DBH_estimates.pdf", onefile = TRUE)
plot(mdl.dbh, pars = c("mu_grand", "muSp"))
plot(mdl.dbh, pars = c("musite"))
plot(mdl.dbh, pars = c("muPhenoSp", "alphaPhenoSp"))
plot(mdl.dbh, pars = c("muForceSp", "alphaForceSp"))
plot(mdl.dbh, pars = c("muChillSp", "alphaChillSp"))
plot(mdl.dbh, pars = c("muPhotoSp", "alphaPhotoSp"))
plot(mdl.dbh, pars = c("betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto"))
plot(mdl.dbh, pars = c("betaTraitxForce","betaForceSp"))
plot(mdl.dbh, pars = c("betaTraitxChill","betaChillSp"))
plot(mdl.dbh, pars = c("betaTraitxPhoto","betaPhotoSp"))
plot(mdl.dbh, pars = c("sigma_traity", "sigma_site", "sigma_sp", "sigmaPhenoSp", "sigmapheno_y"))
dev.off()

saveRDS(object = mdl.dbh, file = "dbh_stanfit.RDS")

##########################################################
# C:N

carbNit <- trtPheno[complete.cases(trtPheno$C.N),]

cn.data <- list(yTraiti = carbNit$C.N,
                 N = nrow(carbNit),
                 n_spec = length(specieslist),
                 trait_species = as.numeric(as.factor(carbNit$species)),
                 n_site = length(sitelist),
                 site = as.numeric(as.factor(carbNit$site)),
                 prior_mu_grand_mu = 20,
                 prior_mu_grand_sigma = 5, #widened
                 prior_sigma_sp_mu = 10,
                 prior_sigma_sp_sigma = 5,
                 prior_sigma_site_mu = 5,
                 prior_sigma_site_sigma = 2,
                 prior_sigma_traity_mu = 5,
                 prior_sigma_traity_sigma = 2,
                 ## Phenology
                 Nph = nrow(pheno.t),
                 phenology_species = as.numeric(as.factor(pheno.t$species)),
                 yPhenoi = pheno.t$bb,
                 forcei = pheno.t$force.z2,
                 chilli = pheno.t$chillport.z2,
                 photoi = pheno.t$photo.z2,
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


mdl.cn <- stan("stan/jointMdl.stan",
                data = cn.data,
                iter = 2000,
                warmup = 1000,
                chains = 4,
                include = FALSE, pars = c("y_hat"))
# bluk ESS is low - run for more iterations

save(mdl.cn, file = "output/cn_raw.Rda")
## N effective?
summary(mdl.cn)$summary[, "n_eff"] # 394.271, 13592.123

### Add species and study names to Stan object
names(mdl.cn)[grep(pattern = "^muSp", x = names(mdl.cn))] <- paste(specieslist, sep = "")
names(mdl.cn)[grep(pattern = "^musite", x = names(mdl.cn))] <- paste(sitelist, sep = "")
##
names(mdl.cn)[grep(pattern = "^alphaForceSp", x = names(mdl.cn))] <- paste(specieslist, sep = "")
names(mdl.cn)[grep(pattern = "^alphaChillSp", x = names(mdl.cn))] <- paste(specieslist, sep = "")
names(mdl.cn)[grep(pattern = "^alphaPhotoSp", x = names(mdl.cn))] <- paste(specieslist, sep = "")
names(mdl.cn)[grep(pattern = "^alphaPhenoSp", x = names(mdl.cn))] <- paste(specieslist, sep = "")
names(mdl.cn)[grep(pattern = "^betaForceSp", x = names(mdl.cn))] <- paste(specieslist, sep = "")
names(mdl.cn)[grep(pattern = "^betaChillSp", x = names(mdl.cn))] <- paste(specieslist, sep = "")
names(mdl.cn)[grep(pattern = "^betaPhotoSp", x = names(mdl.cn))] <- paste(specieslist, sep = "")
names(mdl.cn)[grep(pattern = "^betaPhenoSp", x = names(mdl.cn))] <- paste(specieslist, sep = "")

pdf(file = "CN_estimates.pdf", onefile = TRUE)
plot(mdl.cn, pars = c("mu_grand", "muSp"))
plot(mdl.cn, pars = c("musite"))
plot(mdl.cn, pars = c("muPhenoSp", "alphaPhenoSp"))
plot(mdl.cn, pars = c("muForceSp", "alphaForceSp"))
plot(mdl.cn, pars = c("muChillSp", "alphaChillSp"))
plot(mdl.cn, pars = c("muPhotoSp", "alphaPhotoSp"))
plot(mdl.cn, pars = c("betaTraitxForce", "betaTraitxChill", "betaTraitxPhoto"))
plot(mdl.cn, pars = c("betaTraitxForce","betaForceSp"))
plot(mdl.cn, pars = c("betaTraitxChill","betaChillSp"))
plot(mdl.cn, pars = c("betaTraitxPhoto","betaPhotoSp"))
plot(mdl.cn, pars = c("sigma_traity", "sigma_site", "sigma_sp", "sigmaPhenoSp", "sigmapheno_y"))
dev.off()

saveRDS(object = mdl.cn, file = "cn_stanfit.RDS")


# Plot model fit:
filePathData <- "output/"
traitModelNames <- grep("_stanfit.RDS", list.files(filePathData), value = TRUE) 

#Make a dataframe for saving traiit estimates for results section
traits <- c("Height","CN", "DBH", "SSD")
traitsDF <- data.frame(matrix(NA, 4,18))
names(traitsDF) <- c("Trait", "GrandMean", "GrandMean_upper", "GrandMean_lower", 
                     "SpeciesSigma",  "SpeciesSigma_upper", "SpeciesSigma_lower", 
                     "StudySigma",  "StudySigma_upper", "StudySigma_lower", 
                     "MaxValue",  "MaxValue_upper", "MaxValue_lower", "MaxValueSp", 
                     "MinValue", "MinValueSp", "MinValue_upper", "MinValue_lower")
traitsDF$Trait <- traits

traitPlotList <- list()

#for(traiti in 1:length(traitModelNames)){
  
  
  # 	traiti <- 3
  
  #Load SLA model fit
  traiti <- 1
  cnModel <- readRDS(paste(filePathData,traitModelNames[traiti], sep = "/"))
  traitName <- gsub("_stanfit.RDS", "", traitModelNames[traiti])
  cnModelFit <- rstan::extract(cnModel)
  
  #sensible cue values
  #-------------------------------------
  forcingValue <- 0.85 # 20 degrees C
  chillinValue <- 50 #coudl go up to 2 or 3 
  photoValue <- -0.25 # about 12 Or 0.5(about 16)
  
  #Extracting  postreior values 
  #----------------------------------
  
  #meanInterceptValues
  alphaPhenoSpdf <- data.frame(cnModelFit$alphaPhenoSp)
  alphaPhenoSpMean <- colMeans(alphaPhenoSpdf)
  
  #Forcing slope values 
  betaForceSpdf <- data.frame(cnModelFit$betaForceSp)
  betaForceSpMean <- colMeans(betaForceSpdf)
  
  #Chilling slope values 
  betaChillSpdf <- data.frame(cnModelFit$betaChillSp)
  betaChillSpMean <- colMeans(betaChillSpdf)
  
  
  #Photoperiod slope values 
  betaPhotoSpdf <- data.frame(cnModelFit$betaPhotoSp)
  betaPhotoSpMean <- colMeans(betaPhotoSpdf)
  
  #Overall model 
  sigmapheno_yMean <- mean(cnModelFit$sigmapheno_y)
  
  #Predict DOY based on model estimated parameters, One DOY value per species
  yPhenoi <- vector()
  
  for(ip in 1:length(betaPhotoSpMean)){
    yPhenoi[ip] <- alphaPhenoSpMean[ip] + betaForceSpMean[ip] * forcingValue + betaPhotoSpMean[ip] * photoValue + betaChillSpMean[ip]* chillinValue
  }
  
#  plot(yPhenoi ~ alphaPhenoSpMean)
  betaCombined <- betaPhotoSpMean+betaChillSpMean+betaForceSpMean
  
  mu_grandDf <- data.frame(cnModelFit$mu_grand_sp)
  colnames(mu_grandDf) <- specieslist
  longMeans <- melt(mu_grandDf)
  colnames(longMeans) <- c("species", "traitMean")
  
  mu_grand_mean <- colMeans(mu_grandDf)
  
  meanRealTrait <- aggregate(trtPheno$C.N, by = list(trtPheno$species), FUN = mean, na.rm =T)
  names(meanRealTrait) <- c("species","meanTrait")
  
  color_scheme_set("viridis")
  
  mcmc_intervals(mu_grandDf)+
    theme_classic() + 
    theme(text = element_text(size=20))+
    geom_point(data = trtPheno, aes(y = species, x = C.N), alpha = 0.5)
  
  traitFit <- ggplot(data = trtPheno, aes(y = species, x = C.N, colour = "black"))+
    stat_eye(data = longMeans, aes(y = species, x = traitMean))+
    geom_point( alpha = 0.5, size = 1.2, aes(colour = "red"))+
    theme_classic() +  
    theme(text = element_text(size=16))+
    geom_point(data = meanRealTrait, aes(x = meanTrait,y = species, colour = "purple"), shape = 8, size = 3)+
    labs(title = traitName, y = "Species", x ="Trait Value")+ 
    scale_color_identity(name = "Model fit",
                         breaks = c("black", "red", "purple"),
                         labels = c("Model Posterior", "Raw Data", "Data Mean"),
                         guide = guide_legend(override.aes = list(
                           linetype = c(NA, NA, NA),
                           shape = c(19, 20, 8)))) + 
    theme(legend.title = element_blank())
  
  