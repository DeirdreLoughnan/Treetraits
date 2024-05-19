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
#library(bayesplot)
#library(tidybayes)
library(gridExtra) # for arranging plots 
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

# pheno$site.n <- pheno$population
# pheno$site.n[pheno$site.n == "sm"] <- "1"
# pheno$site.n[pheno$site.n == "mp"] <- "2"
# pheno$site.n[pheno$site.n == "HF"] <- "3"
# pheno$site.n[pheno$site.n == "SH"] <- "4"
# pheno$site.n <- as.numeric(pheno$site.n)

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


##########################################################
# Height

specieslist <- sort(unique(trtPheno$species))
sitelist <- sort(unique(trtPheno$transect))
height <- trtPheno[complete.cases(trtPheno$ht),]

height <- height %>%
  mutate ( site2 = if_else(transect == 2, 1, 0)
  )

ht.data <- list(yTraiti = height$ht, 
                 N = nrow(height),
                 n_spec = length(specieslist),
                 trait_species = as.numeric(as.factor(height$species)),
                n_site = length(sitelist),
                site = as.numeric(as.factor(height$transect)),
                site2 = as.numeric((height$transect)),
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

#ht.data$site
mdl.ht <- stan("stan/jointMdl_Feb13.stan",
                      data = ht.data,
                      iter = 6000,
                      warmup = 3000,
                      chains = 4)

# mdl.htDum <- stan("stan/jointMdl_dummySite.stan",
#                data = ht.data,
#                iter = 6000,
#                warmup = 3000,
#                chains = 4)
                     # include = FALSE, pars = c("y_hat"))

save(mdl.ht, file = "output/ht_Feb13.Rda")
sumer <- summary(mdl.ht)$summary
bsite <- data.frame(sumer[grep("b_site", rownames(sumer)), c("mean","2.5%", "97.5%", "n_eff", "Rhat")])
## N effective?
summary(mdl.ht)$summary[, "n_eff"] # 394.271, 13592.123

### Add species and study names to Stan object
names(mdl.ht)[grep(pattern = "^muSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^b_site", x = names(mdl.ht))] <- paste(sitelist, sep = "")
##
names(mdl.ht)[grep(pattern = "^alphaForceSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^alphaChillSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^alphaPhotoSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^alphaPhenoSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^betaForceSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^betaChillSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^betaPhotoSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")
names(mdl.ht)[grep(pattern = "^betaPhenoSp", x = names(mdl.ht))] <- paste(specieslist, sep = "")

pdf(file = "ht_estimates_37spp.pdf", onefile = TRUE)
plot(mdl.ht, pars = c("mu_grand", "muSp"))
plot(mdl.ht, pars = c("b_site"))
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
sitelist <- sort(unique(trtPheno$transect))
leafMass <- trtPheno[complete.cases(trtPheno$lma),]

lma.data <- list(yTraiti = leafMass$lma,
                N = nrow(leafMass),
                n_spec = length(specieslist),
                trait_species = as.numeric(as.factor(leafMass$species)),
                n_site = length(sitelist),
                site = as.numeric(as.factor(leafMass$transect)),
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
               iter = 4000,
               warmup = 3000,
               chains = 4)
               # include = FALSE, pars = c("y_hat"), 
               # control = list(max_treedepth = 11))

save(mdl.lma, file = "output/lma_raw.Rda")
## N effective?
range(summary(mdl.lma)$summary[, "n_eff"] )# 394.271, 13592.123

### Add species and study names to Stan object
names(mdl.lma)[grep(pattern = "^muSp", x = names(mdl.lma))] <- paste(specieslist, sep = "")
names(mdl.lma)[grep(pattern = "^b_site", x = names(mdl.lma))] <- paste(sitelist, sep = "")
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
plot(mdl.lma, pars = c("b_site"))
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
                 site = as.numeric(as.factor(stemDen$transect)),
                 prior_mu_grand_mu = 1,
                 prior_mu_grand_sigma = 5, #widened
                 prior_sigma_sp_mu = 4, #10
                 prior_sigma_sp_sigma = 5,
                 prior_sigma_site_mu = 2, #5
                 prior_sigma_site_sigma = 5, #2
                 prior_sigma_traity_mu = 2, #5
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


mdl.ssd <- stan("stan/jointMdl.stan",
                data = ssd.data,
                iter = 4000,
                warmup = 2000,
                chains = 4)
                #include = FALSE, pars = c("y_hat"))

save(mdl.ssd, file = "output/ssd_raw_np.Rda")
## N effective?
summary(mdl.ssd)$summary[, "n_eff"] # 1013.934 8423.027

### Add species and study names to Stan object
names(mdl.ssd)[grep(pattern = "^muSp", x = names(mdl.ssd))] <- paste(specieslist, sep = "")
names(mdl.ssd)[grep(pattern = "^b_site", x = names(mdl.ssd))] <- paste(sitelist, sep = "")
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
plot(mdl.ssd, pars = c("b_site"))
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

saveRDS(object = mdl.ssd, file = "SSD_stanfit_np.RDS")

# ##########################################################
# # dbh
# 
diam <- trtPheno[complete.cases(trtPheno$dbh),]

dbh.data <- list(yTraiti = diam$dbh,
                 N = nrow(diam),
                 n_spec = length(specieslist),
                 trait_species = as.numeric(as.factor(diam$species)),
                 n_site = length(sitelist),
                 site = as.numeric(as.factor(diam$transect)),
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
                iter = 6000,
                warmup = 3000,
                chains = 4)
               # include = FALSE, pars = c("y_hat"))
# run for more iterations

save(mdl.dbh, file = "output/dbh_raw.Rda")
## N effective?
range(summary(mdl.dbh)$summary[, "n_eff"]) # 394.271, 13592.123

### Add species and study names to Stan object
names(mdl.dbh)[grep(pattern = "^muSp", x = names(mdl.dbh))] <- paste(specieslist, sep = "")
names(mdl.dbh)[grep(pattern = "^b_site", x = names(mdl.dbh))] <- paste(sitelist, sep = "")
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
plot(mdl.dbh, pars = c("b_site"))
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

# ##########################################################
# # C:N

carbNit <- trtPheno[complete.cases(trtPheno$C.N),]

cn.data <- list(yTraiti = carbNit$C.N,
                 N = nrow(carbNit),
                 n_spec = length(specieslist),
                 trait_species = as.numeric(as.factor(carbNit$species)),
                 n_site = length(sitelist),
                 site = as.numeric(as.factor(carbNit$transect)),
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


mdl.cn <- stan("stan/jointMdl.stan",
                data = cn.data,
                iter = 6000,
                warmup = 3000,
                chains = 4
                # include = FALSE, pars = c("y_hat")
            )
# bluk ESS is low - run for more iterations

save(mdl.cn, file = "output/cn_raw_np.Rda")
## N effective?
range(summary(mdl.cn)$summary[, "n_eff"]) #

### Add species and study names to Stan object
names(mdl.cn)[grep(pattern = "^muSp", x = names(mdl.cn))] <- paste(specieslist, sep = "")
names(mdl.cn)[grep(pattern = "^b_site", x = names(mdl.cn))] <- paste(sitelist, sep = "")
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
plot(mdl.cn, pars = c("b_site"))
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

saveRDS(object = mdl.cn, file = "cn_stanfit_np.RDS")

# plot ypred plots:
y <- carbNit$C.N

postCN <- rstan::extract(mdl.cn)
yrep <-  postCN$y_hat

pdf("figures/cn_ypred_trt.pdf")
ppc_dens_overlay(y, yrep[1:100,])
dev.off()

#### Running the model without site for each of our sites ##########
# height:
# 1. Smithers:
smHt <- subset(height, site == "sm")
smPheno <- subset(pheno.t, population == "sm")

smHt.data <- list(yTraiti = smHt$ht, 
                N = nrow(smHt),
                n_spec = length(unique(smHt$species)),
                trait_species = as.numeric(as.factor(smHt$species)),
                prior_mu_grand_mu = 20,
                prior_mu_grand_sigma = 10,
                prior_sigma_sp_mu = 4,
                prior_sigma_sp_sigma = 5,
                prior_sigma_traity_mu = 3,
                prior_sigma_traity_sigma = 5,
                ## Phenology
                Nph = nrow(smPheno),
                phenology_species = as.numeric(as.factor(smPheno$species)),
                yPhenoi = smPheno$bb,
                forcei = smPheno$force.z2,
                chilli = smPheno$chillport.z2,
                photoi = smPheno$photo.z2,
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
mdl.smHt <- stan("stan/jointMdl_noSite.stan",
               data = smHt.data,
               iter = 6000,
               warmup = 3000,
               chains = 4)

# 1. Manning park:
smHt <- subset(height, site == "sm")
smHt <- subset(smHt, species != "samrac")
smPheno <- subset(pheno.t, population == "sm")

smHt.data <- list(yTraiti = smHt$ht, 
                  N = nrow(smHt),
                  n_spec = length(unique(smHt$species)),
                  trait_species = as.numeric(as.factor(smHt$species)),
                  prior_mu_grand_mu = 20,
                  prior_mu_grand_sigma = 10,
                  prior_sigma_sp_mu = 4,
                  prior_sigma_sp_sigma = 5,
                  prior_sigma_traity_mu = 3,
                  prior_sigma_traity_sigma = 5,
                  ## Phenology
                  Nph = nrow(smPheno),
                  phenology_species = as.numeric(as.factor(smPheno$species)),
                  yPhenoi = smPheno$bb,
                  forcei = smPheno$force.z2,
                  chilli = smPheno$chillport.z2,
                  photoi = smPheno$photo.z2,
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
                  prior_sigmaPhotoSp_mu = 0,
                  prior_sigmaPhotoSp_sigma = 10,
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

smHt.data$trait_species
mdl.smHt <- stan("stan/jointMdl_noSite.stan",
                 data = smHt.data,
                 iter = 6000,
                 warmup = 3000,
                 chains = 4)

smHt.data <- list(yTraiti = smHt$ht, 
                  N = nrow(smHt),
                  n_spec = length(unique(smHt$species)),
                  species = as.numeric(as.factor(smHt$species)),
                  prior_mu_grand = 20,
                  prior_sigma_grand = 10,
                  prior_sigma_sp_mu = 4,
                  prior_sigma_sp_sigma = 5,
                  prior_sigma_traity_mu = 3,
                  prior_sigma_traity_sigma = 5,
                  ## Phenology
                  Nph = nrow(smPheno),
                  phenology_species = as.numeric(as.factor(smPheno$species)),
                  yPhenoi = smPheno$bb,
                  forcei = smPheno$force.z2,
                  chilli = smPheno$chillport.z2,
                  photoi = smPheno$photo.z2,
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
                  prior_sigmaPhotoSp_mu = 0,
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


mdl.smHt <- stan("stan/trait_only_nosite.stan",
                 data = smHt.data,
                 iter = 6000,
                 warmup = 3000,
                 chains = 4)



mdl.smHt <- stan("stan/jointMdl_Feb13_NoPhoto.stan",
                 data = smHt.data,
                 iter = 6000,
                 warmup = 3000,
                 chains = 4)

mdl.smHt <- stan("stan/jointMdl_Feb13_OnlyPhoto.stan",
                 data = smHt.data,
                 iter = 6000,
                 warmup = 3000,
                 chains = 4)

# mdl.htDum <- stan("stan/jointMdl_dummySite.stan",
#                data = ht.data,
#                iter = 6000,
#                warmup = 3000,
#                chains = 4)
# include = FALSE, pars = c("y_hat"))

save(mdl.smHt, file = "output/ht_smNoPhoto.Rda")
sumer <- summary(mdl.smHt)$summary

pairs(mdl.smHt, pars = c("mu_grand", 
                        "muPhenoSp",
                        "muPhotoSp",
                        "betaTraitxPhoto","sigma_sp",
                         "sigma_traity",
                         "sigmaPhenoSp",
                         "sigmaPhotoSp", #sigma_alpha_photo
                         "sigmapheno_y", "lp__")) 

post <- rstan::extract(mdl.smHt)
par(mfrow = c(2,2))
plot(post$muPhotoSp ~ post$sigmaPhotoSp)
plot(post$betaTraitxPhoto ~ post$sigmaPhotoSp)
plot(post$sigma_traity ~post$sigmaPhotoSp)
plot(post$sigmaPhenoSp ~ post$sigmaPhotoSp)
