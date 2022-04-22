# Started April 22, 2022 by Deirdre

# The purpose of this code is to get the traitors joint model working for my trait phenology data! 

rm(list=ls())
options(stringsAsFactors = FALSE)

## Load libraries
library(rstan)
require(shinystan)

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


