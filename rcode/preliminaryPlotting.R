# Started Feb 4, 2022

# The purpose of this code is to plot the trait values for each species against the the slope estimates from the phenology model:

# Get the model output from pheno_bc repo:
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(rstan)
library(shinystan)
#library(reshape2)
library(bayesplot)
library(ggplot2)
#library(RColorBrewer)
library(dplyr)
library(plyr)


if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  else{
  setwd("/home/deirdre/pheno_bc") # for midge
}


load("output/realDummSite.Rda")

setwd("..//Treetraits")
trtData <- read.csv("data/allTrt.csv")
head(trtData)


spp <- c("acegla", "acepen", "acerub", "acesac", "alninc","alnvir", "amealn", "aromel", "betall", "betlen", "betpap",
         "corcor", "corsto", "faggra","franig", "hamvir", "ilemuc", "kalang", "loncan", "loninv", "lyolig", "menfer", "nyssyl", "popbal","popgra", "poptre", "prupen", "quealb", "querub", "quevel", "rhafra", "rhoalb", "rhopri", "riblac", "rubpar", "samrac", "shecan","sorsco", "spialb", "spibet", "spipyr", "symalb", "vacmem", "vacmyr", "vibcas", "vibedu", "viblan")

trtDataSpp <- trtData[trtData$species %in% spp,]
trtDataSpp$species.fact <- as.numeric(as.factor(trtDataSpp$species))

head(trtData)

sumt <- summary(mdl)$summary
bforce <- sumt[grep("b_force", rownames(sumt)), "mean"]; bforce
# bchill <- sumt[grep("b_chill", rownames(sumt)), "mean"]; bchill
# bphoto <- sumt[grep("b_photo", rownames(sumt)), "mean"]; bphoto

for (i in c(1:length(spp.fact))){
  trt.dataset$force[i] <- sumt[grep("b_force", rownames(sumt)),1][spp.fact[i]]
  trt.dataset$chill[i] <- sumt[grep("b_chill\\[", rownames(sumt)),1][spp.fact[i]]
  trt.dataset$photo[i] <- sumt[grep("b_photo\\[", rownames(sumt)),1][spp.fact[i]]

}
##################################################
# now run simple trait model to get trait effects for each species:

lma_data  <- trtDataSpp[complete.cases(trtDataSpp$lma),]
lma_data$species.fact <- as.numeric(as.factor(lma_data$species))

lma_datalist <- list(yTraiti = lma_data$lma, 
                   N = nrow(lma_data), 
                   n_spec = length(unique(lma_data$species)), 
                   species = lma_data$species.fact, 
                   prior_mu_grand_mu = 0.5,
                   prior_mu_grand_sigma = 1,
                   prior_sigma_sp_mu = 1,
                   prior_sigma_sp_sigma = 1,
                   prior_sigma_traity_mu = 1,
                   prior_sigma_traity_sigma = 1
) 


mdl.lma <- stan('stan/bc_trait_only_2.stan',
                  data = lma_datalist,
                  iter = 6000,
                  warmup = 4000,
                  chains = 4,
                  include = FALSE,
                  pars = "mu_y")
save(mdl.lma, file = "output_lma_traitonly.Rda")

sum.lma <- summary(mdl.lma)$summary

########################################################
ht_data  <- trtDataSpp[complete.cases(trtDataSpp$ht),]
ht_datalist <- list(yTraiti = ht_data$ht, 
                     N = nrow(ht_data), 
                     n_spec = length(unique(ht_data$species)), 
                     species = as.numeric(as.factor(ht_data$species)), 
                     prior_mu_grand_mu = 15,
                     prior_mu_grand_sigma = 10,
                     prior_sigma_sp_mu = 5,
                     prior_sigma_sp_sigma = 1,
                     prior_sigma_traity_mu = 5,
                     prior_sigma_traity_sigma = 1
) 

mdl.ht <- stan('stan/bc_trait_only_2.stan',
                data = ht_datalist,
                iter = 6000,
                warmup = 4000,
                chains = 4,
                include = FALSE,
                pars = "mu_y")
save(mdl.ht, file = "output_ht_traitonly.Rda")

sum.ht <- summary(mdl.ht)$summary

######################################################
dbh_data  <- trtDataSpp[complete.cases(trtDataSpp$dbh),]
dbh_data$species.fact <- as.numeric(as.factor(dbh_data$species))

dbh_datalist <- list(yTraiti = dbh_data$dbh, 
                     N = nrow(dbh_data), 
                     n_spec = length(unique(dbh_data$species)), 
                     species = dbh_data$species.fact, 
                     prior_mu_grand_mu = 0.5,
                     prior_mu_grand_sigma = 1,
                     prior_sigma_sp_mu = 1,
                     prior_sigma_sp_sigma = 1,
                     prior_sigma_traity_mu = 1,
                     prior_sigma_traity_sigma = 1
) 


mdl.dbh <- stan('stan/bc_trait_only_2.stan',
                data = dbh_datalist,
                iter = 6000,
                warmup = 4000,
                chains = 4,
                include = FALSE,
                pars = "mu_y")
save(mdl.dbh, file = "output_dbh_traitonly.Rda")

sum.dbh <- summary(mdl.dbh)$summary

########################################################
ssd_data  <- trtDataSpp[complete.cases(trtDataSpp$ssd),]
ssd_data$species.fact <- as.numeric(as.factor(ssd_data$species))

ssd_datalist <- list(yTraiti = ssd_data$ssd, 
                     N = nrow(ssd_data), 
                     n_spec = length(unique(ssd_data$species)), 
                     species = ssd_data$species.fact, 
                     prior_mu_grand_mu = 0.5,
                     prior_mu_grand_sigma = 1,
                     prior_sigma_sp_mu = 1,
                     prior_sigma_sp_sigma = 1,
                     prior_sigma_traity_mu = 1,
                     prior_sigma_traity_sigma = 1
) 


mdl.ssd <- stan('stan/bc_trait_only_2.stan',
                data = ssd_datalist,
                iter = 6000,
                warmup = 4000,
                chains = 4,
                include = FALSE,
                pars = "mu_y")
save(mdl.ssd, file = "output_ssd_traitonly.Rda")

sum.ssd <- summary(mdl.ssd)$summary

########################################################
perN_data  <- trtDataSpp[complete.cases(trtDataSpp$perN),]
perN_data$species.fact <- as.numeric(as.factor(perN_data$species))

perN_datalist <- list(yTraiti = perN_data$perN, 
                     N = nrow(perN_data), 
                     n_spec = length(unique(perN_data$species)), 
                     species = perN_data$species.fact, 
                     prior_mu_grand_mu = 0.5,
                     prior_mu_grand_sigma = 1,
                     prior_sigma_sp_mu = 1,
                     prior_sigma_sp_sigma = 1,
                     prior_sigma_traity_mu = 1,
                     prior_sigma_traity_sigma = 1
) 


mdl.perN <- stan('stan/bc_trait_only_2.stan',
                data = perN_datalist,
                iter = 6000,
                warmup = 4000,
                chains = 4,
                include = FALSE,
                pars = "mu_y")
save(mdl.perN, file = "output_perN_traitonly.Rda")

sum.perN <- summary(mdl.perN)$summary

########################################################
perC_data  <- trtDataSpp[complete.cases(trtDataSpp$perC),]
perC_data$species.fact <- as.numeric(as.factor(perC_data$species))

perC_datalist <- list(yTraiti = perC_data$perC, 
                     N = nrow(perC_data), 
                     n_spec = length(unique(perC_data$species)), 
                     species = perC_data$species.fact, 
                     prior_mu_grand_mu = 0.5,
                     prior_mu_grand_sigma = 1,
                     prior_sigma_sp_mu = 1,
                     prior_sigma_sp_sigma = 1,
                     prior_sigma_traity_mu = 1,
                     prior_sigma_traity_sigma = 1
) 


mdl.perC <- stan('stan/bc_trait_only_2.stan',
                data = perC_datalist,
                iter = 6000,
                warmup = 4000,
                chains = 4,
                include = FALSE,
                pars = "mu_y")
save(mdl.perC, file = "output_perC_traitonly.Rda")

sum.perC <- summary(mdl.perC)$summary

########################################################
plot(force ~ ssd , data = trt.dataset)
plot(chill ~ ssd , data = trt.dataset)
plot(photo ~ ssd , data = trt.dataset)


plot(force ~ per.C , data = trt.dataset)
plot(chill ~ per.C , data = trt.dataset)
plot(photo ~ per.C , data = trt.dataset)

plot(force ~ per.N , data = trt.dataset)
plot(chill ~ per.N , data = trt.dataset)
plot(photo ~ per.N , data = trt.dataset)

plot(force ~ ht , data = trt.dataset)
plot(chill ~ ht , data = trt.dataset)
plot(photo ~ ht , data = trt.dataset)

plot(force ~ lma , data = trt.dataset)
plot(chill ~ lma , data = trt.dataset)
plot(photo ~ lma , data = trt.dataset)

plot(force ~ dbh , data = trt.dataset)
plot(chill ~ dbh , data = trt.dataset)
plot(photo ~ dbh , data = trt.dataset)


