# started Feb 18, 2022
# the purpose of this code is to build test data for my be_trait chapter

# my general idea is to develop a model that builds on the traitors joint model
# except with multiple traits modeled together to get a syndrome effect not just the effect of individual traits
# have a dummy variable for transect perhaps?

# Feb 18:to get started I am just trying to get good test data for the trait, building off of the traitors test data from: trait_pheno_only_models.R

# Started by DL on June 10, 2021.

# There are some persistent issues with the joint model, so we are taking two steps back and going to get the trait only and pheno only model working. Lizzie suggested the following: 
#1. Make a nice, compact version of your R file that runs just your traits model. Make sure you can get nice results back that match your parameters and that the chains look well mixed. Ideally this should not take >4,000 iterations.
#2. Assuming that works, make up a forcing-only phenology model and check via lmer, rstan or such that your code for phenology alone is correct (maybe make a separate R file for this also).
#3. Merge the two models only once 1 and 2 are done.
#4.If all goes well we could add in chilling and photoperiod. To do this, go back and do step 2 (adding chill and photo) then merge.

if(length(grep("deirdreloughnan", getwd()) > 0)) {
  setwd("~/Documents/github/Treetraits")
} else{
  setwd("/home/deirdre/Treetraits")
}

library(rstan)
require(rstanarm)
require(shinystan)
require(bayesplot)
require(truncnorm)
library(ggplot2)

# 50 study, spp, reps, drop sigmatrait_y to 1
# boxplots of all the dots by species and study
#compare sp to study and make sure no correlations
# run test data through rstan arm - how fast? accurate? 
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

Nrep <- 20 # rep per trait
Ntransect <- 2 # number of studies w/ traits (10 seems a little low for early simulation code; remember that you are estimating a distribution of this the same as for species)
Nspp <- 50 # number of species with traits (making this 20 just for speed for now)

# First making a data frame for the test trait data
Ntrt <- Nspp * Ntransect * Nrep # total number of traits observations
Ntrt

#make a dataframe for height
trt.dat <- data.frame(matrix(NA, Ntrt, 1))
names(trt.dat) <- c("rep")
trt.dat$rep <- c(1:Nrep)
trt.dat$species <- rep(1:Nspp)
trt.dat$transect <- rep(c(1:Ntransect), each = Nspp)
# trt.dat$species <- rep(1:Nspp, Nstudy)

# now generating the species trait data, here it is for height
mu.grand <- 10 # the grand mean of the height model
sigma.species <- 1 # we want to keep the variaiton across spp. high

#the alphaTraitSp in Faiths original code:
mu.trtsp <- rnorm(Nspp, 0, sigma.species)
trt.dat$mu.trtsp <- rep(mu.trtsp, Ntransect) #adding ht data for ea. sp

#now generating the effects of study
sigma.transect <- 1
mu.transect <- rnorm(Ntransect, 0, sigma.transect) #intercept for each study
trt.dat$mu.transect <- rep(mu.transect, each = Nspp) # generate data for ea study

# general variance
trt.var <- 1 #sigmaTrait_y in the stan code
trt.dat$trt.er <- rnorm(Ntrt, 0, trt.var)

# generate yhat - heights -  for this first trt model
trt.dat$yTraiti <- mu.grand + trt.dat$mu.trtsp + trt.dat$trt.er
# trt.dat$yTraiti <- mu.grand + trt.dat$mu.trtsp + trt.dat$mu.study + trt.dat$trt.er

## Exploring the test data - boxplots!  ###########################################
# names(trt.dat)
# boxplot(yTraiti ~ study, data = trt.dat)
# boxplot(yTraiti ~ species, data = trt.dat)
# 
# plot(yTraiti ~ study, data = trt.dat)
# # What about running the model with rstanarm?  ####################################
# # library(rstanarm)
# mdl.arm <- stan_lmer(yTraiti ~ mu.trtsp + mu.study + (1 | study) + (1 + |species),
#                      data = trt.dat)
# prior_summary(object = mdl.arm)
# summary(mdl.arm, probs = c(0.025, 0.975),
#         digits = 2)
# mdl.arm
# Stop here and test your work a little ...okay, it's hard to interpret this output but we can check the general variance and intercept and check the relative variance of species and study (maybe using the SD?)

## Trait only stan model ###########################################################
trait_data <- list(yTraiti = trt.dat$yTraiti, 
                   N = Ntrt, 
                   n_spec = Nspp, 
                   trait_species = trt.dat$species, 
                   # study = trt.dat$study, 
                   # n_study = Nstudy,
                   prior_mu_grand_mu = 10,
                   prior_mu_grand_sigma = 1,
                   #prior_mu_sp = 0,
                   prior_sigma_sp_mu = 5,
                   prior_sigma_sp_sigma = 5,
                   # prior_mu_study = 0,
                   # prior_sigma_study_mu = 1,
                   # prior_sigma_study_sigma = 0.5,
                   prior_sigma_traity_mu = 5,
                   prior_sigma_traity_sigma = 2.5
                   ) 

mdl.trait <- stan('stan/bc_trait_only_2.stan',
                  data = trait_data,
                  iter = 6000,
                  warmup = 3000,
                  chains = 4,
                  include = FALSE,
                  pars = "mu_y")

save(mdl.trait, file = "output/output_trait_testdata.Rda")
# bulk ess is low

ssm <-  as.shinystan(mdl.trait)
launch_shinystan(ssm)

sumer <- summary(mdl.trait)$summary
post <- rstan::extract(mdl.trait)

range(sumer[, "n_eff"])
range(sumer[, "Rhat"])

mu_grand <- sumer[grep("mu_grand", rownames(sumer))]; mu_grand
sigma_sp <- sumer[grep("sigma_sp", rownames(sumer))]; sigma_sp
sigma_studyesti <- sumer[grep("sigma_study", rownames(sumer))]; 
sigmaTrait_y <- sumer[grep("sigmaTrait_y", rownames(sumer))];sigmaTrait_y 

