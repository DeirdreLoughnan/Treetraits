rm(list=ls())
options(stringsAsFactors = FALSE)

require(dplyr)
library(stringr)
library(plyr)
library(rstan)

# Set working directory:
# Anyone else working with this code should add their info/path here
if(length(grep("deirdreloughnan", getwd())>0)) {  setwd("~/Documents/github/Treetraits")
} else if
(length(grep("Lizzie", getwd())>0)) {   setwd("~/Documents/git/projects/treegarden/budreview/ospree/analyses/traits")
}

# height:
posterior_ht <- rstan::extract(readRDS(file = "output/height_stanfit.RDS"))

pdf("figures/height_prior_post_dist_traitormdl.pdf", width = 15, height = 25)
par(mfrow = c(4,4))
#plot priors against posterior_hts
h1 <- hist(rnorm(1000, -15,10), col = rgb(1,0,1,1/4), main = "muForceSp")
hist(posterior_ht$muForceSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, -15,10), col = rgb(1,0,1,1/4), main = "muChillSp")
hist(posterior_ht$muChillSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, -15,10), col = rgb(1,0,1,1/4), main = "muPhotoSp")
hist(posterior_ht$muPhotoSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 20,10), col = rgb(1,0,1,1/4), main = "mu_grand")
hist(posterior_ht$mu_grand,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 40,10), col = rgb(1,0,1,1/4), main = "muPhenoSp")
hist(posterior_ht$muPhenoSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 0,1), col = rgb(1,0,1,1/4), main = "betaTraitxForce")
hist(posterior_ht$betaTraitxForce,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 0,1), col = rgb(1,0,1,1/4), main = "betaTraitxChill")
hist(posterior_ht$betaTraitxChill,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 0,1), col = rgb(1,0,1,1/4), main = "betaTraitxPhoto")
hist(posterior_ht$betaTraitxPhoto,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 4,5), col = rgb(1,0,1,1/4), main = "sigma_sp")
hist(posterior_ht$sigma_sp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 2,5), col = rgb(1,0,1,1/4), main = "sigma_study")
hist(posterior_ht$sigma_site,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 3,5), col = rgb(1,0,1,1/4), main = "sigma_traity")
hist(posterior_ht$sigma_traity,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaForceSp")
hist(posterior_ht$sigmaForceSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaChillSp")
hist(posterior_ht$sigmaChillSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaPhotoSp")
hist(posterior_ht$sigmaPhotoSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaPhenoSp")
hist(posterior_ht$sigmaPhenoSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 10,5), col = rgb(1,0,1,1/4), main = "sigmapheno_y")
hist(posterior_ht$sigmapheno_y,add=T,col=rgb(0,0,1,1/4))

dev.off()

#######################################################################
# lma:
posterior_lma <- rstan::extract(readRDS(file = "output/height_stanfit.RDS"))

pdf("figures/height_prior_post_dist_traitormdl.pdf", width = 15, height = 25)
par(mfrow = c(4,4))
#plot priors against posterior_lmas
h1 <- hist(rnorm(1000, -15,10), col = rgb(1,0,1,1/4), main = "muForceSp")
hist(posterior_lma$muForceSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, -15,10), col = rgb(1,0,1,1/4), main = "muChillSp")
hist(posterior_lma$muChillSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, -15,10), col = rgb(1,0,1,1/4), main = "muPhotoSp")
hist(posterior_lma$muPhotoSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 20,10), col = rgb(1,0,1,1/4), main = "mu_grand")
hist(posterior_lma$mu_grand,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 40,10), col = rgb(1,0,1,1/4), main = "muPhenoSp")
hist(posterior_lma$muPhenoSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 0,1), col = rgb(1,0,1,1/4), main = "betaTraitxForce")
hist(posterior_lma$betaTraitxForce,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 0,1), col = rgb(1,0,1,1/4), main = "betaTraitxChill")
hist(posterior_lma$betaTraitxChill,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 0,1), col = rgb(1,0,1,1/4), main = "betaTraitxPhoto")
hist(posterior_lma$betaTraitxPhoto,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 4,5), col = rgb(1,0,1,1/4), main = "sigma_sp")
hist(posterior_lma$sigma_sp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 2,5), col = rgb(1,0,1,1/4), main = "sigma_study")
hist(posterior_lma$sigma_site,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 3,5), col = rgb(1,0,1,1/4), main = "sigma_traity")
hist(posterior_lma$sigma_traity,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaForceSp")
hist(posterior_lma$sigmaForceSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaChillSp")
hist(posterior_lma$sigmaChillSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaPhotoSp")
hist(posterior_lma$sigmaPhotoSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaPhenoSp")
hist(posterior_lma$sigmaPhenoSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 10,5), col = rgb(1,0,1,1/4), main = "sigmapheno_y")
hist(posterior_lma$sigmapheno_y,add=T,col=rgb(0,0,1,1/4))

dev.off()

#######################################################################
# ssd:
posterior_ssd <- rstan::extract(readRDS(file = "output/SSD_stanfit_np.RDS"))

pdf("figures/ssd_prior_post_dist_traitormdl.pdf", width = 15, height = 25)
par(mfrow = c(4,4))
#plot priors against posterior_ssds
h1 <- hist(rnorm(1000, -15,10), col = rgb(1,0,1,1/4), main = "muForceSp")
hist(posterior_ssd$muForceSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, -15,10), col = rgb(1,0,1,1/4), main = "muChillSp")
hist(posterior_ssd$muChillSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, -15,10), col = rgb(1,0,1,1/4), main = "muPhotoSp")
hist(posterior_ssd$muPhotoSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 1,5), col = rgb(1,0,1,1/4), main = "mu_grand")
hist(posterior_ssd$mu_grand,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 40,10), col = rgb(1,0,1,1/4), main = "muPhenoSp")
hist(posterior_ssd$muPhenoSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 0,2), col = rgb(1,0,1,1/4), main = "betaTraitxForce")
hist(posterior_ssd$betaTraitxForce,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 0,2), col = rgb(1,0,1,1/4), main = "betaTraitxChill")
hist(posterior_ssd$betaTraitxChill,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 0,2), col = rgb(1,0,1,1/4), main = "betaTraitxPhoto")
hist(posterior_ssd$betaTraitxPhoto,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 4,5), col = rgb(1,0,1,1/4), main = "sigma_sp")
hist(posterior_ssd$sigma_sp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 2,5), col = rgb(1,0,1,1/4), main = "sigma_study")
hist(posterior_ssd$sigma_site,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 2,5), col = rgb(1,0,1,1/4), main = "sigma_traity")
hist(posterior_ssd$sigma_traity,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaForceSp")
hist(posterior_ssd$sigmaForceSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaChillSp")
hist(posterior_ssd$sigmaChillSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaPhotoSp")
hist(posterior_ssd$sigmaPhotoSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaPhenoSp")
hist(posterior_ssd$sigmaPhenoSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 10,5), col = rgb(1,0,1,1/4), main = "sigmapheno_y")
hist(posterior_ssd$sigmapheno_y,add=T,col=rgb(0,0,1,1/4))

dev.off()


########################################################################
# dbh:
posterior_dbh <- rstan::extract(readRDS(file = "output/dbh_stanfit.RDS"))

pdf("figures/dbh_prior_post_dist_traitormdl.pdf", width = 15, height = 25)
par(mfrow = c(4,4))
#plot priors against posterior_dbhs
h1 <- hist(rnorm(1000, -15,10), col = rgb(1,0,1,1/4), main = "muForceSp")
hist(posterior_dbh$muForceSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, -15,10), col = rgb(1,0,1,1/4), main = "muChillSp")
hist(posterior_dbh$muChillSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, -15,10), col = rgb(1,0,1,1/4), main = "muPhotoSp")
hist(posterior_dbh$muPhotoSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 20,10), col = rgb(1,0,1,1/4), main = "mu_grand")
hist(posterior_dbh$mu_grand,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 40,10), col = rgb(1,0,1,1/4), main = "muPhenoSp")
hist(posterior_dbh$muPhenoSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 0,1), col = rgb(1,0,1,1/4), main = "betaTraitxForce")
hist(posterior_dbh$betaTraitxForce,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 0,1), col = rgb(1,0,1,1/4), main = "betaTraitxChill")
hist(posterior_dbh$betaTraitxChill,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 0,1), col = rgb(1,0,1,1/4), main = "betaTraitxPhoto")
hist(posterior_dbh$betaTraitxPhoto,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 10,5), col = rgb(1,0,1,1/4), main = "sigma_sp")
hist(posterior_dbh$sigma_sp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,2), col = rgb(1,0,1,1/4), main = "sigma_study")
hist(posterior_dbh$sigma_site,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,2), col = rgb(1,0,1,1/4), main = "sigma_traity")
hist(posterior_dbh$sigma_traity,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaForceSp")
hist(posterior_dbh$sigmaForceSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaChillSp")
hist(posterior_dbh$sigmaChillSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaPhotoSp")
hist(posterior_dbh$sigmaPhotoSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaPhenoSp")
hist(posterior_dbh$sigmaPhenoSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 10,5), col = rgb(1,0,1,1/4), main = "sigmapheno_y")
hist(posterior_dbh$sigmapheno_y,add=T,col=rgb(0,0,1,1/4))

dev.off()

################################################################################
# cn:
posterior_cn <- rstan::extract(readRDS(file = "output/cn_stanfit_np.RDS"))

pdf("figures/height_prior_post_dist_traitormdl.pdf", width = 15, height = 25)
par(mfrow = c(4,4))
#plot priors against posterior_cns
h1 <- hist(rnorm(1000, -15,15), col = rgb(1,0,1,1/4), main = "muForceSp")
hist(posterior_cn$muForceSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, -15,15), col = rgb(1,0,1,1/4), main = "muChillSp")
hist(posterior_cn$muChillSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, -15,10), col = rgb(1,0,1,1/4), main = "muPhotoSp")
hist(posterior_cn$muPhotoSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 20,5), col = rgb(1,0,1,1/4), main = "mu_grand")
hist(posterior_cn$mu_grand,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 40,10), col = rgb(1,0,1,1/4), main = "muPhenoSp")
hist(posterior_cn$muPhenoSp,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 0,1), col = rgb(1,0,1,1/4), main = "betaTraitxForce")
hist(posterior_cn$betaTraitxForce,add=T,col=rgb(0,0,1,1/4))

h1 <- hist(rnorm(1000, 0,1), col = rgb(1,0,1,1/4), main = "betaTraitxChill")
hist(posterior_cn$betaTraitxChill,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 0,1), col = rgb(1,0,1,1/4), main = "betaTraitxPhoto")
hist(posterior_cn$betaTraitxPhoto,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 10,5), col = rgb(1,0,1,1/4), main = "sigma_sp")
hist(posterior_cn$sigma_sp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,2), col = rgb(1,0,1,1/4), main = "sigma_study")
hist(posterior_cn$sigma_site,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,2), col = rgb(1,0,1,1/4), main = "sigma_traity")
hist(posterior_cn$sigma_traity,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaForceSp")
hist(posterior_cn$sigmaForceSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaChillSp")
hist(posterior_cn$sigmaChillSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaPhotoSp")
hist(posterior_cn$sigmaPhotoSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 5,5), col = rgb(1,0,1,1/4), main = "sigmaPhenoSp")
hist(posterior_cn$sigmaPhenoSp,col=rgb(0,0,1,1/4),add=T)

h1 <- hist(rnorm(1000, 10,5), col = rgb(1,0,1,1/4), main = "sigmapheno_y")
hist(posterior_cn$sigmapheno_y,add=T,col=rgb(0,0,1,1/4))

dev.off()

#######################################################################

