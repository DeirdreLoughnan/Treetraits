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
  setwd("~/Documents/github/Treetraits") 
}  else{
  setwd("~/Documents/github/Treetraits") # for midge
}

trtPheno <- read.csv("analysis/input/trtPhenoZScore.csv")
pheno.t <- read.csv("analysis/input/bbPhenoZScore.csv")

specieslist <- sort(unique(trtPheno$species))
sitelist <- sort(unique(trtPheno$transect))

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
                forcei = pheno.t$force.z,
                chilli = pheno.t$chillport.z,
                photoi = pheno.t$photo.z
)

mdlHtInt <- stan("analysis/Models/heightDummyIntGrandLat_wInteract.stan",
              data = ht.data,
              iter = 8000, warmup = 7000, chains=4
)

# control = list(adapt_delta = 0.99,
#                max_treedepth = 12))
save(mdlHtInt, file="analysis/output/htContLatHundowLat_interactions.Rdata")

load("analysis/output/htContLatHundowLat_interactions.Rdata")
sumer <- summary(mdlHtInt)$summary
                    
sumer <- data.frame(summary(mdlHtInt)$summary[c(
  "mu_grand", "b_lat",
  "b_tranE","b_tranlat", "muForceSp", "muChillSp", "muPhotoSp",
  "mu_b_inter_wp","mu_b_inter_wc","mu_b_inter_pc",
  "muPhenoSp","betaTraitxForce", "betaTraitxChill","betaTraitxPhoto",
  "sigma_traity" ,"sigma_sp", "sigmaForceSp", "sigmaChillSp", "sigmaPhotoSp","sigmaPhenoSp","sigmapheno_y"),
  c("mean","2.5%","25%","50%", "75%","97.5%")])
sumer

######################################################################
# 2. Leaf mass area
leafMA <- trtPheno[complete.cases(trtPheno$lma),]
leafMA$lma.z <- (leafMA$lma-mean(leafMA$lma,na.rm=TRUE))/(sd(leafMA$lma,na.rm=TRUE))*2

## put it on the same scale as height
leafMA$lma100 <- leafMA$lma*100
leafMA$lmalog <- log10(leafMA$lma)

lma <- list(yTraiti = leafMA$lma, 
            N = nrow(leafMA),
            n_spec = length(specieslist),
            trait_species = as.numeric(as.factor(leafMA$species)),
            n_tran = length(unique(leafMA$transect)),
            lati = leafMA$lat.z,
            tranE = as.numeric(leafMA$transect),
            Nph = nrow(pheno.t),
            phenology_species = as.numeric(as.factor(pheno.t$species)),
            yPhenoi = pheno.t$bb,
            forcei = pheno.t$force.z,
            chilli = pheno.t$chillport.z,
            photoi = pheno.t$photo.z
)

lma.100 <- list(yTraiti = leafMA$lma100, 
                N = nrow(leafMA),
                n_spec = length(specieslist),
                trait_species = as.numeric(as.factor(leafMA$species)),
                n_tran = length(unique(leafMA$transect)),
                lati = leafMA$latitude,
                tranE = as.numeric(leafMA$transect),
                Nph = nrow(pheno.t),
                phenology_species = as.numeric(as.factor(pheno.t$species)),
                yPhenoi = pheno.t$bb,
                forcei = pheno.t$force.z,
                chilli = pheno.t$chillport.z,
                photoi = pheno.t$photo.z
)


mdlLMAInt <- stan("analysis/Models/lmaDummyIntGrandLat_wInteract.stan",
               data = lma.100,
               include = FALSE, pars = c("y_hat"),
               iter = 4000, warmup = 3000, chains= 4,
               control = list(max_treedepth = 12))


save(mdlLMAInt, file="analysis/output/lmaContLatHundowLat_wInteract.Rdata")
# load("analysis/output/lmaContLatHundoLat.Rdata")

sumer <- summary(mdlLMAInt)$summary
sumer <- data.frame(summary(mdlLMAInt)$summary[c("mu_grand",
                                              "b_tranE","b_tranlat","b_lat", "muForceSp", "muChillSp",
                                              "muPhotoSp","muPhenoSp","betaTraitxForce", "betaTraitxChill",
                                              "betaTraitxPhoto","sigma_traity" ,"sigma_sp", "sigmaForceSp", "sigmaChillSp",
                                              "sigmaPhotoSp","sigmaPhenoSp","sigmapheno_y"),c("mean","2.5%","25%","50%", "75%","97.5%")])

# 
postLMA<- data.frame(rstan::extract(mdlLMA))
# 
par(mfrow = c(1,1))
hist(postLMA$betaTraitxForce, main = "betaTraitxForce",  xlim = c(-50, 50))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

####################################################
# 3. diameter at breast height
diam <- trtPheno[complete.cases(trtPheno$dbh),] 
diam$dbh.z2 <- (diam$dbh-mean(diam$dbh,na.rm=TRUE))/(sd(diam$dbh,na.rm=TRUE)*2)

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
                 forcei = pheno.t$force.z,
                 chilli = pheno.t$chillport.z,
                 photoi = pheno.t$photo.z
)

mdlDBHInt <- stan("analysis/Models/diamDummyIntGrandLat_wInteract.stan",
               data = dbh.data,
               iter = 6000, warmup = 5000, chains=4,
               include = FALSE, pars = c("y_hat")
)

save(mdlDBHInt, file="analysis/output/dbhContLatHundowLat_wInteract.Rdata")

######################################################################
# 4. stem specific density
stem <- trtPheno[complete.cases(trtPheno$ssd),]
stem$ssd <- round(stem$ssd, 3)

temp <- aggregate(stem["ssd"], stem["species"], FUN = mean)
stem$ssd.z2 <- (stem$ssd-mean(stem$ssd,na.rm=TRUE))/(sd(stem$ssd,na.rm=TRUE)*2)
stem$ssd100 <- stem$ssd*10
stem$ssdlog <- log10(stem$ssd)

ssd.data <- list(yTraiti = stem$ssd100, 
                 N = nrow(stem),
                 n_spec = length(specieslist),
                 trait_species = as.numeric(as.factor(stem$species)),
                 n_tran = length(unique(stem$transect)),
                 lati = stem$latitude,
                 tranE = as.numeric(stem$transect),
                 Nph = nrow(pheno.t),
                 phenology_species = as.numeric(as.factor(pheno.t$species)),
                 yPhenoi = pheno.t$bb,
                 forcei = pheno.t$force.z,
                 chilli = pheno.t$chillport.z,
                 photoi = pheno.t$photo.z
)


mdlSSDInt <- stan("analysis/Models/justDummyIntWideGrandLat_wInteract.stan",
               data = ssd.data,
               iter = 8000, warmup = 7000, chains=4,
               include = FALSE, pars = c("y_hat"),
               control = list(max_treedepth =12))

save(mdlSSDInt, file="analysis/output/ssdContLatHundowLat10_wInteract.Rdata")

sumer <- summary(mdlSSDInt)$summary
                    
######################################################################
# 5. carbon to nitrogen ratio
carbNit <- trtPheno[complete.cases(trtPheno$C.N),]
carbNit$cn.z2 <- (carbNit$C.N-mean(carbNit$C.N,na.rm=TRUE))/(sd(carbNit$C.N,na.rm=TRUE)*2)
carbNit$cn.port <- (carbNit$C.N/100)
carbNit$cn.port <- (carbNit$C.N/100)

cn.data <- list(yTraiti = carbNit$cn.port, 
                N = nrow(carbNit),
                n_spec = length(specieslist),
                trait_species = as.numeric(as.factor(carbNit$species)),
                n_tran = length(unique(carbNit$transect)),
                lati = carbNit$latitude,
                tranE = as.numeric(carbNit$transect),
                Nph = nrow(pheno.t),
                phenology_species = as.numeric(as.factor(pheno.t$species)),
                yPhenoi = pheno.t$bb,
                forcei = pheno.t$force.z,
                chilli = pheno.t$chillport.z,
                photoi = pheno.t$photo.z
)

mdlCN <- stan("analysis/Models/lmaDummyIntGrandLat_wInteract.stan",
              data = cn.data,
              iter = 4000, warmup = 3000, chains=4,
              include = FALSE, pars = c("y_hat")
              #,control = list(adapt_delta = 0.99, max_treedepth =12)
)

save(mdlCN, file="analysis/output/cnContLat.Rdata")

# mdlCN <- stan("stan/cnDummyIntGrand.stan",
#   data = cnZ.data,
#   iter = 4000, warmup = 3000, chains=4,
#   include = FALSE, pars = c("y_hat")
#   #,control = list(adapt_delta = 0.99, max_treedepth =12)
# )
sumer <- data.frame(summary(mdlCN3)$summary[c( "mu_grand",
                                               "b_tranE","b_tranlat", "muForceSp", "muChillSp", "muPhotoSp",
                                               "muPhenoSp","betaTraitxForce", "betaTraitxChill","betaTraitxPhoto",
                                               "sigma_traity" ,"sigma_sp", "sigmaForceSp", "sigmaChillSp", 
                                               "sigmaPhotoSp","sigmaPhenoSp","sigmapheno_y"),
                                            c("mean","2.5%","25%","50%", "75%","97.5%")])
sumer

postCN<- data.frame(rstan::extract(mdlCN))

par(mfrow = c(1,1))
hist(postCN$mu_grand, main = "betaTraitxForce",  xlim = c(-200, 200))
hist(rnorm(1000, 0.1,10), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postSSD$betaTraitxChill, main = "betaTraitxChill",  xlim = c(-200, 200))
hist(rnorm(1000, 0,90), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postSSD$betaTraitxPhoto, main = "betaTraitxPhoto",  xlim = c(-200, 200))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)
# looks the exact same as if you z-score the lma, so let's not since I am not sure it is a good idea


save(mdlCN, file="output/cnContLat.Rdata")

#####
######################################################################
# 5. leaf nitrogen content
nit <- trtPheno[complete.cases(trtPheno$per.N),]
#carbNit$n.z2 <- (carbNit$-mean(carbNit$C.N,na.rm=TRUE))/(sd(carbNit$C.N,na.rm=TRUE)*2)
# carbNit$n.port <- (carbNit$per.N/100)
# carbNit$nLog <- log10(carbNit$per.N)

nit.data <- list(yTraiti = nit$per.N, 
                 N = nrow(nit),
                 n_spec = length(specieslist),
                 trait_species = as.numeric(as.factor(nit$species)),
                 n_tran = length(unique(nit$transect)),
                 lati = nit$latitude,
                 tranE = as.numeric(nit$transect),
                 Nph = nrow(pheno.t),
                 phenology_species = as.numeric(as.factor(pheno.t$species)),
                 yPhenoi = pheno.t$bb,
                 forcei = pheno.t$force.z,
                 chilli = pheno.t$chillport.z,
                 photoi = pheno.t$photo.z
)

mdlPerN <- stan("analysis/Models/cnDummyIntGrandLat_wInteract.stan",
                data = nit.data,
                iter = 8000, warmup = 4000, chains=4,
                include = FALSE, pars = c("y_hat")
                ,control = list(max_treedepth =12)
)

save(mdlPerN, file="analysis/output/lncContLatHundowLat.Rdata")

# load("output/nitDummyInt.Rdata")
#  ssm <- as.shinystan(mdl)
#  launch_shinystan(ssm)

sumer3 <- data.frame(summary(mdlPerN3)$summary[c( "mu_grand",
                                                  "b_tranE","b_tranlat", "muForceSp", "muChillSp", "muPhotoSp",
                                                  "muPhenoSp","betaTraitxForce", "betaTraitxChill","betaTraitxPhoto",
                                                  "sigma_traity" ,"sigma_sp", "sigmaForceSp", "sigmaChillSp", 
                                                  "sigmaPhotoSp","sigmaPhenoSp","sigmapheno_y"),
                                               c("mean","2.5%","25%","50%", "75%","97.5%")])
sumer3


ModelFit <- rstan::extract(mdlPerN)

muSp <- data.frame(ModelFit$b_muSp)
muSpMean <- colMeans(muSp)

betaForceSp <- data.frame(ModelFit$betaForceSp)
betaForceSpMean <- colMeans(betaForceSp)

quantile2575 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.25, 0.75))
  return(returnQuanilte)
}

bf_quan <- apply(betaForceSp, 2, quantile2575) 
mu_quan <- apply(muSp, 2, quantile2575)

bfs <- rbind(betaForceSpMean, bf_quan)
bfs_t <- t(bfs)
bfs_df <- data.frame(bfs_t)
colnames(bfs_df)[colnames(bfs_df) == "X25."] <- "force25"
colnames(bfs_df)[colnames(bfs_df) == "X75."] <- "force75"
bfs_df$species <- specieslist

mg<- rbind(muSpMean, mu_quan)
mg_t <- t(mg)
mg_df <- data.frame(mg_t)
colnames(mg_df)[colnames(mg_df) == "X25."] <- "trait25"
colnames(mg_df)[colnames(mg_df) == "X75."] <- "trait75"
mg_df$species <- specieslist

muForceSp <- data.frame(ModelFit$muForceSp)
muForceSpMean <- colMeans(muForceSp)

betaTraitxForce <- data.frame(ModelFit$betaTraitxForce)
betaTraitxForceMean <- colMeans(betaTraitxForce)

# mg_df_east <- mg_df[mg_df$species %in% eastSp, ]
# mg_df_west <- mg_df[mg_df$species %in% westSp, ]
# 
# bfs_df_east <- bfs_df[bfs_df$species %in% eastSp, ]
# bfs_df_west <- bfs_df[bfs_df$species %in% westSp, ]

plot( x= mg_df$muSpMean, y = bfs_df$betaForceSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bfs_df$force25), max(bfs_df$force75)), ylab = "Species level forcing slope", xlab = "Trait value", cex.lab = 1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df[,"muSpMean"], # x mean
  bfs_df[,"force25"], # y 25
  mg_df[,"muSpMean"],
  bfs_df[,"force75"],
  length = 0, col= "#218380", lwd = 2
)

arrows(
  mg_df[,"trait25"], # x mean
  bfs_df[,"betaForceSpMean"], # y 25
  mg_df[,"trait75"], # x mean
  bfs_df[,"betaForceSpMean"],
  length = 0, col = "#218380", lwd = 2
)


mtext(side = 3, text = "C:N, Forcing", adj = 0, cex = 1.25)
for(j in 1:length(muForceSp[,1])){
  abline(a = muForceSp[j,], b = betaTraitxForceMean, col=alpha("#73d2de", 0.085))
}
abline(a=muForceSpMean, b=betaTraitxForceMean, col = "black")


#### height #####
post<- data.frame(rstan::extract(mdlHt3))


hist(post$b_tranE, main = "betaTraitxForce",  xlim = c(-200, 200))
hist(rnorm(1000, 0,15), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postSSD$betaTraitxChill, main = "betaTraitxChill",  xlim = c(-200, 200))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)

hist(postSSD$betaTraitxPhoto, main = "betaTraitxPhoto",  xlim = c(-200, 200))
hist(rnorm(1000, 0,25), col=rgb(1,0,1,1/4), add = T)
abline(v = 0, col="red", lwd=3, lty=2)
# looks the exact same as if you z-score the lma, so let's not since I am not sure it is a good idea

