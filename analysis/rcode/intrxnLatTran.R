# February 26, 2025

# aim of this code is to make remake the transect by latitude plots

rm(list=ls())
options(stringsAsFactors = FALSE)

library(stringr)
library(plyr)
library(rstan)
library(reshape2)
library(cowplot)
require(dplyr)
require(plotrix)
require(ggplot2)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/Treetraits") 
}  else{
  setwd("~/Documents/github/Treetraits") # for midge
}

spInfo <- read.csv("analysis/input/species_ring.csv")
pheno <- read.csv("analysis/input/phenoDataWChill.csv")
# trtPheno <- read.csv("input/trtData.csv")
trtPheno <- read.csv("analysis/input/trtPhenoDummy.csv")

# Ht and DBH and CN natural, LMA and SSD rescaled by 100
# cues and latitude still z scored by only by 1 sd
load("analysis/output/htContLatHundowLat.Rdata")
sumerht <- summary(mdlHt)$summary
postHt <- (rstan::extract(mdlHt))

load("analysis/output/lmaContLatHundowLat.Rdata")
postLMA <- rstan::extract(mdlLMA)
sumerlma <- summary(mdlLMA)$summary

load("analysis/output/dbhContLatHundowLat.Rdata")
postDBH <- rstan::extract(mdlDBH)
sumeDBH <- summary(mdlDBH)$summary

load("analysis/output/ssdContLatHundowLat10.Rdata")
postSSD <- rstan::extract(mdlSSD)
sumerSSD <- summary(mdlSSD)$summary

load("analysis/output/lncContLatHundowLat.Rdata")
postLNC <- rstan::extract(mdlPerN)
sumerLNC <- summary(mdlPerN)$summary

tranW <- 0
tranE <- 1

# height:
postHt <- data.frame(postHt)
trtHt <- postHt[2001:4000,c("mu_grand", "b_tranE", "b_tranlat", "b_lat")]

lati <- seq(40,55, length.out = 30)

output <- data.frame(cbind(
  rep(1:(length(lati)*nrow(trtHt))), 
  rep(trtHt$mu_grand, times = length(lati)),
  rep(trtHt$b_tranE, times = length(lati)),
  rep(trtHt$b_tranlat, times = length(lati)),
  rep(trtHt$b_lat, times = length(lati)),
  rep(lati, each = 1000)))

names(output) <- c("iter", "mu_grand","b_tranE", "b_tranLat","b_lat","lat")

for (i in 1:nrow(output)){
  # temp <- output$mu_grand[i] + output$b_tranE[i]* tranE + output$b_lat[i]* output$lat[i] + output$b_tranLat[i]* (tranE*output$lat[i])
  # output$latTrendE[output$iter == i] <- temp
  tempW <- output$mu_grand[i] + output$b_tranE[i]* tranW + output$b_lat[i]* output$lat[i] + output$b_tranLat[i]* (tranW*output$lat[i])
  output$latTrendW[output$iter == i] <- tempW
}

pdf("analysis/figures/traitLatIntJune15.pdf", width = 15, height = 4)
par(mfrow = c(1,5), mar = c(5.1, 4.5, 4.1, 2.1), mgp=c(2.25,1,0))
plot(NA, xlim = c(42,54), ylim = c(-5, 10),
     xlab = "Latitude", ylab = "Height (m)",
     # bty = "n",
     xaxt = "n",
     # yaxt = "n",
     cex.lab = 1.5,
     cex.axis = 1.5)
axis(side = 1, at = seq(42,54, by =1), cex.axis =1.5)
text(42.5, 9.5, 'a)', cex = 2)

polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile(subset(output, lat == 40)$latTrendW, prob =c (0.025)),quantile(subset(output, lat == 55)$latTrendW, prob =c (0.025)), 
              quantile(subset(output, lat == 55)$latTrendW, prob =c (0.975)),quantile(subset(output, lat == 40)$latTrendW, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(0 / 255, 205 / 255, 205 / 255, alpha = 0.3)),border = NA )  

polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile(subset(output, lat == 40)$latTrendE, prob =c (0.025)),quantile(subset(output, lat == 55)$latTrendE, prob =c (0.025)), 
              quantile(subset(output, lat == 55)$latTrendE, prob =c (0.975)),quantile(subset(output, lat == 40)$latTrendE, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(0 / 255, 139 / 255, 139 / 255, alpha = 0.3)),border = NA )  

abline(lm(output$latTrendW ~ output$lat), col = "navy" )
abline(lm(output$latTrendE ~ output$lat ), col = "navy", lty = 2)

legend("bottomleft",legend = c("Eastern transect", 
                               "Western transect"),
col = c("navy", "navy"), lty = c(2,1), bty = "n", lwd = 2, cex = 1.5)



# DBH
postDBH <- data.frame(postDBH)
trtDBH <- postDBH[2001:4000,c("mu_grand", "b_tranE", "b_tranlat", "b_lat")]

output <- data.frame(cbind(
  rep(1:(length(lati)*nrow(trtDBH))), 
  rep(trtDBH$mu_grand, times = length(lati)), 
  rep(trtDBH$b_tranE, times = length(lati)),
  rep(trtDBH$b_tranlat, times = length(lati)),
  rep(trtDBH$b_lat, times = length(lati)),
  rep(lati, each = 1000)))

names(output) <- c("iter", "mu_grand","b_tranE", "b_tranLat","b_lat","lat")
output$latTrend <- NA  

for (i in 1:nrow(output)){
  temp <- output$mu_grand[i] + output$b_tranE[i]* tranE + output$b_lat[i]* output$lat[i] + output$b_tranLat[i]* (tranE*output$lat[i])
  output$latTrendE[output$iter == i] <- temp
  tempW <- output$mu_grand[i] + output$b_tranE[i]* tranW + output$b_lat[i]* output$lat[i] + output$b_tranLat[i]* (tranW*output$lat[i])
  output$latTrendW[output$iter == i] <- tempW
}

plot(NA, xlim = c(42,54), ylim = c(-5, 15),
     xlab = "Latitude", ylab = "Diameter at breast height (cm)",
     xaxt = "n",
     cex.lab = 1.5,
     cex.axis = 1.5)
axis(side = 1, at = seq(42,54, by =1), cex.axis =1.5)
text(42.5, 14, 'b)', cex = 2)

polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile(subset(output, lat == 40)$latTrendW, prob =c (0.025)),quantile(subset(output, lat == 55)$latTrendW, prob =c (0.025)), 
              quantile(subset(output, lat == 55)$latTrendW, prob =c (0.975)),quantile(subset(output, lat == 40)$latTrendW, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(218 / 255, 165 / 255, 32 / 255, alpha = 0.3)),border = NA ) 

polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile(subset(output, lat == 40)$latTrendE, prob =c (0.025)),quantile(subset(output, lat == 55)$latTrendE, prob =c (0.025)), 
              quantile(subset(output, lat == 55)$latTrendE, prob =c (0.975)),quantile(subset(output, lat == 40)$latTrendE, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(139 / 255, 105 / 255, 20 / 255, alpha = 0.3)),border = NA )  

abline(lm(output$latTrendW ~ output$lat), col = "chocolate4" )
abline(lm(output$latTrendE ~ output$lat ), col = "chocolate4", lty = 2)

# SSD:
postSSD <- data.frame(postSSD)
trtSSD <- postSSD[2001:4000,c("mu_grand", "b_tranE", "b_tranlat", "b_lat")]

output <- data.frame(cbind(
  rep(1:(length(lati)*nrow(trtSSD))), 
  rep(trtSSD$mu_grand, times = length(lati)), 
  rep(trtSSD$b_tranE, times = length(lati)),
  rep(trtSSD$b_tranlat, times = length(lati)),
  rep(trtSSD$b_lat, times = length(lati)),
  rep(lati, each = 1000)))

names(output) <- c("iter", "mu_grand","b_tranE", "b_tranLat","b_lat","lat")
output$latTrend <- NA  

for (i in 1:nrow(output)){
  temp <- output$mu_grand[i] + output$b_tranE[i]* tranE + output$b_lat[i]* output$lat[i] + output$b_tranLat[i]* (tranE*output$lat[i])
  output$latTrendE[output$iter == i] <- temp
  tempW <- output$mu_grand[i] + output$b_tranE[i]* tranW + output$b_lat[i]* output$lat[i] + output$b_tranLat[i]* (tranW*output$lat[i])
  output$latTrendW[output$iter == i] <- tempW
}


plot(NA, xlim = c(42,54), ylim = c(0,0.8),
     xlab = "Latitude", ylab = bquote('Wood specific density'~(g/cm^3)),
     xaxt = "n",
     cex.lab = 1.5,
     cex.axis = 1.5)
axis(side = 1, at = seq(42,54, by =1), cex.axis =1.5)
text(42.5, 0.75, 'c)', cex = 2)
polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile((subset(output, lat == 40)$latTrendW)/10, prob =c (0.025)),quantile((subset(output, lat == 55)$latTrendW)/10, prob =c (0.025)), 
              quantile((subset(output, lat == 55)$latTrendW)/10, prob =c (0.975)),quantile((subset(output, lat == 40)$latTrendW)/10, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(176 / 255, 48 / 255, 96 / 255, alpha = 0.3)),border = NA )  

polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile((subset(output, lat == 40)$latTrendE), prob =c (0.025))/10,quantile((subset(output, lat == 55)$latTrendE)/10, prob =c (0.025)), 
              quantile((subset(output, lat == 55)$latTrendE), prob =c (0.975))/10,quantile((subset(output, lat == 40)$latTrendE)/10, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(139 / 255, 28 / 255, 98 / 255, alpha = 0.3)),border = NA )  

abline(lm((output$latTrendW)/10 ~ output$lat ), col = "maroon4" )
abline(lm((output$latTrendE)/10 ~ output$lat ), col = "maroon4", lty = 2)

# LMA:
postLMA <- data.frame(postLMA)
trtLMA <- postLMA[2001:4000,c("mu_grand", "b_tranE", "b_tranlat", "b_lat")]

output <- data.frame(cbind(
  rep(1:(length(lati)*nrow(trtLMA))), 
  rep(trtLMA$mu_grand, times = length(lati)), 
  rep(trtLMA$b_tranE, times = length(lati)),
  rep(trtLMA$b_tranlat, times = length(lati)),
  rep(trtLMA$b_lat, times = length(lati)),
  rep(lati, each = 1000)))

names(output) <- c("iter", "mu_grand","b_tranE", "b_tranLat","b_lat","lat")
output$latTrend <- NA  

for (i in 1:nrow(output)){
  temp <- output$mu_grand[i] + output$b_tranE[i]* tranE + output$b_lat[i]* output$lat[i] + output$b_tranLat[i]* (tranE*output$lat[i])
  output$latTrendE[output$iter == i] <- temp
  tempW <- output$mu_grand[i] + output$b_tranE[i]* tranW + output$b_lat[i]* output$lat[i] + output$b_tranLat[i]* (tranW*output$lat[i])
  output$latTrendW[output$iter == i] <- tempW
}

plot(NA, xlim = c(42,54), ylim = c(00, 0.10),
     xlab = "Latitude", ylab =  bquote('Leaf mass area '~(g/cm^2)),
     xaxt = "n",
     cex.lab = 1.5,
     cex.axis = 1.5)
axis(side = 1, at = seq(42,54, by =1), cex.axis =1.5)
text(42.5, 0.095, 'd)', cex = 2)
polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile(subset(output, lat == 40)$latTrendW/100, prob =c (0.025)),quantile(subset(output, lat == 55)$latTrendW/100, prob =c (0.025)), 
              quantile(subset(output, lat == 55)$latTrendW/100, prob =c (0.975)),quantile(subset(output, lat == 40)$latTrendW/100, prob =c (0.975))),        
        col = (rgb(162 / 255, 205 / 255, 90 / 255, alpha = 0.3)),border = NA )  

polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile(subset(output, lat == 40)$latTrendE/100, prob =c (0.025)),quantile(subset(output, lat == 55)$latTrendE/100, prob =c (0.025)), 
              quantile(subset(output, lat == 55)$latTrendE/100, prob =c (0.975)),quantile(subset(output, lat == 40)$latTrendE/100, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(85 / 255, 107 / 255, 47 / 255, alpha = 0.3)),border = NA )  

abline(lm((output$latTrendW/100) ~ output$lat ), col = "darkgreen" )
abline(lm((output$latTrendE/100) ~ output$lat ), col = "darkgreen", lty = 2)

# LNC:
postLNC <- data.frame(postLNC)
trtLNC <- postLNC[2001:4000,c("mu_grand", "b_tranE", "b_tranlat","b_lat")]

output <- data.frame(cbind(
  rep(1:(length(lati)*nrow(trtLNC))), 
  rep(trtLNC$mu_grand, times = length(lati)), 
  rep(trtLNC$b_tranE, times = length(lati)),
  rep(trtLNC$b_tranlat, times = length(lati)),
  rep(trtLNC$b_lat, times = length(lati)),
  rep(lati, each = 1000)))

names(output) <- c("iter", "mu_grand","b_tranE", "b_tranLat","b_lat","lat")
output$latTrend <- NA  

for (i in 1:nrow(output)){
  # temp <- output$mu_grand[i]  + output$b_lat[i]* output$lat[i] 
  # output$latTrendE[output$iter == i] <- temp
  tempW <- output$mu_grand[i] + output$b_tranE[i]* tranW + output$b_lat[i]* output$lat[i] + output$b_tranLat[i]* (tranW*output$lat[i])
  output$latTrendW[output$iter == i] <- tempW
}

plot(NA, xlim = c(42,54), ylim = c(-5, 10),
     xlab = "Latitude", ylab = "Leaf nitrogen content (%)",
     xaxt = "n",
     cex.lab = 1.5,
     cex.axis = 1.5)
axis(side = 1, at = seq(42,54, by =1), cex.axis =1.5)
text(42.5, 9, 'e)', cex = 2)

polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile(subset(output, lat == 40)$latTrendW, prob =c (0.025)),quantile(subset(output, lat == 55)$latTrendW, prob =c (0.025)), 
              quantile(subset(output, lat == 55)$latTrendW, prob =c (0.975)),quantile(subset(output, lat == 40)$latTrendW, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(155 / 255, 32 / 255, 240 / 255, alpha = 0.3)),border = NA )  


abline(lm(output$latTrendW ~ output$lat ), col = "purple4", lty = 4 )
# abline(lm(output$latTrendE ~ output$lat ), col = "purple4", lty = 2)

dev.off()

sumer <- data.frame(summary(mdlPerN)$summary[c("mu_grand",
                                              "b_tranE","b_tranlat","b_lat", "muForceSp", "muChillSp",
                                              "muPhotoSp","muPhenoSp","betaTraitxForce", "betaTraitxChill",
                                              "betaTraitxPhoto","sigma_traity" ,"sigma_sp", "sigmaForceSp", "sigmaChillSp",
                                              "sigmaPhotoSp","sigmaPhenoSp","sigmapheno_y"),c("mean","2.5%","25%","50%", "75%","97.5%")]);sumer
