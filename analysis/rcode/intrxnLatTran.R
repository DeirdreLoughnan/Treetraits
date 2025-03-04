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
load("analysis/output/htContLatHundoLatFinal.Rdata")
sumerht <- summary(mdlHt)$summary
postHt <- (rstan::extract(mdlHt))

load("analysis/output/lmaContLatHundoLatFinal.Rdata")
postLMA <- rstan::extract(mdlLMA)
sumerlma <- summary(mdlLMA)$summary

load("analysis/output/dbhContLatHundoLatFinal.Rdata")
postDBH <- rstan::extract(mdlDBH)
sumeDBH <- summary(mdlDBH)$summary

load("analysis/output/ssdContLatHundoLatFinal.Rdata")
postSSD <- rstan::extract(mdlSSD)
sumerSSD <- summary(mdlSSD)$summary

load("analysis/output/lncContLatHundoLatFinal.Rdata")
postLNC <- rstan::extract(mdlPerN)
sumerLNC <- summary(mdlPerN)$summary

tranW <- 0
tranE <- 1

# height:
postHt <- data.frame(postHt)
trtHt <- postHt[2001:4000,c("mu_grand", "b_tranE", "b_tranlat")]

lati <- seq(40,55, length.out = 30)

output <- data.frame(cbind(
  rep(1:(length(lati)*nrow(trtHt))), 
  rep(trtHt$mu_grand, times = length(lati)),
  rep(trtHt$b_tranE, times = length(lati)),
  rep(trtHt$b_tranlat, times = length(lati)),
  rep(lati, each = 1000)))

names(output) <- c("iter", "mu_grand","b_tranE", "b_tranLat","lat")
output$latTrend <- NA  

for (i in 1:nrow(output)){
  temp <- output$mu_grand[i] + output$b_tranE[i]* tranE + output$b_tranLat[i]* (tranE*output$lat[i])
  output$latTrend[output$iter == i] <- temp
}

pdf("analysis/figures/traitLatInt.pdf", width = 15, height = 4)
par(mfrow = c(1,5), mar = c(5.1, 4.5, 4.1, 2.1), mgp=c(2.25,1,0))
plot(NA, xlim = c(42,54), ylim = c(-5, 10),
     xlab = "Latitude", ylab = "Height (m)",
     # bty = "n",
     xaxt = "n",
     # yaxt = "n",
     cex.lab = 1.5,
     cex.axis = 1.5)
axis(side = 1, at = seq(42,54, by =1), cex.axis =1.5)
text(42.5, 9, 'a)', cex = 2)

polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile(output$mu_grand, prob =c (0.025)),quantile(output$mu_grand, prob =c (0.025)), 
              quantile(output$mu_grand, prob =c (0.975)),quantile(output$mu_grand, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(0 / 255, 205 / 255, 205 / 255, alpha = 0.3)),border = NA )  

polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile(subset(output, lat == 40)$latTrend, prob =c (0.025)),quantile(subset(output, lat == 55)$latTrend, prob =c (0.025)), 
              quantile(subset(output, lat == 55)$latTrend, prob =c (0.975)),quantile(subset(output, lat == 40)$latTrend, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(0 / 255, 139 / 255, 139 / 255, alpha = 0.3)),border = NA )  

abline(a = mean(output$mu_grand) , b= 0, col = "navy" )
abline(lm(output$latTrend ~ output$lat ), col = "navy", lty = 2)

legend("bottomleft",legend = c("Eastern transect", 
                               "Western transect"),
col = c("navy", "navy"), lty = c(1,2), bty = "n", lwd = 2)



# DBH
postDBH <- data.frame(postDBH)
trtDBH <- postDBH[2001:4000,c("mu_grand", "b_tranE", "b_tranlat")]

output <- data.frame(cbind(
  rep(1:(length(lati)*nrow(trtDBH))), 
  rep(trtDBH$mu_grand, times = length(lati)), 
  rep(trtDBH$b_tranE, times = length(lati)),
  rep(trtDBH$b_tranlat, times = length(lati)),
  rep(lati, each = 1000)))

names(output) <- c("iter", "mu_grand","b_tranE", "b_tranLat","lat")
output$latTrend <- NA  

for (i in 1:nrow(output)){
  temp <- output$mu_grand[i] + output$b_tranE[i]* tranE + output$b_tranLat[i]* (tranE*output$lat[i])
  output$latTrend[output$iter == i] <- temp
}

plot(NA, xlim = c(42,54), ylim = c(-5, 15),
     xlab = "Latitude", ylab = "Diameter at breast height (m)",
     xaxt = "n",
     cex.lab = 1.5,
     cex.axis = 1.5)
axis(side = 1, at = seq(42,54, by =1), cex.axis =1.5)
text(42.5, 14, 'b)', cex = 2)

polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile(output$mu_grand, prob =c (0.025)),quantile(output$mu_grand, prob =c (0.025)), 
              quantile(output$mu_grand, prob =c (0.975)),quantile(output$mu_grand, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(218 / 255, 165 / 255, 32 / 255, alpha = 0.3)),border = NA )  

polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile(subset(output, lat == 40)$latTrend, prob =c (0.025)),quantile(subset(output, lat == 55)$latTrend, prob =c (0.025)), 
              quantile(subset(output, lat == 55)$latTrend, prob =c (0.975)),quantile(subset(output, lat == 40)$latTrend, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(139 / 255, 105 / 255, 20 / 255, alpha = 0.3)),border = NA )  

abline(a = mean(output$mu_grand) , b= 0, col = "chocolate4" )
abline(lm(output$latTrend ~ output$lat ), col = "chocolate4", lty = 2)

# SSD:
postSSD <- data.frame(postSSD)
trtSSD <- postSSD[2001:4000,c("mu_grand", "b_tranE", "b_tranlat")]

output <- data.frame(cbind(
  rep(1:(length(lati)*nrow(trtSSD))), 
  rep(trtSSD$mu_grand, times = length(lati)), 
  rep(trtSSD$b_tranE, times = length(lati)),
  rep(trtSSD$b_tranlat, times = length(lati)),
  rep(lati, each = 1000)))

names(output) <- c("iter", "mu_grand","b_tranE", "b_tranLat","lat")
output$latTrend <- NA  

for (i in 1:nrow(output)){
  temp <- output$mu_grand[i] + output$b_tranE[i]* tranE + output$b_tranLat[i]* (tranE*output$lat[i])
  output$latTrend[output$iter == i] <- temp
}

plot(NA, xlim = c(42,54), ylim = c(0,1),
     xlab = "Latitude", ylab = bquote('Wood specific density'~(g/cm^2)),
     xaxt = "n",
     cex.lab = 1.5,
     cex.axis = 1.5)
axis(side = 1, at = seq(42,54, by =1), cex.axis =1.5)
text(42.5, 0.9, 'c)', cex = 2)
polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile((output$mu_grand/100), prob =c (0.025)),quantile((output$mu_grand/100), prob =c (0.025)), 
              quantile((output$mu_grand/100), prob =c (0.975)),quantile((output$mu_grand/100), prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(176 / 255, 48 / 255, 96 / 255, alpha = 0.3)),border = NA )  

polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile((subset(output, lat == 40)$latTrend/100), prob =c (0.025)),quantile((subset(output, lat == 55)$latTrend/100), prob =c (0.025)), 
              quantile((subset(output, lat == 55)$latTrend/100), prob =c (0.975)),quantile((subset(output, lat == 40)$latTrend/100), prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(139 / 255, 28 / 255, 98 / 255, alpha = 0.3)),border = NA )  

abline(a = (mean(output$mu_grand/100)) , b= 0, col = "maroon4" )
abline(lm((output$latTrend/100) ~ output$lat ), col = "maroon4", lty = 2)

# height:
postLMA <- data.frame(postLMA)
trtLMA <- postLMA[2001:4000,c("mu_grand", "b_tranE", "b_tranlat")]

output <- data.frame(cbind(
  rep(1:(length(lati)*nrow(trtLMA))), 
  rep(trtLMA$mu_grand, times = length(lati)), 
  rep(trtLMA$b_tranE, times = length(lati)),
  rep(trtLMA$b_tranlat, times = length(lati)),
  rep(lati, each = 1000)))

names(output) <- c("iter", "mu_grand","b_tranE", "b_tranLat","lat")
output$latTrend <- NA  

for (i in 1:nrow(output)){
  temp <- output$mu_grand[i] + output$b_tranE[i]* tranE + output$b_tranLat[i]* (tranE*output$lat[i])
  output$latTrend[output$iter == i] <- temp
}

plot(NA, xlim = c(42,54), ylim = c(00, 0.10),
     xlab = "Latitude", ylab =  bquote('Leaf mass area '~(g/m^2)),
     xaxt = "n",
     cex.lab = 1.5,
     cex.axis = 1.5)
axis(side = 1, at = seq(42,54, by =1), cex.axis =1.5)
text(42.5, 0.09, 'd)', cex = 2)
polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile(output$mu_grand/100, prob =c (0.025)),quantile(output$mu_grand/100, prob =c (0.025)), 
              quantile(output$mu_grand/100, prob =c (0.975)),quantile(output$mu_grand/100, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(162 / 255, 205 / 255, 90 / 255, alpha = 0.3)),border = NA )  

polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile(subset(output, lat == 40)$latTrend/100, prob =c (0.025)),quantile(subset(output, lat == 55)$latTrend/100, prob =c (0.025)), 
              quantile(subset(output, lat == 55)$latTrend/100, prob =c (0.975)),quantile(subset(output, lat == 40)$latTrend/100, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(85 / 255, 107 / 255, 47 / 255, alpha = 0.3)),border = NA )  

abline(a = (mean(output$mu_grand/100)) , b= 0, col = "darkgreen" )
abline(lm((output$latTrend/100) ~ output$lat ), col = "darkgreen", lty = 2)

# LNC:
postLNC <- data.frame(postLNC)
trtLNC <- postLNC[2001:4000,c("mu_grand", "b_tranE", "b_tranlat")]

output <- data.frame(cbind(
  rep(1:(length(lati)*nrow(trtLNC))), 
  rep(trtLNC$mu_grand, times = length(lati)), 
  rep(trtLNC$b_tranE, times = length(lati)),
  rep(trtLNC$b_tranlat, times = length(lati)),
  rep(lati, each = 1000)))

names(output) <- c("iter", "mu_grand","b_tranE", "b_tranLat","lat")
output$latTrend <- NA  

for (i in 1:nrow(output)){
  temp <- output$mu_grand[i] + output$b_tranE[i]* tranE + output$b_tranLat[i]* (tranE*output$lat[i])
  output$latTrend[output$iter == i] <- temp
}

plot(NA, xlim = c(42,54), ylim = c(-5, 10),
     xlab = "Latitude", ylab = "Leaf nitrogen content (%)",
     xaxt = "n",
     cex.lab = 1.5,
     cex.axis = 1.5)
axis(side = 1, at = seq(42,54, by =1), cex.axis =1.5)
text(42.5, 9, 'e)', cex = 2)

polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile(output$mu_grand, prob =c (0.025)),quantile(output$mu_grand, prob =c (0.025)), 
              quantile(output$mu_grand, prob =c (0.975)),quantile(output$mu_grand, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(155 / 255, 32 / 255, 240 / 255, alpha = 0.3)),border = NA )  

polygon(x = c(min(output$lat), max(output$lat), max(output$lat), min(output$lat)),    # X-Coordinates of polygon
        y = c(quantile(subset(output, lat == 40)$latTrend, prob =c (0.025)),quantile(subset(output, lat == 55)$latTrend, prob =c (0.025)), 
              quantile(subset(output, lat == 55)$latTrend, prob =c (0.975)),quantile(subset(output, lat == 40)$latTrend, prob =c (0.975))),                             # Y-Coordinates of polygon
        col = (rgb(85 / 255, 26 / 255, 139 / 255, alpha = 0.3)),border = NA )  

abline(a = mean(output$mu_grand) , b= 0, col = "purple" )
abline(lm(output$latTrend ~ output$lat ), col = "purple4", lty = 2)

dev.off()
