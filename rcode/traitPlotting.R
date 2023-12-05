# started March 18, 2023 by Deirdre

# aim of this code is to make the plots for BC traits projects

#1. Figure comapring the site differences for each trait
#2. Figure comparing the betatraitCues for each figure

rm(list=ls())
options(stringsAsFactors = FALSE)

library(stringr)
library(plyr)
library(rstan)
library(reshape2)
library(cowplot)
require(dplyr)
require(plotrix)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/Treetraits") 
}  else{
  setwd("/home/deirdre/Treetraits") # for midge
}

spInfo <- read.csv("input/species_ring.csv")
 pheno <- read.csv("input/phenoDataWChill.csv")
# trtPheno <- read.csv("input/trtData.csv")
trtPheno <- read.csv("input/trtPhenoDummy.csv")

load("output/heightDummyIntGrandZ25.Rdata")
sumerht <- summary(mdlHt)$summary
postHt <- rstan::extract(mdlHt)

load("output/lmaDummyIntGrandZ25.Rdata")
postLMA <- rstan::extract(mdlLMA)
sm.sum <- summary(mdlLMA)$summary

load("output/dbhDummyIntGrandZ25.Rdata")
postDBH <- rstan::extract(mdlDBH)

load("output/ssdDummyIntGrandZ25.Rdata")
postCN <- rstan::extract(mdl)

load("output/cnDummyIntGrandZ25.Rdata")
postSSD <- rstan::extract(mdlSSD)

###### Compare the spp and site level effects across traits ###########
# col1 <- rgb(204 / 255, 102 / 255, 119 / 255, alpha = 0.8)
# col2 <- rgb(68 / 255, 170 / 255, 153 / 255, alpha = 0.6)

col1.sp <-c( rgb(204 / 255, 105 / 255, 112 / 255, alpha = 0.5)) # red
col2.sp <- c( rgb(205 / 255, 122 / 255, 0 / 255, alpha = 0.5)) # yellow
col3.sp <-c( rgb(9/ 255, 168 / 255, 82 / 255, alpha = 0.6)) # green
col4.sp <- c( rgb(34 / 255, 166 / 255, 167 / 255, alpha = 0.5)) # blue
col5.sp <- c( rgb(141 / 255, 34 / 255, 171 / 255, alpha = 0.5)) # purple

hist(postHt$b_tranE, col = col2.sp, main = "", xlim = c(-10,10))
hist(postDBH$b_tranE, col = col3.sp, main = "", add = T)
hist(postCN$b_tranE, col = col1.sp, main = "", add = T)
hist(postLMA$b_tranE, col = col4.sp, main = "", add = T)
hist(postSSD$b_tranE, col = col5.sp, main = "", add = T)

hist(postHt$muForceSp, col = col2.sp, main = "", xlim = c(-20,0))
hist(postDBH$muForceSp, col = col3.sp, main = "", add = T)
hist(postCN$muForceSp, col = col1.sp, main = "", add = T)
hist(postLMA$muForceSp, col = col4.sp, main = "", add = T)
hist(postSSD$muForceSp, col = col5.sp, main = "", add = T)

hist(postDBH$muChillSp, col = col2.sp, main = "", xlim = c(-25,0))
hist(postHt$muChillSp, col = col3.sp, main = "", add = T)
hist(postCN$muChillSp, col = col1.sp, main = "", add = T)
hist(postLMA$muChillSp, col = col4.sp, main = "", add = T)
hist(postSSD$muChillSp, col = col5.sp, main = "", add = T)

hist(postHt$muPhotoSp, col = col2.sp, main = "", xlim = c(-10, 10), ylim = c(0,1800))
hist(postDBH$muPhotoSp, col = col3.sp, main = "", add = T)
hist(postCN$muPhotoSp, col = col1.sp, main = "", add = T)
hist(postLMA$muPhotoSp, col = col4.sp, main = "", add = T)
hist(postSSD$muPhotoSp, col = col5.sp, main = "", add = T)

legend("topright",legend = c(expression("Height"),
                            expression("LMA"),
                            expression("DBH"),
                            expression("SSD"),
                            expression("C:N")
                            ),
       col = c(col2.sp, col4.sp, col3.sp, col5.sp, col1.sp),
       lty = "solid", lwd = 7, cex= 1.5, bty = "n")


## But how do the betatraitcue differ?

hist(postCN$betaTraitxForce, col = col2.sp, main = "", xlim = c(-10,10), ylim = c(0,1000))
hist(postLMA$betaTraitxForce, col = col3.sp, main = "", add = T)
hist(postSSD$betaTraitxForce, col = col5.sp, main = "", add = T)
hist(postHt$betaTraitxForce, col = col1.sp, main = "", add = T)
hist(postDBH$betaTraitxForce, col = col4.sp, main = "", add = T)

hist(postCN$betaTraitxChill, col = col2.sp, main = "", xlim = c(-5,5), ylim = c(0,1200))
hist(postLMA$betaTraitxChill, col = col3.sp, main = "", add = T)
hist(postSSD$betaTraitxChill, col = col5.sp, main = "", add = T)
hist(postDBH$betaTraitxChill, col = col1.sp, main = "", add = T)
hist(postHt$betaTraitxChill, col = col4.sp, main = "", add = T)

hist(postCN$betaTraitxPhoto, col = col2.sp, main = "", xlim = c(-5,5), ylim = c(0,1500))
hist(postLMA$betaTraitxPhoto, col = col3.sp, main = "", add = T)
hist(postSSD$betaTraitxPhoto, col = col5.sp, main = "", add = T)
hist(postDBH$betaTraitxPhoto, col = col1.sp, main = "", add = T)
hist(postHt$betaTraitxPhoto, col = col4.sp, main = "", add = T)

legend("topright",legend = c(expression("Height"),
  expression("LMA"),
  expression("DBH"),
  expression("SSD"),
  expression("C:N")
),
  col = c(col4.sp, col3.sp, col1.sp, col5.sp, col2.sp),
  lty = "solid", lwd = 7, cex= 1.5, bty = "n")

##########################################################################
# ggplot() + 
#   stat_eye(data = longPhotoSiteInter, aes(x = site, y = photoSiteInter, fill = "cyan4"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
#   theme_classic() +   
#   theme(legend.position = "none") +
#   labs( x = "Site", y = "Photoperiod response", main = NA)+
#   scale_fill_manual(values = c("cyan4"))
# 
# +
#   geom_text(aes(label=species),hjust= 0.5, vjust= 1.5, show.legend = F) +
#   geom_errorbar(aes(ymin= bChill25, ymax = bChill75), width= 0) +
#   geom_errorbar(aes(xmin= bPhoto25, xmax = bPhoto75), width= 0) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),
#     legend.key=element_rect(fill="white")) # removed grey boxes around legends


# How do transect effects differ?
sumHt <- summary(mdlHt)$summary
muGrand = (sumHt[grep("mu_grand", rownames(sumHt)), 1])
b_trtSpHt = (sumHt[grep("b_muSp", rownames(sumHt)), 1])
a_trtSpHt = mean((sumHt[grep("mu_grand_sp", rownames(sumHt)), 1]))
b_tranEHt = sumHt[grep("b_tranE", rownames(sumHt)), 1]
b_tranlatHt = sumHt[grep("b_tranlat", rownames(sumHt)), 1]

b_phenoSpHt = (sumHt[grep("alphaPhenoSp", rownames(sumHt)), 1])
a_phenoSpHt = (sumHt[grep("muPhenoSp", rownames(sumHt)), 1])

a_chillSpHt = sumHt[grep("alphaChillSp", rownames(sumHt)), 1]
a_forceSpHt = sumHt[grep("alphaForceSp", rownames(sumHt)), 1]
a_photoSpHt = sumHt[grep("alphaPhotoSp", rownames(sumHt)), 1]

b_photoSpHt = sumHt[grep("muPhotoSp", rownames(sumHt)), 1]
b_forceSpHt = sumHt[grep("muForceSp", rownames(sumHt)), 1]
b_chillSpHt = sumHt[grep("muChillSp", rownames(sumHt)), 1]

bTrtChillHt = sumHt[grep("betaTraitxChill", rownames(sumHt)), 1]
bTrtForceHt = sumHt[grep("betaTraitxForce", rownames(sumHt)), 1]
bTrtPhotoHt = sumHt[grep("betaTraitxPhoto", rownames(sumHt)), 1]

a_trtsp5Ht <- vector()
for(i in 1:ncol(postHt$muSp)){  
  quantU <- round(quantile(postHt$muSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  a_trtsp5Ht <- rbind(a_trtsp5Ht, quantU)
}
colnames(a_trtsp5Ht) <- c("Int5","Int95","Int25","Int75")

b_tran5Ht <- round(quantile(postHt$b_tranE, c(0.05, 0.95, 0.25, 0.75)),1)
b_tranlat5Ht <- round(quantile(postHt$b_tranlat, c(0.05, 0.95, 0.25, 0.75)),1)

b_chill5 <- vector()
for(i in 1:ncol(postHt$betaChillSp)){
  quantU <- round(quantile(postHt$betaChillSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_chill5 <- rbind(b_chill5, quantU)
}
colnames(b_chill5) <- c("chill5","chill95","chill25","chill75")

b_force5 <- vector()
for(i in 1:ncol(postHt$betaForceSp)){
  quantU <- round(quantile(postHt$betaForceSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_force5 <- rbind(b_force5, quantU)
}
colnames(b_force5) <- c("force5","force95","force25","force75")

b_photo5 <- vector()
for(i in 1:ncol(postHt$betaPhotoSp)){
  quantU <- round(quantile(postHt$betaPhotoSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_photo5 <- rbind(b_photo5, quantU)
}
colnames(b_photo5) <- c("photo5","photo95","photo25","photo75")

## Simulate interaction with transect and latitude:

eData <- subset(trtPheno, transect == "1" )
wData <- subset(trtPheno, transect == "0" )

# Make the other parameters constant

lati <- seq(40, 60, by = 0.5)
latZ <- (lati-mean(lati,na.rm=TRUE))/(sd(lati,na.rm=TRUE)*2)
tranW <- -0.4406491
tranE <- 0.5669498

# plot first for west coast
ht_w = a_trtSpHt + b_tranEHt * tranW + b_tranlatHt * (tranW*latZ)
ht_e = a_trtSpHt + b_tranEHt * tranE + b_tranlatHt * (tranE*latZ)

par(mfrow = c(1,1))
plot(0, type = "n", xlim = c(25,60), ylim = c(-1,1),
  xlab = "Latitude",
  ylab = "Trait")
abline(lm(ht_w ~ lati), col = "darkslategray", lwd = 3, lty = 2)
abline(lm(ht_e ~lati), col = "darkslategray", lwd = 3, lty =1)

################################################
##lma

sumLMA <- summary(mdlLMA)$summary
muGrand = (sumLMA[grep("mu_grand", rownames(sumLMA)), 1])
b_trtSpLMA = (sumLMA[grep("b_muSp", rownames(sumLMA)), 1])
a_trtSpLMA = mean((sumLMA[grep("mu_grand_sp", rownames(sumLMA)), 1]))
b_tranELMA = sumLMA[grep("b_tranE", rownames(sumLMA)), 1]
b_tranlatLMA = sumLMA[grep("b_tranlat", rownames(sumLMA)), 1]

b_phenoSpLMA = (sumLMA[grep("alphaPhenoSp", rownames(sumLMA)), 1])
a_phenoSpLMA = (sumLMA[grep("muPhenoSp", rownames(sumLMA)), 1])

a_chillSpLMA = sumLMA[grep("alphaChillSp", rownames(sumLMA)), 1]
a_forceSpLMA = sumLMA[grep("alphaForceSp", rownames(sumLMA)), 1]
a_photoSpLMA = sumLMA[grep("alphaPhotoSp", rownames(sumLMA)), 1]

b_photoSpLMA = sumLMA[grep("muPhotoSp", rownames(sumLMA)), 1]
b_forceSpLMA = sumLMA[grep("muForceSp", rownames(sumLMA)), 1]
b_chillSpLMA = sumLMA[grep("muChillSp", rownames(sumLMA)), 1]

bTrtChillLMA = sumLMA[grep("betaTraitxChill", rownames(sumLMA)), 1]
bTrtForceLMA = sumLMA[grep("betaTraitxForce", rownames(sumLMA)), 1]
bTrtPhotoLMA = sumLMA[grep("betaTraitxPhoto", rownames(sumLMA)), 1]

a_trtsp5LMA <- vector()
for(i in 1:ncol(postLMA$muSp)){  
  quantU <- round(quantile(postLMA$muSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  a_trtsp5LMA <- rbind(a_trtsp5LMA, quantU)
}
colnames(a_trtsp5LMA) <- c("Int5","Int95","Int25","Int75")

b_tran5LMA <- round(quantile(postLMA$b_tranE, c(0.05, 0.95, 0.25, 0.75)),1)
b_tranlat5LMA <- round(quantile(postLMA$b_tranlat, c(0.05, 0.95, 0.25, 0.75)),1)

b_chill5 <- vector()
for(i in 1:ncol(postLMA$betaChillSp)){
  quantU <- round(quantile(postLMA$betaChillSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_chill5 <- rbind(b_chill5, quantU)
}
colnames(b_chill5) <- c("chill5","chill95","chill25","chill75")

b_force5 <- vector()
for(i in 1:ncol(postLMA$betaForceSp)){
  quantU <- round(quantile(postLMA$betaForceSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_force5 <- rbind(b_force5, quantU)
}
colnames(b_force5) <- c("force5","force95","force25","force75")

b_photo5 <- vector()
for(i in 1:ncol(postLMA$betaPhotoSp)){
  quantU <- round(quantile(postLMA$betaPhotoSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_photo5 <- rbind(b_photo5, quantU)
}
colnames(b_photo5) <- c("photo5","photo95","photo25","photo75")

## Simulate interaction with transect and latitude:

eData <- subset(trtPheno, transect == "1" )
wData <- subset(trtPheno, transect == "0" )

# Make the other parameters constant

# plot first for west coast
LMA_w = a_trtSpLMA + b_tranELMA * tranW + b_tranlatLMA * (tranW*latZ)
LMA_e = a_trtSpLMA + b_tranELMA * tranE + b_tranlatLMA * (tranE*latZ)

par(mfrow = c(1,1))
plot(0, type = "n", xlim = c(25,60), ylim = c(-1,1),
     xlab = "Latitude",
     ylab = "Trait")
abline(lm(LMA_w ~ lati), col = "darkslategray", lwd = 3, lty = 2)
abline(lm(LMA_e ~lati), col = "darkslategray", lwd = 3, lty =1)
#############################################
##dbh

sumDBH <- summary(mdlDBH)$summary
muGrand = (sumDBH[grep("mu_grand", rownames(sumDBH)), 1])
b_trtSpDBH = (sumDBH[grep("b_muSp", rownames(sumDBH)), 1])
a_trtSpDBH = mean((sumDBH[grep("mu_grand_sp", rownames(sumDBH)), 1]))
b_tranEDBH = sumDBH[grep("b_tranE", rownames(sumDBH)), 1]
b_tranlatDBH = sumDBH[grep("b_tranlat", rownames(sumDBH)), 1]

b_phenoSpDBH = (sumDBH[grep("alphaPhenoSp", rownames(sumDBH)), 1])
a_phenoSpDBH = (sumDBH[grep("muPhenoSp", rownames(sumDBH)), 1])

a_chillSpDBH = sumDBH[grep("alphaChillSp", rownames(sumDBH)), 1]
a_forceSpDBH = sumDBH[grep("alphaForceSp", rownames(sumDBH)), 1]
a_photoSpDBH = sumDBH[grep("alphaPhotoSp", rownames(sumDBH)), 1]

b_photoSpDBH = sumDBH[grep("muPhotoSp", rownames(sumDBH)), 1]
b_forceSpDBH = sumDBH[grep("muForceSp", rownames(sumDBH)), 1]
b_chillSpDBH = sumDBH[grep("muChillSp", rownames(sumDBH)), 1]

bTrtChillDBH = sumDBH[grep("betaTraitxChill", rownames(sumDBH)), 1]
bTrtForceDBH = sumDBH[grep("betaTraitxForce", rownames(sumDBH)), 1]
bTrtPhotoDBH = sumDBH[grep("betaTraitxPhoto", rownames(sumDBH)), 1]

a_trtsp5DBH <- vector()
for(i in 1:ncol(postDBH$muSp)){  
  quantU <- round(quantile(postDBH$muSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  a_trtsp5DBH <- rbind(a_trtsp5DBH, quantU)
}
colnames(a_trtsp5DBH) <- c("Int5","Int95","Int25","Int75")

b_tran5DBH <- round(quantile(postDBH$b_tranE, c(0.05, 0.95, 0.25, 0.75)),1)
b_tranlat5DBH <- round(quantile(postDBH$b_tranlat, c(0.05, 0.95, 0.25, 0.75)),1)

b_chill5 <- vector()
for(i in 1:ncol(postDBH$betaChillSp)){
  quantU <- round(quantile(postDBH$betaChillSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_chill5 <- rbind(b_chill5, quantU)
}
colnames(b_chill5) <- c("chill5","chill95","chill25","chill75")

b_force5 <- vector()
for(i in 1:ncol(postDBH$betaForceSp)){
  quantU <- round(quantile(postDBH$betaForceSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_force5 <- rbind(b_force5, quantU)
}
colnames(b_force5) <- c("force5","force95","force25","force75")

b_photo5 <- vector()
for(i in 1:ncol(postDBH$betaPhotoSp)){
  quantU <- round(quantile(postDBH$betaPhotoSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_photo5 <- rbind(b_photo5, quantU)
}
colnames(b_photo5) <- c("photo5","photo95","photo25","photo75")

## Simulate interaction with transect and latitude:

eData <- subset(trtPheno, transect == "1" )
wData <- subset(trtPheno, transect == "0" )


# plot first for west coast
DBH_w = a_trtSpDBH + b_tranEDBH * tranW + b_tranlatDBH * (tranW*latZ)
DBH_e = a_trtSpDBH + b_tranEDBH * tranE + b_tranlatDBH * (tranE*latZ)

par(mfrow = c(1,1))
plot(0, type = "n", xlim = c(25,60), ylim = c(-10,10),
     xlab = "Latitude",
     ylab = "Trait")
abline(lm(DBH_w ~ lati), col = "darkslategray", lwd = 3, lty = 2)
abline(lm(DBH_e ~lati), col = "darkslategray", lwd = 3, lty =1)


#############################################
sumSSD <- summary(mdlSSD)$summary
muGrand = (sumSSD[grep("mu_grand", rownames(sumSSD)), 1])
b_trtSpSSD = (sumSSD[grep("b_muSp", rownames(sumSSD)), 1])
a_trtSpSSD = mean((sumSSD[grep("mu_grand_sp", rownames(sumSSD)), 1]))
b_tranESSD = sumSSD[grep("b_tranE", rownames(sumSSD)), 1]
b_tranlatSSD = sumSSD[grep("b_tranlat", rownames(sumSSD)), 1]

b_phenoSpSSD = (sumSSD[grep("alphaPhenoSp", rownames(sumSSD)), 1])
a_phenoSpSSD = (sumSSD[grep("muPhenoSp", rownames(sumSSD)), 1])

a_chillSpSSD = sumSSD[grep("alphaChillSp", rownames(sumSSD)), 1]
a_forceSpSSD = sumSSD[grep("alphaForceSp", rownames(sumSSD)), 1]
a_photoSpSSD = sumSSD[grep("alphaPhotoSp", rownames(sumSSD)), 1]

b_photoSpSSD = sumSSD[grep("muPhotoSp", rownames(sumSSD)), 1]
b_forceSpSSD = sumSSD[grep("muForceSp", rownames(sumSSD)), 1]
b_chillSpSSD = sumSSD[grep("muChillSp", rownames(sumSSD)), 1]

bTrtChillSSD = sumSSD[grep("betaTraitxChill", rownames(sumSSD)), 1]
bTrtForceSSD = sumSSD[grep("betaTraitxForce", rownames(sumSSD)), 1]
bTrtPhotoSSD = sumSSD[grep("betaTraitxPhoto", rownames(sumSSD)), 1]

a_trtsp5SSD <- vector()
for(i in 1:ncol(postSSD$muSp)){  
  quantU <- round(quantile(postSSD$muSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  a_trtsp5SSD <- rbind(a_trtsp5SSD, quantU)
}
colnames(a_trtsp5SSD) <- c("Int5","Int95","Int25","Int75")

b_tran5SSD <- round(quantile(postSSD$b_tranE, c(0.05, 0.95, 0.25, 0.75)),1)
b_tranlat5SSD <- round(quantile(postSSD$b_tranlat, c(0.05, 0.95, 0.25, 0.75)),1)

b_chill5 <- vector()
for(i in 1:ncol(postSSD$betaChillSp)){
  quantU <- round(quantile(postSSD$betaChillSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_chill5 <- rbind(b_chill5, quantU)
}
colnames(b_chill5) <- c("chill5","chill95","chill25","chill75")

b_force5 <- vector()
for(i in 1:ncol(postSSD$betaForceSp)){
  quantU <- round(quantile(postSSD$betaForceSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_force5 <- rbind(b_force5, quantU)
}
colnames(b_force5) <- c("force5","force95","force25","force75")

b_photo5 <- vector()
for(i in 1:ncol(postSSD$betaPhotoSp)){
  quantU <- round(quantile(postSSD$betaPhotoSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_photo5 <- rbind(b_photo5, quantU)
}
colnames(b_photo5) <- c("photo5","photo95","photo25","photo75")

## Simulate interaction with transect and latitude:

eData <- subset(trtPheno, transect == "1" )
wData <- subset(trtPheno, transect == "0" )


# plot first for west coast
SSD_w = a_trtSpSSD + b_tranESSD * tranW + b_tranlatSSD * (tranW*latZ)
SSD_e = a_trtSpSSD + b_tranESSD * tranE + b_tranlatSSD * (tranE*latZ)

par(mfrow = c(1,1))
plot(0, type = "n", xlim = c(25,60), ylim = c(-100,100),
     xlab = "Latitude",
     ylab = "Trait")
abline(lm(SSD_w ~ lati), col = "darkslategray", lwd = 3, lty = 2)
abline(lm(SSD_e ~lati), col = "darkslategray", lwd = 3, lty =1)

#############################################
sumCN <- summary(mdl)$summary
muGrand = (sumCN[grep("mu_grand", rownames(sumCN)), 1])
b_trtSpCN = (sumCN[grep("b_muSp", rownames(sumCN)), 1])
a_trtSpCN = mean((sumCN[grep("mu_grand_sp", rownames(sumCN)), 1]))
b_tranECN = sumCN[grep("b_tranE", rownames(sumCN)), 1]
b_tranlatCN = sumCN[grep("b_tranlat", rownames(sumCN)), 1]

b_phenoSpCN = (sumCN[grep("alphaPhenoSp", rownames(sumCN)), 1])
a_phenoSpCN = (sumCN[grep("muPhenoSp", rownames(sumCN)), 1])

a_chillSpCN = sumCN[grep("alphaChillSp", rownames(sumCN)), 1]
a_forceSpCN = sumCN[grep("alphaForceSp", rownames(sumCN)), 1]
a_photoSpCN = sumCN[grep("alphaPhotoSp", rownames(sumCN)), 1]

b_photoSpCN = sumCN[grep("muPhotoSp", rownames(sumCN)), 1]
b_forceSpCN = sumCN[grep("muForceSp", rownames(sumCN)), 1]
b_chillSpCN = sumCN[grep("muChillSp", rownames(sumCN)), 1]

bTrtChillCN = sumCN[grep("betaTraitxChill", rownames(sumCN)), 1]
bTrtForceCN = sumCN[grep("betaTraitxForce", rownames(sumCN)), 1]
bTrtPhotoCN = sumCN[grep("betaTraitxPhoto", rownames(sumCN)), 1]

a_trtsp5CN <- vector()
for(i in 1:ncol(postCN$muSp)){  
  quantU <- round(quantile(postCN$muSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  a_trtsp5CN <- rbind(a_trtsp5CN, quantU)
}
colnames(a_trtsp5CN) <- c("Int5","Int95","Int25","Int75")

b_tran5CN <- round(quantile(postCN$b_tranE, c(0.05, 0.95, 0.25, 0.75)),1)
b_tranlat5CN <- round(quantile(postCN$b_tranlat, c(0.05, 0.95, 0.25, 0.75)),1)

b_chill5 <- vector()
for(i in 1:ncol(postCN$betaChillSp)){
  quantU <- round(quantile(postCN$betaChillSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_chill5 <- rbind(b_chill5, quantU)
}
colnames(b_chill5) <- c("chill5","chill95","chill25","chill75")

b_force5 <- vector()
for(i in 1:ncol(postCN$betaForceSp)){
  quantU <- round(quantile(postCN$betaForceSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_force5 <- rbind(b_force5, quantU)
}
colnames(b_force5) <- c("force5","force95","force25","force75")

b_photo5 <- vector()
for(i in 1:ncol(postCN$betaPhotoSp)){
  quantU <- round(quantile(postCN$betaPhotoSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_photo5 <- rbind(b_photo5, quantU)
}
colnames(b_photo5) <- c("photo5","photo95","photo25","photo75")

## Simulate interaction with transect and latitude:

eData <- subset(trtPheno, transect == "1" )
wData <- subset(trtPheno, transect == "0" )


# plot first for west coast
CN_w = a_trtSpCN + b_tranECN * tranW + b_tranlatCN * (tranW*latZ)
CN_e = a_trtSpCN + b_tranECN * tranE + b_tranlatCN * (tranE*latZ)

par(mfrow = c(1,1))
plot(0, type = "n", xlim = c(25,60), ylim = c(-10,10),
     xlab = "Latitude",
     ylab = "Trait")
abline(lm(CN_w ~ lati), col = "darkslategray", lwd = 3, lty = 2)
abline(lm(CN_e ~lati), col = "darkslategray", lwd = 3, lty =1)
#################################
pdf("figures/transectIntrxnZ25.pdf", width =13, height =3)
par(mfrow = c(1,5))
plot(0, type = "n", xlim = c(40,55), ylim = c(-2,2),
     xlab = "Latitude",
     ylab = "Height", cex.lab = 1.3)

#abline(lm(ht_w ~ lati), col = "darkslategray", lwd = 3, lty = 2)
abline(lm(ht_e ~lati), col = "darkslategray4", lwd = 3, lty =1)
abline(lm(ht_w ~lati), col = "darkslategray4", lwd = 3, lty =2)
text(40.8, 2, label = "a)", cex = 1.25)

plot(0, type = "n", xlim = c(40,55), ylim = c(-2,2),
     xlab = "Latitude",
     ylab = "Diameter at breast height", cex.lab = 1.3)
abline(lm(DBH_e ~lati), col = "goldenrod", lwd = 3, lty =1)
abline(lm(DBH_w ~lati), col = "goldenrod", lwd = 3, lty =2)
text(40.8, 2, label = "b)", cex = 1.25)

plot(0, type = "n", xlim = c(40,55), ylim = c(-2,2),
     xlab = "Latitude",
     ylab = "Leaf mass area", cex.lab = 1.3)
abline(lm(LMA_e ~lati), col = "darkolivegreen", lwd = 3, lty =1)
abline(lm(LMA_w ~lati), col = "darkolivegreen", lwd = 3, lty =2)
text(40.8,2, label = "c)", cex = 1.25)


plot(0, type = "n", xlim = c(40,55), ylim = c(-2,2),
     xlab = "Latitude",
     ylab = "Stem specific density", cex.lab = 1.3)
abline(lm(SSD_e ~lati), col = "maroon", lwd = 3, lty =1)
abline(lm(SSD_w ~lati), col = "maroon", lwd = 3, lty =2)
text(40.8,2, label = "d)", cex = 1.25)

plot(0, type = "n", xlim = c(40,55), ylim = c(-2,2),
     xlab = "Latitude",
     ylab = "Carbon:Nitrogen", cex.lab = 1.3)
abline(lm(CN_e ~lati), col = "purple4", lwd = 3, lty =1)
abline(lm(CN_w ~lati), col = "purple4", lwd = 3, lty =2)
text(40.8,2, label = "e)", cex = 1.25)

legend("topright",legend = c(expression("Western"),
                             expression("Eastern")
),
col = c("black", "black"),
lty = c(2,1), lwd = 1, cex= 1.25, bty = "n")

dev.off()
# #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# I think what we want is a loop that goes through each iteration of the posteriors and calculates the bb, but using 20 for forcing, 12 for photoperiod, 75 (75/10 when rescaled), and smithers to start
# 

#If we are using the old model, we will use the z-scored values for the parameters
photo <- -0.5033863 #8 h photo
siteSM <- 0
force <- -0.3568628 #5/15 C trt
chill <- -0.3546922 # low chill

m <- matrix(nrow = 1000, ncol = 47)

for(sp in 1:47){
  for (it in 1:nrow(m)){
    m[it,sp] <- post$a_sp[it,sp]+ post$b_site2[it] * siteSM + post$b_site3[it] * siteSM + post$b_site4[it] * siteSM + 
      post$b_warm[it,sp] * force + post$b_photo[it, sp] * photo + post$b_chill[it,sp] * chill +
      post$b_inter_wp[it,sp] * (force*photo) + post$b_inter_wc1[it,sp] * (force*chill) + post$b_inter_pc1[it,sp] * (photo*chill) +
      post$b_inter_s2c1[it,sp] * (chill*siteSM) + post$b_inter_ws2[it,sp] * (force*siteSM) + post$b_inter_ps2[it,sp] * (photo*siteSM) +
      post$b_inter_s3c1[it,sp] * (chill*siteSM) + post$b_inter_ws3[it,sp] * (force*siteSM) + post$b_inter_ps3[it,sp] * (photo*siteSM) +
      post$b_inter_s4c1[it,sp] * (chill*siteSM) + post$b_inter_ws4[it,sp] * (force*siteSM) + post$b_inter_ps4[it,sp] * (photo*siteSM)
  }
}

############ SHRUB VS TREE ##############################
spInfo <- read.csv("..//input/species_ring.csv")
spInfo <- spInfo[, 1:5]
head(spInfo)

## Start with height:
fit <- rstan::extract(mdlHt)

chillB <- data.frame(fit$betaChillSp)

colnames(chillB) <- sort(spInfo$species.name)

longChill <- melt(chillB)
colnames(longChill) <- c("species.name","betaCueSp")
longChill$cue <- "Chilling"

head(longChill)

###################################################################
# Photoperiod:
photoB <- data.frame(fit$betaPhotoSp)

colnames(photoB) <- sort(spInfo$species.name)

longPhoto <- melt(photoB)
colnames(longPhoto) <- c("species.name","betaCueSp")
longPhoto$cue <- "Photoperiod"
head(longPhoto)

###################################################################
# Forcing:
forceB <- data.frame(fit$betaForceSp)

colnames(forceB) <- sort(spInfo$species.name)

longForce <- melt(forceB)
colnames(longForce) <- c("species.name","betaCueSp")
longForce$cue <- "Forcing"

head(longForce)

longCues <- rbind(longForce, longChill, longPhoto)
longCues <- merge(longCues, spInfo, by = "species.name")

ring <- subset(longCues, ringType == "Ring ")
diffuse <- subset(longCues, ringType == "Diffuse")
difRing <- subset(longCues, ringType == "Diffuse/semi-ring")
semi <- subset(longCues, ringType == "Semi-ring" )

longCuesRing <- subset(longCues, ringType != "")

ggplot() + 
  geom_violin(data = longCuesRing, aes(x = as.factor(cue), y = betaCueSp, col = factor(cue))) +
  facet_grid(col = vars(ringType), scales = "free_y") + theme(strip.background = element_blank(), strip.text.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.title=element_blank()) +
   ylab ("Cue response") + xlab ("Cue")


ggplot() + 
  geom_violin(data = longCues, aes(x = as.factor(cue), y = betaCueSp, col = factor(cue))) +
  facet_grid(col = vars(type), scales = "free_y") + theme(strip.background = element_blank(), strip.text.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.title=element_blank()) +
  ylab ("Cue response") + xlab ("Cue")

  facet_grid(col = vars(type))
              
              +
  stat_pointinterval(data = longest, aes(x = as.factor(cue), y = value, fill = factor(site, level = siteOrder)), .width = c(.5, .95) ,position = position_dodge(0.9)) +
  theme_classic() +   
  theme(legend.position = "right", 
        legend.title = element_blank(),
        axis.text.x = element_text( size= 16),
        axis.text.y = element_text( size= 12),
        axis.title=element_text(size = 14)) +
  labs( x = "Treatment cue", y = "Cue response (days/standardized unit)", main = NA) +
  scale_color_manual(values = c("Smithers" = "deepskyblue3",
                                "Manning Park" = "palegreen4", 
                                "St. Hippolyte"="darkorchid3", 
                                "Harvard Forest" = "tomato3"))+
  scale_fill_manual(values = c("Smithers" = "deepskyblue3",
                               "Manning Park" = "palegreen4", 
                               "St. Hippolyte"="darkorchid3", 
                               "Harvard Forest" = "tomato3"))
  
  
  # compare points:
  # chilling

longChill <- merge(longChill, spInfo, by = "species.name")
  
meanChill <- aggregate(longChill[c("betaCueSp")], longChill[c("ringType")], FUN = mean)
names(meanChill) <- c( "ringType", "betaCueSp")

meanError75 <- aggregate(longChill[c("betaCueSp")], longChill[c( "ringType")], FUN = function(i) quantile(i, probs = 0.75, na.rm = T))
names(meanError75) <- c( "ringType", "error75")

meanError25 <- aggregate(longChill[c("betaCueSp")], longChill[c( "ringType")], FUN = function(i) quantile(i, probs = 0.25, na.rm = T))
names(meanError25) <- c( "ringType", "error25")

meanError95 <- aggregate(longChill[c("betaCueSp")], longChill[c("ringType")], FUN = function(i) quantile(i, probs = 0.95, na.rm = T))
names(meanError95) <- c( "ringType", "error95")

meanError05 <- aggregate(longChill[c("betaCueSp")], longChill[c( "ringType")], FUN = function(i) quantile(i, probs = 0.05, na.rm = T))
names(meanError05) <- c( "ringType", "error05")


meanChill2 <- merge(meanChill, meanError05, by = c( "ringType"))
meanChill2 <- merge(meanChill2, meanError25, by = c("ringType"))
meanChill2 <- merge(meanChill2, meanError75, by = c("ringType"))
meanChill2 <- merge(meanChill2, meanError95, by = c("ringType"))

# require(dplyr)
# longCues %>% group_by(species.name) %>%
#   summarise(p90 = quantile(betaCueSp, probs=0.9, na.rm=TRUE))

ringChill <- ggplot(meanChill2,aes(y= betaCueSp, x = ringType), size = 7) +
  geom_point(size = 7,  color = "cyan4") +
  geom_errorbar(aes(ymin= error05, ymax = error95,xmin= ringType, xmax = ringType), width= 0, size = 0.5, color = "cyan4") +
  geom_errorbar(aes(ymin= error25, ymax = error75,xmin= ringType, xmax = ringType), width= 0, size = 1.5, color = "cyan4") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + labs( x = "Ring Type", y = "Chilling response (days/standardized unit)", main = NA) +
   theme(legend.title = element_blank()) 
# theme(axis.text.x = element_text( size=17,angle = 78,  hjust=1),
#       axis.text.y=element_text(size = 15),
#       axis.title=element_text(size=  17),
#       legend.position = "none") +
#   labs( x = "Ring Type", y = "Chilling response (days/standardized unit)", main = NA) +
#   theme(legend.title = element_blank()) 
  
longforce <- merge(longForce, spInfo, by = "species.name")

meanforce <- aggregate(longforce[c("betaCueSp")], longforce[c("ringType")], FUN = mean)
names(meanforce) <- c( "ringType", "betaCueSp")

meanError75 <- aggregate(longforce[c("betaCueSp")], longforce[c( "ringType")], FUN = function(i) quantile(i, probs = 0.75, na.rm = T))
names(meanError75) <- c( "ringType", "error75")

meanError25 <- aggregate(longforce[c("betaCueSp")], longforce[c( "ringType")], FUN = function(i) quantile(i, probs = 0.25, na.rm = T))
names(meanError25) <- c( "ringType", "error25")

meanError95 <- aggregate(longforce[c("betaCueSp")], longforce[c("ringType")], FUN = function(i) quantile(i, probs = 0.95, na.rm = T))
names(meanError95) <- c( "ringType", "error95")

meanError05 <- aggregate(longforce[c("betaCueSp")], longforce[c( "ringType")], FUN = function(i) quantile(i, probs = 0.05, na.rm = T))
names(meanError05) <- c( "ringType", "error05")


meanforce2 <- merge(meanforce, meanError05, by = c( "ringType"))
meanforce2 <- merge(meanforce2, meanError25, by = c("ringType"))
meanforce2 <- merge(meanforce2, meanError75, by = c("ringType"))
meanforce2 <- merge(meanforce2, meanError95, by = c("ringType"))

# require(dplyr)
# longCues %>% group_by(species.name) %>%
#   summarise(p90 = quantile(betaCueSp, probs=0.9, na.rm=TRUE))

ringForce <- ggplot(meanforce2,aes(y= betaCueSp, x = ringType), size = 7) +
  geom_point(size = 7, color = "goldenrod") +
  geom_errorbar(aes(ymin= error05, ymax = error95,xmin= ringType, xmax = ringType), width= 0, size = 0.5, color = "goldenrod") +
  geom_errorbar(aes(ymin= error25, ymax = error75,xmin= ringType, xmax = ringType), width= 0, size = 1.5, color = "goldenrod") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  # theme(axis.text.x = element_text( size=17,angle = 78,  hjust=1),
  #       axis.text.y=element_text(size = 15),
  #       axis.title=element_text(size=  17),
  #       legend.position = "none") +
  labs( x = "Ring Type", y = "Forceing response (days/standardized unit)", main = NA) +
  theme(legend.title = element_blank())

longphoto <- merge(longPhoto, spInfo, by = "species.name")

meanphoto <- aggregate(longphoto[c("betaCueSp")], longphoto[c("ringType")], FUN = mean)
names(meanphoto) <- c( "ringType", "betaCueSp")

meanError75 <- aggregate(longphoto[c("betaCueSp")], longphoto[c( "ringType")], FUN = function(i) quantile(i, probs = 0.75, na.rm = T))
names(meanError75) <- c( "ringType", "error75")

meanError25 <- aggregate(longphoto[c("betaCueSp")], longphoto[c( "ringType")], FUN = function(i) quantile(i, probs = 0.25, na.rm = T))
names(meanError25) <- c( "ringType", "error25")

meanError95 <- aggregate(longphoto[c("betaCueSp")], longphoto[c("ringType")], FUN = function(i) quantile(i, probs = 0.95, na.rm = T))
names(meanError95) <- c( "ringType", "error95")

meanError05 <- aggregate(longphoto[c("betaCueSp")], longphoto[c( "ringType")], FUN = function(i) quantile(i, probs = 0.05, na.rm = T))
names(meanError05) <- c( "ringType", "error05")


meanphoto2 <- merge(meanphoto, meanError05, by = c( "ringType"))
meanphoto2 <- merge(meanphoto2, meanError25, by = c("ringType"))
meanphoto2 <- merge(meanphoto2, meanError75, by = c("ringType"))
meanphoto2 <- merge(meanphoto2, meanError95, by = c("ringType"))

# require(dplyr)
# longCues %>% group_by(species.name) %>%
#   summarise(p90 = quantile(betaCueSp, probs=0.9, na.rm=TRUE))

ringPhoto <- ggplot(meanphoto2,aes(y= betaCueSp, x = ringType), size = 7) +
  geom_point(size = 7, color = "maroon") +
  geom_errorbar(aes(ymin= error05, ymax = error95,xmin= ringType, xmax = ringType), width= 0, size = 0.5, color = "maroon") +
  geom_errorbar(aes(ymin= error25, ymax = error75,xmin= ringType, xmax = ringType), width= 0, size = 1.5, color = "maroon") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  # theme(axis.text.x = element_text( size=17,angle = 78,  hjust=1),
  #       axis.text.y=element_text(size = 15),
  #       axis.title=element_text(size=  17),
  #       legend.position = "none") +
  labs( x = "Ring Type", y = "photoing response (days/standardized unit)", main = NA) +
  theme(legend.title = element_blank()) 

pdf("figures/ringPorosity.pdf", width = 12, height = 4)
plot_grid(ringChill, ringForce, ringPhoto, nrow = 1, ncol = 3, align = "v")
dev.off()
