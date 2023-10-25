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

spInfo <- read.csv("input/species_list.csv")
pheno <- read.csv("input/phenoData.csv")
trtPheno <- read.csv("input/trtData.csv")

load("output/heightDummyIntGrand.Rdata")
sumerht <- summary(mdl)$summary
mdlHt <- mdl
postHt <- rstan::extract(mdlHt)

load("output/lmaDummyInt.Rdata")
postLMA <- rstan::extract(mdlLMA)

load("output/dbhDummyInt.Rdata")
mdlDBH <- mdl
postDBH <- rstan::extract(mdlDBH)

load("output/ssdDummyInt.Rdata")
mdlCN <- mdl
postCN <- rstan::extract(mdlCN)

load("output/cnDummyInt.Rdata")
mdlSSD <- mdl
postSSD <- rstan::extract(mdlSSD)

###### Compare the spp and site level effects across traits ###########
# col1 <- rgb(204 / 255, 102 / 255, 119 / 255, alpha = 0.8)
# col2 <- rgb(68 / 255, 170 / 255, 153 / 255, alpha = 0.6)

col1.sp <-c( rgb(204 / 255, 105 / 255, 112 / 255, alpha = 0.5)) # red
col2.sp <- c( rgb(205 / 255, 122 / 255, 0 / 255, alpha = 0.5)) # yellow
col3.sp <-c( rgb(9/ 255, 168 / 255, 82 / 255, alpha = 0.6)) # green
col4.sp <- c( rgb(34 / 255, 166 / 255, 167 / 255, alpha = 0.5)) # blue
col5.sp <- c( rgb(141 / 255, 34 / 255, 171 / 255, alpha = 0.5)) # purple


hist(postHt$b_tranE, col = col2.sp, main = "")
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
for(i in 1:ncol(postHt$b_muSp)){
  quantU <- round(quantile(postHt$b_muSp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  a_trtsp5Ht <- rbind(a_trtsp5Ht, quantU)
}
colnames(a_trtsp5Ht) <- c("Int5","Int95","Int25","Int75")

b_tran5Ht <- round(quantile(postHt$b_tranE, c(0.05, 0.95, 0.25, 0.75)),1)
b_tranlat5Ht <- round(quantile(postHt$b_tranlat, c(0.05, 0.95, 0.25, 0.75)),1)

b_chill5 <- vector()
for(i in 1:ncol(post$b_chill1)){
  quantU <- round(quantile(post$b_chill1[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_chill5 <- rbind(b_chill5, quantU)
}
colnames(b_chill5) <- c("chill5","chill95","chill25","chill75")

b_force5 <- vector()
for(i in 1:ncol(post$b_warm)){
  quantU <- round(quantile(post$b_warm[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_force5 <- rbind(b_force5, quantU)
}
colnames(b_force5) <- c("force5","force95","force25","force75")

b_photo5 <- vector()
for(i in 1:ncol(post$b_photo)){
  quantU <- round(quantile(post$b_photo[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_photo5 <- rbind(b_photo5, quantU)
}
colnames(b_photo5) <- c("photo5","photo95","photo25","photo75")

## Simulate interaction with transect and latitude:

eData <- subset(trtPheno, transect == "1" )
wData <- subset(trtPheno, transect == "0" )

# Make the other parameters constant

lati <- seq(42, 55, by = 0.5)
tranW <- 0
tranE <- 1

# plot first for west coast
ht_w = a_trtSpHt + b_tranEHt * tranW + b_tranlatHt * (tranW*lati)
ht_e = a_trtSpHt + b_tranEHt * tranE + b_tranlatHt * (tranE*lati)

par(mfrow = c(1,2))
plot(0, type = "n", xlim = c(40,60), ylim = c(0,10),
  xlab = "Latitude",
  ylab = "Trait")
abline(lm(ht_w ~ lati), col = "darkslategray", lwd = 3, lty = 2)
abline(lm(ht_e ~lati), col = "darkslategray", lwd = 3, lty =1)

##dbh
sumDBH <- summary(mdlDBH)$summary
b_trtSpDBH = (sumDBH [grep("b_muSp", rownames(sumDBH )), 1])
a_trtSpDBH = (sumDBH [grep("muSp", rownames(sumDBH )), 1]); a_trtSpDBH <- a_trtSpDBH[48]
b_tranEDBH = sumDBH [grep("b_tranE", rownames(sumDBH )), 1]
b_tranlatDBH = sumDBH [grep("b_tranlat", rownames(sumDBH )), 1]

DBH_w = a_trtSpDBH + b_tranEDBH * tranW + b_tranlatDBH * (tranW*lati)
DBH_e = a_trtSpDBH + b_tranEDBH * tranE + b_tranlatDBH * (tranE*lati)


abline(lm(DBH_w ~ lati), col = "goldenrod", lwd = 3, lty = 2)
abline(lm(DBH_e ~lati), col = "goldenrod", lwd = 3, lty = 1)

legend("topright",legend = c(expression("Eastern"),
  expression("Western"),
  "Height",
  "DBH"),
  col = c("black", "black","darkslategray","goldenrod"),
  lty = c(1,2, 1,1), bty = "n", lwd =2)

##lma
sumLMA <- summary(mdlLMA)$summary
b_trtSpLMA = (sumLMA [grep("b_muSp", rownames(sumLMA )), 1])
a_trtSpLMA = (sumLMA [grep("muSp", rownames(sumLMA )), 1]); a_trtSpLMA <- a_trtSpLMA[48]
b_tranELMA = sumLMA [grep("b_tranE", rownames(sumLMA )), 1]
b_tranlatLMA = sumLMA [grep("b_tranlat", rownames(sumLMA )), 1]

lma_w = a_trtSpLMA + b_tranELMA * tranW + b_tranlatLMA * (tranW*lati)
lma_e = a_trtSpLMA + b_tranELMA * tranE + b_tranlatLMA * (tranE*lati)


plot(0, type = "n", xlim = c(40,60), ylim = c(0,1),
  xlab = "Latitude",
  ylab = "Trait")
abline(lm(lma_w ~ lati), col = "darkred", lwd = 3, lty = 2)
abline(lm(lma_e ~lati), col = "darkred", lwd = 3, lty = 1)


# SSD
sumSSD <- summary(mdlSSD)$summary
b_trtSpSSD = (sumSSD [grep("b_muSp", rownames(sumSSD )), 1])
a_trtSpSSD = (sumSSD [grep("muSp", rownames(sumSSD )), 1]); a_trtSpSSD <- a_trtSpSSD[48]
b_tranESSD = sumSSD [grep("b_tranE", rownames(sumSSD )), 1]
b_tranlatSSD = sumSSD [grep("b_tranlat", rownames(sumSSD )), 1]

SSD_w = a_trtSpSSD + b_tranESSD * tranW + b_tranlatSSD * (tranW*lati)
SSD_e = a_trtSpSSD + b_tranESSD * tranE + b_tranlatSSD * (tranE*lati)


abline(lm(SSD_w ~ lati), col = "orange2", lwd = 3, lty = 2)
abline(lm(SSD_e ~lati), col = "orange2", lwd = 3, lty = 1)

# CN
sumCN <- summary(mdlCN)$summary
b_trtSpCN = (sumCN [grep("b_muSp", rownames(sumCN )), 1])
a_trtSpCN = (sumCN [grep("muSp", rownames(sumCN )), 1]); a_trtSpCN <- a_trtSpCN[48]
b_tranECN = sumCN [grep("b_tranE", rownames(sumCN )), 1]
b_tranlatCN = sumCN [grep("b_tranlat", rownames(sumCN )), 1]

CN_w = a_trtSpCN + b_tranECN * tranW + b_tranlatCN * (tranW*lati)
CN_e = a_trtSpCN + b_tranECN * tranE + b_tranlatCN * (tranE*lati)

abline(lm(CN_w ~ lati), col = "darkolivegreen", lwd = 3, lty = 2)
abline(lm(CN_e ~lati), col = "darkolivegreen", lwd = 3, lty = 1)

legend("topright",legend = c(expression("eastern"),
  expression("western"),
  "LMA","SSD","CN"),
  col = c("black", "black","maroon","orange2","darkolivegreen"),
  lty = c(1,2, 1,1,1), bty = "n", lwd =2)

###########################################################
pdf("figures/latTranInteraction.pdf", width = 10, height =5)
par(mfrow = c(1,2))
plot(0, type = "n", xlim = c(40,60), ylim = c(0,10),
  xlab = "Latitude",
  ylab = "Trait")
abline(lm(ht_w ~ lati), col = "darkslategray", lwd = 3, lty = 2)
abline(lm(ht_e ~lati), col = "darkslategray", lwd = 3, lty =1)

abline(lm(DBH_w ~ lati), col = "goldenrod", lwd = 3, lty = 2)
abline(lm(DBH_e ~lati), col = "goldenrod", lwd = 3, lty = 1)

legend("topright",legend = c(expression("Eastern"),
  expression("Western"),
  "Height",
  "DBH"),
  col = c("black", "black","darkslategray","goldenrod"),
  lty = c(1,2, 1,1), bty = "n", lwd =2)

plot(0, type = "n", xlim = c(40,60), ylim = c(0,1),
  xlab = "Latitude",
  ylab = "Trait")
abline(lm(lma_w ~ lati), col = "darkred", lwd = 3, lty = 2)
abline(lm(lma_e ~lati), col = "darkred", lwd = 3, lty = 1)

abline(lm(SSD_w ~ lati), col = "orange2", lwd = 3, lty = 2)
abline(lm(SSD_e ~lati), col = "orange2", lwd = 3, lty = 1)

abline(lm(CN_w ~ lati), col = "darkolivegreen", lwd = 3, lty = 2)
abline(lm(CN_e ~lati), col = "darkolivegreen", lwd = 3, lty = 1)

legend("topright",legend = c(expression("eastern"),
  expression("western"),
  "LMA","SSD","CN"),
  col = c("black", "black","maroon","orange2","darkolivegreen"),
  lty = c(1,2, 1,1,1), bty = "n", lwd =2)
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


# # Compare the musites for the 4 sites
# muSiteHtW <- data.frame(postHtW$musite)
# site <- c("Alex Fraser", "Kamloops","Manning Park", "Smithers")
# names(muSiteHtW) <- site
# 
# longSiteHtW <- melt(muSiteHtW)
# names(longSiteHtW) <- c("Population", "value")
# 
# htSiteW <- ggplot() + 
#   stat_eye(data = longSiteHtW, aes(x = Population, y = value, fill = Population), .width = c(.90, .5), cex = 0.75) +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(legend.position = "none",
#         axis.text.x = element_text(size=12, angle = 78,hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Populations", y = "Population effect", main = "Height") +
#     scale_fill_manual(values = c("cyan4","#cc6a70ff", "#f9b641ff", "mediumpurple2")) + 
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 5, label = "a)    Height", cex = 5) 
# 
# # LMA
# muSiteLMAW <- data.frame(postLMAW$musite)
# site <- c("Alex Fraser", "Kamloops","Manning Park", "Smithers")
# names(muSiteLMAW) <- site
# 
# longSiteLMAW <- melt(muSiteLMAW)
# names(longSiteLMAW) <- c("Population", "value")
# 
# LMASiteW <- ggplot() + 
#   stat_eye(data = longSiteLMAW, aes(x = Population, y = value, fill = Population), .width = c(.90, .5), cex = 0.75) +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(legend.position = "none",
#         axis.text.x = element_text(size=12, angle = 78,hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Populations", y = "Population effect", main = "LMA") +
#   scale_fill_manual(values = c("cyan4","#cc6a70ff", "#f9b641ff", "mediumpurple2")) + 
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 0.1, label = "b)    LMA", cex =5) 
# 
# # DBH
# muSiteDBHW <- data.frame(postDbhW$musite)
# site <- c("Alex Fraser", "Kamloops","Manning Park", "Smithers")
# names(muSiteDBHW) <- site
# 
# longSiteDBHW <- melt(muSiteDBHW)
# names(longSiteDBHW) <- c("Population", "value")
# 
# DBHSiteW <- ggplot() + 
#   stat_eye(data = longSiteDBHW, aes(x = Population, y = value, fill = Population), .width = c(.90, .5), cex = 0.75) +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(legend.position = "none",
#         axis.text.x = element_text(size=12, angle = 78,hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Populations", y = "Population effect", main = "DBH") +
#   scale_fill_manual(values = c("cyan4","#cc6a70ff", "#f9b641ff", "mediumpurple2")) + 
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 7, label = "c)    DBH", cex =5) 
# 
# # C:N
# muSiteCNW <- data.frame(postCNW$musite)
# site <- c("Alex Fraser", "Kamloops","Manning Park", "Smithers")
# names(muSiteCNW) <- site
# 
# longSiteCNW <- melt(muSiteCNW)
# names(longSiteCNW) <- c("Population", "value")
# 
# CNSiteW <- ggplot() + 
#   stat_eye(data = longSiteCNW, aes(x = Population, y = value, fill = Population), .width = c(.90, .5), cex = 0.75) +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(axis.text.x = element_text(size=12, angle = 78, hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Populations", y = "Population effect", main = "Carbon:Nitrogen") +
#   scale_fill_manual(values = c("cyan4","#cc6a70ff", "#f9b641ff", "mediumpurple2")) + 
#   theme(legend.title = element_blank()) +  annotate("text", x = 2, y = 10, label = "d)    Carbon:Nitrogen", cex =5) 
# 
# pdf("figures/traitMuSiteW.pdf", width = 20, height =5)
# plot_grid(htSiteW,LMASiteW,DBHSiteW,CNSiteW, ncol = 4, align = "v")
# dev.off()
# 
# pdf("figures/sigmaSpSite.pdf", width = 12, height = 4)
# par(mfrow = c(1, 5))
# hist(postHtW$sigma_site, main = "Height", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 30, ylim = c(0, 4000))
# hist(postHtW$sigma_sp, add = T, col = col2)
# text(0.2,4000, label = "a)", cex = 2)
# 
# hist(postLMAW$sigma_site, main = "Leaf Mass Area", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 50, xlim = c(0,0.2), ylim = c(0, 6000))
# hist(postLMAW$sigma_sp, add = T, col = col2)
# text(0.002,6000, label = "b)", cex = 2)
# 
# hist(postDbhW$sigma_site, main = "Diameter at Breast Height", xlab = "Variation", ylab = "Frequency", col = col1, ylim = c(0, 6000))
# hist(postDbhW$sigma_sp, add = T, col = col2)
# text(0.2,6000, label = "c)", cex = 2)
# 
# hist(postCNW$sigma_site, main = "Carbon:Nitrogen", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 30, ylim = c(0, 4000), xlim = c(0,12))
# hist(postCNW$sigma_sp, add = T, col = col2)
# text(0.2,4000, label = "d)", cex = 2)
# 
# hist(postSsdW$sigma_site, main = "Stem Specific Density", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 30, ylim = c(0, 4000), xlim = c(0,6))
# hist(postCNW$sigma_sp, add = T, col = col2, breaks = 50)
# text(0.1,4000, label = "e)", cex = 2)
# 
# legend("topright",legend = c(expression("Species variation"),
#                              expression("Study variation")),
#        col = c(col2,col1),
#        inset = 0.02, cex = 1, bty = "n", lwd = 3)
# 
# dev.off()
# 
# 
# 
# load("output/rda/ht_eastern_Feb21.Rda")
# Model <- mdl.ht4
# 
# load("output/rda/lma_eastern_Feb21.Rda")
# Model <- mdl.lma4
# 
# load("output/rda/dbh_eastern_Feb21.Rda")
# Model <- mdl.dbh4
# 
# load("output/rda/cn_eastern_Feb21.Rda")
# Model <- mdl.cn4
# 
# load("output/rda/ssd_eastern_Feb21.Rda")
# Model <- mdl.ssd1

# #pdf("figures/florum_sigmaSpSite.pdf", width = 12, height = 4)
# par(mfrow = c(1, 4))
# hist(postHtW$sigma_site, main = "Height", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 30, ylim = c(0, 4000))
# hist(postHtW$sigma_sp, add = T, col = col2)
# text(0.2,4000, label = "a)", cex = 2)
#      
# hist(postLMAW$sigma_site, main = "Leaf Mass Area", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 50, xlim = c(0,0.2), ylim = c(0, 6000))
# hist(postLMAW$sigma_sp, add = T, col = col2)
# text(0.002,6000, label = "b)", cex = 2)
# 
# hist(postDbhW$sigma_site, main = "Diameter at Breast Height", xlab = "Variation", ylab = "Frequency", col = col1, ylim = c(0, 6000))
# hist(postDbhW$sigma_sp, add = T, col = col2)
# text(0.2,6000, label = "c)", cex = 2)
# 
# hist(postCNW$sigma_site, main = "Carbon:Nitrogen", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 30, ylim = c(0, 4000), xlim = c(0,12))
# hist(postCNW$sigma_sp, add = T, col = col2)
# text(0.2,4000, label = "d)", cex = 2)
# 
# hist(postSsdW$sigma_site, main = "Stem Specific Density", xlab = "Variation", ylab = "Frequency", col = col1, breaks = 30, ylim = c(0, 4000), xlim = c(0,6))
# hist(postCNW$sigma_sp, add = T, col = col2, breaks = 50)
# text(0.1,4000, label = "e)", cex = 2)
# 
# legend("topright",legend = c(expression("Species variation"),
#                              expression("Study variation")),
#        col = c(col2,col1),
#        inset = 0.02, cex = 1, bty = "n", lwd = 3)
# 
# dev.off()