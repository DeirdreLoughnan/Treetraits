# aim of this code is to generate the values referenced in the sweave file:

library(rstan)
require(shinystan)
library(tidybayes)
library(reshape2)


########## General values ###############################
spInfo <- read.csv("..//input/species_ring.csv")

spInfoE <- subset(spInfo, transect != "west")
spInfoW <- subset(spInfo, transect != "east")

spInfoEShrub <- subset(spInfoE, type == "shrub")
spInfoWShrub <- subset(spInfoW, type == "shrub")

spInfoETree <- subset(spInfoE, type == "tree")
spInfoWTree <- subset(spInfoW, type == "tree")

nSpp <- length(unique(spInfo$species.name))
nSppE <- length(unique(spInfoE$species.name))
nSppW <- length(unique(spInfoW$species.name))
nSppES <- length(unique(spInfoEShrub$species.name))
nSppWS <- length(unique(spInfoWShrub$species.name))
nSppET <- length(unique(spInfoETree$species.name))
nSppWT <- length(unique(spInfoWTree$species.name))


pheno <- read.csv("..//input/phenoDataWChill.csv")
# trtPheno <- read.csv("input/trtData.csv")
trtPheno <- read.csv("..//input/trtPhenoDummy.csv")

height <- trtPheno[complete.cases(trtPheno$ht),]
nHt <- nrow(height)

leafMA <- trtPheno[complete.cases(trtPheno$lma),]
nLMA <- nrow(leafMA)

diam <- trtPheno[complete.cases(trtPheno$dbh),] 
nDBH <- nrow(diam)

stem <- trtPheno[complete.cases(trtPheno$ssd),]
nSSD <- nrow(stem)

carbNit <- trtPheno[complete.cases(trtPheno$C.N),]
nCNit <- nrow(carbNit)

nit <- trtPheno[complete.cases(trtPheno$per.N),]
nLNC <- nrow(nit)

temp <- trtPheno[,c("sample","species","site")]
temp$count <- 1
noObsTrait <- aggregate(temp["count"], trtPheno[,c("species","site")], FUN = sum)

noObsTrtMin <- min(noObsTrait$count)
noObsTrtMin <- min(noObsTrait$count)

eTrt <- subset(trtPheno, transect == "1")
wTrt <- subset(trtPheno, transect == "0")
noIndivE <- length(unique(eTrt$sample))
noIndivW<- length(unique(wTrt$sample))
noIndiv <- length(unique(trtPheno$sample))

pheno$count <- 1
noObsPheno <- aggregate(pheno["count"], pheno[,c("species","population")], FUN = sum)
noIndivPheno <- nrow(pheno)


# col4fig <- c("mean","sd","25%","50%","75%","Rhat")
# col4table <- c("mean","sd","2.5%","50%","97.5%","Rhat")
# 
# mu_params_4 <- c( 
#   #"a_z",
#   # "lam_interceptsa",
#   "mu_b_warm",
#   "mu_b_photo",
#   "mu_b_chill1",
#   "b_site2",
#   "b_site3",
#   "b_site4",
#   "mu_b_inter_wp",
#   "mu_b_inter_wc1",
#   "mu_b_inter_pc1",
#   "mu_b_inter_ws2",
#   "mu_b_inter_ps2",
#   "mu_b_inter_s2c1",
#   "mu_b_inter_ws3",
#   "mu_b_inter_ps3",
#   "mu_b_inter_s3c1",
#   "mu_b_inter_ws4",
#   "mu_b_inter_ps4",
#   "mu_b_inter_s4c1")
# 
# meanz4 <- sum[mu_params_4, col4table]
# 
# rownames(meanz4) = c( 
#   #"Root trait intercept", "Lambda",
#   "Forcing",
#   "Photoperiod",
#   "Chilling",
#   "Manning Park",
#   "Harvard Forest",
#   "St. Hippolyte",
#   "Forcing x photoperiod",
#   "Forcing x chilling",
#   "Photoperiod x chilling",
#   "Forcing x Manning Park",
#   "Photoperiod x Manning Park",
#   "Chilling x Manning Park",
#   "Forcing x Harvard Forest",
#   "Photoperiod x Harvard Forest",
#   "Chilling x Harvard Forest",
#   "Forcing x St. Hippolyte",
#   "Photoperiod x St. Hippolyte",
#   "Chilling x St. Hippolyte"            
# )
# 
# meanz4.table <- sum[mu_params_4, col4table]
# row.names(meanz4.table) <- row.names(meanz4)
# head(meanz4.table)
# #write.table(meanzew.table , "output/term.mdl.esti.dldf.csv", sep = ",", row.names = FALSE)
# 
# # call table values:
# 
# chillCue <- round(meanz4["Chilling",1],1)
# chillCueU <- round(quantile(fit$mu_b_chill1, c(0.95)),1)
# chillCueL <- round(quantile(fit$mu_b_chill1, c(0.05)),1)


load("..//output/htContLatHundoLat.Rdata")
htModelFit <- rstan::extract(mdlHt)

htTran <- as.numeric(round(mean(htModelFit$b_tranE),1))
lower_htTran <- format(round(quantile(htModelFit$b_tranE, prob = 0.05),1), nsmall =1)
upper_htTran <- round(quantile(htModelFit$b_tranE, prob = 0.95),1)

htLatTran <- as.numeric(round(mean(htModelFit$b_tranlat),1))
lower_htLatTran <- format(round(quantile(htModelFit$b_tranlat, prob = 0.05),1), nsmall =1)
upper_htLatTran <- round(quantile(htModelFit$b_tranlat, prob = 0.95),1)

htBFSpMean <- as.numeric(round(mean(htModelFit$betaTraitxForce),1))
lower_htBFSpMean <- format(round(quantile(htModelFit$betaTraitxForce, prob = 0.05),1), nsmall =1)
upper_htBFSpMean <- round(quantile(htModelFit$betaTraitxForce, prob = 0.95),1)
htBFSpMean; lower_htBFSpMean; upper_htBFSpMean

htBCSpMean <- as.numeric(round(mean(htModelFit$betaTraitxChill),1))
lower_htBCSpMean <- round(quantile(htModelFit$betaTraitxChill, prob = 0.05),1)
upper_htBCSpMean <- round(quantile(htModelFit$betaTraitxChill, prob = 0.95),1)
htBCSpMean; lower_htBCSpMean; upper_htBCSpMean

htBPSpMean <- as.numeric(round(mean(htModelFit$betaTraitxPhoto),1))
lower_htBPSpMean <- round(quantile(htModelFit$betaTraitxPhoto, prob = 0.05),1)
upper_htBPSpMean <- round(quantile(htModelFit$betaTraitxPhoto, prob = 0.95),1)
htBPSpMean; lower_htBPSpMean; upper_htBPSpMean

# muForceSp
htMuFSpMean <- as.numeric(round(mean(htModelFit$muForceSp),1))
lower_htMuFSpMean <- format(round(quantile(htModelFit$muForceSp, prob = 0.05),1), nsmall =1)
upper_htMuFSpMean <- round(quantile(htModelFit$muForceSp, prob = 0.95),1)
htMuFSpMean; lower_htMuFSpMean; upper_htMuFSpMean

#muChillSP
htMuCSpMean <- as.numeric(round(mean(htModelFit$muChillSp),1))
lower_htMuCSpMean <- format(round(quantile(htModelFit$muChillSp, prob = 0.05),1), nsmall =1)
upper_htMuCSpMean <- round(quantile(htModelFit$muChillSp, prob = 0.95),1)
htMuCSpMean; lower_htMuCSpMean; upper_htMuCSpMean

#muPhotoSp
htMuPSpMean <- as.numeric(round(mean(htModelFit$muPhotoSp),1))
lower_htMuPSpMean <- format(round(quantile(htModelFit$muPhotoSp, prob = 0.05),1), nsmall =1)
upper_htMuPSpMean <- round(quantile(htModelFit$muPhotoSp, prob = 0.95),1)
htMuPSpMean; lower_htMuPSpMean; upper_htMuPSpMean

htchill <- data.frame(htModelFit$alphaChillSp)
htChillSpMean <- colMeans(htchill)
htChillMax <- max(htChillSpMean)

htsigmaSp <- as.numeric(round(mean(htModelFit$sigma_sp),1))
lower_htsigmaSp <- format(round(quantile(htModelFit$sigma_sp, prob = 0.05),1), nsmall =1)
upper_htsigmaSp <- round(quantile(htModelFit$sigma_sp, prob = 0.95),1)
######## Leaf mass area###############################

load("..//output/lmaContLatHundoLat.Rdata")
lmaModelFit <- rstan::extract(mdlLMA)

lmaBFSpMean <- as.numeric(round(mean(lmaModelFit$betaTraitxForce),1))
lower_lmaBFSpMean <- format(round(quantile(lmaModelFit$betaTraitxForce, prob = 0.05),1), nsmall =1)
upper_lmaBFSpMean <- round(quantile(lmaModelFit$betaTraitxForce, prob = 0.95),1)
lmaBFSpMean; lower_lmaBFSpMean; upper_lmaBFSpMean

lmaBCSpMean <- as.numeric(round(mean(lmaModelFit$betaTraitxChill),1))
lower_lmaBCSpMean <- round(quantile(lmaModelFit$betaTraitxChill, prob = 0.05),1)
upper_lmaBCSpMean <- round(quantile(lmaModelFit$betaTraitxChill, prob = 0.95),1)
lmaBCSpMean; lower_lmaBCSpMean; upper_lmaBCSpMean

lmaBPSpMean <- as.numeric(round(mean(lmaModelFit$betaTraitxPhoto),1))
lower_lmaBPSpMean <- round(quantile(lmaModelFit$betaTraitxPhoto, prob = 0.05),1)
upper_lmaBPSpMean <- format(round(quantile(lmaModelFit$betaTraitxPhoto, prob = 0.95),1), nsmall =1)
lmaBPSpMean; lower_lmaBPSpMean; upper_lmaBPSpMean

lmachill <- data.frame(lmaModelFit$alphaChillSp)
lmaChillSpMean <- colMeans(lmachill)
lmaChillMax <- max(lmaChillSpMean)

# muForceSp
lmaMuFSpMean <- as.numeric(round(mean(lmaModelFit$muForceSp),1))
lower_lmaMuFSpMean <- format(round(quantile(lmaModelFit$muForceSp, prob = 0.05),1), nsmall =1)
upper_lmaMuFSpMean <- round(quantile(lmaModelFit$muForceSp, prob = 0.95),1)
lmaMuFSpMean; lower_lmaMuFSpMean; upper_lmaMuFSpMean

#muChillSP
lmaMuCSpMean <- as.numeric(round(mean(lmaModelFit$muChillSp),1))
lower_lmaMuCSpMean <- format(round(quantile(lmaModelFit$muChillSp, prob = 0.05),1), nsmall =1)
upper_lmaMuCSpMean <- round(quantile(lmaModelFit$muChillSp, prob = 0.95),1)
lmaMuCSpMean; lower_lmaMuCSpMean; upper_lmaMuCSpMean

#muPhotoSp
lmaMuPSpMean <- as.numeric(round(mean(lmaModelFit$muPhotoSp),1))
lower_lmaMuPSpMean <- format(round(quantile(lmaModelFit$muPhotoSp, prob = 0.05),1), nsmall =1)
upper_lmaMuPSpMean <- round(quantile(lmaModelFit$muPhotoSp, prob = 0.95),1)
lmaMuPSpMean; lower_lmaMuPSpMean; upper_lmaMuPSpMean

lmachill <- data.frame(lmaModelFit$alphaChillSp)
lmaChillSpMean <- colMeans(lmachill)
lmaChillMax <- max(lmaChillSpMean)

lmaLatTran <- as.numeric(round(mean(lmaModelFit$b_tranlat),1))
lower_lmaLatTran <- format(round(quantile(lmaModelFit$b_tranlat, prob = 0.05),1), nsmall =1)
upper_lmaLatTran <- round(quantile(lmaModelFit$b_tranlat, prob = 0.95),1)

lmasigmaSp <- as.numeric(round(mean(lmaModelFit$sigma_sp),1))
lower_lmasigmaSp <- format(round(quantile(lmaModelFit$sigma_sp, prob = 0.05),1), nsmall =1)
upper_lmasigmaSp <- round(quantile(lmaModelFit$sigma_sp, prob = 0.95),1)
##### C:N ###############################
load("..//output/lncContLatHundoLat.Rdata")
lncModelFit <- rstan::extract(mdlPerN)

lncBFSpMean <- as.numeric(round(mean(lncModelFit$betaTraitxForce),1))
lower_lncBFSpMean <- format(round(quantile(lncModelFit$betaTraitxForce, prob = 0.05),1), nsmall =1)
upper_lncBFSpMean <- round(quantile(lncModelFit$betaTraitxForce, prob = 0.95),1)
lncBFSpMean; lower_lncBFSpMean; upper_lncBFSpMean

lncBCSpMean <- as.numeric(round(mean(lncModelFit$betaTraitxChill),1))
lower_lncBCSpMean <- round(quantile(lncModelFit$betaTraitxChill, prob = 0.05),1)
upper_lncBCSpMean <- round(quantile(lncModelFit$betaTraitxChill, prob = 0.95),1)
lncBCSpMean; lower_lncBCSpMean; upper_lncBCSpMean

lncBPSpMean <- as.numeric(round(mean(lncModelFit$betaTraitxPhoto),1))
lower_lncBPSpMean <- round(quantile(lncModelFit$betaTraitxPhoto, prob = 0.05),1)
upper_lncBPSpMean <- format(round(quantile(lncModelFit$betaTraitxPhoto, prob = 0.95),1), nsmall =1)
lncBPSpMean; lower_lncBPSpMean; upper_lncBPSpMean

lncchill <- data.frame(lncModelFit$alphaChillSp)
lncChillSpMean <- colMeans(lncchill)
lncChillMax <- max(lncChillSpMean)

# muForceSp
lncMuFSpMean <- as.numeric(round(mean(lncModelFit$muForceSp),1))
lower_lncMuFSpMean <- format(round(quantile(lncModelFit$muForceSp, prob = 0.05),1), nsmall =1)
upper_lncMuFSpMean <- round(quantile(lncModelFit$muForceSp, prob = 0.95),1)
lncMuFSpMean; lower_lncMuFSpMean; upper_lncMuFSpMean

#muChillSP
lncMuCSpMean <- as.numeric(round(mean(lncModelFit$muChillSp),1))
lower_lncMuCSpMean <- format(round(quantile(lncModelFit$muChillSp, prob = 0.05),1), nsmall =1)
upper_lncMuCSpMean <- round(quantile(lncModelFit$muChillSp, prob = 0.95),1)
lncMuCSpMean; lower_lncMuCSpMean; upper_lncMuCSpMean

#muPhotoSp
lncMuPSpMean <- as.numeric(round(mean(lncModelFit$muPhotoSp),1))
lower_lncMuPSpMean <- format(round(quantile(lncModelFit$muPhotoSp, prob = 0.05),1), nsmall =1)
upper_lncMuPSpMean <- round(quantile(lncModelFit$muPhotoSp, prob = 0.95),1)
lncMuPSpMean; lower_lncMuPSpMean; upper_lncMuPSpMean

lncchill <- data.frame(lncModelFit$alphaChillSp)
lncChillSpMean <- colMeans(lncchill)
lncChillMax <- max(lncChillSpMean)

lncLatTran <- as.numeric(round(mean(lncModelFit$b_tranlat),1))
lower_lncLatTran <- format(round(quantile(lncModelFit$b_tranlat, prob = 0.05),1), nsmall =1)
upper_lncLatTran <- round(quantile(lncModelFit$b_tranlat, prob = 0.95),1)
#lncmuSp <- apply(posterior_ssd$muSp, MARGIN = 2, FUN = mean)

lncTran <- as.numeric(round(mean(lncModelFit$b_tranE),1))
lower_lncTran <- format(round(quantile(lncModelFit$b_tranE, prob = 0.05),1), nsmall =1)
upper_lncTran <- round(quantile(lncModelFit$b_tranE, prob = 0.95),1)
#lncmuSp <- apply(posterior_ssd$muSp, MARGIN = 2, FUN = mean)

lncsigmaSp <- as.numeric(round(mean(lncModelFit$sigma_sp),1))
lower_lncsigmaSp <- format(round(quantile(lncModelFit$sigma_sp, prob = 0.05),1), nsmall =1)
upper_lncsigmaSp <- round(quantile(lncModelFit$sigma_sp, prob = 0.95),1)
##### SSD  ###############################
load("..//output/ssdContLatHundoLat.Rdata")
ssdModelFit <- rstan::extract(mdlSSD)

ssdBFSpMean <- as.numeric(round(mean(ssdModelFit$betaTraitxForce),1))
lower_ssdBFSpMean <- format(round(quantile(ssdModelFit$betaTraitxForce, prob = 0.05),1), nsmall =1)
upper_ssdBFSpMean <- round(quantile(ssdModelFit$betaTraitxForce, prob = 0.95),1)
ssdBFSpMean; lower_ssdBFSpMean; upper_ssdBFSpMean

ssdBCSpMean <- as.numeric(round(mean(ssdModelFit$betaTraitxChill),1))
lower_ssdBCSpMean <- round(quantile(ssdModelFit$betaTraitxChill, prob = 0.05),1)
upper_ssdBCSpMean <- round(quantile(ssdModelFit$betaTraitxChill, prob = 0.95),1)
ssdBCSpMean; lower_ssdBCSpMean; upper_ssdBCSpMean

ssdBPSpMean <- as.numeric(round(mean(ssdModelFit$betaTraitxPhoto),1))
lower_ssdBPSpMean <- round(quantile(ssdModelFit$betaTraitxPhoto, prob = 0.05),1)
upper_ssdBPSpMean <- format(round(quantile(ssdModelFit$betaTraitxPhoto, prob = 0.95),1), nsmall =1)
ssdBPSpMean; lower_ssdBPSpMean; upper_ssdBPSpMean

ssdchill <- data.frame(ssdModelFit$alphaChillSp)
ssdChillSpMean <- colMeans(ssdchill)
ssdChillMax <- max(ssdChillSpMean)

# muForceSp
ssdMuFSpMean <- as.numeric(round(mean(ssdModelFit$muForceSp),1))
lower_ssdMuFSpMean <- format(round(quantile(ssdModelFit$muForceSp, prob = 0.05),1), nsmall =1)
upper_ssdMuFSpMean <- round(quantile(ssdModelFit$muForceSp, prob = 0.95),1)
ssdMuFSpMean; lower_ssdMuFSpMean; upper_ssdMuFSpMean

#muChillSP
ssdMuCSpMean <- as.numeric(round(mean(ssdModelFit$muChillSp),1))
lower_ssdMuCSpMean <- format(round(quantile(ssdModelFit$muChillSp, prob = 0.05),1), nsmall =1)
upper_ssdMuCSpMean <- round(quantile(ssdModelFit$muChillSp, prob = 0.95),1)
ssdMuCSpMean; lower_ssdMuCSpMean; upper_ssdMuCSpMean

#muPhotoSp
ssdMuPSpMean <- as.numeric(round(mean(ssdModelFit$muPhotoSp),1))
lower_ssdMuPSpMean <- format(round(quantile(ssdModelFit$muPhotoSp, prob = 0.05),1), nsmall =1)
upper_ssdMuPSpMean <- round(quantile(ssdModelFit$muPhotoSp, prob = 0.95),1)
ssdMuPSpMean; lower_ssdMuPSpMean; upper_ssdMuPSpMean

ssdchill <- data.frame(ssdModelFit$alphaChillSp)
ssdChillSpMean <- colMeans(ssdchill)
ssdChillMax <- max(ssdChillSpMean)

ssdLatTran <- as.numeric(round(mean(ssdModelFit$b_tranlat),1))
lower_ssdLatTran <- format(round(quantile(ssdModelFit$b_tranlat, prob = 0.05),1), nsmall =1)
upper_ssdLatTran <- round(quantile(ssdModelFit$b_tranlat, prob = 0.95),1)

ssdsigmaSp <- as.numeric(round(mean(ssdModelFit$sigma_sp),1))
lower_ssdsigmaSp <- format(round(quantile(ssdModelFit$sigma_sp, prob = 0.05),1), nsmall =1)
upper_ssdsigmaSp <- round(quantile(ssdModelFit$sigma_sp, prob = 0.95),1)
########## DBH ###############################
load("..//output/dbhContLatHundoLat.Rdata")
dbhModelFit <- rstan::extract(mdlDBH)

dbhBFSpMean <- as.numeric(round(mean(dbhModelFit$betaTraitxForce),1))
lower_dbhBFSpMean <- format(round(quantile(dbhModelFit$betaTraitxForce, prob = 0.05),1), nsmall =1)
upper_dbhBFSpMean <- round(quantile(dbhModelFit$betaTraitxForce, prob = 0.95),1)
dbhBFSpMean; lower_dbhBFSpMean; upper_dbhBFSpMean

dbhBCSpMean <- as.numeric(round(mean(dbhModelFit$betaTraitxChill),1))
lower_dbhBCSpMean <- round(quantile(dbhModelFit$betaTraitxChill, prob = 0.05),1)
upper_dbhBCSpMean <- round(quantile(dbhModelFit$betaTraitxChill, prob = 0.95),1)
dbhBCSpMean; lower_dbhBCSpMean; upper_dbhBCSpMean

dbhBPSpMean <- as.numeric(round(mean(dbhModelFit$betaTraitxPhoto),1))
lower_dbhBPSpMean <- round(quantile(dbhModelFit$betaTraitxPhoto, prob = 0.05),1)
upper_dbhBPSpMean <- format(round(quantile(dbhModelFit$betaTraitxPhoto, prob = 0.95),1), nsmall =1)
dbhBPSpMean; lower_dbhBPSpMean; upper_dbhBPSpMean

dbhchill <- data.frame(dbhModelFit$alphaChillSp)
dbhChillSpMean <- colMeans(dbhchill)
dbhChillMax <- max(dbhChillSpMean)

# muForceSp
dbhMuFSpMean <- as.numeric(round(mean(dbhModelFit$muForceSp),1))
lower_dbhMuFSpMean <- format(round(quantile(dbhModelFit$muForceSp, prob = 0.05),1), nsmall =1)
upper_dbhMuFSpMean <- round(quantile(dbhModelFit$muForceSp, prob = 0.95),1)
dbhMuFSpMean; lower_dbhMuFSpMean; upper_dbhMuFSpMean

#muChillSP
dbhMuCSpMean <- as.numeric(round(mean(dbhModelFit$muChillSp),1))
lower_dbhMuCSpMean <- format(round(quantile(dbhModelFit$muChillSp, prob = 0.05),1), nsmall =1)
upper_dbhMuCSpMean <- round(quantile(dbhModelFit$muChillSp, prob = 0.95),1)
dbhMuCSpMean; lower_dbhMuCSpMean; upper_dbhMuCSpMean

#muPhotoSp
dbhMuPSpMean <- as.numeric(round(mean(dbhModelFit$muPhotoSp),1))
lower_dbhMuPSpMean <- format(round(quantile(dbhModelFit$muPhotoSp, prob = 0.05),1), nsmall =1)
upper_dbhMuPSpMean <- round(quantile(dbhModelFit$muPhotoSp, prob = 0.95),1)
dbhMuPSpMean; lower_dbhMuPSpMean; upper_dbhMuPSpMean

dbhchill <- data.frame(dbhModelFit$alphaChillSp)
dbhChillSpMean <- colMeans(dbhchill)
dbhChillMax <- max(dbhChillSpMean)

dbhLatTran <- as.numeric(round(mean(dbhModelFit$b_tranlat),1))
lower_dbhLatTran <- format(round(quantile(dbhModelFit$b_tranlat, prob = 0.05),1), nsmall =1)
upper_dbhLatTran <- round(quantile(dbhModelFit$b_tranlat, prob = 0.95),1)

dbhsigmaSp <- as.numeric(round(mean(dbhModelFit$sigma_sp),1))
lower_dbhsigmaSp <- format(round(quantile(dbhModelFit$sigma_sp, prob = 0.05),1), nsmall =1)
upper_dbhsigmaSp <- round(quantile(dbhModelFit$sigma_sp, prob = 0.95),1)

tot <- read.csv("..//..//pheno_bc/input/phenoMini.csv")

#Only want the total number that reached 7:
below <- subset(tot, bbch.t < 8)
totalObs <- nrow(below)


totalDays <- max(tot$day)