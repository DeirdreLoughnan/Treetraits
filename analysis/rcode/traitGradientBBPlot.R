# Aim of this code is to see whether we still get the gradients in bb in the joint model

rm(list=ls())
options(stringsAsFactors = FALSE)

# the set cues will be: 12h photoperiod, 20C, high chill - 75/10
require(rstan)
library(forcats)
library(ggdist)
library(reshape2)
require(cowplot)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/Treetraits") 
}  else{
  setwd("~/Documents/github/Treetraits") # for midge
}
spInfo <- read.csv("analysis/input/species_ring.csv")
colnames(spInfo)[colnames(spInfo) == "X"] <- "ringType"
spInfo <- spInfo[,1:5]

trtPheno <- read.csv("analysis/input/trtPhenoDummy.csv")
trtMeans <- aggregate(trtPheno[c("ssd","ht","lma","dbh","C.N","per.N")], trtPheno[c("species")], FUN = mean, na.rm = T)

spInfo <- merge(spInfo, trtMeans, by = "species")

# load("output/heightDummyIntGrandZ25.Rdata")
# sumHt <- summary(mdlHt)$summary
# postHt <- rstan::extract(mdlHt)

load("analysis/output/htContLatHundoLatFinal.Rdata")
sumHt <- summary(mdlHt)$summary
postHt <- rstan::extract(mdlHt)

# load("output/lmaContLat.Rdata")
# postLMA <- rstan::extract(mdlLMA)
# sumLMA - summary(mdlLMA)$summary
# 
# load("output/dbhContLat.Rdata")
# postDBH <- rstan::extract(mdlDBH)
# sumDBH<- summary(mdlDBH)$summary
# 
# load("output/ssdContLat.Rdata")
# postSSD <- rstan::extract(mdlSSD)
# sumSSD<- summary(mdlSSD)$summary

load("analysis/output/lncContLatHundoLatFinal.Rdata")
postCN <- rstan::extract(mdlPerN)
sumCN<- summary(mdlPerN)$summary


a_sp = (sumHt[grep("mu_grand_sp", rownames(sumHt)), 1])
b_photo = sumHt[grep("betaPhotoSp\\[", rownames(sumHt)), 1]
b_chill = sumHt[grep("betaChillSp\\[", rownames(sumHt)), 1]
b_force = sumHt[grep("betaForceSp\\[", rownames(sumHt)), 1]

a_sp5 <- vector()
for(i in 1:ncol(postHt$mu_grand_sp)){
  quantU <- round(quantile(postHt$mu_grand_sp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  a_sp5 <- rbind(a_sp5, quantU)
}
colnames(a_sp5) <- c("Int5","Int95","Int25","Int75")
#a_sp5 <- a_sp5/100

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

# #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#

#If we are using the old model, we will use the z-scored values for the parameters
photo <- -0.5033863 #8 h photo
force <- -0.3568628 #5/15 C trt
chill <- -0.3546922 # low chill

m <- matrix(nrow = 1000, ncol = 47)

for(sp in 1:47){
  for (it in 1:nrow(m)){
    m[it,sp] <- postHt$mu_grand_sp[it,sp]+  
      postHt$betaForceSp[it,sp] * force + 
      postHt$betaPhotoSp[it, sp] * photo + 
      postHt$betaChillSp[it,sp] * chill 
  }
}

photoHigh <- 0.4965051 #8 h photo
siteSM <- 0
forceHigh <- 0.5877121 
chillHigh <- 0.3660412 # high chill for Smithers

mHigh <- matrix(nrow = 1000, ncol = 47)

for(sp in 1:47){
  for (it in 1:nrow(mHigh)){
    mHigh[it,sp] <- postHt$mu_grand_sp[it,sp]+  
      postHt$betaForceSp[it,sp] * forceHigh + 
      postHt$betaPhotoSp[it, sp] * photoHigh + 
      postHt$betaChillSp[it,sp] * chillHigh 
  }
}

spInfo <- spInfo[order(spInfo$species),]
head(spInfo)
spInfo$meanBB <- colMeans(m)
colnames(m) <- spInfo$species.name

spInfo$meanBBHigh <- colMeans(mHigh)
colnames(mHigh) <- spInfo$species.name

spInfo$Int <- a_sp
spInfo <- cbind(spInfo, a_sp5,b_force5, b_chill5,b_photo5)

spInfo$force <- b_force
spInfo$chill <- b_chill
spInfo$photo <- b_photo

quantile595 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.05, 0.95, 0.25,0.75))
  return(returnQuanilte)
}

bb_quan <- apply(m, 2, quantile595)
bb_t <- t(bb_quan)
bb_df <- data.frame(bb_t)
colnames(bb_df)[colnames(bb_df) == "X5."] <- "bb5"
colnames(bb_df)[colnames(bb_df) == "X95."] <- "bb95"
colnames(bb_df)[colnames(bb_df) == "X25."] <- "bb25"
colnames(bb_df)[colnames(bb_df) == "X75."] <- "bb75"


bb_quanHigh <- apply(mHigh, 2, quantile595)
bb_tHigh <- t(bb_quanHigh)
bb_dfHigh <- data.frame(bb_tHigh)
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X5."] <- "bb5High"
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X95."] <- "bb95High"
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X25."] <- "bb25High"
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X75."] <- "bb75High"



spInfo <- cbind(spInfo, bb_df)
spInfo$value <- spInfo$meanBB

spInfo <- cbind(spInfo, bb_dfHigh)
spInfo$valueHigh <- spInfo$meanBBHigh


m <- data.frame(m)

long <- reshape2::melt(m)
names(long) <- c("species.name", "valueLow")

mHigh <- data.frame(mHigh)

longHigh <- reshape2::melt(mHigh)
names(longHigh) <- c("species.name", "valueHigh")

long <- cbind(long, longHigh[,2])

long <- merge(long,spInfo, by = "species.name")

spOrderData <- spInfo[order(spInfo$meanBB),]
spOrder <- as.factor(spOrderData$species.name)

spOrderHtData <- spInfo[order(spInfo$ht),]
spOrderHt <- as.factor(spOrderData$species.name)

long <- long[order(long$species),]

# longPhotoInfo$mean <- rowMeans(longPhotoInfo[,c("Site1","Site2","Site3","Site4")], na.rm=TRUE)

bChill <- data.frame(postHt$betaChillSp[1:1000,])
colnames(bChill) <- (spInfo$species.name)
longChill <- reshape2::melt(bChill)
names(longChill) <- c("species.name", "chill")

long <- cbind(long, longChill$chill)

# Add forcing
bForce <- data.frame(postHt$betaForceSp[1:1000,])
colnames(bForce) <- (spInfo$species.name)
longForce <- reshape2::melt(bForce)
names(longForce) <- c("species.name", "force")

long <- cbind(long, longForce$force)

# photoperiod
bPhoto <- data.frame(postHt$betaPhotoSp[1:1000,])
colnames(bPhoto) <- (spInfo$species.name)
longPhoto <- reshape2::melt(bPhoto)
names(longPhoto) <- c("species.name", "photo")

long <- cbind(long, longPhoto$photo)

# intercept
aSp <- data.frame(postHt$mu_grand_sp[1:1000,])
colnames(aSp) <- (spInfo$species.name)
longInt <- reshape2::melt(aSp)
names(longInt) <- c("species.name", "int")

long <- cbind(long, longInt$int)

data <- long[order(long$meanBB),]
# 
# data$species.name <- factor(data$species.name, levels=unique(data$species.name) )
#data <- transform(data, variable=reorder(species.name, -meanBB) ) 

names(data) <- c("species.name","valueLow", "valueHigh","species","type","transect","ringType","ssd","ht","lma","dbh","C.N","per.N","meanBB","meanBBHigh", "Int","Int5","Int95","Int25","Int75","force5","force95","force25","force75","chill5", "chill95", "chill25", "chill75","photo5", "photo95", "photo25", "photo75","spMeanForce", "spMeanChill", "spMeanPhoto","bb5","bb95","bb25","bb75", "spacing","bb5High","bb95High","bb75High","bb25High", "valueHigh","chill","force","photo","intercept")

east <- subset(spInfo, transect != "west")
eastSp <- unique(east$species.name)

dataEast <- data[data$species.name %in% eastSp,]

meanPtE <- aggregate(dataEast[c("meanBB", "meanBBHigh","Int","ht","per.N")], dataEast[c("species.name","type","transect")], FUN = mean)
names(meanPtE) <- c("species.name","type","transect","Budburst", "BudburstHigh","Intercept","ht","lnc")


 htE <- ggplot(meanPtE) +
  geom_point(aes(y= Budburst, x = Budburst, shape = "Budburst", col=type ), size = 5) +
  geom_point(aes(y= BudburstHigh, x = Budburst, shape = "BudburstHigh", col=type), size = 5) +
  geom_point(aes(y= Intercept, x = Budburst, shape = "Intercept", col=type), size = 5) +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Budburst), data = meanPtE, col = "black") +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = BudburstHigh), data = meanPtE, col = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=20)) +
   ylim (-20,35) +
  scale_x_continuous( breaks = east$meanBB, labels = east$species,limits = c(7,28)) +
  labs( x = "", y = "Day of budburst", main = NA) +
  theme(legend.title = element_blank()) + annotate("text", x = 15, y = 35, label = "a) Eastern transect - height", cex =8) +
  scale_color_manual(values = c("maroon","cyan4")) +
  scale_shape_discrete( labels = c("low cue",
                                   "high cue",
                                   "intercept" ),
                        breaks = c("Budburst","BudburstHigh", "Intercept"))
htE
# htOrderE <- ggplot(meanPtE) +
#    geom_point(aes(y= Budburst, x = ht, shape = "Budburst", col=type ), size = 5) +
#    geom_point(aes(y= BudburstHigh, x = ht, shape = "BudburstHigh", col=type), size = 5) +
#    geom_point(aes(y= Intercept, x = ht, shape = "Intercept", col=type), size = 5) +
#    geom_segment(aes(x = ht, y = Intercept, xend = ht, yend = Budburst), data = meanPtE, col = "black") +
#    geom_segment(aes(x = ht, y = Intercept, xend = ht, yend = BudburstHigh), data = meanPtE, col = "black") +
#    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#      panel.background = element_blank(), axis.line = element_line(colour = "black"),
#      legend.key=element_rect(fill="white")) +
#    theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
#      axis.text.y=element_text(size = 15),
#      axis.title=element_text(size=20)) +
#    scale_x_continuous( breaks = east$ht, labels = east$species,limits = c(0.5,18)) +
#    labs( x = "", y = "Day of budburst (days/standardized units)", main = NA) +
#    theme(legend.title = element_blank()) +  annotate("text", x = 4, y = 20, label = "a) Eastern transect - height", cex =8) +
#    scale_color_manual(values = c("maroon","cyan4")) +
#    scale_shape_discrete( labels = c("low cue",
#      "high cue",
#      "intercept" ),
#      breaks = c("Budburst","BudburstHigh", "Intercept"))
 
# ggplot(meanPtE) +
#   geom_point(aes(y= Budburst, x = ht, shape = "Budburst", col=type ), size = 5) +
#   geom_point(aes(y= BudburstHigh, x = Budburst, shape = "BudburstHigh", col=type), size = 5) +
#   geom_point(aes(y= Intercept, x = Budburst, shape = "Intercept", col=type), size = 5) +
#   geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Budburst), data = meanPtE, col = "black") +
#   geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = BudburstHigh), data = meanPtE, col = "black") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),
#     legend.key=element_rect(fill="white")) +
#   theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
#     axis.text.y=element_text(size = 15),
#     axis.title=element_text(size=20)) +
#   scale_x_continuous( breaks = east$meanBB, labels = east$species,limits = c(2,22)) +
#   labs( x = "", y = "Day of budburst (days/standardized units)", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 12, y = 20, label = "a) Eastern transect - height", cex =8) +
#   scale_color_manual(values = c("maroon","cyan4")) +
#   scale_shape_discrete( labels = c("low cue",
#     "high cue",
#     "intercept" ),
#     breaks = c("Budburst","BudburstHigh", "Intercept"))
# 

###################
west <- subset(spInfo, transect != "east")
westSp <- unique(west$species.name)

datawest <- data[data$species.name %in% westSp,]

meanPtW <- aggregate(datawest[c("meanBB", "meanBBHigh","Int","ht","per.N")], datawest[c("species.name","type","transect")], FUN = mean)
names(meanPtW) <- c("species.name","type","transect","Budburst", "BudburstHigh","Intercept","ht","lnc")

htW <- ggplot(meanPtW) +
  geom_point(aes(y= Budburst, x = Budburst, shape = "Budburst", col=type ), size = 5) +
  geom_point(aes(y= BudburstHigh, x = Budburst, shape = "BudburstHigh", col=type), size = 5) +
  geom_point(aes(y= Intercept, x = Budburst, shape = "Intercept", col=type), size = 5) +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Budburst), data = meanPtW, col = "black") +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = BudburstHigh), data = meanPtW, col = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=20)) +
  ylim (-20,35) +
  scale_x_continuous( breaks = west$meanBB, labels = west$species,limits = c(7,30)) +
  labs( x = "", y = "Day of budburst", main = NA) +
  theme(legend.title = element_blank()) + annotate("text", x = 15, y = 35, label = "b) Western transect - height", cex =8) +
  scale_color_manual(values = c("maroon","cyan4")) +
  scale_shape_discrete( labels = c("low cue",
                                   "high cue",
                                   "intercept" ),
                        breaks = c("Budburst","BudburstHigh", "Intercept"))
htW

# htOrderW <- ggplot(meanPtW) +
#   geom_point(aes(y= Budburst, x = ht, shape = "Budburst", col=type ), size = 5) +
#   geom_point(aes(y= BudburstHigh, x = ht, shape = "BudburstHigh", col=type), size = 5) +
#   geom_point(aes(y= Intercept, x = ht, shape = "Intercept", col=type), size = 5) +
#   geom_segment(aes(x = ht, y = Intercept, xend = ht, yend = Budburst), data = meanPtW, col = "black") +
#   geom_segment(aes(x = ht, y = Intercept, xend = ht, yend = BudburstHigh), data = meanPtW, col = "black") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),
#     legend.key=element_rect(fill="white")) +
#   theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
#     axis.text.y=element_text(size = 15),
#     axis.title=element_text(size=20)) +
#   scale_x_continuous( breaks = west$ht, labels = west$species,limits = c(0.5,18)) +
#   labs( x = "", y = "Day of budburst (days/standardized units)", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 4, y = 20, label = "c) Western transect - height", cex =8) +
#   scale_color_manual(values = c("maroon","cyan4")) +
#   scale_shape_discrete( labels = c("low cue",
#     "high cue",
#     "intercept" ),
#     breaks = c("Budburst","BudburstHigh", "Intercept"))


rank <- spInfo[,c("species.name","species","type","transect","meanBB","meanBBHigh","Int")]

rankE <- subset(rank, transect == "east")

rankE <- rankE[order(rankE$Int),]
rankE$rankInt <- seq(1:nrow(rankE))

rankE <- rankE[order(rankE$meanBB),]
rankE$rankLowC <- seq(1:nrow(rankE))

rankE <- rankE[order(rankE$meanBBHigh),]
rankE$rankHighC <- seq(1:nrow(rankE))

rankW <- subset(rank, transect == "west")

rankW <- rankW[order(rankW$Int),]
rankW$rankInt <- seq(1:nrow(rankW))

rankW <- rankW[order(rankW$meanBB),]
rankW$rankLowC <- seq(1:nrow(rankW))

rankW <- rankW[order(rankW$meanBBHigh),]
rankW$rankHighC <- seq(1:nrow(rankW))

########### CN ##############################

a_sp = (sumCN[grep("mu_grand_sp", rownames(sumCN)), 1])
b_photo = sumCN[grep("betaPhotoSp\\[", rownames(sumCN)), 1]
b_chill = sumCN[grep("betaChillSp\\[", rownames(sumCN)), 1]
b_force = sumCN[grep("betaForceSp\\[", rownames(sumCN)), 1]

a_sp5 <- vector()
for(i in 1:ncol(postCN$mu_grand_sp)){
  quantU <- round(quantile(postCN$a_sp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  a_sp5 <- rbind(a_sp5, quantU)
}
colnames(a_sp5) <- c("Int5","Int95","Int25","Int75")

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

# #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#


mCN <- matrix(nrow = 1000, ncol = 47)

for(sp in 1:47){
  for (it in 1:nrow(mCN)){
    mCN[it,sp] <- postCN$mu_grand_sp[it,sp]+  
      postCN$betaForceSp[it,sp] * force + 
      postCN$betaPhotoSp[it, sp] * photo + 
      postCN$betaChillSp[it,sp] * chill 
  }
}

photoHigh <- 0.4965051 #8 h photo
siteSM <- 0
forceHigh <- 0.5877121 
chillHigh <- 0.3660412 # high chill for Smithers

mCNHigh <- matrix(nrow = 1000, ncol = 47)

for(sp in 1:47){
  for (it in 1:nrow(mCNHigh)){
    mCNHigh[it,sp] <- postCN$mu_grand_sp[it,sp]+  
      postCN$betaForceSp[it,sp] * forceHigh + 
      postCN$betaPhotoSp[it, sp] * photoHigh + 
      postCN$betaChillSp[it,sp] * chillHigh 
  }
}

spInfo <- spInfo[order(spInfo$species),]
head(spInfo)
spInfo$meanBB <- colMeans(mCN)
colnames(mCN) <- spInfo$species.name

spInfo$meanBBHigh <- colMeans(mCNHigh)
colnames(mCNHigh) <- spInfo$species.name

spInfo$Int <- a_sp
spInfo <- cbind(spInfo, a_sp5,b_force5, b_chill5,b_photo5)

spInfo$force <- b_force
spInfo$chill <- b_chill
spInfo$photo <- b_photo




quantile595 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.05, 0.95))
  return(returnQuanilte)
}

quantile75.25 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.75, 0.25))
  return(returnQuanilte)
}

bb_quan <- apply(mCN, 2, quantile595)
bb_t <- t(bb_quan)
bb_df <- data.frame(bb_t)
colnames(bb_df)[colnames(bb_df) == "X5."] <- "bb5"
colnames(bb_df)[colnames(bb_df) == "X95."] <- "bb95"

bb_quan75.25 <- apply(mCN, 2, quantile75.25)
bb_t75.25 <- t(bb_quan75.25)
bb_df75.25 <- data.frame(bb_t75.25)
colnames(bb_df75.25)[colnames(bb_df75.25) == "X75."] <- "bb75"
colnames(bb_df75.25)[colnames(bb_df75.25) == "X25."] <- "bb25"
colnames(bb_df)[colnames(bb_df) == "X25."] <- "bb25"
colnames(bb_df)[colnames(bb_df) == "X75."] <- "bb75"

bb_quanHigh <- apply(mCNHigh, 2, quantile595)
bb_tHigh <- t(bb_quanHigh)
bb_dfHigh <- data.frame(bb_tHigh)
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X5."] <- "bb5High"
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X95."] <- "bb95High"

bb_quan75.25High <- apply(mCNHigh, 2, quantile75.25)
bb_t75.25High <- t(bb_quan75.25High)
bb_df75.25High <- data.frame(bb_t75.25High)
colnames(bb_df75.25High)[colnames(bb_df75.25High) == "X75."] <- "bb75High"
colnames(bb_df75.25High)[colnames(bb_df75.25High) == "X25."] <- "bb25High"
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X25."] <- "bb25High"
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X75."] <- "bb75High"


spInfo <- cbind(spInfo, bb_df)
spInfo <- cbind(spInfo, bb_df75.25)
spInfo$value <- spInfo$meanBB

spInfo <- cbind(spInfo, bb_dfHigh)
spInfo <- cbind(spInfo, bb_df75.25High)
spInfo$valueHigh <- spInfo$meanBBHigh


mCN2 <- data.frame(mCN)

long <- reshape2::melt(mCN2)
names(long) <- c("species.name", "valueLow")

mHigh <- data.frame(mCNHigh)

longHigh <- melt(mHigh)
names(longHigh) <- c("species.name", "valueHigh")

long <- cbind(long, longHigh[,2])

long <- merge(long,spInfo, by = "species.name")

spOrderData <- spInfo[order(spInfo$meanBB),]
spOrder <- as.factor(spOrderData$species.name)

long <- long[order(long$species),]

# longPhotoInfo$mean <- rowMeans(longPhotoInfo[,c("Site1","Site2","Site3","Site4")], na.rm=TRUE)

bChill <- data.frame(postCN$betaChillSp[1:1000,])
colnames(bChill) <- (spInfo$species.name)
longChill <- melt(bChill)
names(longChill) <- c("species.name", "chill")

long <- cbind(long, longChill$chill)

# Add forcing
bForce <- data.frame(postCN$betaForceSp[1:1000,])
colnames(bForce) <- (spInfo$species.name)
longForce <- melt(bForce)
names(longForce) <- c("species.name", "force")

long <- cbind(long, longForce$force)

# photoperiod
bPhoto <- data.frame(postCN$betaPhotoSp[1:1000,])
colnames(bPhoto) <- (spInfo$species.name)
longPhoto <- melt(bPhoto)
names(longPhoto) <- c("species.name", "photo")

long <- cbind(long, longPhoto$photo)

# intercept
aSp <- data.frame(postCN$mu_grand_sp[1:1000,])
colnames(aSp) <- (spInfo$species.name)
longInt <- melt(aSp)
names(longInt) <- c("species.name", "int")

long <- cbind(long, longInt$int)

data <- long[order(long$meanBB),]
# 
# data$species.name <- factor(data$species.name, levels=unique(data$species.name) )
#data <- transform(data, variable=reorder(species.name, -meanBB) ) 

names(data) <- c("species.name","valueLow", "valueHigh","species","type","transect","ringType","ssd","ht","lma", "dbh","C.N","per.N","meanBB","meanBBHigh", "Int","Int5","Int95","Int25","Int75","force5","force95","force25","force75","chill5", "chill95", "chill25", "chill75","photo5", "photo95", "photo25", "photo75","spMeanForce", "spMeanChill", "spMeanPhoto","bb5","bb95","bb25","bb75", "spacing","bb5High","bb95High","bb75High","bb25High", "valueHigh","chill","force","photo","intercept")

east <- subset(spInfo, transect != "west")
eastSp <- unique(east$species.name)

dataEast <- data[data$species.name %in% eastSp,]

meanPtE <- aggregate(dataEast[c("meanBB", "meanBBHigh","Int","ht","per.N")], dataEast[c("species.name","type","transect")], FUN = mean)
names(meanPtE) <- c("species.name","type","transect","Budburst", "BudburstHigh","Intercept","ht","lnc")

###################

meanPtW <- aggregate(datawest[c("meanBB", "meanBBHigh","Int","ht","per.N")], datawest[c("species.name","type","transect")], FUN = mean)
names(meanPtW) <- c("species.name","type","transect","Budburst", "BudburstHigh","Intercept","ht","lnc")


CNE <- ggplot(meanPtE) +
  geom_point(aes(y= Budburst, x = Budburst, shape = "Budburst", col=type ), size = 5) +
  geom_point(aes(y= BudburstHigh, x = Budburst, shape = "BudburstHigh", col=type), size = 5) +
  geom_point(aes(y= Intercept, x = Budburst, shape = "Intercept", col=type), size = 5) +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Budburst), data = meanPtE, col = "black") +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = BudburstHigh), data = meanPtE, col = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=20))+
  ylim (-20,35) +
  scale_x_continuous( breaks = east$meanBB, labels = east$species,limits = c(5,21)) +
  labs( x = "", y = "Day of budburst", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 11, y = 35, label = "c) Eastern transect - LNC", cex =8) +
  scale_color_manual(values = c("maroon","cyan4")) +
  scale_shape_discrete( labels = c("low cue",
                                   "high cue",
                                   "intercept" ),
                        breaks = c("Budburst","BudburstHigh", "Intercept"))
CNE

# CNOrderE <- ggplot(meanPtE) +
#   geom_point(aes(y= Budburst, x = CN, shape = "Budburst", col=type ), size = 5) +
#   geom_point(aes(y= BudburstHigh, x = CN, shape = "BudburstHigh", col=type), size = 5) +
#   geom_point(aes(y= Intercept, x = CN, shape = "Intercept", col=type), size = 5) +
#   geom_segment(aes(x = CN, y = Intercept, xend = CN, yend = Budburst), data = meanPtE, col = "black") +
#   geom_segment(aes(x = CN, y = Intercept, xend = CN, yend = BudburstHigh), data = meanPtE, col = "black") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),
#     legend.key=element_rect(fill="white")) +
#   theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
#     axis.text.y=element_text(size = 15),
#     axis.title=element_text(size=20)) +
#   scale_x_continuous( breaks = east$CN, labels = east$species,limits = c(0.021,0.063)) +
#   labs( x = "", y = "Day of budburst (days/standardized units)", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 0.038, y = 17, label = "b) Eastern transect - CN", cex =8) +
#   scale_color_manual(values = c("maroon","cyan4")) +
#   scale_shape_discrete( labels = c("low cue",
#     "high cue",
#     "intercept" ),
#     breaks = c("Budburst","BudburstHigh", "Intercept"))

###################
west <- subset(spInfo, transect != "east")
westSp <- unique(west$species.name)

datawest <- data[data$species.name %in% westSp,]

meanPtW <- aggregate(datawest[c("meanBB", "meanBBHigh","Int")], datawest[c("species.name","type","transect")], FUN = mean)
names(meanPtW) <- c("species.name","type","transect","Budburst", "BudburstHigh","Intercept")

CNW <- ggplot(meanPtW) +
  geom_point(aes(y= Budburst, x = Budburst, shape = "Budburst", col=type ), size = 5) +
  geom_point(aes(y= BudburstHigh, x = Budburst, shape = "BudburstHigh", col=type), size = 5) +
  geom_point(aes(y= Intercept, x = Budburst, shape = "Intercept", col=type), size = 5) +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Budburst), data = meanPtW, col = "black") +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = BudburstHigh), data = meanPtW, col = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=20)) +
  scale_x_continuous( breaks = west$meanBB, labels = west$species,limits = c(8,24)) +
  labs( x = "", y = "Day of budburst", main = NA) +
  theme(legend.title = element_blank()) +#  annotate("text", x = 14, y = 35, label = "d) western transect - LNC", cex =8) +
  scale_color_manual(values = c("maroon","cyan4")) +
  scale_shape_discrete( labels = c("low cue",
                                   "high cue",
                                   "intercept" ),
                        breaks = c("Budburst","BudburstHigh", "Intercept"))

# CNOrderW <- ggplot(meanPtW) +
#   geom_point(aes(y= Budburst, x = CN, shape = "Budburst", col=type ), size = 5) +
#   geom_point(aes(y= BudburstHigh, x = CN, shape = "BudburstHigh", col=type), size = 5) +
#   geom_point(aes(y= Intercept, x = CN, shape = "Intercept", col=type), size = 5) +
#   geom_segment(aes(x = CN, y = Intercept, xend = CN, yend = Budburst), data = meanPtW, col = "black") +
#   geom_segment(aes(x = CN, y = Intercept, xend = CN, yend = BudburstHigh), data = meanPtW, col = "black") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),
#     legend.key=element_rect(fill="white")) +
#   theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
#     axis.text.y=element_text(size = 15),
#     axis.title=element_text(size=20)) +
#   scale_x_continuous( breaks = west$CN, labels = west$species,limits = c(0.021,0.065)) +
#   labs( x = "", y = "Day of budburst (days/standardized units)", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 0.04, y = 20, label = "d) Western transect - CN", cex =8) +
#   scale_color_manual(values = c("maroon","cyan4")) +
#   scale_shape_discrete( labels = c("low cue",
#     "high cue",
#     "intercept" ),
#     breaks = c("Budburst","BudburstHigh", "Intercept"))
# pdf("figures/dotShrubTreeHtCNHundoa.pdf", width = 7, height = 5)
# htE
# dev.off()
# 
# pdf("figures/dotShrubTreeHtCNHundob.pdf", width = 7, height = 5)
# htW
# dev.off()
# 
# pdf("figures/dotShrubTreeHtCNHundoc.pdf", width = 7, height = 5)
# CNE
# dev.off()
# 
# pdf("figures/dotShrubTreeHtCNHundod.pdf", width = 7, height = 5)
# CNW
# dev.off()

pdf("analysis/figures/dotShrubTreeHtCNHundo.pdf", width = 15, height = 10)
plot_grid(htE, htW, CNE,CNW, nrow = 2, ncol = 2, align = "v")
dev.off()
# 
# CNOrderE <- ggplot(meanPtE) +
#   geom_point(aes(y= Budburst, x = CN, shape = "Budburst", col=type ), size = 5) +
#   geom_point(aes(y= BudburstHigh, x = CN, shape = "BudburstHigh", col=type), size = 5) +
#   geom_point(aes(y= Intercept, x = CN, shape = "Intercept", col=type), size = 5) +
#   geom_segment(aes(x = CN, y = Intercept, xend = CN, yend = Budburst), data = meanPtE, col = "black") +
#   geom_segment(aes(x = CN, y = Intercept, xend = CN, yend = BudburstHigh), data = meanPtE, col = "black") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),
#     legend.key=element_rect(fill="white")) +
#   theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
#     axis.text.y=element_text(size = 15),
#     axis.title=element_text(size=20)) +
#   ylim (-20,35) +
#   scale_x_continuous( breaks = east$CN, labels = east$species,limits = c(0.021,0.063)) +
#   labs( x = "", y = "Day of budburst (days/standardized units)", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 0.035, y = 17, label = "c) Eastern transect - CN", cex =8) +
#   scale_color_manual(values = c("maroon","cyan4")) +
#   scale_shape_discrete( labels = c("low cue",
#     "high cue",
#     "intercept" ),
#     breaks = c("Budburst","BudburstHigh", "Intercept"))

# CNOrderW <- ggplot(meanPtW) +
#   geom_point(aes(y= Budburst, x = CN, shape = "Budburst", col=type ), size = 5) +
#   geom_point(aes(y= BudburstHigh, x = CN, shape = "BudburstHigh", col=type), size = 5) +
#   geom_point(aes(y= Intercept, x = CN, shape = "Intercept", col=type), size = 5) +
#   geom_segment(aes(x = CN, y = Intercept, xend = CN, yend = Budburst), data = meanPtW, col = "black") +
#   geom_segment(aes(x = CN, y = Intercept, xend = CN, yend = BudburstHigh), data = meanPtW, col = "black") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),
#     legend.key=element_rect(fill="white")) +
#   theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
#     axis.text.y=element_text(size = 15),
#     axis.title=element_text(size=20)) +
#   scale_x_continuous( breaks = west$CN, labels = west$species,limits = c(0.021,0.065)) +
#   labs( x = "", y = "Day of budburst (days/standardized units)", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 0.037, y = 20, label = "d) Western transect - CN", cex =8) +
#   scale_color_manual(values = c("maroon","cyan4")) +
#   scale_shape_discrete( labels = c("low cue",
#     "high cue",
#     "intercept" ),
#     breaks = c("Budburst","BudburstHigh", "Intercept"))
# 
# htOrderW <- ggplot(meanPtW) +
#   geom_point(aes(y= Budburst, x = ht, shape = "Budburst", col=type ), size = 5) +
#   geom_point(aes(y= BudburstHigh, x = ht, shape = "BudburstHigh", col=type), size = 5) +
#   geom_point(aes(y= Intercept, x = ht, shape = "Intercept", col=type), size = 5) +
#   geom_segment(aes(x = ht, y = Intercept, xend = ht, yend = Budburst), data = meanPtW, col = "black") +
#   geom_segment(aes(x = ht, y = Intercept, xend = ht, yend = BudburstHigh), data = meanPtW, col = "black") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),
#     legend.key=element_rect(fill="white")) +
#   theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
#     axis.text.y=element_text(size = 15),
#     axis.title=element_text(size=20)) +
#   scale_x_continuous( breaks = west$ht, labels = west$species,limits = c(0.5,18)) +
#   labs( x = "", y = "Day of budburst (days/standardized units)", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 7, y = 20, label = "b) Western transect - height", cex =8) +
#   scale_color_manual(values = c("maroon","cyan4")) +
#   scale_shape_discrete( labels = c("low cue",
#     "high cue",
#     "intercept" ),
#     breaks = c("Budburst","BudburstHigh", "Intercept"))
# 
# htOrderE <- ggplot(meanPtE) +
#   geom_point(aes(y= Budburst, x = ht, shape = "Budburst", col=type ), size = 5) +
#   geom_point(aes(y= BudburstHigh, x = ht, shape = "BudburstHigh", col=type), size = 5) +
#   geom_point(aes(y= Intercept, x = ht, shape = "Intercept", col=type), size = 5) +
#   geom_segment(aes(x = ht, y = Intercept, xend = ht, yend = Budburst), data = meanPtE, col = "black") +
#   geom_segment(aes(x = ht, y = Intercept, xend = ht, yend = BudburstHigh), data = meanPtE, col = "black") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     panel.background = element_blank(), axis.line = element_line(colour = "black"),
#     legend.key=element_rect(fill="white")) +
#   theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
#     axis.text.y=element_text(size = 15),
#     axis.title=element_text(size=20)) +
#   scale_x_continuous( breaks = east$ht, labels = east$species,limits = c(0.5,18)) +
#   labs( x = "", y = "Day of budburst (days/standardized units)", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 7, y = 20, label = "a) Eastern transect - height", cex =8) +
#   scale_color_manual(values = c("maroon","cyan4")) +
#   scale_shape_discrete( labels = c("low cue",
#     "high cue",
#     "intercept" ),
#     breaks = c("Budburst","BudburstHigh", "Intercept"))
# 
# pdf("figures/dotShrubTreeHtCNTraitOrder.pdf", width = 15, height = 10)
# plot_grid(htOrderE, htOrderW, CNOrderE, CNOrderW, nrow = 2, ncol = 2, align = "v")
# dev.off()

# 
# htE <- ggplot(meanPtE) +
#   geom_point(aes(y= Budburst, x = Budburst, shape = "Budburst", col=type ), size = 5) +
#   geom_point(aes(y= BudburstHigh, x = Budburst, shape = "BudburstHigh", col=type), size = 5) +
#   geom_point(aes(y= Intercept, x = Budburst, shape = "Intercept", col=type), size = 5) +
#   geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Budburst), data = meanPtE, col = "black") +
#   geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = BudburstHigh), data = meanPtE, col = "black") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         legend.key=element_rect(fill="white")) +
#   theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
#         axis.text.y=element_text(size = 15),
#         axis.title=element_text(size=20)) +
#   scale_x_continuous( breaks = east$meanBB, labels = east$species,limits = c(5,22)) +
#   labs( x = "", y = "Day of budburst (days/standardized units)", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 11, y = 20, label = "a) Eastern transect - height", cex =8) +
#   scale_color_manual(values = c("maroon","cyan4")) +
#   scale_shape_discrete( labels = c("low cue",
#                                    "high cue",
#                                    "intercept" ),
#                         breaks = c("Budburst","BudburstHigh", "Intercept"))
# 
# 
# htW <- ggplot(meanPtW) +
#   geom_point(aes(y= Budburst, x = Budburst, shape = "Budburst", col=type ), size = 5) +
#   geom_point(aes(y= BudburstHigh, x = Budburst, shape = "BudburstHigh", col=type), size = 5) +
#   geom_point(aes(y= Intercept, x = Budburst, shape = "Intercept", col=type), size = 5) +
#   geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Budburst), data = meanPtW, col = "black") +
#   geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = BudburstHigh), data = meanPtW, col = "black") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         legend.key=element_rect(fill="white")) +
#   theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
#         axis.text.y=element_text(size = 15),
#         axis.title=element_text(size=20)) +
#   scale_x_continuous( breaks = west$meanBB, labels = west$species,limits = c(5,20)) +
#   labs( x = "", y = "Day of budburst (days/standardized units)", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 10.5, y = 30, label = "c) western transect - height", cex =8) +
#   scale_color_manual(values = c("maroon","cyan4")) +
#   scale_shape_discrete( labels = c("low cue",
#                                    "high cue",
#                                    "intercept" ),
#                         breaks = c("Budburst","BudburstHigh", "Intercept"))
# 
# CNE <- ggplot(meanPtE) +
#   geom_point(aes(y= Budburst, x = Budburst, shape = "Budburst", col=type ), size = 5) +
#   geom_point(aes(y= BudburstHigh, x = Budburst, shape = "BudburstHigh", col=type), size = 5) +
#   geom_point(aes(y= Intercept, x = Budburst, shape = "Intercept", col=type), size = 5) +
#   geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Budburst), data = meanPtE, col = "black") +
#   geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = BudburstHigh), data = meanPtE, col = "black") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         legend.key=element_rect(fill="white")) +
#   theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
#         axis.text.y=element_text(size = 15),
#         axis.title=element_text(size=20)) +
#   scale_x_continuous( breaks = east$meanBB, labels = east$species,limits = c(4,22)) +
#   labs( x = "", y = "Day of budburst (days/standardized units)", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 10.5, y = 30, label = "b) Eastern transect - CN", cex =8) +
#   scale_color_manual(values = c("maroon","cyan4")) +
#   scale_shape_discrete( labels = c("low cue",
#                                    "high cue",
#                                    "intercept" ),
#                         breaks = c("Budburst","BudburstHigh", "Intercept"))
# 
# CNW <- ggplot(meanPtW) +
#   geom_point(aes(y= Budburst, x = Budburst, shape = "Budburst", col=type ), size = 5) +
#   geom_point(aes(y= BudburstHigh, x = Budburst, shape = "BudburstHigh", col=type), size = 5) +
#   geom_point(aes(y= Intercept, x = Budburst, shape = "Intercept", col=type), size = 5) +
#   geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Budburst), data = meanPtW, col = "black") +
#   geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = BudburstHigh), data = meanPtW, col = "black") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         legend.key=element_rect(fill="white")) +
#   theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
#         axis.text.y=element_text(size = 15),
#         axis.title=element_text(size=20)) +
#   scale_x_continuous( breaks = west$meanBB, labels = west$species,limits = c(4,22)) +
#   labs( x = "", y = "Day of budburst (days/standardized units)", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 14.5, y = 30, label = "d) western transect - CN", cex =8) +
#   scale_color_manual(values = c("maroon","cyan4")) +
#   scale_shape_discrete( labels = c("low cue",
#                                    "high cue",
#                                    "intercept" ),
#                         breaks = c("Budburst","BudburstHigh", "Intercept"))
# 
# pdf("figures/dotShrubTree.pdf", width = 15, height = 10)
# plot_grid(htE,CNE,htW, CNW, nrow = 2, ncol = 2, align = "v")
# dev.off()
# 
# ##################################################
# # do different traits alter the timing of spp in the community?
# 
# rankCN <- spInfo[,c("species.name","species","type","transect","meanBB","meanBBHigh","Int")]
# rankCN <- rankCN[order(rankCN$Int),]
# 
# rankCNE <- subset(rankCN, transect == "east")
# rankCNE <- rankCNE[order(rankCNE$Int),]
# rankCNE$rankCNInt <- seq(1:nrow(rankCNE))
# 
# rankCNE <- rankCNE[order(rankCNE$meanBB),]
# rankCNE$rankCNLowC <- seq(1:nrow(rankCNE))
# 
# rankCNE <- rankCNE[order(rankCNE$meanBBHigh),]
# rankCNE$rankCNHighC <- seq(1:nrow(rankCNE))
# 
# rankCNW <- subset(rankCN, transect == "west")
# rankCNW <- rankCNW[order(rankCNW$Int),]
# rankCNW$rankCNInt <- seq(1:nrow(rankCNW))
# 
# rankCNW <- rankCNW[order(rankCNW$meanBB),]
# rankCNW$rankCNLowC <- seq(1:nrow(rankCNW))
# 
# rankCNW <- rankCNW[order(rankCNW$meanBBHigh),]
# rankCNW$rankCNHighC <- seq(1:nrow(rankCNW))
# 
# eastRank <- merge(rankE, rankCNE, by = c("species.name", "species", "type","transect"))
# 
# pdf("figures/rankEstiBB.pdf", width = 9, height =3)
# colTran <- c("maroon","navy","forestgreen")
# par(mfrow = c(1,3), mar = c(5.1, 4.8, 4.1, 2.1))
# plot(eastRank$rankLowC~eastRank$rankCNLowC, 
#      col = "darkolivegreen",
#      pch = 19,
#      xlab = "Height low cue rank",
#      ylab = "CN low cue rank",
#      cex.lab =1.5,
#      cex =1.5)
# abline(0,1)
# 
# plot(eastRank$rankHighC~eastRank$rankCNHighC, 
#      col = "darkolivegreen",
#      pch = 19,
#      xlab = "Height high cue rank",
#      ylab = "CN high cue rank",
#      cex.lab =1.5,
#      cex =1.5)
# abline(0,1)
# 
# plot(eastRank$rankInt~eastRank$rankCNInt, 
#      col = "darkolivegreen",
#      pch = 19,
#      xlab = "Height high cue rank",
#      ylab = "CN high cue rank",
#      cex.lab =1.5,
#      cex =1.5)
# abline(0,1)
# 
# ###### Western ########
# westRank <- merge(rankW, rankCNW, by = c("species.name", "species", "type","transect"))
# 
# pdf("figures/rankEstiBB.pdf", width = 9, height =3)
# colTran <- c("maroon","navy","forestgreen")
# par(mfrow = c(1,3), mar = c(5.1, 4.8, 4.1, 2.1))
# plot(westRank$rankLowC~westRank$rankCNLowC, 
#      col = "darkolivegreen",
#      pch = 19,
#      xlab = "Height low cue rank",
#      ylab = "CN low cue rank",
#      cex.lab =1.5,
#      cex =1.5)
# abline(0,1)
# 
# plot(westRank$rankHighC~westRank$rankCNHighC, 
#      col = "darkolivegreen",
#      pch = 19,
#      xlab = "Height high cue rank",
#      ylab = "CN high cue rank",
#      cex.lab =1.5,
#      cex =1.5)
# abline(0,1)
# 
# plot(westRank$rankInt~westRank$rankCNInt, 
#      col = "darkolivegreen",
#      pch = 19,
#      xlab = "Height high cue rank",
#      ylab = "CN high cue rank",
#      cex.lab =1.5,
#      cex =1.5)
# abline(0,1)