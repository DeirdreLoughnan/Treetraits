# want to get a sense of how much trait-cue rel contribute to spp overall bb

sumHt <- summary(mdlHt)$summary
muGrand = (sumHt[grep("mu_grand", rownames(sumHt)), 1])
b_trtSpHt = (sumHt[grep("b_muSp", rownames(sumHt)), 1])
a_trtSpHt = mean((sumHt[grep("mu_grand_sp", rownames(sumHt)), 1]))
a_grandSpHt = sumHt[grep("mu_grand_sp", rownames(sumHt)), 1]
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
for(i in 1:ncol(postHt$b_chill1)){
  quantU <- round(quantile(postHt$b_chill1[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_chill5 <- rbind(b_chill5, quantU)
}
colnames(b_chill5) <- c("chill5","chill95","chill25","chill75")

b_force5 <- vector()
for(i in 1:ncol(postHt$b_warm)){
  quantU <- round(quantile(postHt$b_warm[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_force5 <- rbind(b_force5, quantU)
}
colnames(b_force5) <- c("force5","force95","force25","force75")

b_photo5 <- vector()
for(i in 1:ncol(postHt$b_photo)){
  quantU <- round(quantile(postHt$b_photo[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_photo5 <- rbind(b_photo5, quantU)
}
colnames(b_photo5) <- c("photo5","photo95","photo25","photo75")

# start with low cues:
photo <- -0.5033863 #8 h photo
tranE <- 0
lat <- 42.53150 # Harvard Forest
force <- -0.3568628 #5/15 C trt
chill <- -0.3546922 # low chill

m <- matrix(nrow = 1000, ncol = 47)


for(sp in 1:47){
  for (it in 1:nrow(m)){
    m[it,sp] <- postHt$alphaPhenoSp[it,sp]+ (postHt$mu_grand_sp[it,sp] * postHt$betaTraitxForce[it] + postHt$alphaForceSp[it,sp]) * force +
      (postHt$mu_grand_sp[it,sp] * postHt$betaTraitxChill[it] + postHt$alphaChillSp[it,sp]) * chill +
       (postHt$mu_grand_sp[it,sp] * postHt$betaTraitxPhoto[it] + postHt$alphaPhotoSp[it,sp])*photo
  }
}

mNoTrt <- matrix(nrow = 1000, ncol = 47)


for(sp in 1:47){
  for (it in 1:nrow(mNoTrt)){
    mNoTrt[it,sp] <- postHt$alphaPhenoSp[it,sp]+ (postHt$mu_grand_sp[it,sp] * 0 + postHt$alphaForceSp[it,sp]) * force +
      (postHt$mu_grand_sp[it,sp] * 0+ postHt$alphaChillSp[it,sp]) * chill +
      (postHt$mu_grand_sp[it,sp] * 0 + postHt$alphaPhotoSp[it,sp])*photo
  }
}

#################################################
spInfo <- read.csv("input/species_list.csv")

spInfo <- spInfo[order(spInfo$species),]
head(spInfo)
spInfo$meanBB <- colMeans(m)
colnames(m) <- spInfo$species.name

spInfo$meanBBNoTrt <- colMeans(mNoTrt)
colnames(mNoTrt) <- spInfo$species.name

spInfo$Int <- b_phenoSpHt

ggplot(spInfo) +
  geom_point(aes(y= meanBB, x = meanBB, colour = "Budburst"), size = 5) +
  geom_point(aes(y= Int, x = meanBB, colour = "Intercept"), size = 5) +
  geom_point(aes(y= meanBBNoTrt, x = meanBB, colour = "BudburstNoTrt"), size = 5) +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Budburst), data = meanPt, col = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.key=element_rect(fill="white")) + ylim(-1.5,75) +
  theme(axis.text.x = element_text( size=15,angle = 78,  hjust=1),
    axis.text.y=element_text(size = 15),
    axis.title=element_text(size=20), legend.position = "none") + 
  # scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,62)) +
  labs( x = "Species ordered by predicted budburst date", y = "Estimated parameter (days/standardized units)", main = NA) +
  theme(legend.title = element_blank(), legend.text = element_text(size =25), legend.position = "top") +  annotate("text", x = 25, y = 70, label = "a)      Western transect", cex =8) +
  scale_color_manual(values = c("cyan4", "firebrick4", "goldenrod"))


#spInfo <- cbind(spInfo, a_sp5,b_force5, b_chill5,b_photo5)

# spInfo$bforceSp <- b_forceSpHt
# spInfo$bchillSp <- b_chillSpHt
# spInfo$bphotoSp <- b_photoSpHt
# 
# spInfo$aforceSp <- a_forceSpHt
# spInfo$achillSp <- a_chillSpHt
# spInfo$aphotoSp <- a_photoSpHt
# 
# spInfo$trtxforce <- bTrtChillHt
# spInfo$trtxchillSp <- bTrtChillHt
# spInfo$trtxphotoSp <- bTrtPhotoHt
# 
# quantile595 <- function(x){
#   returnQuanilte <- quantile(x, prob = c(0.05, 0.95))
#   return(returnQuanilte)
# }
# 
# quantile75.25 <- function(x){
#   returnQuanilte <- quantile(x, prob = c(0.75, 0.25))
#   return(returnQuanilte)
# }
# 
# bb_quan <- apply(m, 2, quantile595)
# bb_t <- t(bb_quan)
# bb_df <- data.frame(bb_t)
# colnames(bb_df)[colnames(bb_df) == "X5."] <- "bb5"
# colnames(bb_df)[colnames(bb_df) == "X95."] <- "bb95"
# 
# bb_quan75.25 <- apply(m, 2, quantile75.25)
# bb_t75.25 <- t(bb_quan75.25)
# bb_df75.25 <- data.frame(bb_t75.25)
# colnames(bb_df75.25)[colnames(bb_df75.25) == "X75."] <- "bb75"
# colnames(bb_df75.25)[colnames(bb_df75.25) == "X25."] <- "bb25"
# colnames(bb_df)[colnames(bb_df) == "X25."] <- "bb25"
# colnames(bb_df)[colnames(bb_df) == "X75."] <- "bb75"

# bb_quanHigh <- apply(mHigh, 2, quantile595)
# bb_tHigh <- t(bb_quanHigh)
# bb_dfHigh <- data.frame(bb_tHigh)
# colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X5."] <- "bb5High"
# colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X95."] <- "bb95High"
# 
# bb_quan75.25High <- apply(mHigh, 2, quantile75.25)
# bb_t75.25High <- t(bb_quan75.25High)
# bb_df75.25High <- data.frame(bb_t75.25High)
# colnames(bb_df75.25High)[colnames(bb_df75.25High) == "X75."] <- "bb75High"
# colnames(bb_df75.25High)[colnames(bb_df75.25High) == "X25."] <- "bb25High"
# colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X25."] <- "bb25High"
# colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X75."] <- "bb75High"


# spInfo <- cbind(spInfo, bb_df)
# spInfo <- cbind(spInfo, bb_df75.25)
# spInfo$value <- spInfo$meanBB

# spInfo <- cbind(spInfo, bb_dfHigh)
# spInfo <- cbind(spInfo, bb_df75.25High)
# spInfo$valueHigh <- spInfo$meanBBHigh
spInfo <- read.csv("input/species_list.csv")

spInfo <- spInfo[order(spInfo$species),]
head(spInfo)

bTrtChillHt <- data.frame(postHt$betaTraitxChill)
bTrtForceHt <- data.frame(postHt$betaTraitxForce)
bTrtPhotoHt <- data.frame(postHt$betaTraitxPhoto)

bHtCue <- cbind(bTrtChillHt,bTrtForceHt,bTrtPhotoHt)

ggplot() +
  stat_eye(data = bHtCue, aes(x = site, y = photoSiteInter, fill = "cyan4"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
  theme_classic() +
  theme(legend.position = "none") +
  labs( x = "Site", y = "Photoperiod response", main = NA)+
  scale_fill_manual(values = c("cyan4"))

+
  geom_text(aes(label=species),hjust= 0.5, vjust= 1.5, show.legend = F) +
  geom_errorbar(aes(ymin= bChill25, ymax = bChill75), width= 0) +
  geom_errorbar(aes(xmin= bPhoto25, xmax = bPhoto75), width= 0) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.key=element_rect(fill="white")) # removed grey boxes around legends

bChill <- data.frame(postHt$alphaChillSp[1:1000,])
colnames(bChill) <- (spInfo$species.name)
longChill <- melt(bChill)
names(longChill) <- c("species.name", "chill")

long <- cbind(long, longChill$chill)

# Add forcing
bForce <- data.frame(postHt$alphaForceSp[1:1000,])
colnames(bForce) <- (spInfo$species.name)
longForce <- melt(bForce)
names(longForce) <- c("species.name", "force")

long <- cbind(long, longForce$force)

# photoperiod
bPhoto <- data.frame(postHt$alphaPhotoSp[1:1000,])
colnames(bPhoto) <- (spInfo$species.name)
longPhoto <- melt(bPhoto)
names(longPhoto) <- c("species.name", "photo")

long <- cbind(long, longPhoto$photo)

# intercept
aSp <- data.frame(postHt$alphaPhenoSp[1:1000,])
colnames(aSp) <- (spInfo$species.name)
longInt <- melt(aSp)
names(longInt) <- c("species.name", "int")

long <- cbind(long, longInt$int)

data <- long[order(long$meanBB),]
# 
# data$species.name <- factor(data$species.name, levels=unique(data$species.name) )
#data <- transform(data, variable=reorder(species.name, -meanBB) ) 

names(data) <- c("species.name","valueLow", "valueNoTrt","species","type","transect","meanBB", "Int",#"Int5","Int95","Int25","Int75","force5","force95","force25","force75","chill5", "chill95", "chill25", "chill75","photo5", "photo95", "photo25", "photo75",
  "spMeanForce", "spMeanChill", "spMeanPhoto","bb5","bb95","spacing","bb75","bb25","chill", "force","photo","intercept")


#####################################
meanPt <- aggregate(data[c("valueLow","valueNoTrt","Int")], data[c("species.name","type","transect")], FUN = mean)
names(meanPt) <- c("species.name","type","transect","Budburst","BudburstNoTrt","Intercept")

ggplot(meanPt) +
  geom_point(aes(y= Budburst, x = Budburst, colour = "Budburst"), size = 5) +
  geom_point(aes(y= Intercept, x = Budburst, colour = "Intercept"), size = 5) +
  geom_point(aes(y= BudburstNoTrt, x = Budburst, colour = "BudburstNoTrt"), size = 5) +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Budburst), data = meanPt, col = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    legend.key=element_rect(fill="white")) + ylim(-1.5,75) +
  theme(axis.text.x = element_text( size=15,angle = 78,  hjust=1),
    axis.text.y=element_text(size = 15),
    axis.title=element_text(size=20), legend.position = "none") + 
 # scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,62)) +
  labs( x = "Species ordered by predicted budburst date", y = "Estimated parameter (days/standardized units)", main = NA) +
  theme(legend.title = element_blank(), legend.text = element_text(size =25), legend.position = "top") +  annotate("text", x = 18, y = 60, label = "b)      Western transect", cex =8) +
  scale_color_manual(values = c("cyan4", "firebrick4", "goldenrod"))

