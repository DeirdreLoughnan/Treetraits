# Started January 2025 by D. Loughnan
# aim of this codeput.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, 
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

put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, 
                           offset=c(0, 0), ...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.05,0.98),
                     topcenter = c(0.5525,0.98),
                     topright = c(0.985, 0.98),
                     bottomleft = c(0.015, 0.02), 
                     bottomcenter = c(0.5525, 0.02), 
                     bottomright = c(0.985, 0.02),
                     c(0.015, 0.98) )
  } else {
    coords <- c(x,y)
  }
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  text(labels=label[1], x=this.x, y=this.y, xpd=T, cex = 2, ...)
}

spInfo <- read.csv("analysis/input/species_ring.csv")
trtPheno <- read.csv("analysis/input/trtPhenoDummy.csv")
specieslist <- sort(unique(trtPheno$species))

load("analysis/output/ssdContLatHundoLatFinal.Rdata")

ModelFit <- rstan::extract(mdlSSD)

muSp <- data.frame(ModelFit$mu_grand_sp)
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
a5 <- quantile(muForceSp$ModelFit.muForceSp, prob = c(0.25))
a95 <- quantile(muForceSp$ModelFit.muForceSp, prob = c(0.75))
apoly <- subset(muForceSp, ModelFit.muForceSp > a5 & ModelFit.muForceSp < a95)

betaTraitxForce<- data.frame(ModelFit$betaTraitxForce)
b5 <- quantile(betaTraitxForce$ModelFit.betaTraitxForce, prob = c(0.25))
b95 <- quantile(betaTraitxForce$ModelFit.betaTraitxForce, prob = c(0.75))
bpolly <- subset(betaTraitxForce, ModelFit.betaTraitxForce > b5 & ModelFit.betaTraitxForce < b95)
betaTraitxForceMean <- colMeans(betaTraitxForce)

diffuse <- unique(subset(spInfo, ring.type == "Diffuse")$species)
unknown <- unique(subset(spInfo, is.na(ring.type))$species)
diffRing <- unique(subset(spInfo, ring.type == "Diffuse/semi-ring")$species)
semiRing <- unique(subset(spInfo, ring.type == "Semi-ring")$species)
ring <- unique(subset(spInfo, ring.type == "Ring")$species)


mg_df_diffuse <- mg_df[mg_df$species %in% diffuse, ]
mg_df_unknown <- mg_df[mg_df$species %in% unknown, ]
mg_df_diffRing <- mg_df[mg_df$species %in% diffRing, ]
mg_df_semiRing <- mg_df[mg_df$species %in% semiRing, ]
mg_df_ring <- mg_df[mg_df$species %in% ring, ]

bfs_df_diffuse <- bfs_df[bfs_df$species %in% diffuse, ]
bfs_df_unknown <- bfs_df[bfs_df$species %in% unknown, ]
bfs_df_diffRing <- bfs_df[bfs_df$species %in% diffRing, ]
bfs_df_semiRing <- bfs_df[bfs_df$species %in% semiRing, ]
bfs_df_ring <- bfs_df[bfs_df$species %in% ring, ]

pdf("analysis/figures/cuetraitRingType_Known.pdf", height = 4, width = 15)
par(mar = c(5, 5, 2, 2), mfrow = c(1,3))
plot( x= mg_df$muSpMean, y = bfs_df$betaForceSpMean, 
      type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bfs_df$force25), max(bfs_df$force75)), 
      ylab = "Species level forcing slope", xlab = bquote('Wood specific density'~(g/cm^2)), cex.lab = 2.25, cex.axis = 2) # blank plot with x range 
text ( 25,-0.5, "a)", cex =2)
# 3 columns, mean, quantile
# min and max defined by quantiles
#mtext(side = 3, text = "Forcing", adj = 0, cex = 1.5)

for(j in 1:length(apoly[,1])){   
  abline(a = apoly[j,], b = bpolly[j,], col= "#E0E0E0") 
}

# for(j in 1:length(muForceSp[,1])){
#   abline(a = muForceSp[j,], b = betaTraitxForceMean, col=alpha("#73d2de", 0.085))
# }
abline(a=muForceSpMean, b=betaTraitxForceMean, col = "black")


arrows(
  mg_df_diffuse[,"muSpMean"], # x mean
  bfs_df_diffuse[,"force25"], # y 25
  mg_df_diffuse[,"muSpMean"],
  bfs_df_diffuse[,"force75"],
  length = 0, col= "#218380", lwd = 4
)

arrows(
  mg_df_diffuse[,"trait25"], # x mean
  bfs_df_diffuse[,"betaForceSpMean"], # y 25
  mg_df_diffuse[,"trait75"], # x mean
  bfs_df_diffuse[,"betaForceSpMean"],
  length = 0, col = "#218380", lwd = 4
)

arrows(
  mg_df_diffRing[,"muSpMean"], # x mean
  bfs_df_diffRing[,"force25"], # y 25
  mg_df_diffRing[,"muSpMean"],
  bfs_df_diffRing[,"force75"],
  length = 0, col= "#8f2d56", lwd = 4 
)

arrows(
  mg_df_diffRing[,"trait25"], # x mean
  bfs_df_diffRing[,"betaForceSpMean"], # y 25
  mg_df_diffRing[,"trait75"], # x mean
  bfs_df_diffRing[,"betaForceSpMean"],
  length = 0, col = "#8f2d56", lwd = 4
)

# arrows(
#   mg_df_unknown[,"muSpMean"], # x mean
#   bfs_df_unknown[,"force25"], # y 25
#   mg_df_unknown[,"muSpMean"],
#   bfs_df_unknown[,"force75"],
#   length = 0, col= "green4", lwd = 4
# )
# 
# arrows(
#   mg_df_unknown[,"trait25"], # x mean
#   bfs_df_unknown[,"betaForceSpMean"], # y 25
#   mg_df_unknown[,"trait75"], # x mean
#   bfs_df_unknown[,"betaForceSpMean"],
#   length = 0, col = "green4", lwd = 4
# )

arrows(
  mg_df_semiRing[,"muSpMean"], # x mean
  bfs_df_semiRing[,"force25"], # y 25
  mg_df_semiRing[,"muSpMean"],
  bfs_df_semiRing[,"force75"],
  length = 0, col= "purple4", lwd = 4 
)

arrows(
  mg_df_semiRing[,"trait25"], # x mean
  bfs_df_semiRing[,"betaForceSpMean"], # y 25
  mg_df_semiRing[,"trait75"], # x mean
  bfs_df_semiRing[,"betaForceSpMean"],
  length = 0, col = "purple4", lwd = 4
)

arrows(
  mg_df_ring[,"muSpMean"], # x mean
  bfs_df_ring[,"force25"], # y 25
  mg_df_ring[,"muSpMean"],
  bfs_df_ring[,"force75"],
  length = 0, col= "goldenrod", lwd = 4
)

arrows(
  mg_df_ring[,"trait25"], # x mean
  bfs_df_ring[,"betaForceSpMean"], # y 25
  mg_df_ring[,"trait75"], # x mean
  bfs_df_ring[,"betaForceSpMean"],
  length = 0, col = "goldenrod", lwd = 4
)


# my.label <- paste("j", ".", sep="")
# put.fig.letter(label=my.label, location= "topleft", font=2)
#dev.off()
###############################################################
betaChillSp <- data.frame(ModelFit$betaChillSp)
betaChillSpMean <- colMeans(betaChillSp)
bc_quan <- apply(betaChillSp, 2, quantile2575)

bcs <- rbind(betaChillSpMean, bc_quan)
bcs_t <- t(bcs)
bcs_df <- data.frame(bcs_t)
colnames(bcs_df)[colnames(bcs_df) == "X25."] <- "chill25"
colnames(bcs_df)[colnames(bcs_df) == "X75."] <- "chill75"
bcs_df$species <- specieslist

muChillSp <- data.frame(ModelFit$muChillSp)
muChillSpMean <- colMeans(muChillSp)
a5 <- quantile(muChillSp$ModelFit.muChillSp, prob = c(0.25))
a95 <- quantile(muChillSp$ModelFit.muChillSp, prob = c(0.75))
apoly <- subset(muChillSp, ModelFit.muChillSp > a5 & ModelFit.muChillSp < a95)

betaTraitxChill<- data.frame(ModelFit$betaTraitxChill)
b5 <- quantile(betaTraitxChill$ModelFit.betaTraitxChill, prob = c(0.25))
b95 <- quantile(betaTraitxChill$ModelFit.betaTraitxChill, prob = c(0.75))
bpolly <- subset(betaTraitxChill, ModelFit.betaTraitxChill > b5 & ModelFit.betaTraitxChill < b95)
betaTraitxChillMean <- colMeans(betaTraitxChill)

mg_df_diffuse <- mg_df[mg_df$species %in% diffuse, ]
mg_df_unknown <- mg_df[mg_df$species %in% unknown, ]
mg_df_diffRing <- mg_df[mg_df$species %in% diffRing, ]
mg_df_semiRing <- mg_df[mg_df$species %in% semiRing, ]
mg_df_ring <- mg_df[mg_df$species %in% ring, ]

bcs_df_diffuse <- bcs_df[bcs_df$species %in% diffuse, ]
bcs_df_unknown <- bcs_df[bcs_df$species %in% unknown, ]
bcs_df_diffRing <- bcs_df[bcs_df$species %in% diffRing, ]
bcs_df_semiRing <- bcs_df[bcs_df$species %in% semiRing, ]
bcs_df_ring <- bcs_df[bcs_df$species %in% ring, ]

#pdf("figures/cuetraitHundok.pdf", height = 4, width = 5)
#par(mar = c(5, 5, 2, 2), mfrow = c(1,1))
plot( x= mg_df$muSpMean, y = bcs_df$betaChillSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), 
      ylim = c(min(bcs_df$chill25), max(bcs_df$chill75)), 
      ylab = "Species level chilling slope", 
      xlab = bquote('Wood specific density'~(g/cm^2)),
      cex.lab = 2.25, cex.axis = 2) # blank plot with x range 
text ( 25,8, "b)", cex =2)
# 3 columns, mean, quantile
# min and max defined by quantiles
#mtext(side = 3, text = "Chilling", adj = 0, cex = 1.5)
for(j in 1:length(apoly[,1])){   
  abline(a = apoly[j,], b = bpolly[j,], col= "#E0E0E0") 
}

# for(j in 1:length(muChillSp[,1])){
#   abline(a = muChillSp[j,], b = betaTraitxChillMean, col=alpha("#73d2de", 0.085))
# }
abline(a=muChillSpMean, b=betaTraitxChillMean, col = "black")


arrows(
  mg_df_diffuse[,"muSpMean"], # x mean
  bcs_df_diffuse[,"chill25"], # y 25
  mg_df_diffuse[,"muSpMean"],
  bcs_df_diffuse[,"chill75"],
  length = 0, col= "#218380", lwd = 4
)

arrows(
  mg_df_diffuse[,"trait25"], # x mean
  bcs_df_diffuse[,"betaChillSpMean"], # y 25
  mg_df_diffuse[,"trait75"], # x mean
  bcs_df_diffuse[,"betaChillSpMean"],
  length = 0, col = "#218380", lwd = 4
)

arrows(
  mg_df_diffRing[,"muSpMean"], # x mean
  bcs_df_diffRing[,"chill25"], # y 25
  mg_df_diffRing[,"muSpMean"],
  bcs_df_diffRing[,"chill75"],
  length = 0, col= "#8f2d56", lwd = 4 
)

arrows(
  mg_df_diffRing[,"trait25"], # x mean
  bcs_df_diffRing[,"betaChillSpMean"], # y 25
  mg_df_diffRing[,"trait75"], # x mean
  bcs_df_diffRing[,"betaChillSpMean"],
  length = 0, col = "#8f2d56", lwd = 4
)

# arrows(
#   mg_df_unknown[,"muSpMean"], # x mean
#   bcs_df_unknown[,"chill25"], # y 25
#   mg_df_unknown[,"muSpMean"],
#   bcs_df_unknown[,"chill75"],
#   length = 0, col= "green4", lwd = 4
# )
# 
# arrows(
#   mg_df_unknown[,"trait25"], # x mean
#   bcs_df_unknown[,"betaChillSpMean"], # y 25
#   mg_df_unknown[,"trait75"], # x mean
#   bcs_df_unknown[,"betaChillSpMean"],
#   length = 0, col = "green4", lwd = 4
# )

arrows(
  mg_df_semiRing[,"muSpMean"], # x mean
  bcs_df_semiRing[,"chill25"], # y 25
  mg_df_semiRing[,"muSpMean"],
  bcs_df_semiRing[,"chill75"],
  length = 0, col= "purple4", lwd = 4 
)

arrows(
  mg_df_semiRing[,"trait25"], # x mean
  bcs_df_semiRing[,"betaChillSpMean"], # y 25
  mg_df_semiRing[,"trait75"], # x mean
  bcs_df_semiRing[,"betaChillSpMean"],
  length = 0, col = "purple4", lwd = 4
)

arrows(
  mg_df_ring[,"muSpMean"], # x mean
  bcs_df_ring[,"chill25"], # y 25
  mg_df_ring[,"muSpMean"],
  bcs_df_ring[,"chill75"],
  length = 0, col= "goldenrod", lwd = 4
)

arrows(
  mg_df_ring[,"trait25"], # x mean
  bcs_df_ring[,"betaChillSpMean"], # y 25
  mg_df_ring[,"trait75"], # x mean
  bcs_df_ring[,"betaChillSpMean"],
  length = 0, col = "goldenrod", lwd = 4
)
# my.label <- paste("k", ".", sep="")
# put.fig.letter(label=my.label, location= "topleft", font=2)
#dev.off()
###############################################################
betaPhotoSp <- data.frame(ModelFit$betaPhotoSp)
betaPhotoSpMean <- colMeans(betaPhotoSp)
bp_quan <- apply(betaPhotoSp, 2, quantile2575)

bps <- rbind(betaPhotoSpMean, bp_quan)
bps_t <- t(bps)
bps_df <- data.frame(bps_t)
colnames(bps_df)[colnames(bps_df) == "X25."] <- "photo25"
colnames(bps_df)[colnames(bps_df) == "X75."] <- "photo75"
bps_df$species <- specieslist

muPhotoSp <- data.frame(ModelFit$muPhotoSp)
muPhotoSpMean <- colMeans(muPhotoSp)
a5 <- quantile(muPhotoSp$ModelFit.muPhotoSp, prob = c(0.25))
a95 <- quantile(muPhotoSp$ModelFit.muPhotoSp, prob = c(0.75))
apoly <- subset(muPhotoSp, ModelFit.muPhotoSp > a5 & ModelFit.muPhotoSp < a95)

betaTraitxPhoto<- data.frame(ModelFit$betaTraitxPhoto)
b5 <- quantile(betaTraitxPhoto$ModelFit.betaTraitxPhoto, prob = c(0.25))
b95 <- quantile(betaTraitxPhoto$ModelFit.betaTraitxPhoto, prob = c(0.75))
bpolly <- subset(betaTraitxPhoto, ModelFit.betaTraitxPhoto > b5 & ModelFit.betaTraitxPhoto < b95)
betaTraitxPhotoMean <- colMeans(betaTraitxPhoto)

mg_df_diffuse <- mg_df[mg_df$species %in% diffuse, ]
mg_df_unknown <- mg_df[mg_df$species %in% unknown, ]
mg_df_diffRing <- mg_df[mg_df$species %in% diffRing, ]
mg_df_semiRing <- mg_df[mg_df$species %in% semiRing, ]
mg_df_ring <- mg_df[mg_df$species %in% ring, ]

bps_df_diffuse <- bps_df[bps_df$species %in% diffuse, ]
bps_df_unknown <- bps_df[bps_df$species %in% unknown, ]
bps_df_diffRing <- bps_df[bps_df$species %in% diffRing, ]
bps_df_semiRing <- bps_df[bps_df$species %in% semiRing, ]
bps_df_ring <- bps_df[bps_df$species %in% ring, ]

#pdf("figures/cuetraitHundol.pdf", height = 4, width = 5)
#par(mar = c(5, 5, 2, 2), mfrow = c(1,1))
plot( x= mg_df$muSpMean, y = bps_df$betaPhotoSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bps_df$photo25), 5), ylab = "Species level photoperiod slope", xlab = bquote('Wood specific density'~(g/cm^2)), cex.lab = 2.25, cex.axis = 2) # blank plot with x range 
text ( 25,4.5, "c)", cex =2)# min and max defined by quantiles
#mtext(side = 3, text = "Photoperiod", adj = 0, cex = 1.5)
for(j in 1:length(apoly[,1])){   
  abline(a = apoly[j,], b = bpolly[j,], col = "#E0E0E0") 
}

# for(j in 1:length(muPhotoSp[,1])){
#   abline(a = muPhotoSp[j,], b = betaTraitxPhotoMean, col=alpha("#73d2de", 0.085))
# }
abline(a=muPhotoSpMean, b=betaTraitxPhotoMean, col = "black")


arrows(
  mg_df_diffuse[,"muSpMean"], # x mean
  bps_df_diffuse[,"photo25"], # y 25
  mg_df_diffuse[,"muSpMean"],
  bps_df_diffuse[,"photo75"],
  length = 0, col= "#218380", lwd = 4
)

arrows(
  mg_df_diffuse[,"trait25"], # x mean
  bps_df_diffuse[,"betaPhotoSpMean"], # y 25
  mg_df_diffuse[,"trait75"], # x mean
  bps_df_diffuse[,"betaPhotoSpMean"],
  length = 0, col = "#218380", lwd = 4
)

arrows(
  mg_df_diffRing[,"muSpMean"], # x mean
  bps_df_diffRing[,"photo25"], # y 25
  mg_df_diffRing[,"muSpMean"],
  bps_df_diffRing[,"photo75"],
  length = 0, col= "#8f2d56", lwd = 4 
)

arrows(
  mg_df_diffRing[,"trait25"], # x mean
  bps_df_diffRing[,"betaPhotoSpMean"], # y 25
  mg_df_diffRing[,"trait75"], # x mean
  bps_df_diffRing[,"betaPhotoSpMean"],
  length = 0, col = "#8f2d56", lwd = 4
)

# arrows(
#   mg_df_unknown[,"muSpMean"], # x mean
#   bps_df_unknown[,"photo25"], # y 25
#   mg_df_unknown[,"muSpMean"],
#   bps_df_unknown[,"photo75"],
#   length = 0, col= "green4", lwd = 4
# )
# 
# arrows(
#   mg_df_unknown[,"trait25"], # x mean
#   bps_df_unknown[,"betaPhotoSpMean"], # y 25
#   mg_df_unknown[,"trait75"], # x mean
#   bps_df_unknown[,"betaPhotoSpMean"],
#   length = 0, col = "green4", lwd = 4
# )

arrows(
  mg_df_semiRing[,"muSpMean"], # x mean
  bps_df_semiRing[,"photo25"], # y 25
  mg_df_semiRing[,"muSpMean"],
  bps_df_semiRing[,"photo75"],
  length = 0, col= "purple4", lwd = 4 
)

arrows(
  mg_df_semiRing[,"trait25"], # x mean
  bps_df_semiRing[,"betaPhotoSpMean"], # y 25
  mg_df_semiRing[,"trait75"], # x mean
  bps_df_semiRing[,"betaPhotoSpMean"],
  length = 0, col = "purple4", lwd = 4
)

arrows(
  mg_df_ring[,"muSpMean"], # x mean
  bps_df_ring[,"photo25"], # y 25
  mg_df_ring[,"muSpMean"],
  bps_df_ring[,"photo75"],
  length = 0, col= "goldenrod", lwd = 4
)

arrows(
  mg_df_ring[,"trait25"], # x mean
  bps_df_ring[,"betaPhotoSpMean"], # y 25
  mg_df_ring[,"trait75"], # x mean
  bps_df_ring[,"betaPhotoSpMean"],
  length = 0, col = "goldenrod", lwd = 4
)

legend(53, 5.5,legend = c("Ring", "Semi-ring", "Diffuse", "Diffuse-ring"),
       col  = c( "goldenrod", "purple4", "#218380","#8f2d56"), lwd =4, bty = "n", cex = 1.5)

dev.off()
