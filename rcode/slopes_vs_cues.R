rm(list=ls())
options(stringsAsFactors = FALSE)

require(dplyr)
library(stringr)
library(plyr)
library(rstan)

# Set working directory:
# Anyone else working with this code should add their info/path here
if(length(grep("deirdreloughnan", getwd())>0)) {  setwd("~/Documents/github/pheno_bc")
} else if
(length(grep("Lizzie", getwd())>0)) {   setwd("~/Documents/git/projects/treegarden/budreview/ospree/analyses/traits")
}

#Function for labels:
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

# Get the data
dl <- read.csv("input/dl_allbb.csv")

temp <- str_split_fixed(dl$trt, "_", 3); head(temp)
dl$chill<- temp[,1]
dl$photo <- temp[,2]
dl$force <- temp[,3]

dl.chill <- read.csv("input/chilling_values_Hope_Smithers.csv")

dl.wchill <- merge(dl, dl.chill, by = c("population","chill"))
dl.wchill$lab3 <- dl.wchill$lab2
dl.wchill$lab2 <- paste(dl.wchill$species, dl.wchill$population, dl.wchill$rep, sep = "_")

df <- read.csv("input/df_dxb_prepped_data.csv")
df.chill <- read.csv("input/chilling_values_eastern.csv")
df.wchill <- merge(df, df.chill, by =c("population","chill"))
df.wchill <- df.wchill[, c("population", "chill","force","photo","lab2", "bb","species", "treatment","Chill_portions","Utah_Model")]

# mergeing the my data with DF
pheno <- rbind.fill(dl.wchill, df.wchill)
pheno$force.n <- pheno$force
pheno$force.n[pheno$force.n == "HF"] <- "1"
pheno$force.n[pheno$force.n == "LF"] <- "0"
pheno$force.n <- as.numeric(pheno$force.n)

pheno$photo.n <- pheno$photo
pheno$photo.n[pheno$photo.n == "HP"] <- "1"
pheno$photo.n[pheno$photo.n == "LP"] <- "0"
pheno$photo.n <- as.numeric(pheno$photo.n)

pheno$site.n <- pheno$population
pheno$site.n[pheno$site.n == "sm"] <- "1"
pheno$site.n[pheno$site.n == "mp"] <- "2"
pheno$site.n[pheno$site.n == "HF"] <- "3"
pheno$site.n[pheno$site.n == "SH"] <- "4"
pheno$site.n <- as.numeric(pheno$site.n)

head(pheno)
#add dummy/ site level effects:
pheno <- pheno %>%
  mutate ( site2 = if_else(site.n == 2, 1, 0),
           site3 = if_else(site.n == 3, 1, 0),
           site4 = if_else(site.n == 4, 1, 0))

# standardize the 0/1 and standardize sites? 
pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)

pheno$site2.z2 <- (pheno$site2-mean(pheno$site2,na.rm=TRUE))/(sd(pheno$site2,na.rm=TRUE)*2)
pheno$site3.z2 <- (pheno$site3-mean(pheno$site3,na.rm=TRUE))/(sd(pheno$site3,na.rm=TRUE)*2)
pheno$site4.z2 <- (pheno$site4-mean(pheno$site4,na.rm=TRUE))/(sd(pheno$site4,na.rm=TRUE)*2)

#going to split it into analysis of terminal bb and lateral bb
# Starting with the terminal buds:
#pheno.term <- pheno[,c("tbb", "chill.n", "force.n", "photo.n", "site.n", "species", "lab2")]
pheno.term <- pheno[,c("bb", "force.z2", "photo.z2", "population", "species", "lab2","Utah_Model","Chill_portions","chillport.z2", "site2.z2", "site3.z2","site4.z2")]
pheno.t <- pheno.term[complete.cases(pheno.term$bb), ] # none,great!
length(unique(pheno.t$species))
#pheno.t <- pheno.term[complete.cases(pheno.term$bb), ] # 1780 rows data 
pheno.t$species <- tolower(pheno.t$species)
pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 47 bc two species occur in both transects

# Now get the trait data and subset to only include spp we have pheno data for:
setwd("..//Treetraits")
trtData <- read.csv("data/allTrt.csv", stringsAsFactors = FALSE)
head(trtData)

phenoSp <- sort(unique(pheno.t$species))

trtPheno <- trtData[trtData$species %in% phenoSp, ]
length(unique(trtPheno$species))

head(trtPheno)

trtPheno$transect <- trtPheno$site
trtPheno$transect[trtPheno$transect == "GR"] <- "east"
trtPheno$transect[trtPheno$transect == "WM"] <- "east"
trtPheno$transect[trtPheno$transect == "SH"] <- "east"
trtPheno$transect[trtPheno$transect == "HF"] <- "east"

trtPheno$transect[trtPheno$transect == "af"] <- "west"
trtPheno$transect[trtPheno$transect == "kl"] <- "west"
trtPheno$transect[trtPheno$transect == "sm"] <- "west"
trtPheno$transect[trtPheno$transect == "mp"] <- "west"


#####################################################################
# # Comparisons of trees vs shrubs:
shrubs = c("VIBLAN","RHAFRA","RHOPRI","SPIALB","VACMYR","VIBCAS", "AROMEL","ILEMUC", "KALANG", "LONCAN", "LYOLIG", "alninc","alnvir","amelan", "corsto","loninv", "menfer","rhoalb", "riblac","rubpar","samrac","shecan","sorsco","spibet","spipyr","symalb","vacmem","vibedu")
trees = tolower(c("ACEPEN", "ACERUB", "ACESAC", "ALNINC", "BETALL", "BETLEN", "BETPAP", "CORCOR", "FAGGRA", "FRANIG", "HAMVIR", "NYSSYL", "POPGRA", "PRUPEN", "QUEALB" , "QUERUB", "QUEVEL", "acegla","betpap", "poptre", "popbal"))

east <- subset(trtPheno, transect == "east")
eastSp <- unique(east$species)

west <- subset(trtPheno, transect == "west")
westSp <- unique(west$species)
#####################################################################

cnMean <- aggregate(trtPheno["C.N"], trtPheno[c("species")], FUN = mean)
htMean <- aggregate(trtPheno["ht"], trtPheno[c("species")], FUN = mean)
ssdMean <- aggregate(trtPheno["ssd"], trtPheno[c("species")], FUN = mean)
dbhMean <- aggregate(trtPheno["dbh"], trtPheno[c("species")], FUN = mean)
lmaMean <- aggregate(trtPheno["lma"], trtPheno[c("species")], FUN = mean)

# Sorted species and study list
specieslist <- sort(unique(trtPheno$species))

# Write a loop to run all the different traits plots:
################################
files <- list.files(path = "output", pattern ="_stanfit.RDS" )
files

################################
# C:N
Model <- readRDS(paste("output/", files[1], sep = ""))

ModelFit <- rstan::extract(Model)

muGrandSp <- data.frame(ModelFit$mu_grand_sp)
muGrandSpMean <- colMeans(muGrandSp)

betaForceSp <- data.frame(ModelFit$betaForceSp)
betaForceSpMean <- colMeans(betaForceSp)

quantile2575 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.25, 0.75))
  return(returnQuanilte)
}

bf_quan <- apply(betaForceSp, 2, quantile2575) 
mugrand_quan <- apply(muGrandSp, 2, quantile2575)

bfs <- rbind(betaForceSpMean, bf_quan)
bfs_t <- t(bfs)
bfs_df <- data.frame(bfs_t)
colnames(bfs_df)[colnames(bfs_df) == "X25."] <- "force25"
colnames(bfs_df)[colnames(bfs_df) == "X75."] <- "force75"
bfs_df$species <- specieslist

mg<- rbind(muGrandSpMean, mugrand_quan)
mg_t <- t(mg)
mg_df <- data.frame(mg_t)
colnames(mg_df)[colnames(mg_df) == "X25."] <- "trait25"
colnames(mg_df)[colnames(mg_df) == "X75."] <- "trait75"
mg_df$species <- specieslist

muForceSp <- data.frame(ModelFit$muForceSp)
muForceSpMean <- colMeans(muForceSp)

betaTraitxForce<- data.frame(ModelFit$betaTraitxForce)
betaTraitxForceMean <- colMeans(betaTraitxForce)

mg_df_east <- mg_df[mg_df$species %in% eastSp, ]
mg_df_west <- mg_df[mg_df$species %in% westSp, ]

bfs_df_east <- bfs_df[bfs_df$species %in% eastSp, ]
bfs_df_west <- bfs_df[bfs_df$species %in% westSp, ]

# # Add coloured ones for the two species:"Corylus_avellana" = 11, "Acer_pseudoplatanus" = 2
# col.sp <- c( rgb(149 / 255, 216 / 255, 64 / 255, alpha = 0.9), rgb(72 / 255, 38 / 255, 119 / 255, alpha = 0.8))
# col1.sp <- c( rgb(149 / 255, 216 / 255, 64 / 255, alpha = 0.2), rgb(72 / 255, 38 / 255, 119 / 255, alpha = 0.14))
# col2.sp <- c( rgb(149 / 255, 216 / 255, 64 / 255, alpha = 0.5), rgb(72 / 255, 38 / 255, 119 / 255, alpha = 0.4))

#pdf(paste("figures/force", files[i], ".pdf", sep = ""))
pdf(paste("figures/cue", "trait_wtrend", ".pdf", sep = ""), height = 16, width = 12)
par(mar = c(5, 5, 2, 2), mfrow = c(5,3))
plot( x= mg_df$muGrandSpMean, y = bfs_df$betaForceSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bfs_df$force25), max(bfs_df$force75)), ylab = "Species level forcing slope", xlab = "Trait value", cex.lab = 1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df_east[,"muGrandSpMean"], # x mean
  bfs_df_east[,"force25"], # y 25
  mg_df_east[,"muGrandSpMean"],
  bfs_df_east[,"force75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df_east[,"trait25"], # x mean
  bfs_df_east[,"betaForceSpMean"], # y 25
  mg_df_east[,"trait75"], # x mean
  bfs_df_east[,"betaForceSpMean"],
  length = 0, col = "maroon", lwd = 2
)

arrows(
  mg_df_west[,"muGrandSpMean"], # x mean
  bfs_df_west[,"force25"], # y 25
  mg_df_west[,"muGrandSpMean"],
  bfs_df_west[,"force75"],
  length = 0, col= "darkslategray4", lwd = 2
)

arrows(
  mg_df_west[,"trait25"], # x mean
  bfs_df_west[,"betaForceSpMean"], # y 25
  mg_df_west[,"trait75"], # x mean
  bfs_df_west[,"betaForceSpMean"],
  length = 0, col = "darkslategray4", lwd = 2
)

mtext(side = 3, text = "CN, Forcing", adj = 0, cex = 1.25)
for(j in 1:length(muForceSp[,1])){
  abline(a = muForceSp[j,], b = betaTraitxForceMean, col=alpha("darkslategray4", 0.015))
}
abline(a=muForceSpMean, b=betaTraitxForceMean, col = "black")

my.label <- paste("a", ".", sep="")
put.fig.letter(label=my.label, location= "topleft", font=2)

#------------------------------------------------------------------------------#
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

betaTraitxChill<- data.frame(ModelFit$betaTraitxChill)
betaTraitxChillMean <- colMeans(betaTraitxChill)

mg_df_east <- mg_df[mg_df$species %in% eastSp, ]
mg_df_west <- mg_df[mg_df$species %in% westSp, ]

bcs_df_east <- bcs_df[bcs_df$species %in% eastSp, ]
bcs_df_west <- bcs_df[bcs_df$species %in% westSp, ]

#pdf(paste("figures/chill", files[i], ".pdf", sep = ""))
#pdf(paste("figures/chill", "lnc", ".pdf", sep = ""), height = 5, width = 5)
plot( x= mg_df$muGrandSpMean, y = bcs_df$betaChillSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bcs_df$chill25), max(bcs_df$chill75)), ylab = "Species level chilling slope", xlab = "Trait value", cex.lab = 1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles

arrows(
  mg_df_east[,"muGrandSpMean"], # x mean
  bcs_df_east[,"chill25"], # y 25
  mg_df_east[,"muGrandSpMean"],
  bcs_df_east[,"chill75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df_east[,"trait25"], # x mean
  bcs_df_east[,"betaChillSpMean"], # y 25
  mg_df_east[,"trait75"], # x mean
  bcs_df_east[,"betaChillSpMean"],
  length = 0, col = "maroon", lwd = 2
)

arrows(
  mg_df_west[,"muGrandSpMean"], # x mean
  bcs_df_west[,"chill25"], # y 25
  mg_df_west[,"muGrandSpMean"],
  bcs_df_west[,"chill75"],
  length = 0, col= "darkslategray4", lwd = 2
)

arrows(
  mg_df_west[,"trait25"], # x mean
  bcs_df_west[,"betaChillSpMean"], # y 25
  mg_df_west[,"trait75"], # x mean
  bcs_df_west[,"betaChillSpMean"],
  length = 0, col = "darkslategray4", lwd = 2
)

mtext(side = 3, text = "C:N, Chilling", adj = 0, cex = 1.25)
for(j in 1:length(muChillSp[,1])){
  abline(a = muChillSp[j,], b = betaTraitxChillMean, col=alpha("darkslategray4", 0.015))
}
abline(a=muChillSpMean, b=betaTraitxChillMean, col = "black")

my.label <- paste("b", ".", sep="")
put.fig.letter(label=my.label, location= "topleft", font=2)
#dev.off()
#------------------------------------------------------------------------------#
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

betaTraitxPhoto<- data.frame(ModelFit$betaTraitxPhoto)
betaTraitxPhotoMean <- colMeans(betaTraitxPhoto)

bps_df_east <- bps_df[bps_df$species %in% eastSp, ]
bps_df_west <- bps_df[bps_df$species %in% westSp, ]

#pdf(paste("figures/photo", files[i], ".pdf", sep = ""))
#pdf(paste("figures/photo", "height", ".pdf", sep = ""), height = 5, width = 5)
plot( x= mg_df$muGrandSpMean, y = bps_df$betaPhotoSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bps_df$photo25), max(bps_df$photo75)), ylab = "Species level photoperiod slope", xlab = "Trait value", cex.lab = 1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df_east[,"muGrandSpMean"], # x mean
  bps_df_east[,"photo25"], # y 25
  mg_df_east[,"muGrandSpMean"],
  bps_df_east[,"photo75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df_east[,"trait25"], # x mean
  bps_df_east[,"betaPhotoSpMean"], # y 25
  mg_df_east[,"trait75"], # x mean
  bps_df_east[,"betaPhotoSpMean"],
  length = 0, col = "maroon", lwd = 2
)

arrows(
  mg_df_west[,"muGrandSpMean"], # x mean
  bps_df_west[,"photo25"], # y 25
  mg_df_west[,"muGrandSpMean"],
  bps_df_west[,"photo75"],
  length = 0, col= "darkslategray4", lwd = 2
)

arrows(
  mg_df_west[,"trait25"], # x mean
  bps_df_west[,"betaPhotoSpMean"], # y 25
  mg_df_west[,"trait75"], # x mean
  bps_df_west[,"betaPhotoSpMean"],
  length = 0, col = "darkslategray4", lwd = 2
)

mtext(side = 3, text = "CN, Photoperiod", adj = 0, cex = 1.25)
for(j in 1:length(muPhotoSp[,1])){
  abline(a = muPhotoSp[j,], b = betaTraitxPhotoMean, col=alpha("darkslategray4", 0.015))
}
abline(a=muPhotoSpMean, b=betaTraitxPhotoMean, col = "black")


my.label <- paste("c", ".", sep="")
put.fig.letter(label=my.label, location= "topleft", font=2)
###############################################
# DBH
Model <- readRDS(paste("output/", files[2], sep = ""))

ModelFit <- rstan::extract(Model)

muGrandSp <- data.frame(ModelFit$mu_grand_sp)
muGrandSpMean <- colMeans(muGrandSp)

betaForceSp <- data.frame(ModelFit$betaForceSp)
betaForceSpMean <- colMeans(betaForceSp)

quantile2575 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.25, 0.75))
  return(returnQuanilte)
}

bf_quan <- apply(betaForceSp, 2, quantile2575) 
mugrand_quan <- apply(muGrandSp, 2, quantile2575)

bfs <- rbind(betaForceSpMean, bf_quan)
bfs_t <- t(bfs)
bfs_df <- data.frame(bfs_t)
colnames(bfs_df)[colnames(bfs_df) == "X25."] <- "force25"
colnames(bfs_df)[colnames(bfs_df) == "X75."] <- "force75"
bfs_df$species <- specieslist

mg<- rbind(muGrandSpMean, mugrand_quan)
mg_t <- t(mg)
mg_df <- data.frame(mg_t)
colnames(mg_df)[colnames(mg_df) == "X25."] <- "trait25"
colnames(mg_df)[colnames(mg_df) == "X75."] <- "trait75"
mg_df$species <- specieslist

muForceSp <- data.frame(ModelFit$muForceSp)
muForceSpMean <- colMeans(muForceSp)

betaTraitxForce<- data.frame(ModelFit$betaTraitxForce)
betaTraitxForceMean <- colMeans(betaTraitxForce)

mg_df_east <- mg_df[mg_df$species %in% eastSp, ]
mg_df_west <- mg_df[mg_df$species %in% westSp, ]

bfs_df_east <- bfs_df[bfs_df$species %in% eastSp, ]
bfs_df_west <- bfs_df[bfs_df$species %in% westSp, ]

plot( x= mg_df$muGrandSpMean, y = bfs_df$betaForceSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bfs_df$force25), max(bfs_df$force75)), ylab = "Species level forcing slope", xlab = "Trait value", cex.lab = 1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df_east[,"muGrandSpMean"], # x mean
  bfs_df_east[,"force25"], # y 25
  mg_df_east[,"muGrandSpMean"],
  bfs_df_east[,"force75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df_east[,"trait25"], # x mean
  bfs_df_east[,"betaForceSpMean"], # y 25
  mg_df_east[,"trait75"], # x mean
  bfs_df_east[,"betaForceSpMean"],
  length = 0, col = "maroon", lwd = 2
)

arrows(
  mg_df_west[,"muGrandSpMean"], # x mean
  bfs_df_west[,"force25"], # y 25
  mg_df_west[,"muGrandSpMean"],
  bfs_df_west[,"force75"],
  length = 0, col= "darkslategray4", lwd = 2
)

arrows(
  mg_df_west[,"trait25"], # x mean
  bfs_df_west[,"betaForceSpMean"], # y 25
  mg_df_west[,"trait75"], # x mean
  bfs_df_west[,"betaForceSpMean"],
  length = 0, col = "darkslategray4", lwd = 2
)

mtext(side = 3, text = "DBH, Forcing", adj = 0, cex = 1.25)
for(j in 1:length(muForceSp[,1])){
  abline(a = muForceSp[j,], b = betaTraitxForceMean, col=alpha("darkslategray4", 0.015))
}
abline(a=muForceSpMean, b=betaTraitxForceMean, col = "black")

my.label <- paste("d", ".", sep="")
put.fig.letter(label=my.label, location= "topleft", font=2)

######################################################
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

betaTraitxChill<- data.frame(ModelFit$betaTraitxChill)
betaTraitxChillMean <- colMeans(betaTraitxChill)

mg_df_east <- mg_df[mg_df$species %in% eastSp, ]
mg_df_west <- mg_df[mg_df$species %in% westSp, ]

bcs_df_east <- bcs_df[bcs_df$species %in% eastSp, ]
bcs_df_west <- bcs_df[bcs_df$species %in% westSp, ]


plot( x= mg_df$muGrandSpMean, y = bcs_df$betaChillSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bcs_df$chill25), max(bcs_df$chill75)), ylab = "Species level chilling slope", xlab = "Trait value", cex.lab =1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df_east[,"muGrandSpMean"], # x mean
  bcs_df_east[,"chill25"], # y 25
  mg_df_east[,"muGrandSpMean"],
  bcs_df_east[,"chill75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df_east[,"trait25"], # x mean
  bcs_df_east[,"betaChillSpMean"], # y 25
  mg_df_east[,"trait75"], # x mean
  bcs_df_east[,"betaChillSpMean"],
  length = 0, col = "maroon", lwd = 2
)

arrows(
  mg_df_west[,"muGrandSpMean"], # x mean
  bcs_df_west[,"chill25"], # y 25
  mg_df_west[,"muGrandSpMean"],
  bcs_df_west[,"chill75"],
  length = 0, col= "darkslategray4", lwd = 2
)

arrows(
  mg_df_west[,"trait25"], # x mean
  bcs_df_west[,"betaChillSpMean"], # y 25
  mg_df_west[,"trait75"], # x mean
  bcs_df_west[,"betaChillSpMean"],
  length = 0, col = "darkslategray4", lwd = 2
)

mtext(side = 3, text = "DBH, Chilling", adj = 0, cex = 1.25)
for(j in 1:length(muChillSp[,1])){
  abline(a = muChillSp[j,], b = betaTraitxChillMean, col=alpha("darkslategray4", 0.015))
}
abline(a=muChillSpMean, b=betaTraitxChillMean, col = "black")

my.label <- paste("e", ".", sep="")
put.fig.letter(label=my.label, location= "topleft", font=2)
#######################################################################
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

betaTraitxPhoto<- data.frame(ModelFit$betaTraitxPhoto)
betaTraitxPhotoMean <- colMeans(betaTraitxPhoto)

betaTraitxPhoto<- data.frame(ModelFit$betaTraitxPhoto)
betaTraitxPhotoMean <- colMeans(betaTraitxPhoto)

bps_df_east <- bps_df[bps_df$species %in% eastSp, ]
bps_df_west <- bps_df[bps_df$species %in% westSp, ]

plot( x= mg_df$muGrandSpMean, y = bps_df$betaPhotoSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bps_df$photo25), max(bps_df$photo75)), ylab = "Species level photoperiod slope", xlab = "Trait value",cex.lab = 1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df_east[,"muGrandSpMean"], # x mean
  bps_df_east[,"photo25"], # y 25
  mg_df_east[,"muGrandSpMean"],
  bps_df_east[,"photo75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df_east[,"trait25"], # x mean
  bps_df_east[,"betaPhotoSpMean"], # y 25
  mg_df_east[,"trait75"], # x mean
  bps_df_east[,"betaPhotoSpMean"],
  length = 0, col = "maroon", lwd = 2
)

arrows(
  mg_df_west[,"muGrandSpMean"], # x mean
  bps_df_west[,"photo25"], # y 25
  mg_df_west[,"muGrandSpMean"],
  bps_df_west[,"photo75"],
  length = 0, col= "darkslategray4", lwd = 2
)

arrows(
  mg_df_west[,"trait25"], # x mean
  bps_df_west[,"betaPhotoSpMean"], # y 25
  mg_df_west[,"trait75"], # x mean
  bps_df_west[,"betaPhotoSpMean"],
  length = 0, col = "darkslategray4", lwd = 2
)

mtext(side = 3, text = "DBH, Photoperiod", adj = 0, cex = 1.25)
for(j in 1:length(muPhotoSp[,1])){
  abline(a = muPhotoSp[j,], b = betaTraitxPhotoMean, col=alpha("darkslategray4", 0.015))
}
abline(a=muPhotoSpMean, b=betaTraitxPhotoMean, col = "black")

my.label <- paste("f", ".", sep="")
put.fig.letter(label=my.label, location= "topleft", font=2)
#dev.off()
####################
# height
#pdf(paste("figures/cue", "trait_wtrend_supp", ".pdf", sep = ""), height = 8, width = 12)

Model <- readRDS(paste("output/", files[3], sep = ""))

ModelFit <- rstan::extract(Model)

muGrandSp <- data.frame(ModelFit$mu_grand_sp)
muGrandSpMean <- colMeans(muGrandSp)

betaForceSp <- data.frame(ModelFit$betaForceSp)
betaForceSpMean <- colMeans(betaForceSp)

quantile2575 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.25, 0.75))
  return(returnQuanilte)
}

bf_quan <- apply(betaForceSp, 2, quantile2575) 
mugrand_quan <- apply(muGrandSp, 2, quantile2575)

bfs <- rbind(betaForceSpMean, bf_quan)
bfs_t <- t(bfs)
bfs_df <- data.frame(bfs_t)
colnames(bfs_df)[colnames(bfs_df) == "X25."] <- "force25"
colnames(bfs_df)[colnames(bfs_df) == "X75."] <- "force75"
bfs_df$species <- specieslist

mg<- rbind(muGrandSpMean, mugrand_quan)
mg_t <- t(mg)
mg_df <- data.frame(mg_t)
colnames(mg_df)[colnames(mg_df) == "X25."] <- "trait25"
colnames(mg_df)[colnames(mg_df) == "X75."] <- "trait75"
mg_df$species <- specieslist


muForceSp <- data.frame(ModelFit$muForceSp)
muForceSpMean <- colMeans(muForceSp)

betaTraitxForce<- data.frame(ModelFit$betaTraitxForce)
betaTraitxForceMean <- colMeans(betaTraitxForce)

mg_df_east <- mg_df[mg_df$species %in% eastSp, ]
mg_df_west <- mg_df[mg_df$species %in% westSp, ]

bfs_df_east <- bfs_df[bfs_df$species %in% eastSp, ]
bfs_df_west <- bfs_df[bfs_df$species %in% westSp, ]

plot( x= mg_df$muGrandSpMean, y = bfs_df$betaForceSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bfs_df$force25), max(bfs_df$force75)), ylab = "Species level forcing slope", xlab = "Trait value",  cex.lab =1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df_east[,"muGrandSpMean"], # x mean
  bfs_df_east[,"force25"], # y 25
  mg_df_east[,"muGrandSpMean"],
  bfs_df_east[,"force75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df_east[,"trait25"], # x mean
  bfs_df_east[,"betaForceSpMean"], # y 25
  mg_df_east[,"trait75"], # x mean
  bfs_df_east[,"betaForceSpMean"],
  length = 0, col = "maroon", lwd = 2
)

arrows(
  mg_df_west[,"muGrandSpMean"], # x mean
  bfs_df_west[,"force25"], # y 25
  mg_df_west[,"muGrandSpMean"],
  bfs_df_west[,"force75"],
  length = 0, col= "darkslategray4", lwd = 2
)

arrows(
  mg_df_west[,"trait25"], # x mean
  bfs_df_west[,"betaForceSpMean"], # y 25
  mg_df_west[,"trait75"], # x mean
  bfs_df_west[,"betaForceSpMean"],
  length = 0, col = "darkslategray4", lwd = 2
)

mtext(side = 3, text = "Height, Forcing", adj = 0, cex = 1.25)
for(j in 1:length(muForceSp[,1])){
  abline(a = muForceSp[j,], b = betaTraitxForceMean, col=alpha("darkslategray4", 0.015))
}
abline(a=muForceSpMean, b=betaTraitxForceMean, col = "black")

my.label <- paste("a", ".", sep="")
put.fig.letter(label=my.label, location= "topleft", font=2)
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

betaTraitxChill<- data.frame(ModelFit$betaTraitxChill)
betaTraitxChillMean <- colMeans(betaTraitxChill)

mg_df_east <- mg_df[mg_df$species %in% eastSp, ]
mg_df_west <- mg_df[mg_df$species %in% westSp, ]

bcs_df_east <- bcs_df[bcs_df$species %in% eastSp, ]
bcs_df_west <- bcs_df[bcs_df$species %in% westSp, ]

plot( x= mg_df$muGrandSpMean, y = bcs_df$betaChillSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bcs_df$chill25), max(bcs_df$chill75)), ylab = "Species level chilling slope", xlab = "Trait value", cex.lab = 1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df_east[,"muGrandSpMean"], # x mean
  bcs_df_east[,"chill25"], # y 25
  mg_df_east[,"muGrandSpMean"],
  bcs_df_east[,"chill75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df_east[,"trait25"], # x mean
  bcs_df_east[,"betaChillSpMean"], # y 25
  mg_df_east[,"trait75"], # x mean
  bcs_df_east[,"betaChillSpMean"],
  length = 0, col = "maroon", lwd = 2
)

arrows(
  mg_df_west[,"muGrandSpMean"], # x mean
  bcs_df_west[,"chill25"], # y 25
  mg_df_west[,"muGrandSpMean"],
  bcs_df_west[,"chill75"],
  length = 0, col= "darkslategray4", lwd = 2
)

arrows(
  mg_df_west[,"trait25"], # x mean
  bcs_df_west[,"betaChillSpMean"], # y 25
  mg_df_west[,"trait75"], # x mean
  bcs_df_west[,"betaChillSpMean"],
  length = 0, col = "darkslategray4", lwd = 2
)
mtext(side = 3, text = "Height, Chilling", adj = 0, cex = 1.25)
for(j in 1:length(muChillSp[,1])){
  abline(a = muChillSp[j,], b = betaTraitxChillMean, col=alpha("darkslategray4", 0.015))
}
abline(a=muChillSpMean, b=betaTraitxChillMean, col = "black")


my.label <- paste("b", ".", sep="")
put.fig.letter(label=my.label, location= "topleft", font=2)
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

betaTraitxPhoto<- data.frame(ModelFit$betaTraitxPhoto)
betaTraitxPhotoMean <- colMeans(betaTraitxPhoto)

betaTraitxPhoto<- data.frame(ModelFit$betaTraitxPhoto)
betaTraitxPhotoMean <- colMeans(betaTraitxPhoto)

bps_df_east <- bps_df[bps_df$species %in% eastSp, ]
bps_df_west <- bps_df[bps_df$species %in% westSp, ]

plot( x= mg_df$muGrandSpMean, y = bps_df$betaPhotoSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bps_df$photo25), max(bps_df$photo75)), ylab = "Species level photoperiod slope", xlab = "Trait value", cex.lab =1.5 ) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df_east[,"muGrandSpMean"], # x mean
  bps_df_east[,"photo25"], # y 25
  mg_df_east[,"muGrandSpMean"],
  bps_df_east[,"photo75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df_east[,"trait25"], # x mean
  bps_df_east[,"betaPhotoSpMean"], # y 25
  mg_df_east[,"trait75"], # x mean
  bps_df_east[,"betaPhotoSpMean"],
  length = 0, col = "maroon", lwd = 2
)

arrows(
  mg_df_west[,"muGrandSpMean"], # x mean
  bps_df_west[,"photo25"], # y 25
  mg_df_west[,"muGrandSpMean"],
  bps_df_west[,"photo75"],
  length = 0, col= "darkslategray4", lwd = 2
)

arrows(
  mg_df_west[,"trait25"], # x mean
  bps_df_west[,"betaPhotoSpMean"], # y 25
  mg_df_west[,"trait75"], # x mean
  bps_df_west[,"betaPhotoSpMean"],
  length = 0, col = "darkslategray4", lwd = 2
)

mtext(side = 3, text = "Height, Photoperiod", adj = 0, cex = 1.25)  
for(j in 1:length(muPhotoSp[,1])){
  abline(a = muPhotoSp[j,], b = betaTraitxPhotoMean, col=alpha("darkslategray4", 0.015))
}
abline(a=muPhotoSpMean, b=betaTraitxPhotoMean, col = "black")

my.label <- paste("c", ".", sep="")
put.fig.letter(label=my.label, location= "topleft", font=2)
###############################################################
####################
# SSD
Model <- readRDS(paste("output/", files[5], sep = ""))

ModelFit <- rstan::extract(Model)

muGrandSp <- data.frame(ModelFit$mu_grand_sp)
muGrandSpMean <- colMeans(muGrandSp)

betaForceSp <- data.frame(ModelFit$betaForceSp)
betaForceSpMean <- colMeans(betaForceSp)

quantile2575 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.25, 0.75))
  return(returnQuanilte)
}

bf_quan <- apply(betaForceSp, 2, quantile2575) 
mugrand_quan <- apply(muGrandSp, 2, quantile2575)

bfs <- rbind(betaForceSpMean, bf_quan)
bfs_t <- t(bfs)
bfs_df <- data.frame(bfs_t)
colnames(bfs_df)[colnames(bfs_df) == "X25."] <- "force25"
colnames(bfs_df)[colnames(bfs_df) == "X75."] <- "force75"
bfs_df$species <- specieslist

mg<- rbind(muGrandSpMean, mugrand_quan)
mg_t <- t(mg)
mg_df <- data.frame(mg_t)
colnames(mg_df)[colnames(mg_df) == "X25."] <- "trait25"
colnames(mg_df)[colnames(mg_df) == "X75."] <- "trait75"
mg_df$species <- specieslist

muForceSp <- data.frame(ModelFit$muForceSp)
muForceSpMean <- colMeans(muForceSp)

betaTraitxForce <- data.frame(ModelFit$betaTraitxForce)
betaTraitxForceMean <- colMeans(betaTraitxForce)

mg_df_east <- mg_df[mg_df$species %in% eastSp, ]
mg_df_west <- mg_df[mg_df$species %in% westSp, ]

bfs_df_east <- bfs_df[bfs_df$species %in% eastSp, ]
bfs_df_west <- bfs_df[bfs_df$species %in% westSp, ]

plot( x= mg_df$muGrandSpMean, y = bfs_df$betaForceSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bfs_df$force25), max(bfs_df$force75)), ylab = "Species level forcing slope", xlab = "Trait value", cex.lab = 1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df_east[,"muGrandSpMean"], # x mean
  bfs_df_east[,"force25"], # y 25
  mg_df_east[,"muGrandSpMean"],
  bfs_df_east[,"force75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df_east[,"trait25"], # x mean
  bfs_df_east[,"betaForceSpMean"], # y 25
  mg_df_east[,"trait75"], # x mean
  bfs_df_east[,"betaForceSpMean"],
  length = 0, col = "maroon", lwd = 2
)

arrows(
  mg_df_west[,"muGrandSpMean"], # x mean
  bfs_df_west[,"force25"], # y 25
  mg_df_west[,"muGrandSpMean"],
  bfs_df_west[,"force75"],
  length = 0, col= "darkslategray4", lwd = 2
)

arrows(
  mg_df_west[,"trait25"], # x mean
  bfs_df_west[,"betaForceSpMean"], # y 25
  mg_df_west[,"trait75"], # x mean
  bfs_df_west[,"betaForceSpMean"],
  length = 0, col = "darkslategray4", lwd = 2
)

mtext(side = 3, text = "SSD, Forcing", adj = 0, cex = 1.25)
for(j in 1:length(muForceSp[,1])){
  abline(a = muForceSp[j,], b = betaTraitxForceMean, col=alpha("darkslategray4", 0.015))
}
abline(a=muForceSpMean, b=betaTraitxForceMean, col = "black")


my.label <- paste("d", ".", sep="")
put.fig.letter(label=my.label, location= "topleft", font=2)
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

betaTraitxChill<- data.frame(ModelFit$betaTraitxChill)
betaTraitxChillMean <- colMeans(betaTraitxChill)

mg_df_east <- mg_df[mg_df$species %in% eastSp, ]
mg_df_west <- mg_df[mg_df$species %in% westSp, ]

bcs_df_east <- bcs_df[bcs_df$species %in% eastSp, ]
bcs_df_west <- bcs_df[bcs_df$species %in% westSp, ]

plot( x= mg_df$muGrandSpMean, y = bcs_df$betaChillSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bcs_df$chill25), max(bcs_df$chill75)), ylab = "Species level chilling slope", xlab = "Trait value", cex.lab =1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df_east[,"muGrandSpMean"], # x mean
  bcs_df_east[,"chill25"], # y 25
  mg_df_east[,"muGrandSpMean"],
  bcs_df_east[,"chill75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df_east[,"trait25"], # x mean
  bcs_df_east[,"betaChillSpMean"], # y 25
  mg_df_east[,"trait75"], # x mean
  bcs_df_east[,"betaChillSpMean"],
  length = 0, col = "maroon", lwd = 2
)

arrows(
  mg_df_west[,"muGrandSpMean"], # x mean
  bcs_df_west[,"chill25"], # y 25
  mg_df_west[,"muGrandSpMean"],
  bcs_df_west[,"chill75"],
  length = 0, col= "darkslategray4", lwd = 2
)

arrows(
  mg_df_west[,"trait25"], # x mean
  bcs_df_west[,"betaChillSpMean"], # y 25
  mg_df_west[,"trait75"], # x mean
  bcs_df_west[,"betaChillSpMean"],
  length = 0, col = "darkslategray4", lwd = 2
)
mtext(side = 3, text = "SSD, Chilling", adj = 0, cex = 1.25)
for(j in 1:length(muChillSp[,1])){
  abline(a = muChillSp[j,], b = betaTraitxChillMean, col=alpha("darkslategray4", 0.015))
}
abline(a=muChillSpMean, b=betaTraitxChillMean, col = "black")

my.label <- paste("e", ".", sep="")
put.fig.letter(label=my.label, location= "topleft", font=2)
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

betaTraitxPhoto<- data.frame(ModelFit$betaTraitxPhoto)
betaTraitxPhotoMean <- colMeans(betaTraitxPhoto)

betaTraitxPhoto<- data.frame(ModelFit$betaTraitxPhoto)
betaTraitxPhotoMean <- colMeans(betaTraitxPhoto)

bps_df_east <- bps_df[bps_df$species %in% eastSp, ]
bps_df_west <- bps_df[bps_df$species %in% westSp, ]

plot( x= mg_df$muGrandSpMean, y = bps_df$betaPhotoSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bps_df$photo25), max(bps_df$photo75)), ylab = "Species level photoperiod slope", xlab = "Trait value", cex.lab = 1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df_east[,"muGrandSpMean"], # x mean
  bps_df_east[,"photo25"], # y 25
  mg_df_east[,"muGrandSpMean"],
  bps_df_east[,"photo75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df_east[,"trait25"], # x mean
  bps_df_east[,"betaPhotoSpMean"], # y 25
  mg_df_east[,"trait75"], # x mean
  bps_df_east[,"betaPhotoSpMean"],
  length = 0, col = "maroon", lwd = 2
)

arrows(
  mg_df_west[,"muGrandSpMean"], # x mean
  bps_df_west[,"photo25"], # y 25
  mg_df_west[,"muGrandSpMean"],
  bps_df_west[,"photo75"],
  length = 0, col= "darkslategray4", lwd = 2
)

arrows(
  mg_df_west[,"trait25"], # x mean
  bps_df_west[,"betaPhotoSpMean"], # y 25
  mg_df_west[,"trait75"], # x mean
  bps_df_west[,"betaPhotoSpMean"],
  length = 0, col = "darkslategray4", lwd = 2
)

mtext(side = 3, text = "SSD, Photoperiod", adj = 0, cex = 1.25)
for(j in 1:length(muPhotoSp[,1])){
  abline(a = muPhotoSp[j,], b = betaTraitxPhotoMean, col=alpha("darkslategray4", 0.015))
}
abline(a=muPhotoSpMean, b=betaTraitxPhotoMean, col = "black")


my.label <- paste("f", ".", sep="")
put.fig.letter(label=my.label, location= "topleft", font=2)
##############################################################################################################################
####################
# LMA
Model <- readRDS(paste("output/", files[4], sep = ""))

ModelFit <- rstan::extract(Model)

muGrandSp <- data.frame(ModelFit$mu_grand_sp)
muGrandSpMean <- colMeans(muGrandSp)

betaForceSp <- data.frame(ModelFit$betaForceSp)
betaForceSpMean <- colMeans(betaForceSp)

quantile2575 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.25, 0.75))
  return(returnQuanilte)
}

bf_quan <- apply(betaForceSp, 2, quantile2575) 
mugrand_quan <- apply(muGrandSp, 2, quantile2575)

bfs <- rbind(betaForceSpMean, bf_quan)
bfs_t <- t(bfs)
bfs_df <- data.frame(bfs_t)
colnames(bfs_df)[colnames(bfs_df) == "X25."] <- "force25"
colnames(bfs_df)[colnames(bfs_df) == "X75."] <- "force75"
bfs_df$species <- specieslist

mg<- rbind(muGrandSpMean, mugrand_quan)
mg_t <- t(mg)
mg_df <- data.frame(mg_t)
colnames(mg_df)[colnames(mg_df) == "X25."] <- "trait25"
colnames(mg_df)[colnames(mg_df) == "X75."] <- "trait75"
mg_df$species <- specieslist

muForceSp <- data.frame(ModelFit$muForceSp)
muForceSpMean <- colMeans(muForceSp)

betaTraitxForce <- data.frame(ModelFit$betaTraitxForce)
betaTraitxForceMean <- colMeans(betaTraitxForce)

mg_df_east <- mg_df[mg_df$species %in% eastSp, ]
mg_df_west <- mg_df[mg_df$species %in% westSp, ]

bfs_df_east <- bfs_df[bfs_df$species %in% eastSp, ]
bfs_df_west <- bfs_df[bfs_df$species %in% westSp, ]

plot( x= mg_df$muGrandSpMean, y = bfs_df$betaForceSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bfs_df$force25), max(bfs_df$force75)), ylab = "Species level forcing slope", xlab = "Trait value", cex.lab = 1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df_east[,"muGrandSpMean"], # x mean
  bfs_df_east[,"force25"], # y 25
  mg_df_east[,"muGrandSpMean"],
  bfs_df_east[,"force75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df_east[,"trait25"], # x mean
  bfs_df_east[,"betaForceSpMean"], # y 25
  mg_df_east[,"trait75"], # x mean
  bfs_df_east[,"betaForceSpMean"],
  length = 0, col = "maroon", lwd = 2
)

arrows(
  mg_df_west[,"muGrandSpMean"], # x mean
  bfs_df_west[,"force25"], # y 25
  mg_df_west[,"muGrandSpMean"],
  bfs_df_west[,"force75"],
  length = 0, col= "darkslategray4", lwd = 2
)

arrows(
  mg_df_west[,"trait25"], # x mean
  bfs_df_west[,"betaForceSpMean"], # y 25
  mg_df_west[,"trait75"], # x mean
  bfs_df_west[,"betaForceSpMean"],
  length = 0, col = "darkslategray4", lwd = 2
)

mtext(side = 3, text = "LMA, Forcing", adj = 0, cex = 1.25)
for(j in 1:length(muForceSp[,1])){
  abline(a = muForceSp[j,], b = betaTraitxForceMean, col=alpha("darkslategray4", 0.015))
}
abline(a=muForceSpMean, b=betaTraitxForceMean, col = "black")


my.label <- paste("d", ".", sep="")
put.fig.letter(label=my.label, location= "topleft", font=2)
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

betaTraitxChill<- data.frame(ModelFit$betaTraitxChill)
betaTraitxChillMean <- colMeans(betaTraitxChill)

mg_df_east <- mg_df[mg_df$species %in% eastSp, ]
mg_df_west <- mg_df[mg_df$species %in% westSp, ]

bcs_df_east <- bcs_df[bcs_df$species %in% eastSp, ]
bcs_df_west <- bcs_df[bcs_df$species %in% westSp, ]

plot( x= mg_df$muGrandSpMean, y = bcs_df$betaChillSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bcs_df$chill25), max(bcs_df$chill75)), ylab = "Species level chilling slope", xlab = "Trait value", cex.lab =1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df_east[,"muGrandSpMean"], # x mean
  bcs_df_east[,"chill25"], # y 25
  mg_df_east[,"muGrandSpMean"],
  bcs_df_east[,"chill75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df_east[,"trait25"], # x mean
  bcs_df_east[,"betaChillSpMean"], # y 25
  mg_df_east[,"trait75"], # x mean
  bcs_df_east[,"betaChillSpMean"],
  length = 0, col = "maroon", lwd = 2
)

arrows(
  mg_df_west[,"muGrandSpMean"], # x mean
  bcs_df_west[,"chill25"], # y 25
  mg_df_west[,"muGrandSpMean"],
  bcs_df_west[,"chill75"],
  length = 0, col= "darkslategray4", lwd = 2
)

arrows(
  mg_df_west[,"trait25"], # x mean
  bcs_df_west[,"betaChillSpMean"], # y 25
  mg_df_west[,"trait75"], # x mean
  bcs_df_west[,"betaChillSpMean"],
  length = 0, col = "darkslategray4", lwd = 2
)
mtext(side = 3, text = "LMA, Chilling", adj = 0, cex = 1.25)
for(j in 1:length(muChillSp[,1])){
  abline(a = muChillSp[j,], b = betaTraitxChillMean, col=alpha("darkslategray4", 0.015))
}
abline(a=muChillSpMean, b=betaTraitxChillMean, col = "black")

my.label <- paste("e", ".", sep="")
put.fig.letter(label=my.label, location= "topleft", font=2)
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

betaTraitxPhoto<- data.frame(ModelFit$betaTraitxPhoto)
betaTraitxPhotoMean <- colMeans(betaTraitxPhoto)

betaTraitxPhoto<- data.frame(ModelFit$betaTraitxPhoto)
betaTraitxPhotoMean <- colMeans(betaTraitxPhoto)

bps_df_east <- bps_df[bps_df$species %in% eastSp, ]
bps_df_west <- bps_df[bps_df$species %in% westSp, ]

plot( x= mg_df$muGrandSpMean, y = bps_df$betaPhotoSpMean, type="n", xlim = c(min(mg_df$trait25), max(mg_df$trait75)), ylim = c(min(bps_df$photo25), max(bps_df$photo75)), ylab = "Species level photoperiod slope", xlab = "Trait value", cex.lab = 1.5) # blank plot with x range 
# 3 columns, mean, quantile
# min and max defined by quantiles
arrows(
  mg_df_east[,"muGrandSpMean"], # x mean
  bps_df_east[,"photo25"], # y 25
  mg_df_east[,"muGrandSpMean"],
  bps_df_east[,"photo75"],
  length = 0, col= "maroon", lwd = 2
)

arrows(
  mg_df_east[,"trait25"], # x mean
  bps_df_east[,"betaPhotoSpMean"], # y 25
  mg_df_east[,"trait75"], # x mean
  bps_df_east[,"betaPhotoSpMean"],
  length = 0, col = "maroon", lwd = 2
)

arrows(
  mg_df_west[,"muGrandSpMean"], # x mean
  bps_df_west[,"photo25"], # y 25
  mg_df_west[,"muGrandSpMean"],
  bps_df_west[,"photo75"],
  length = 0, col= "darkslategray4", lwd = 2
)

arrows(
  mg_df_west[,"trait25"], # x mean
  bps_df_west[,"betaPhotoSpMean"], # y 25
  mg_df_west[,"trait75"], # x mean
  bps_df_west[,"betaPhotoSpMean"],
  length = 0, col = "darkslategray4", lwd = 2
)

mtext(side = 3, text = "LMA, Photoperiod", adj = 0, cex = 1.25)
for(j in 1:length(muPhotoSp[,1])){
  abline(a = muPhotoSp[j,], b = betaTraitxPhotoMean, col=alpha("darkslategray4", 0.015))
}
abline(a=muPhotoSpMean, b=betaTraitxPhotoMean, col = "black")


my.label <- paste("f", ".", sep="")
put.fig.letter(label=my.label, location= "topleft", font=2)
dev.off()
