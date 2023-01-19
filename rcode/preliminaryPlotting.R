# Started Feb 4, 2022

# The purpose of this code is to plot the trait values for each species against the the slope estimates from the phenology model:

# Get the model output from pheno_bc repo:
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(rstan)
library(shinystan)
#library(reshape2)
library(bayesplot)
library(ggplot2)
#library(RColorBrewer)
library(dplyr)
library(plyr)


if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  else{
  setwd("/home/deirdre/pheno_bc") # for midge
}

# Get the raw pheno data:
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

head(pheno)

# get the model output from my data
load("output/final/bb_4sites_phylo_mini.Rda")

sumt <- summary(mdl.4phyloMini)$summary
bforce <- sumt[grep("b_warm", rownames(sumt)), "mean"]; bforce <- bforce[2:48]
bchill <- sumt[grep("b_chill", rownames(sumt)), "mean"]; bchill <- bchill[2:48]
bphoto <- sumt[grep("b_photo", rownames(sumt)), "mean"]; bphoto <- bphoto[50:96]

spList <- read.csv("input/species_list.csv")
spList$spFact <- as.numeric(as.factor(spList$species.name))

spList <- spList[order(spList$spFact),]

spList$bforce <- bforce
spList$bchill <- bchill
spList$bphoto <- bphoto

##################################################
setwd("..//Treetraits")
trtDataSpp <- read.csv("data/allTrt.csv")
head(trtDataSpp)

spp <- spList$species

dlspp <- subset(spList, transect != "east")

dfspp <- subset(spList, transect != "west")

trtDataSpp$spFact <- as.numeric(as.factor(trtDataSpp$species))

trtDataSpp <- merge(trtDataSpp, spList, by = "species")
# Make a ton of exploratory plots!

# 1. Get species means for traits and for day bb and cues:
library(plyr)

trtMean <- ddply(trtDataSpp, c("species"), summarise,
      N = length(lma),
      meanLMA = mean(lma, na.rm = T),
      sdLMA   = sd(lma, na.rm = T),
      meanSSD = mean(ssd, na.rm = T),
      sdSSD   = sd(ssd, na.rm = T),
      meanCN = mean(C.N, na.rm = T),
      sdCN  = sd(C.N, na.rm = T),
      meanht = mean(ht, na.rm = T),
      sdht   = sd(ht, na.rm = T)
      #seLMA   = sd / sqrt(N)
)
trtMean$seLMA <- trtMean$sdLMA/sqrt(trtMean$N)
trtMean$seSSD <- trtMean$sdSSD/sqrt(trtMean$N)
trtMean$seCN <- trtMean$sdCN/sqrt(trtMean$N)
trtMean$seht <- trtMean$sdht/sqrt(trtMean$N)

# Get mean day bb per spp:
pheno$species <- tolower(pheno$species)
meanBB <- aggregate(pheno["bb"], pheno[c("species")], FUN = mean, na.rm = T)

#merge the cue, trait and bb data:

bbTrt <- merge(meanBB, trtMean, by = "species")
bbTrt <- merge(bbTrt, spList, by = "species")

## LOTS AND LOTS OF PLOTS!!!

# species by trait:
par(mfrow = c(1,4))

library(scales)
library(viridis)
q_colors = 47 # for no particular reason
v_colors =  viridis(q_colors)
bbTrt$color <- v_colors

temp <- bbTrt[order(bbTrt$bb),]
temp$color <- v_colors

pdf("figures/exploratory/bb_vs_trait.pdf", width = 18, height = 5)
par(mfrow = c(1,4))
plot(bb ~ meanLMA, data = temp, bg = color, pch = 21, cex = 3, cex.lab = 2, ylab = "Mean budburst day", xlab = "Mean LMA")
abline(lm(bb~meanLMA, data= temp), lwd = 2)

arrows( # x mean
  temp[,"meanLMA"] + temp[,"seLMA"] , # y 25
  temp[,"bb"],
  temp[,"meanLMA"] -temp[,"seLMA"],
  temp[,"bb"],
  length = 0, 
  col = temp$color
)


plot(bb ~ meanSSD, data =temp, bg = color, pch = 21, cex = 3, cex.lab = 2, ylab = "Mean budburst day", xlab = "Mean SSD")
abline(lm(bb~meanSSD, data= temp), lwd = 2)

arrows( # x mean
  temp[,"meanSSD"] + temp[,"seSSD"] , # y 25
  temp[,"bb"],
  temp[,"meanSSD"] -temp[,"seSSD"],
  temp[,"bb"],
  length = 0, 
  col = temp$color
)

plot(bb ~ meanCN, data =temp, bg = color, pch = 21, cex = 3, cex.lab = 2, ylab = "Mean budburst day", xlab = "Mean C:N")
abline(lm(bb~meanCN, data= temp), lwd = 2)

arrows( # x mean
  temp[,"meanCN"] + temp[,"seCN"] , # y 25
  temp[,"bb"],
  temp[,"meanCN"] -temp[,"seCN"],
  temp[,"bb"],
  length = 0, 
  col = temp$color
)

plot(bb ~ meanht, data =temp, bg = color, pch = 21, cex = 3, ylab = "Mean budburst day", xlab = "Mean Height")
abline(lm(bb~meanht, data= temp), lwd = 2)

arrows( # x mean
  temp[,"meanht"] + temp[,"seht"] , # y 25
  temp[,"bb"],
  temp[,"meanht"] -temp[,"seht"],
  temp[,"bb"],
  length = 0, 
  col = temp$color
)

dev.off()

# now separating the total from e vs w data:
pdf("figures/exploratory/bb_vs_trait_EW.pdf", width = 20, height = 20)
par(mfrow = c(3,4))
plot(bb ~ meanLMA, data = temp, bg = color, pch = 21, cex = 3, cex.lab = 2, ylab = "Mean budburst day", xlab = "Mean LMA", main = "Both transects", cex.main = 3, xlim = c(0,0.08), ylim = c(0,60))
abline(lm(bb~meanLMA, data= temp), lwd = 2)

arrows( # x mean
  temp[,"meanLMA"] + temp[,"seLMA"] , # y 25
  temp[,"bb"],
  temp[,"meanLMA"] -temp[,"seLMA"],
  temp[,"bb"],
  length = 0, 
  col = temp$color
)


plot(bb ~ meanSSD, data =temp, bg = color, pch = 21, cex = 3, cex.lab = 2, ylab = "Mean budburst day", xlab = "Mean SSD", xlim = c(0.1,0.8), ylim = c(0,60))
abline(lm(bb~meanSSD, data= temp), lwd = 2)

arrows( # x mean
  temp[,"meanSSD"] + temp[,"seSSD"] , # y 25
  temp[,"bb"],
  temp[,"meanSSD"] -temp[,"seSSD"],
  temp[,"bb"],
  length = 0, 
  col = temp$color
)

plot(bb ~ meanCN, data =temp, bg = color, pch = 21, cex = 3, cex.lab = 2, ylab = "Mean budburst day", xlab = "Mean C:N", xlim = c(5,40), ylim = c(0,60))
abline(lm(bb~meanCN, data= temp), lwd = 2)

arrows( # x mean
  temp[,"meanCN"] + temp[,"seCN"] , # y 25
  temp[,"bb"],
  temp[,"meanCN"] -temp[,"seCN"],
  temp[,"bb"],
  length = 0, 
  col = temp$color
)

plot(bb ~ meanht, data =temp, bg = color, pch = 21, cex = 3, cex.lab = 2, ylab = "Mean budburst day", xlab = "Mean Height", xlim = c(0,20), ylim = c(0,60))
abline(lm(bb~meanht, data= temp), lwd = 2)

arrows( # x mean
  temp[,"meanht"] + temp[,"seht"] , # y 25
  temp[,"bb"],
  temp[,"meanht"] -temp[,"seht"],
  temp[,"bb"],
  length = 0, 
  col = temp$color
)

## Eastern sites only:
tempE <- subset(temp, transect != "west")
plot(bb ~ meanLMA, data = tempE, bg = color, pch = 21, cex = 3, cex.lab = 2, ylab = "Mean budburst day", xlab = "Mean LMA", main = "Eastern Transect", cex.main = 3, xlim = c(0,0.08), ylim = c(0,60))
abline(lm(bb~meanLMA, data= tempE), lwd = 2)

arrows( # x mean
  tempE[,"meanLMA"] + tempE[,"seLMA"] , # y 25
  tempE[,"bb"],
  tempE[,"meanLMA"] -tempE[,"seLMA"],
  tempE[,"bb"],
  length = 0, 
  col = tempE$color
)


plot(bb ~ meanSSD, data =tempE, bg = color, pch = 21, cex = 3, cex.lab = 2, ylab = "Mean budburst day", xlab = "Mean SSD", xlim = c(0.1,0.8), ylim = c(0,60))
abline(lm(bb~meanSSD, data= tempE), lwd = 2)

arrows( # x mean
  tempE[,"meanSSD"] + tempE[,"seSSD"] , # y 25
  tempE[,"bb"],
  tempE[,"meanSSD"] -tempE[,"seSSD"],
  tempE[,"bb"],
  length = 0, 
  col = tempE$color
)

plot(bb ~ meanCN, data =tempE, bg = color, pch = 21, cex = 3, cex.lab = 2, ylab = "Mean budburst day", xlab = "Mean C:N", xlim = c(5,40), ylim = c(0,60))
abline(lm(bb~meanCN, data= tempE), lwd = 2)

arrows( # x mean
  tempE[,"meanCN"] + tempE[,"seCN"] , # y 25
  tempE[,"bb"],
  tempE[,"meanCN"] -tempE[,"seCN"],
  tempE[,"bb"],
  length = 0, 
  col = tempE$color
)

plot(bb ~ meanht, data =tempE, bg = color, pch = 21, cex = 3, cex.lab = 2, ylab = "Mean budburst day", xlab = "Mean Height", xlim = c(0,20), ylim = c(0,60))
abline(lm(bb~meanht, data= tempE), lwd = 2)

arrows( # x mean
  tempE[,"meanht"] + tempE[,"seht"] , # y 25
  tempE[,"bb"],
  tempE[,"meanht"] -tempE[,"seht"],
  tempE[,"bb"],
  length = 0, 
  col = tempE$color
)

## Western sites only:
tempW <- subset(temp, transect != "east")
plot(bb ~ meanLMA, data = tempW, bg = color, pch = 21, cex = 3, cex.lab = 2, ylab = "Mean budburst day", xlab = "Mean LMA", main = "Western Transect", cex.main = 3, xlim = c(0,0.08), ylim = c(0,60))
abline(lm(bb~meanLMA, data= tempW), lwd = 2)

arrows( # x mean
  tempW[,"meanLMA"] + tempW[,"seLMA"] , # y 25
  tempW[,"bb"],
  tempW[,"meanLMA"] -tempW[,"seLMA"],
  tempW[,"bb"],
  length = 0, 
  col = tempW$color
)


plot(bb ~ meanSSD, data =tempW, bg = color, pch = 21, cex = 3, cex.lab = 2, ylab = "Mean budburst day", xlab = "Mean SSD", xlim = c(0.1,0.8), ylim = c(0,60))
abline(lm(bb~meanSSD, data= tempW), lwd = 2)

arrows( # x mean
  tempW[,"meanSSD"] + tempW[,"seSSD"] , # y 25
  tempW[,"bb"],
  tempW[,"meanSSD"] -tempW[,"seSSD"],
  tempW[,"bb"],
  length = 0, 
  col = tempW$color
)

plot(bb ~ meanCN, data =tempW, bg = color, pch = 21, cex = 3, cex.lab = 2, ylab = "Mean budburst day", xlab = "Mean C:N", xlim = c(5,40), ylim = c(0,60))
abline(lm(bb~meanCN, data= tempW), lwd = 2)

arrows( # x mean
  tempW[,"meanCN"] + tempW[,"seCN"] , # y 25
  tempW[,"bb"],
  tempW[,"meanCN"] -tempW[,"seCN"],
  tempW[,"bb"],
  length = 0, 
  col = tempW$color
)

plot(bb ~ meanht, data =tempW, bg = color, pch = 21, cex = 3, cex.lab = 2, ylab = "Mean budburst day", xlab = "Mean Height", xlim = c(0,20), ylim = c(0,60))
abline(lm(bb~meanht, data= tempW), lwd = 2)

arrows( # x mean
  tempW[,"meanht"] + tempW[,"seht"] , # y 25
  tempW[,"bb"],
  tempW[,"meanht"] -tempW[,"seht"],
  tempW[,"bb"],
  length = 0, 
  col = tempW$color
)

dev.off()

# plot bb vs trait, but with colours by cue responses:
chill <- bbTrt[order(bbTrt$bchill),]
chill$color <- v_colors

pdf("figures/exploratory/bb_vs_trait_chillColor_EW.pdf", width = 18, height = 15)
par(mfrow = c(3,4))
plot(bb ~ meanLMA, data = chill, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean LMA", main = "Chilling", cex.main = 2)

plot(bb ~ meanSSD, data = chill, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean SSD")

plot(bb ~ meanCN, data =chill, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean C:N")

plot(bb ~ meanht, data = chill, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean Height")

# Eastern transect
chillE <- subset(chill, transect != "west")
plot(bb ~ meanLMA, data = chillE, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean LMA", main = "Chilling Eastern transects", cex.main = 2)

plot(bb ~ meanSSD, data = chillE, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean SSD")

plot(bb ~ meanCN, data =chillE, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean C:N")

plot(bb ~ meanht, data = chillE, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean Height")

# Western transect
chillW <- subset(chill, transect != "east")
plot(bb ~ meanLMA, data = chillW, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean LMA", main = "Western transects", cex.main = 2)

plot(bb ~ meanSSD, data = chillW, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean SSD")

plot(bb ~ meanCN, data =chillW, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean C:N")

plot(bb ~ meanht, data = chillW, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean Height")
dev.off()

#<><><><><><>
forcing <- bbTrt[order(bbTrt$bforce),]
forcing$color <- v_colors

pdf("figures/exploratory/bb_vs_trait_forceColor_EW.pdf", width = 18, height = 15)
par(mfrow = c(3,4))
plot(bb ~ meanLMA, data = forcing, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean LMA", main = "Forcing both transects", cex.main = 2)

plot(bb ~ meanSSD, data = forcing, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean SSD")

plot(bb ~ meanCN, data = forcing, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean C:N")

plot(bb ~ meanht, data = forcing, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean Height")

# Eastern
forcingE <- subset(forcing, transect != "west")
plot(bb ~ meanLMA, data = forcingE, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean LMA", main = "Eastern transect", cex.main = 2)

plot(bb ~ meanSSD, data = forcingE, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean SSD")

plot(bb ~ meanCN, data = forcingE, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean C:N")

plot(bb ~ meanht, data = forcingE, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean Height")

# western
forcingW <- subset(forcing, transect != "east")

plot(bb ~ meanLMA, data = forcingW, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean LMA", main = "Western transect", cex.main = 2)

plot(bb ~ meanSSD, data = forcingW, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean SSD")

plot(bb ~ meanCN, data = forcingW, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean C:N")

plot(bb ~ meanht, data = forcingW, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean Height")

dev.off()
#<><><><><><>
photo <- bbTrt[order(bbTrt$bphoto),]
photo$color <- v_colors

pdf("figures/exploratory/bb_vs_trait_photoColor_EW.pdf", width = 18, height = 15)
par(mfrow = c(3,4))
plot(bb ~ meanLMA, data = photo, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean LMA", main = "Photoperiod both transects", cex.main = 2)

plot(bb ~ meanSSD, data = photo, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean SSD")

plot(bb ~ meanCN, data = photo, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean C:N")

plot(bb ~ meanht, data = photo, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean Height")

# eastern transect
photoE <- subset(photo, transect != "west")
plot(bb ~ meanLMA, data = photoE, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean LMA", main = "Eastern transect", cex.main = 2)

plot(bb ~ meanSSD, data = photoE, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean SSD")

plot(bb ~ meanCN, data = photoE, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean C:N")

plot(bb ~ meanht, data = photoE, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean Height")

# western 
photoW <- subset(photo, transect != "east")
plot(bb ~ meanLMA, data = photoW, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean LMA", main = "Western transect", cex.main = 2)

plot(bb ~ meanSSD, data = photoW, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean SSD")

plot(bb ~ meanCN, data = photoW, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean C:N")

plot(bb ~ meanht, data = photoW, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean Height")
dev.off()

# plot with cue responses on y axis, colour by species 
pdf("..//figures/exploratory/cues_vs_trait_spColor.pdf", width = 18, height = 15)
par(mfrow = c(3,4))
plot(bchill ~ meanLMA, data = temp, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Chilling response", xlab = "Mean LMA")

plot(bchill ~ meanSSD, data =temp, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Chilling response", xlab = "Mean SSD")

plot(bchill ~ meanCN, data =temp, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Chilling response", xlab = "Mean C:N")

plot(bchill ~ meanht, data =temp, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Chilling response", xlab = "Mean Height")

## Forcing
plot(bforce ~ meanLMA, data = temp, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Forcing response", xlab = "Mean LMA")

plot(bforce ~ meanSSD, data =temp, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Forcing response", xlab = "Mean SSD")

plot(bforce ~ meanCN, data =temp, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Forcing response", xlab = "Mean C:N")

plot(bforce ~ meanht, data =temp, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Forcing response", xlab = "Mean Height")

## Photoperiod
plot(bphoto ~ meanLMA, data = temp, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Photoperiod response", xlab = "Mean LMA")

plot(bphoto ~ meanSSD, data =temp, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Photoperiod response", xlab = "Mean SSD")

plot(bphoto ~ meanCN, data =temp, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Photoperiod response", xlab = "Mean C:N")

plot(bphoto ~ meanht, data =temp, bg = color, pch = 21, cex = 3, cex.lab = 3, ylab = "Photoperiod response", xlab = "Mean Height")
dev.off()


################################################
# repeat the above plots but colour by east vs west:
temp <- bbTrt[order(bbTrt$bb),]
temp$ewcolor <- temp$transect
temp$ewcolor[temp$ewcolor == "east"] <- "cyan4"
temp$ewcolor[temp$ewcolor == "west"] <- "maroon"
temp$ewcolor[temp$ewcolor == "east/west"] <- "#FDE725FF"


pdf("..//figures/exploratory/bb_vs_trait_eastwest.pdf", width = 18, height = 5)
par(mfrow = c(1,4))
plot(bb ~ meanLMA, data =temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean LMA")

plot(bb ~ meanSSD, data =temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean SSD")

plot(bb ~ meanCN, data =temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean C:N")

plot(bb ~ meanht, data =temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Mean budburst day", xlab = "Mean Height")
dev.off()

# plot bb vs trait, but with colours by cue responses:
chill <- bbTrt[order(bbTrt$bchill),]
chill$color <- v_colors

# plot with cue responses on y axis, colour by species 
pdf("..//figures/exploratory/cues_vs_trait_eastwestColor.pdf", width = 18, height = 15)
par(mfrow = c(3,4))
plot(bchill ~ meanLMA, data = temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Chilling response", xlab = "Mean LMA")

plot(bchill ~ meanSSD, data =temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Chilling response", xlab = "Mean SSD")

plot(bchill ~ meanCN, data =temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Chilling response", xlab = "Mean C:N")

plot(bchill ~ meanht, data =temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Chilling response", xlab = "Mean Height")

## Forcing
plot(bforce ~ meanLMA, data = temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Forcing response", xlab = "Mean LMA")

plot(bforce ~ meanSSD, data =temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Forcing response", xlab = "Mean SSD")

plot(bforce ~ meanCN, data =temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Forcing response", xlab = "Mean C:N")

plot(bforce ~ meanht, data =temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Forcing response", xlab = "Mean Height")

## Photoperiod
plot(bphoto ~ meanLMA, data = temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Photoperiod response", xlab = "Mean LMA")

plot(bphoto ~ meanSSD, data =temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Photoperiod response", xlab = "Mean SSD")

plot(bphoto ~ meanCN, data =temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Photoperiod response", xlab = "Mean C:N")

plot(bphoto ~ meanht, data =temp, bg = ewcolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Photoperiod response", xlab = "Mean Height")
dev.off()


################################################
# repeat the above plots but colour by tree vs shrub
temp <- bbTrt[order(bbTrt$bb),]
temp$tscolor <- temp$type
temp$tscolor[temp$tscolor == "shrub"] <- "cyan4"
temp$tscolor[temp$tscolor == "tree"] <- "maroon"


pdf("..//figures/exploratory/bb_vs_trait_treeShrub.pdf", width = 18, height = 5)
par(mfrow = c(1,4))
plot(bb ~ meanLMA, data =temp, bg = tscolor, pch = 21, cex = 3, ylab = "Mean budburst day", xlab = "Mean LMA")

plot(bb ~ meanSSD, data =temp, bg = tscolor, pch = 21, cex = 3, ylab = "Mean budburst day", xlab = "Mean SSD")

plot(bb ~ meanCN, data =temp, bg = tscolor, pch = 21, cex = 3, ylab = "Mean budburst day", xlab = "Mean C:N")

plot(bb ~ meanht, data =temp, bg = tscolor, pch = 21, cex = 3, ylab = "Mean budburst day", xlab = "Mean Height")
dev.off()

# plot bb vs trait, but with colours by cue responses:
chill <- bbTrt[order(bbTrt$bchill),]
chill$color <- v_colors

# plot with cue responses on y axis, colour by species 
pdf("..//figures/exploratory/cues_vs_trait_treeShrubColor.pdf", width = 18, height = 15)
par(mfrow = c(3,4))
plot(bchill ~ meanLMA, data = temp, bg = tscolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Chilling response", xlab = "Mean LMA")

plot(bchill ~ meanSSD, data =temp, bg = tscolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Chilling response", xlab = "Mean SSD")

plot(bchill ~ meanCN, data =temp, bg = tscolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Chilling response", xlab = "Mean C:N")

plot(bchill ~ meanht, data =temp, bg = tscolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Chilling response", xlab = "Mean Height")

## Forcing
plot(bforce ~ meanLMA, data = temp, bg = tscolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Forcing response", xlab = "Mean LMA")

plot(bforce ~ meanSSD, data =temp, bg = tscolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Forcing response", xlab = "Mean SSD")

plot(bforce ~ meanCN, data =temp, bg = tscolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Forcing response", xlab = "Mean C:N")

plot(bforce ~ meanht, data =temp, bg = tscolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Forcing response", xlab = "Mean Height")

## Photoperiod
plot(bphoto ~ meanLMA, data = temp, bg = tscolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Photoperiod response", xlab = "Mean LMA")

plot(bphoto ~ meanSSD, data =temp, bg = tscolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Photoperiod response", xlab = "Mean SSD")

plot(bphoto ~ meanCN, data =temp, bg = tscolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Photoperiod response", xlab = "Mean C:N")

plot(bphoto ~ meanht, data =temp, bg = tscolor, pch = 21, cex = 3, cex.lab = 3, ylab = "Photoperiod response", xlab = "Mean Height")
dev.off()

### Plot just the shrubs and trees
shrub <- subset(bbTrt, type == "shrub")
tree <- subset(bbTrt, type == "tree")

shrub$transNum <- as.numeric(as.factor(shrub$transect))

pdf("figures/exploratory/bb_shrubs_transects.pdf", width = 18, height = 5)
par(mfrow = c(1, 4))
plot(bb ~ meanLMA, data = shrub, col = "maroon", cex =4, cex.lab =2 ,pch = c(16, 17, 21)[shrub$transNum], xlab = "Mean LMA", ylab = "Mean budburst day")
plot(bb ~ meanSSD, data = shrub, col = "maroon", cex =4, cex.lab =2, pch = c(16, 17, 21)[shrub$transNum], xlab = "Mean SSD", ylab = "Mean budburst day")
plot(bb ~ meanCN, data = shrub, col = "maroon", cex =4, cex.lab =2, pch = c(16, 17, 21)[shrub$transNum], xlab = "Mean C:N", ylab = "Mean budburst day")
plot(bb ~ meanht, data = shrub, col = "maroon", cex =4, cex.lab =2, pch = c(16, 17, 21)[shrub$transNum], xlab = "Mean height", ylab = "Mean budburst day")

legend("topright", legend = c(expression("Western"),
                           expression("Eastern"),
                             expression("Both")),
        col = c("maroon", "maroon", "maroon"),
        #pt.bg = c("maroon", "maroon", "maroon"),
        inset = 0.02, pch = c(21,16,17 ), cex = 2.5, bty = "n")
dev.off()

tree$transNum <- as.numeric(as.factor(tree$transect))


pdf("figures/exploratory/bb_trees_transects.pdf", width = 18, height = 5)
par(mfrow = c(1, 4))
plot(bb ~ meanLMA, data = tree, col = "maroon", cex =4, cex.lab =2, pch = c(16, 17, 21)[tree$transNum], xlab = "Mean LMA", ylab = "Mean budburst day")
plot(bb ~ meanSSD, data = tree, col = "maroon", cex =4, cex.lab =2, pch = c(16, 17, 21)[tree$transNum], xlab = "Mean SSD", ylab = "Mean budburst day")
plot(bb ~ meanCN, data = tree, col = "maroon", cex =4, cex.lab =2, pch = c(16, 17, 21)[tree$transNum], xlab = "Mean C:N", ylab = "Mean budburst day")
plot(bb ~ meanht, data = tree, col = "maroon", cex =4, cex.lab =2, pch = c(16, 17, 21)[tree$transNum], xlab = "Mean height", ylab = "Mean budburst day")

legend("topright", legend = c(expression("Western"),
                              expression("Eastern"),
                              expression("Both")),
       col = c("maroon", "maroon", "maroon"),
       #pt.bg = c("maroon", "maroon", "maroon"),
       inset = 0.02, pch = c(21,16,17 ), cex = 2.5, bty = "n")
dev.off()
# Plots of the raw trait data look 


# now run simple trait model to get trait effects for each species:

lma_data  <- trtDataSpp[complete.cases(trtDataSpp$lma),]
lma_data$species.fact <- as.numeric(as.factor(lma_data$species))

lma_datalist <- list(yTraiti = lma_data$lma, 
                   N = nrow(lma_data), 
                   n_spec = length(unique(lma_data$species)), 
                   species = lma_data$species.fact, 
                   prior_mu_grand_mu = 0.5,
                   prior_mu_grand_sigma = 1,
                   prior_sigma_sp_mu = 1,
                   prior_sigma_sp_sigma = 1,
                   prior_sigma_traity_mu = 1,
                   prior_sigma_traity_sigma = 1
) 


mdl.lma <- stan('stan/bc_trait_only_2.stan',
                  data = lma_datalist,
                  iter = 6000,
                  warmup = 4000,
                  chains = 4,
                  include = FALSE,
                  pars = c("mu_y","y_hat"))
save(mdl.lma, file = "output_lma_traitonly.Rda")

sum.lma <- summary(mdl.lma)$summary

lmaTrt_df <- sum.lma[grep("muSp", rownames(sum.lma)), "mean"]; lmaTrt_df
########################################################
ht_data  <- trtDataSpp[complete.cases(trtDataSpp$ht),]
ht_datalist <- list(yTraiti = ht_data$ht, 
                     N = nrow(ht_data), 
                     n_spec = length(unique(ht_data$species)), 
                     species = as.numeric(as.factor(ht_data$species)), 
                     prior_mu_grand_mu = 15,
                     prior_mu_grand_sigma = 10,
                     prior_sigma_sp_mu = 5,
                     prior_sigma_sp_sigma = 1,
                     prior_sigma_traity_mu = 5,
                     prior_sigma_traity_sigma = 1
) 

mdl.ht <- stan('stan/bc_trait_only_2.stan',
                data = ht_datalist,
                iter = 6000,
                warmup = 4000,
                chains = 4,
                include = FALSE,
                pars = "mu_y")
save(mdl.ht, file = "output_ht_traitonly_df.Rda")

sum.ht <- summary(mdl.ht)$summary

htTrt_df <- sum.ht[grep("muSp", rownames(sum.ht)), "mean"]; htTrt_df
######################################################
dbh_data  <- trtDataSpp[complete.cases(trtDataSpp$dbh),]
dbh_data$species.fact <- as.numeric(as.factor(dbh_data$species))

dbh_datalist <- list(yTraiti = dbh_data$dbh, 
                     N = nrow(dbh_data), 
                     n_spec = length(unique(dbh_data$species)), 
                     species = dbh_data$species.fact, 
                     prior_mu_grand_mu = 0.5,
                     prior_mu_grand_sigma = 1,
                     prior_sigma_sp_mu = 1,
                     prior_sigma_sp_sigma = 1,
                     prior_sigma_traity_mu = 1,
                     prior_sigma_traity_sigma = 1
) 


mdl.dbh <- stan('stan/bc_trait_only_2.stan',
                data = dbh_datalist,
                iter = 6000,
                warmup = 4000,
                chains = 4,
                include = FALSE,
                pars = "mu_y")
save(mdl.dbh, file = "output_dbh_traitonly_df.Rda")

sum.dbh <- summary(mdl.dbh)$summary

dbhTrt_df <- sum.dbh[grep("muSp", rownames(sum.dbh)), "mean"]; dbhTrt_df
########################################################
ssd_data  <- trtDataSpp[complete.cases(trtDataSpp$ssd),]
ssd_data$species.fact <- as.numeric(as.factor(ssd_data$species))

ssd_datalist <- list(yTraiti = ssd_data$ssd, 
                     N = nrow(ssd_data), 
                     n_spec = length(unique(ssd_data$species)), 
                     species = ssd_data$species.fact, 
                     prior_mu_grand_mu = 0.5,
                     prior_mu_grand_sigma = 1,
                     prior_sigma_sp_mu = 1,
                     prior_sigma_sp_sigma = 1,
                     prior_sigma_traity_mu = 1,
                     prior_sigma_traity_sigma = 1
) 


mdl.ssd <- stan('stan/bc_trait_only_2.stan',
                data = ssd_datalist,
                iter = 6000,
                warmup = 4000,
                chains = 4,
                include = FALSE,
                pars = "mu_y")
save(mdl.ssd, file = "output_ssd_traitonly_df.Rda")

sum.ssd <- summary(mdl.ssd)$summary

ssdTrt_df <- sum.ssd[grep("muSp", rownames(sum.ssd)), "mean"]; ssdTrt_df
########################################################
perCN_data  <- trtDataSpp[complete.cases(trtDataSpp$C.N),]
perCN_data$species.fact <- as.numeric(as.factor(perCN_data$species))

perCN_datalist <- list(yTraiti = perCN_data$C.N, 
                     N = nrow(perCN_data), 
                     n_spec = length(unique(perCN_data$species)), 
                     species = perCN_data$species.fact, 
                     prior_mu_grand_mu = 0.5,
                     prior_mu_grand_sigma = 1,
                     prior_sigma_sp_mu = 1,
                     prior_sigma_sp_sigma = 1,
                     prior_sigma_traity_mu = 1,
                     prior_sigma_traity_sigma = 1
) 


mdl.perCN <- stan('stan/bc_trait_only_2.stan',
                data = perCN_datalist,
                iter = 6000,
                warmup = 4000,
                chains = 4,
                include = FALSE,
                pars = "mu_y")
save(mdl.perCN, file = "output_perCN_traitonly_df.Rda")

sum.perCN <- summary(mdl.perCN)$summary
perCNTrt_df <- sum.perCN[grep("muSp", rownames(sum.perCN)), "mean"]; perCNTrt_df
########################################################

# get the phenology data:
load("output/dl_df_allbb_4sites.Rda")

sum.pheno <- summary(mdl.t)$summary

spp <- sort(unique(trtData$species))

sppPheno <- tolower(sppPheno)
bforce <- data.frame(sum.pheno[grep("^b_warm", rownames(sum.pheno)), "mean"])
bforce$spp <- sppPheno
colnames(bforce) <- c("muSp","species")

bchill <- data.frame(sum.pheno[grep("^b_chill", rownames(sum.pheno)), "mean"])
bchill$spp <- sppPheno
colnames(bchill) <- c("muSp","species")

bphoto <- data.frame(sum.pheno[grep("^b_photo", rownames(sum.pheno)), "mean"])
bphoto$spp <- sppPheno
colnames(bphoto) <- c("muSp","species")
bphoto<- bphoto[50:98,]

perCNTrt_df <- data.frame(sum.perCN[grep("muSp", rownames(sum.perCN)), "mean"])
perCNTrt_df$spp <- dfspp
colnames(perCNTrt_df) <- c("muSp","species")

perCNTrt_dl <- data.frame(sum.perCN.dl[grep("muSp", rownames(sum.perCN.dl)), "mean"])
perCNTrt_dl$spp <- dlspp
colnames(perCNTrt_dl) <- c("muSp","species")

trtOut <- rbind(perCNTrt_df, perCNTrt_dl)
trtOut <- aggregate(trtOut["muSp"], trtOut[c("species")], FUN = mean)

head(trtOut)
head(bphoto)

bforce <- subset(bforce, species != "menfer")
bforce <- subset(bforce, species != "rhoalb")

plot(bforce$muSp ~ trtOut$muSp, col = "maroon", pch =19)
dfPt <- bforce[bforce$species %in% dfspp, ]
dftrt <- trtOut[trtOut$species %in% dfspp, ]
points(dfPt$muSp ~ dftrt$muSp, col = "darkslategray4", pch =19)

bchill <- subset(bchill, species != "menfer")
bchill <- subset(bchill, species != "rhoalb")

plot(bchill$muSp ~ trtOut$muSp, col = "maroon", pch =19)
dfPt <- bchill[bchill$species %in% dfspp, ]
dftrt <- trtOut[trtOut$species %in% dfspp, ]
points(dfPt$muSp ~ dftrt$muSp, col = "darkslategray4", pch =19)

bphoto <- subset(bphoto, species != "menfer")
bphoto <- subset(bphoto, species != "rhoalb")

plot(bphoto$muSp ~ trtOut$muSp, col = "maroon", pch =19)
dfPt <- bphoto[bphoto$species %in% dfspp, ]
dftrt <- trtOut[trtOut$species %in% dfspp, ]
points(dfPt$muSp ~ dftrt$muSp, col = "darkslategray4", pch =19)

########################################################
# get dl mdl output
load("output/output_lma_traitonly.Rda")
sum.lma <- summary(mdl.lma)$summary
lmaTrt <- sum.lma[grep("muSp", rownames(sum.lma)), "mean"]; 

load("output/output_ht_traitonly.Rda")
sum.ht <- summary(mdl.ht)$summary
htTrt <- sum.ht[grep("muSp", rownames(sum.ht)), "mean"];

load("output/output_ssd_traitonly.Rda")
sum.ssd <- summary(mdl.ssd)$summary
ssdTrt <- sum.ssd[grep("muSp", rownames(sum.ssd)), "mean"];

load("output/output_dbh_traitonly.Rda")
sum.dbh <- summary(mdl.dbh)$summary
dbhTrt <- sum.dbh[grep("muSp", rownames(sum.dbh)), "mean"]

load("output/output_perCN_traitonly_dl.Rda")
sum.perCN.dl<- summary(mdl.perCN)$summary
perCNTrt_dl <- sum.perCN.dl[grep("muSp", rownames(sum.perCN.dl)), "mean"]; 

perCNTrt_df

pdf("figures/muSpvscue.pdf", width = 15, height = 25)
par(mfrow = c(6, 3))
plot(bforce_dl ~ ssdTrt, col = "maroon", pch =19, ylim = c(-15, 5), xlim = c(-0.3,0.3))
points(bforce_df ~ssdTrt_df, col = "darkslategray4", pch =19)

plot(bchill_dl ~ ssdTrt, col = "maroon", pch =19, ylim = c(-10, 30), xlim = c(-0.1, 0.2))
points(bchill_df ~ ssdTrt_df, col = "darkslategray4", pch =19)

plot(bphoto_dl ~ ssdTrt, col = "maroon", pch =19, ylim = c(-15, 0), xlim = c(-0.25, 0.2))
points(bphoto_df ~ ssdTrt_df, col = "darkslategray4", pch =19)

plot(bforce_dl ~ perCNTrt, col = "maroon", pch =19, ylim = c(-15, 4), xlim = c(37,47))
points(bforce_df ~ perCTrt_df, col = "darkslategray4", pch =19)

plot(bchill_dl ~ perCTrt, col = "maroon", pch =19, ylim = c(-15, 30), xlim = c(38, 53))
points(bchill_df ~ perCTrt_df, col = "darkslategray4", pch =19)

plot(bphoto_dl ~ perCTrt, col = "maroon", pch =19, ylim = c(-15, 0), xlim = c(38, 53))
points(bphoto_df ~ perCTrt_df, col = "darkslategray4", pch =19)

plot(bforce_dl ~ perNTrt, col = "maroon", pch =19, ylim = c(-15, 0), xlim = c(4,-1))
points(bforce_df ~ perNTrt_df, col = "darkslategray4", pch =19)

plot(bchill_dl ~ perNTrt, col = "maroon", pch =19, ylim = c(-15, 30), xlim = c(4,-1))
points(bchill_df ~ perNTrt_df, col = "darkslategray4", pch =19)

plot(bphoto_dl ~ perNTrt, col = "maroon", pch =19, ylim = c(-15, 1), xlim = c(4,-1))
points(bphoto_df ~ perNTrt_df, col = "darkslategray4", pch =19)

plot(bforce_dl ~ htTrt, col = "maroon", pch =19, ylim = c(-15, 5), xlim = c(-5,15))
points(bforce_df ~ htTrt_df, col = "darkslategray4", pch =19)

plot(bchill_dl ~ htTrt, col = "maroon", pch =19, ylim = c(-15, 35), xlim = c(-5,15))
points(bchill_df ~ htTrt_df, col = "darkslategray4", pch =19)

plot(bphoto_dl ~ htTrt, col = "maroon", pch =19, ylim = c(-15, 1), xlim = c(-5,15))
points(bphoto_df ~ htTrt_df, col = "darkslategray4", pch =19)

plot(bforce_dl ~ lmaTrt, col = "maroon", pch =19, ylim = c(-15, 5), xlim = c(-0.05,0.05))
points(bforce_df ~ lmaTrt_df, col = "darkslategray4", pch =19)

plot(bchill_dl ~ lmaTrt, col = "maroon", pch =19, ylim = c(-15, 30), xlim = c(-0.05,0.05))
points(bchill_df ~ lmaTrt_df, col = "darkslategray4", pch =19)

plot(bphoto_dl ~ lmaTrt, col = "maroon", pch =19)
plot(bphoto_df ~ lmaTrt_df, col = "darkslategray4", pch =19)


plot(bforce_dl ~ dbhTrt, col = "maroon", pch =19)
plot(bforce_df ~ dbhTrt_df, col = "darkslategray4", pch =19)

plot(bchill_dl ~ dbhTrt, col = "maroon", pch =19)
plot(bchill_df ~ dbhTrt_df, col = "darkslategray4", pch =19)

plot(bphoto_dl ~ dbhTrt, col = "maroon", pch =19)
plot(bphoto_df ~ dbhTrt_df, col = "darkslategray4", pch =19)

dev.off()

