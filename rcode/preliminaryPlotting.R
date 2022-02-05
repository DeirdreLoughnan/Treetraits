# Started Feb 4, 2022

# The purpose of this code is to plot the trait values for each species against the the slope estimates from the phenology model:

# Get the model output from pheno_bc repo:

library(rstan)
library(shinystan)
#library(reshape2)
library(bayesplot)
library(ggplot2)
#library(RColorBrewer)
library(dplyr)
library(plyr)

rm(list=ls()) 
options(stringsAsFactors = FALSE)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  else{
  setwd("/home/deirdre/pheno_bc") # for midge
}


load("output/realDummSite.Rda")

setwd("..//Treetraits")
trtData <- read.csv("data/allTrt.csv")
head(trtData)


spp <- c("acegla", "acepen", "acerub", "acesac", "alninc","alnvir", "amealn", "aromel", "betall", "betlen", "betpap",
         "corcor", "corsto", "faggra","franig", "hamvir", "ilemuc", "kalang", "loncan", "loninv", "lyolig", "menfer", "nyssyl", "popbal","popgra", "poptre", "prupen", "quealb", "querub", "quevel", "rhafra", "rhoalb", "rhopri", "riblac", "rubpar", "samrac", "shecan","sorsco", "spialb", "spibet", "spipyr", "symalb", "vacmem", "vacmyr", "vibcas", "vibedu", "viblan")

trtDataSpp <- trtData[trtData$species %in% spp,]
trtDataSpp$species.fact <- as.numeric(as.factor(trtDataSpp$species))

head(trtData)


trt.dataset <- trtDataSpp %>%
  group_by(species,species.fact, site) %>%
  summarise_at(vars("ssd","per.N","per.C","ht","dbh","lma"), mean, na.rm = T)


spp.fact <- sort(unique(trt.dataset$species.fact))

sumt <- summary(mdl)$summary
bforce <- sumt[grep("b_force", rownames(sumt)), "mean"]; bforce
# bchill <- sumt[grep("b_chill", rownames(sumt)), "mean"]; bchill
# bphoto <- sumt[grep("b_photo", rownames(sumt)), "mean"]; bphoto

for (i in c(1:length(spp.fact))){
  trt.dataset$force[i] <- sumt[grep("b_force", rownames(sumt)),1][spp.fact[i]]
  trt.dataset$chill[i] <- sumt[grep("b_chill\\[", rownames(sumt)),1][spp.fact[i]]
  trt.dataset$photo[i] <- sumt[grep("b_photo\\[", rownames(sumt)),1][spp.fact[i]]

}

plot(force ~ ssd , data = trt.dataset)
plot(chill ~ ssd , data = trt.dataset)
plot(photo ~ ssd , data = trt.dataset)


plot(force ~ per.C , data = trt.dataset)
plot(chill ~ per.C , data = trt.dataset)
plot(photo ~ per.C , data = trt.dataset)

plot(force ~ per.N , data = trt.dataset)
plot(chill ~ per.N , data = trt.dataset)
plot(photo ~ per.N , data = trt.dataset)

plot(force ~ ht , data = trt.dataset)
plot(chill ~ ht , data = trt.dataset)
plot(photo ~ ht , data = trt.dataset)

plot(force ~ lma , data = trt.dataset)
plot(chill ~ lma , data = trt.dataset)
plot(photo ~ lma , data = trt.dataset)

plot(force ~ dbh , data = trt.dataset)
plot(chill ~ dbh , data = trt.dataset)
plot(photo ~ dbh , data = trt.dataset)


