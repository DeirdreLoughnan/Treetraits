# Started January 6, 2022 by Deirdre 

#cleaning trait data for western data:
# Traits:
# LMA - mass/area
#stem specific density - mass/ volume
# LCC/LNN
# height 
# DBH - of primary, secondary, tertiary stem

# I think I will need two datasets: 
# 1. Individual level metrics - height, DBH, SSD, LCC/LNN, 
# 2. Leaf level - SLA - measured five leaves
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(stringr)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/Treetraits/data") 
}  else{
  setwd("/home/deirdre/treetrait") # for midge
}

 species <- c( "acegla","betpap", "poptre", "popbal", 
           "alninc","alnvir","amealn", "corsto",
           "loninv", "menfer","rhoalb", "riblac",
           "rubpar","samrac","shecan","sorsco",
           "spibet","spipyr","symalb","vacmem","vibedu")
 
 type <- c("tree", "tree", "tree","tree", 
           "shrub", "shrub", "shrub","shrub", 
           "shrub", "shrub", "shrub","shrub", 
           "shrub", "shrub", "shrub","shrub", 
           "shrub", "shrub", "shrub","shrub", "shrub")
 
 sp.typ <- data.frame(species, type)
 
 #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
 
 # Start with height:
ht <- read.csv("western/uncleaned/traits_ht_calc_dbh_final.csv")
head(ht) 
sort(unique(ht$species)) # has extra species
ht$site[ht$site == "afrf"] <- "af"

ht.spp <- ht[ht$species %in% species, ]; sort(unique(ht.spp$species))
ht.merge <- merge(ht.spp, sp.typ, by = "species")

ht.tree <- subset(ht.merge, type == "tree")
# round to the nearest cm and convert to m
ht.tree$ht <- (round(ht.tree$ht, 0)/100)
head(ht.tree)

ht.shrub <- subset(ht.merge, type == "shrub")
ht.shrub$ht <- (round(ht.shrub$ht, 0)/100)
head(ht.shrub)

ht.treeshrub <- rbind(ht.tree, ht.shrub)
ht.treeshrub$sample <- paste(ht.treeshrub$site,ht.treeshrub$species, ht.treeshrub$no, sep = "_")
ht.dbh <- ht.treeshrub[, c("species","sample","type","date","site","no","ht","dbh","dbh2","dbh3")]
head(ht.dbh)

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#

#LCC/LNN 
cn <- read.csv("western/uncleaned/Traits_cn_2019.csv")

head(cn)
sort(unique(ht.dbh$sample))
sort(unique(cn$sample))

cn$temp <- cn$sample
cn.spp <- cn %>% separate(temp, c("site", "species","no"))

sort(unique(cn.spp$species))

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# Stem specific density

vol <- read.csv("western/uncleaned/traits_stemvol.csv")
wt <- read.csv("western/uncleaned/stem.wt.total.csv")

vol$site[vol$site == "afrf"] <- "af"

vol.spp <- vol[vol$species %in% species, ]; sort(unique(vol.spp$species))
vol.spp$sample <- paste(vol.spp$site,vol.spp$species, vol.spp$no, sep = "_")

wt.spp <- wt[wt$species %in% species, ]; sort(unique(wt.spp$species))
wt.spp$sample <- paste(wt.spp$site,wt.spp$species, wt.spp$no, sep = "_")


ssd <- merge(vol.spp, wt.spp, by = c("species","no","site","sample"), all = TRUE)

temp <- ssd
temp$count <- 1

trt <- temp %>%
  group_by(sample) %>%
  summarise(no_rows = sum(count), .groups = 'drop')

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#

lma <- read.csv("western/uncleaned/mergedareamass_Jan2022.csv")
lma$species[lma$species == "symalb "] <- "symalb"

lma.spp <- lma[lma$species %in% species, ]; sort(unique(lma.spp$species))
lma.spp$sample <- paste(lma.spp$site,lma.spp$species, lma.spp$indiv.no, sep = "_")

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
