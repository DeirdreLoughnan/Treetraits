# Started January 6, 2022 by Deirdre 

#cleaning trait data for western data:
# Traits:
# LMA - mass/area
#stem specific density - mass/ volume
# LCC/LNN
# height 
# DBH - of primary, secondary, tertiary stem

# should have about 780 - but with missing indiv = 769

# I think I will need two datasets: 
# 1. Individual level metrics - height, DBH, SSD, LCC/LNN, 
# 2. Leaf level - SLA - measured five leaves
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(stringr)
library(plyr)

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
ht$species[ht$species == "amelan"] <- "amealn"
ht$species[ht$species == "betpap "] <- "betpap"
ht$species[ht$species == "spibet "] <- "spibet"

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

length(unique(ht.dbh$sample))

# ht.dbh$count <- 1
# 
# trt <- ht.dbh %>%
#   group_by(site, species,) %>%
#   summarise(no_rows = sum(count), .groups = 'drop')

head(ht.dbh)
# remove date column
ht.dbh <- ht.dbh[, c("sample","species","site","no","type","ht","dbh","dbh2","dbh3")]
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#

#LCC/LNN 
cn <- read.csv("western/uncleaned/Traits_cn_2019.csv")

head(cn)
sort(unique(ht.dbh$sample))
sort(unique(cn$sample))

cn$temp <- cn$sample
cn.spp <- cn %>% separate(temp, c("site", "species","no"))

sort(unique(cn.spp$species))

# remove date column:
cn.spp <- cn.spp[, c("sample","species","site","no","Weight..mg.","per.N","per.C","C.N")]
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

ssd$ssd <- ssd$stem.weight/ssd$vol
ssd.temp <- ssd
ssd.temp$count <- 1

sdd.tempy <-  aggregate(ssd.temp["count"], ssd.temp[c("sample")], FUN = sum)

ssd <- ssd[, c("sample","species","site","no","vol","stem.weight","ssd")]
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#

lma <- read.csv("western/uncleaned/mergedareamass_final.csv")
lma$species[lma$species == "symalb "] <- "symalb"

# Smithers had some miss identified Rhoalb 

lma.spp <- lma[lma$species %in% species, ]; sort(unique(lma.spp$species))
lma.spp$sample <- paste(lma.spp$site,lma.spp$species, lma.spp$indiv.no, sep = "_")

lma.spp$lma <- lma.spp$leaf.mass/lma.spp$leaf.area

lma.temp <-  aggregate(lma.spp["lma"], lma.spp[c("species", "site","indiv.no","sample")], FUN = mean)
colnames(lma.temp)[colnames(lma.temp) == "indiv.no"] <- "no"

#remove collected by:
lma.spp <- lma.spp[, c("sample","species","site","indiv.no","leaf.area","leaf.mass","lma")]
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
incom.smpl <- c("mp_acegla10_02","mp_loninv_4","mp_menfer_8","mp_rubpar_21","mp_rubpar_22", "mp_rubpar_7_2","sm_symalb_9_02","sm_alnvir_7", "sm_alnvir_1_June3","sm_alnvir_1", "sm_alnvir_16", "sm_betpap_7_2", "sm_poptre_10","sm_poptre_14","sm_riblac_8", "af_betpap_1","af_poptre_11","af_poptre12_02","af_rubpar_19","af_rubpar9","af_vibedu_3_02","kl_betpap_9_02") #24

length(unique(ssd$sample))
length(unique(cn.spp$sample))
length(unique(ht.dbh$sample))
length(unique(lma.spp$sample)) # 769
length(unique(lma.temp$sample)) # 769

# Let's get a list of the samples I have complete data for:

temp <- merge(ssd,cn.spp, by = c("sample", "site", "no", "species"), all =T)
temp2 <- merge(temp, ht.dbh, by = c("sample", "site", "no", "species"), all =T)
temp3 <- merge(temp2, lma.temp, by = c("sample", "site", "no", "species"), all =T)

ht.temp <- temp3[complete.cases(temp3$ht),]

ht.temp2 <- ht.temp[!ht.temp$sample %in% incom.smpl, ]

ht.temp2$count <- 1
tempy <- aggregate(ht.temp2$count, ht.temp2[c("species","site")], FUN = sum)

write.csv(ht.temp, "temp_traitdata.csv", row.names = F)
# smp <- unique(ssd$sample)
# 
# t <- lma.spp[!lma.spp$sample %in% smp, ]
# sort(unique(t$sample))
# 
# lma.spp$count <- 1
# lma.spp$site <- as.factor(lma.spp$site)
# lma.spp$species <- as.factor(lma.spp$species)
# 
# aggregate(lma.spp$count, lma.spp[c("species", "site")], FUN = sum)

###########################################
# subset the other traits to match the cleaned leaf data:

smp <- unique(lma.spp$sample)

ssd.smp <- ssd[ssd$sample %in% smp, ]
length(unique(ssd.smp$sample))

ssd.smp$count <- 1
temp <- aggregate(ssd.smp$count, ssd.smp[c("sample")], FUN = sum)

#####
cn.smp <- cn.spp[cn.spp$sample %in% smp, ]
length(unique(cn.smp$sample))

cn.smp$count <- 1
temp <- aggregate(cn.smp$count, cn.smp[c("sample")], FUN = sum)

#####
ht.dbh.smp <- ht.dbh[ht.dbh$sample %in% smp, ]
length(unique(ht.dbh.smp$sample))

ht.dbh.smp$count <- 1
temp <- aggregate(ht.dbh.smp$count, ht.dbh.smp[c("sample")], FUN = sum)
