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
# rm(list=ls()) 
# options(stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(stringr)
library(plyr)

# if(length(grep("deirdreloughnan", getwd()) > 0)) {
#   setwd("~/Documents/github/Treetraits/data")
# }  else{
#   setwd("/home/deirdre/treetrait") # for midge
# }

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
 
 incom.smpl <- c("mp_acegla10_02","mp_loninv_4","mp_menfer_8","mp_rubpar_21","mp_rubpar_22", "mp_rubpar_7_2","sm_symalb_9_02","sm_alnvir_7", "sm_alnvir_1_June3","sm_alnvir_1", "sm_alnvir_16", "sm_betpap_7_2", "sm_poptre_10","sm_poptre_14","sm_riblac_8", "af_betpap_1","af_poptre_11","af_poptre12_02","af_rubpar_19","af_rubpar9","af_vibedu_3_02","kl_betpap_9_02") #22
 #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
 
 # Start with height:
ht <- read.csv("western/uncleaned/traits_ht_calc_dbh_final.csv")

head(ht) 
sort(unique(ht$species)) # has extra species
ht$site <- as.character(ht$site)
ht$species <- as.character(ht$species)
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
ht.dbh <- ht.dbh[!ht.dbh$sample %in% incom.smpl, ]
length(unique(ht.dbh$sample)) # 789

ht.dbh$count <- 1
tempy <- aggregate(ht.dbh$count, ht.dbh[c("species","site")], FUN = sum)

# remove date column
ht.dbh <- ht.dbh[, c("sample","species","site","no","type","ht","dbh","dbh2","dbh3")]
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#

#LCC/LNN 
# cn <- read.csv("western/uncleaned/Traits_cn_2019.csv")
# 
# dup <- c("af_acegla_1_2", "af_poptre_3_2","af_poptre_8_2","af_rubpar_13_2","af_rubpar_11_2","af_rubpar_10_2","af_rubpar_3_2","af_rubpar_2_2","af_rubpar_4_2","af_rubpar_5_2","af_rubpar_6_8", "af_vacmem_6_2","mp_acegla_10_2","mp_rubpar_7_2","mp_alninc_8_2","mp_spipyr_10_2","kl_riblac_2_2", "kl_riblac_10_2", "kl_loninv_8_2","kl_popbal_7_2","mp_alninc_8_2","mp_vacmem_5_2","mp_sorsco_10_2","sm_alnvir_10_2","sm_alnvir_2_2","sm_alnvir_7_ai")
# 
# cn <- cn[!cn$sample %in% dup, ]
# 
# head(cn)
# # sort(unique(ht.dbh$sample))
# # sort(unique(cn$sample))
# 
# cn$temp <- cn$sample
# cn.spp <- cn %>% separate(temp, c("site", "species","no"))
# 
# sort(unique(cn.spp$species))
# 
# length(unique(cn.spp$sample))
# cn.spp <- cn.spp[!cn.spp$sample %in% incom.smpl, ]
# length(unique(cn.spp$sample)) # 800 - dups of some samples
# 
# cn.spp$count <- 1
# tempy <- aggregate(cn.spp$count, cn.spp[c("species","site")], FUN = sum)

# remove date column:
cn.cl <- read.csv("western/cleaned/cn_cleaned.csv")
dup <- c("af_acegla_1_2", "af_poptre_3_2","af_poptre_8_2","af_rubpar_13_2","af_rubpar_11_2","af_rubpar_10_2","af_rubpar_3_2","af_rubpar_2_2","af_rubpar_4_2","af_rubpar_5_2","af_rubpar_6_8", "af_vacmem_6_2","mp_acegla_10_2","mp_rubpar_7_2","mp_alninc_8_2","mp_spipyr_10_2","kl_riblac_2_2", "kl_riblac_10_2", "kl_loninv_8_2","kl_popbal_7_2","mp_alninc_8_2","mp_vacmem_5_2","mp_sorsco_10_2","sm_alnvir_10_2","sm_alnvir_2_2","sm_alnvir_7_ai","af_popbal_2_2")

cn.cl <- cn.cl[!cn.cl$sample %in% dup, ]
cn.spp <- cn.cl[, c("sample","species","site","no","Weight..mg.","per.N","per.C","C.N")]
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# Stem specific density

vol <- read.csv("western/uncleaned/traits_stemvol.csv")
wt <- read.csv("western/uncleaned/stem.wt.total.csv")

vol$site <- as.character(vol$site)
vol$site[vol$site == "afrf"] <- "af"

vol.spp <- vol[vol$species %in% species, ]; sort(unique(vol.spp$species))
vol.spp$sample <- paste(vol.spp$site,vol.spp$species, vol.spp$no, sep = "_")

wt.spp <- wt[wt$species %in% species, ]; sort(unique(wt.spp$species))
wt.spp$sample <- paste(wt.spp$site,wt.spp$species, wt.spp$no, sep = "_")

ssd <- merge(vol.spp, wt.spp, by = c("species","no","site","sample"), all = TRUE)

ssd$ssd <- ssd$stem.weight/ssd$vol

length(unique(ssd$sample))
ssd <- ssd[!ssd$sample %in% incom.smpl, ]
length(unique(ssd$sample)) # 790

ssd$count <- 1

tempy <-  aggregate(ssd["count"], ssd[c("site","species")], FUN = sum)

ssd <- ssd[, c("sample","species","site","no","vol","stem.weight","ssd")]
ssd.comp <- ssd[complete.cases(ssd),]
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#

lma <- read.csv("western/uncleaned/mergedareamass_final.csv")
lma$species[lma$species == "symalb "] <- "symalb"

# Smithers had some miss identified Rhoalb 

lma.spp <- lma[lma$species %in% species, ]; sort(unique(lma.spp$species))
lma.spp$sample <- paste(lma.spp$site,lma.spp$species, lma.spp$indiv.no, sep = "_")

# last minute corrections:
# lma.spp$sample <- as.character(lma.spp$sample)
lma.spp$sample[lma.spp$sample == "kl_symalb_02_02"] <- "kl_symalb_2_02"
lma.spp$sample[lma.spp$sample == "mp_corsto_10_2"] <- "mp_corsto_10_02"
lma.spp$sample[lma.spp$sample == "mp_menfer_10_2"] <- "mp_menfer_10_02"

lma.spp$lma <- lma.spp$leaf.mass/lma.spp$leaf.area

length(unique(lma.spp$sample))
lma.spp <- lma.spp[!lma.spp$sample %in% incom.smpl, ]
length(unique(lma.spp$sample)) # 800 - dups of some samples

# lma.spp$count <- 1
# tempy <-  aggregate(lma.spp["count"], lma.spp[c("site","species")], FUN = sum)

lma.mean <-  aggregate(lma.spp["lma"], lma.spp[c("species", "site","indiv.no","sample")], FUN = mean)
colnames(lma.mean)[colnames(lma.mean) == "indiv.no"] <- "no"

#remove collected by:
lma.mean <- lma.mean[, c("sample","species","site","no","lma")]
lma.comp <- lma.mean[complete.cases(lma.mean),]

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
dim(ht.spp)

# write.csv(ssd, "western/cleaned/ssd_cleaned.csv", row.names = F)
# write.csv(cn.spp, "western/cleaned/cn_cleaned.csv", row.names = F)
# write.csv(ht.dbh, "western/cleaned/ht_dbh_cleaned.csv", row.names = F)
# write.csv(lma.temp, "western/cleaned/lma_cleaned.csv", row.names = F)

length(unique(ssd$sample))
length(unique(cn.spp$sample))
length(unique(ht.dbh$sample))
length(unique(lma.spp$sample)) # 769
length(unique(lma.comp$sample)) # 769


# # Let's get a list of the samples I have complete data for:
dim(ssd)
dim(cn.spp)
ssd.cn <- merge(ssd,cn.spp, by = c("sample", "site", "no", "species"), all =T)
ssd.cn.ht <- merge(ssd.cn, ht.dbh, by = c("sample", "no", "site", "species"), all =T)
west <- merge(ssd.cn.ht, lma.comp, by = c("sample", "site", "species"), all =T)

# 2 rows with no trait data:

length(unique(west$sample))
rm <- c("af_riblac_9","sm_menfer_2")

west <- west[!west$sample %in% rm, ]

#write.csv(west, "western/cleaned/westTrait.csv", row.names = F)
