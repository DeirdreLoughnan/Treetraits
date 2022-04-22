# Started Feb 4, 2022

# Aim of this code is to combine the eastern and western trait data into a single file
rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("~/Documents/github/Treetraits/data")

source("cleaning_western_traits.R")

source("cleaning_eastern_traits.R")

head(trt)

head(west)
#Start by changing column names in the datasets so they are least match
colnames(trt)[colnames(trt)=="Individual"] <- "sample"
colnames(trt)[colnames(trt)=="Site"] <- "site"
colnames(trt)[colnames(trt)=="Species"] <- "species"
colnames(trt)[colnames(trt)=="X.N"] <- "per.N"
colnames(trt)[colnames(trt)=="X.C"] <- "per.C"
colnames(trt)[colnames(trt)=="Height"] <- "ht"
colnames(trt)[colnames(trt)=="DBH"] <- "dbh"
colnames(trt)[colnames(trt)=="DBH.2"] <- "dbh2"
colnames(trt)[colnames(trt)=="DBH.3"] <- "dbh3"
colnames(trt)[colnames(trt)=="Stomatal.Density"] <- "stomatal.density"

trt$species <- tolower(trt$species)
# fix the names of species on the 
#Many of the columns are character, but should be numeric with NA
trt.sub <- trt[,c("sample","site","species","ssd","per.N","per.C","ht","dbh","dbh2","dbh3","lma","C.N")]

west.sub <- west[, c("sample","site","species","type","ssd","per.N","per.C","C.N","ht","dbh","dbh2","dbh3","lma")]

allTrt <- rbind.fill(trt.sub, west.sub)
write.csv(allTrt, "allTrt.csv", row.names = F)

hist(west$lma)
hist(trt$lma)
