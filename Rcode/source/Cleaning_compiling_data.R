### Started January 29 2019 ###

## DL writing merging data for phenology and functional traits ##

#The aim of this code is to combine the phenology and functional trait data collected by DF in 2015
rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("~/Documents/github/Treetraits")

phen<-read.csv("input/Budburst.csv", header=T, na.strings=c("","NA"))

trt<-read.csv("input/Tree_Traits_2015.csv",header=T, na.strings=c("","NA"))

str(trt)
#remove columns of logistics (23-25)
trt<-trt[,c(1:12,18:22,26:29)]

#Start by changing column names in the datasets so they are least match
colnames(trt)[colnames(trt)=="Individual"] <- "ind"
colnames(trt)[colnames(trt)=="Site"] <- "site"
colnames(trt)[colnames(trt)=="Species"] <- "sp"
colnames(trt)[colnames(trt)=="X.N"] <- "per.N"
colnames(trt)[colnames(trt)=="X.C"] <- "per.C"

#Many of the columns are character, but should be numeric with NA
str(trt)
trt$DBH<-as.numeric(as.character(trt$DBH)); unique(trt$DBH)
trt$DBH.2<-as.numeric(as.character(trt$DBH.2)); unique(trt$DBH.2)
trt$DBH.3<-as.numeric(as.character(trt$DBH.3)); unique(trt$DBH.2)
trt$DBH.4<-as.numeric(as.character(trt$DBH.4)); unique(trt$DBH.2)
trt$DBH.5<-as.numeric(as.character(trt$DBH.5)); unique(trt$DBH.2)

str(trt)
trt$per.C<-as.numeric(as.character(trt$per.C)); unique(trt$per.C)
trt$per.N<-as.numeric(as.character(trt$per.N)); unique(trt$per.N)

names(phen)
# head(phen)

names(trt)
#head(trt)

#what Col should be in the new dataset used for Project 2 in 507?
# id, sp, rep, site, ind, treatcode, warm, photo, chill, Term.fl, Lat.fl, Term.lf, Lat.lf,tlea,leaf

#Individual, Site, Species, Latitude, Longitude, Elevation, Leaf.area, Fresh.mass, Dry.mass, Stem.vol, Stem.mass, Height, DBH, X.N, X.C, Stomatal. Length, Stomatal.Density


#Note that the individual ID are different between the two datasets; I am going to try and break them apart for sorting

# library(reshape2)
# 
# #I don't fully understand what this code is doing, but it does work to split the character from the numeric parts of the species ID
# #https://stackoverflow.com/questions/9756360/split-character-data-into-numbers-and-letters
# temp<-colsplit(trt$Individual, "(?<=\\p{L})(?=[\\d+$])", c("char", "digit"))
# 
# #A bit hacky,but the order appears to be the same, so this shouldn't be an issue. 
# trt2<-cbind(trt,temp)
# 
# #now doing the same for the phenology data
# tempphen<-colsplit(phen$id, "(?<=\\p{L})(?=[\\d+$])", c("char", "digit"))
# 
# 
# head(trt2)
# 
# #Now paste them together to get a species ID that matches the one in the phenology dataset, need to further sep out replicates 
# trt2$code<-paste(trt2$char, trt2$digit, sep = '_')
# 
# head(trt2)

#Subset both for only the two populations both phenology and trait data were collected:
unique(phen$site)
unique(trt$site)

trtsites<-subset(trt, site==c("SH","HF"))

#Now I want to put together the phenology dataset and the trait data. How inconsistent are these two datasets?
#Very different in length, what about the number of trees and 

#Note id= twig level id
length(unique(phen$ind)) #274
length(unique(trtsites$ind)) #2196

#I will think on this more, but here I am taking the mean of the phenology for a given combination of 
require(plyr)
avg.phen<-ddply(phen, c("ind", "sp","site","treatcode","warm","photo","chill"), summarise, 
                Term.flm=mean(Term.fl, na.rm=TRUE),
                Lat.flm=mean(Lat.fl, na.rm=TRUE),
                # Term.lfm=mean(Term.lf, na.rm=TRUE),
                # Lat.lfm=mean(Lat.lf, na.rm=TRUE),
                tleafm=mean(tleaf, na.rm=TRUE),
                lleafm=mean(lleaf, na.rm=TRUE))


new<-left_join(avg.phen, trtsites, by= c("ind","sp","site"))

comb<-subset(new, Latitude>0 & Leaf.area>0)


#Now calculating the required traits: sla, wood density, 

names(comb)
comb$sla<-comb$Leaf.area/comb$Dry.mass
comb$wood.den<-comb$Stem.volume/comb$Stem.mass
#comb$m.dbh<-colMeans(comb[,26:30], na.rm=TRUE) #this doesn't do what it should            
comb$cn<-(comb$per.C/comb$per.N)
str(comb)
#Now I can prune the dataset to just the values I will be working with for this project
