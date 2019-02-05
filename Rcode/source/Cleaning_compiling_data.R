### Started January 29 2019 ###

## DL writing merging data for phenologyology and functional traits ##

#The aim of this code is to combine the phenologyology and functional trait data collected by DF in 2015
rm(list=ls()) 
options(stringsAsFactors = FALSE)

forlatex = TRUE # set to FALSE if just trying new figures, TRUE if outputting for final
runstan = FALSE # set to TRUE to actually run stan models. FALSE if loading from previous runs

# Analysis of bud burst experiment 2015. 

library(memisc) # for getSummary 
library(xtable)
library(scales) # for alpha
library(ggplot2)
library(caper) # for pgls
library(png) # readPNG for Fig 1

setwd("~/Documents/github/Treetraits")

#phenology<-read.csv("input/Budburst.csv", header=T, na.strings=c("","NA"))
(toload <- sort(dir("./input")[grep("Budburst Data", dir('./input'))], T)[1])

load(file.path("input", toload))

if(forlatex) figpath = "../docs/ms/images" else figpath = "graphs"

phenology<- dx

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

names(phenology)
# head(phenology)

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
# #now doing the same for the phenologyology data
# tempphenology<-colsplit(phenology$id, "(?<=\\p{L})(?=[\\d+$])", c("char", "digit"))
# 
# 
# head(trt2)
# 
# #Now paste them together to get a species ID that matches the one in the phenologyology dataset, need to further sep out replicates 
# trt2$code<-paste(trt2$char, trt2$digit, sep = '_')
# 
# head(trt2)

#Subset both for only the two populations both phenologyology and trait data were collected:
unique(phenology$site)
unique(trt$site)

trtsites<-subset(trt, site==c("SH","HF"))

#Now I want to put together the phenologyology dataset and the trait data. How inconsistent are these two datasets?
#Very different in length, what about the number of trees and 

#Note id= twig level id
length(unique(phenology$ind)) #274
length(unique(trtsites$ind)) #2196

#I will think on this more, but here I am taking the mean of the phenologyology for a given combination of 
require(plyr)
require(dplyr)
avg.phenology<-ddply(phenology, c("ind", "sp","site","treatcode","warm","photo","chill"), summarise, 
                Term.flm=mean(Term.fl, na.rm=TRUE),
                Lat.flm=mean(Lat.fl, na.rm=TRUE),
                # Term.lfm=mean(Term.lf, na.rm=TRUE),
                # Lat.lfm=mean(Lat.lf, na.rm=TRUE),
                tleafm=mean(tleaf, na.rm=TRUE),
                lleafm=mean(lleaf, na.rm=TRUE))


new<-left_join(avg.phenology, trtsites, by= c("ind","sp","site"))

comb<-subset(new, Latitude>0 & Leaf.area>0)


#Now calculating the required traits: sla, wood density, 

#want sla to be in mm^2/g
comb$Leaf.area.sq<-comb$Leaf.area*100
comb$Dry.mass.g<-comb$Dry.mass*1000
names(comb)
comb$sla<-comb$Leaf.area.sq/comb$Dry.mass.g
comb$wood.den<-comb$Stem.volume/comb$Stem.mass
#comb$m.dbh<-colMeans(comb[,26:30], na.rm=TRUE) #this doesn't do what it should            
comb$cn<-(comb$per.C/comb$per.N)
str(comb)
#Now I can prune the dataset to just the values I will be working with for this project
comb<-comb[,c(1:7,10:11,20:21,29,32:34) ]; head(comb)
