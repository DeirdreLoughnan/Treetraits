### Started January 29 2019 ###

## DL writing dummy data for 507 Project ##

#The aim of this code is to build test data for modeling the relationship between phenology and functional traits


rm(list=ls()) 
options(stringsAsFactors = FALSE)

#For this early model, I will have the following variables that need estimates:
#Traits: SLA, ht, m.dbh, wood density, C:N
#Phenology: Chill, warm, photo

unique(comb)
setwd("~/Documents/github/Treetraits")
source('Rcode/source/Cleaning_compiling_data.R')

head(comb)
str(comb)

#Starting by looking for high correlations between traits
library(psych)
pairs.panels(comb[,c(20:21, 28:32)])

# Correlation panel
#there are many ways to do this, but this provides a busy but I think useful visual
library("PerformanceAnalytics")
chart.Correlation(comb[,c(20:21, 28:32)], histogram=TRUE, pch=19)

length(unique(comb$sp))

#To start with the simplest model possible, I am going to develop a linear model with the five traits and have the growth chamber treatment as a categorical variable, I am also using the indiviudal mean per treatment but this needs to be better thought through

nsite = 2 #number of sites
nsp = 28 #number of species

nwarm = 2 
nphoto = 2
nchill = 3
ntrt<-nwarm*nphoto*nchill #there are 12 treatments
length(unique(comb$treatcode)) #perfect, this matches, so i can just use the treatment code as a dummy variable for now
rep = 10 # I think the maximum number is 10, but for many it is closer to 5
ntot<-ntrt*rep*nsite; ntot
#240
#############################################
site.coef
#############################################
#Here I have species as the random effect, ignoring for now population 
inter=20

#allow intercepts to differ by species 
spint<-inter +c(1:nsp)-mean(1:nsp) #spint = species intercept

head(comb)
range(comb$Stomatal.Density, na.rm=TRUE)
unique(comb$Height)
(115.4456-547.2000)/2

fake<-vector()
for(i in 1:nsp){
  cfp=rnorm(ntot[i], 5,2) #gc treatment effect coeff
  sla=rnorm(ntot[i], 250,50)
  ht=rnorm(ntot[i], 10,5) #Based on personal experience, mean will be 10, but since there are  woody shrubs, the SD will need to be high
  cn=rnorm(ntot[i], 20,10) #This one I have less expert knowledge but usally C >>N so it should be positive
  wood=rnorm(ntot[i], 2,1)
  stom=rnorm(ntot[i], 500,50)
  
#This is for calculating the random species intercepts?
  
  cfp.c= #gc treatment effect coeff
  sla.c
  ht.c=
  cn.c=
  wood.c=
  stom.c=
    
  cfp.sd= #gc treatment effect coeff
  sla.sd=
  ht.sd=
  cn.sd=
  wood.sd=
  stom.sd=
    
    
  coeff<-c(spint[i],
           rnorm(1,2,0.5)
           )
}





#for the growth chamber model the ntot would be:
ntot.gc<-nsite*nwarm*nphoto*nchill*rep# 240
ntot

# Build up the data frame
site = gl(nsite, rep, length = ntot)
warm = gl(nwarm, rep*nsite, length = ntot); warm
photo = gl(nphoto, rep*nsite*nwarm, length = ntot); photo
chill = gl(nchill, rep*nsite*nwarm*nphoto, length = ntot); chill

#making dummy variables for the other chill treatment
chill1 = ifelse(chill == 2, 1, 0); chill1 
chill2 = ifelse(chill == 3, 1, 0); chill2 
