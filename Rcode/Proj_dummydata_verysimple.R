### Started January 29 2019 ###

## DL writing dummy data for 507 Project ##

#The aim of this code is to build test data for modeling the relationship between phenology and functional traits


rm(list=ls()) 
options(stringsAsFactors = FALSE)

#For this early model, I will have the following variables that need estimates:
#Traits: SLA, ht, m.dbh, wood density, C:N

unique(comb)
setwd("~/Documents/github/Treetraits")
source('Rcode/source/Cleaning_compiling_data.R')


#Starting by looking for high correlations between traits

# Correlation panel
#there are many ways to do this, but this provides a busy but I think useful visual
library("PerformanceAnalytics")
chart.Correlation(comb[,c(20:21, 28:32)], histogram=TRUE, pch=19)


#To start with the simpler, I am going to develop a linear model with the five traits
#in the source file the mean phenology is taken for each indiviudal tree of for a given treatment factor
int<-11
nsite = 2 #number of sites
nsp = 28 #number of species
rep = 10 # I think the maximum number is 10, but for many it is closer to 5
ntot<-nsite*nsp*rep #560

#Building the required dataframe for all of the replication and traits
site=gl(nsite, rep, length=ntot); site

#mean values for each trait for 100 data points
# sla.v= rnorm(1e2, 250, 50)
# ht.v=rnorm(1e2, 10, 5) 
# cn.v=rnorm(1e2, 20, 5)
# wood.v=rnorm(1e2, 2, 1)
# stom.v=rnorm(1e2, 500, 50)

sla.v= rnorm(ntot, 250, 50)
ht.v=rnorm(ntot, 10, 5) #sigma should be really high bc it includes trees and woody shrubs
cn.v=rnorm(ntot, 20, 5)
wood.v=rnorm(ntot, 2, 1)
stom.v=rnorm(ntot, 500, 50)

#I am uncertain how to calculate dummy data for y
#the y is phenology, so it has to be between 0 and the duration of the experiment
x<-seq(1:45)
sigma<-0.5

#effect sizes
sla.diff= -5
ht.diff=3
cn.diff=-2
wood.diff=3
stom.diff=0.5

for(i in 1:ntot){
  sla=rnorm(ntot[i], 250,50)
  ht=rnorm(ntot[i], 10,5) #Based on personal experience, mean will be 10, but since there are  woody shrubs, the SD will need to be high
  cn=rnorm(ntot[i], 20,10) #This one I have less expert knowledge but usally C >>N so it should be positive
  wood=rnorm(ntot[i], 2,1)
  stom=rnorm(ntot[i], 500,50)
}

# #sigmas
# sla.sd=50
# ht.sd=5
# cn.sd=5
# wood.sd=1
# stom.sd=50

df<-data.frame(site,sla.v,ht.v,cn.v,wood.v,stom.v)

summary(lm(bb~(site+sla+ht+cn+wood+stom),data=df))

# #allow intercepts to differ by species 
# spint<-inter +c(1:nsp)-mean(1:nsp) #spint = species intercept
d<-data.frame(site)
for(i in 1:length(nsp)){
  cfp=rnorm(ntot[i], 5,2) #gc treatment effect coeff
  sla=rnorm(ntot[i], 250,50)
  ht=rnorm(ntot[i], 10,5) #Based on personal experience, mean will be 10, but since there are  woody shrubs, the SD will need to be high
  cn=rnorm(ntot[i], 20,10) #This one I have less expert knowledge but usally C >>N so it should be positive
  wood=rnorm(ntot[i], 2,1)
  stom=rnorm(ntot[i], 500,50)
}


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
