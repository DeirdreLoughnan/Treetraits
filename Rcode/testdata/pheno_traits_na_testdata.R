## Started Sept 8, 2019 ##

## DL writing code for fake data to test for relationships between functional traits and phenology under variable climate conditions

rm(list=ls()) 
options(stringsAsFactors = FALSE)

#For this model, I will use a single trait with C, F, and P
#partial pooling

setwd("~/Documents/github/Treetraits/Rcode/testdata")

int <-11
sigma <-0.1

nsite = 4 #number of sites with trait data
nsp = 28 #number of species
rep = 50 #
ntot<-nsite*nsp*rep #5600
