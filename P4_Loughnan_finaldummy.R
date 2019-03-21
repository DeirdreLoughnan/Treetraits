### Started January 29 2019 ###

## DL writing dummy data for 507 Project ##

#The aim of this code is to build test data for modeling the relationship between phenology and functional traits
#Update: Added partial pooling across species & also changed data to not exclude sites with no phenology data
# March 19 adding effect of chilling, forcing and photoperiod

rm(list=ls()) 
options(stringsAsFactors = FALSE)

#For this early model, I will have the following variables:
#Traits: SLA, height (ht), m.dbh, wood density, C:N, & stomatal density
#Note: After looking for collinearity in the x varibales, I decided to exclude DBH and just use ht, it is biologically more relevant as well

library(rethinking)

setwd("~/Documents/github/Treetraits")
#source('Rcode/source/Cleaning_compiling_data.R') #Here I combined the budburst and trait data, averaging over individual trees for a given treatment
##################################################
# Making test data
###################################################
#Now developing a model with all traits and partial pooling across species
#This code is devleoped from existing code: #FakeBudburst_Generate_ind.R
int<-11 # days into the experiment of bb 
sigma<-2

nsite = 2 #number of sites
nsp = 28 #number of species
rep = 1 # individuals will be looped over below

#Adding growth chamber component
nwarm=2
nphoto=2
nchill=3

ntot<-nsite*nwarm*nphoto*nchill*rep #FakeBudburst_Generate_ind.R, nsp was not included here

#Building the required dataframe for all of the replication and traits
#Assigning a variable to site **This will be returned to later when partially pooling across sites
site=gl(nsite, rep, length=ntot)
warm=gl(nwarm, rep*nsite, length=ntot)
photo=gl(nphoto, rep*nsite*nwarm, length=ntot)
chill=gl(nchill, rep*nsite*nwarm*nphoto, length=ntot)

treatcombo= paste(warm,photo,chill, sep="_")

d<-data.frame(site,warm,photo, chill, treatcombo); head(d); unique(d$treatcombo)
#randomly generating numbers for each trait based on what I think a small tree would be
#sitev=rnorm(ntot,5,2) #Ultimatley I will be pooling across sites, so this is not necessary at this point in time

# Trait values (trait name + v for value)
slav= rnorm(ntot, 5, 1); 
htv=rnorm(ntot, 11, 3) #sigma should be really high (3m) bc it includes trees and woody shrubs
cnv=rnorm(ntot, 10, 2)
woodv=rnorm(ntot, 0.5, 0.05)
stomv=rnorm(ntot, 200, 50)

#effect sizes
sitediff=2 # I think site 2 will be delayed by 2 days due to the 5 degree diff in lat
sladiff= -0.5
htdiff=0.5
cndiff=-0.5
wooddiff=0.3
stomdiff=1

warmdiff = -20 # days earlier from 1 to 2
photodiff = -14
chilldiff= -10

#making a matrix of the x dummy variables
mm <- model.matrix(~(site+slav+htv+cnv+woodv+stomv+warm+photo+chill), data.frame(site,slav,htv,cnv,woodv,stomv,warm,photo,chill))
mm

#Creating centered data:
# centering the x dummy data
mm.cent<-scale(mm[,3:7]);head(mm.cent)
#pairs(mm.cent) # These plots are just clouds of points and not helpful
mm.cent<-cbind(mm[,1:2], mm.cent)
head(mm.cent)

#making a simple data frame of the x dummy variables
df<-data.frame(site,slav,htv,cnv,woodv,stomv,warm,photo,chill)
range(df$stomv)
str(df)
#####################################################################
# Now adding separate intercepts for each species ** Starting very simple with just SLA

fake_spint<-vector()
sp_mean <- rnorm(nsp, mean=int, sd=5)
for (j in c(1:length(sp_mean))){
  #phen_spint<- int+sladiff*slav+rnorm(rep*nsite, sp_mean[j], sigma)
  phen_spint<- sladiff*slav+htdiff*htv+cndiff*cnv+wooddiff*woodv+stomdiff*stomv+warm*warmdiff+photo*photodiff+chill*chilldiff+rnorm(rep*nsite, sp_mean[j], sigma)
  temp<-data.frame(phen_spint)
  phen_spint_cent<- rbind(fake_spint, temp)
}
# fake_spint<-vector()
# sp_mean <- rnorm(nsp, mean=int, sd=sigma)
# for (j in c(1:length(sp_mean))){
#   phen_spint<- int+sladiff*slav+htdiff*htv+cndiff*cnv+wooddiff*woodv+stomdiff*stomv+rnorm(rep*nsite, sp_mean[j], sigma)
#   temp<-data.frame(phen_spint)
#   phen_spint_cent<- rbind(fake_spint, temp)
# }

### Using matrix mathematics 

mm <- model.matrix(~(site+warm+photo+chill), data.frame(site, warm, photo,chill))
colnames(mm)

coeff <- c(1, sitediff, warmdiff, photodiff,chilldiff
          , sladiff,htdiff,cndiff,wooddiff,stomdiff
          )

bb <- rnorm(n = length(warm), mean = mm %*% coeff, sd = 1) # should be able to do sd = mm %*% sd.coeff as well, with a different sd for each parameter.

(fake <- data_frame(bb, site, warm, photo))





head(phen_spint)
range(phen_spint)

#Add species - there is probably a better way to do this
sp=gl(nsp, rep*nsite, length= ntot)


#Combine the phenology data with the species identity 
phen_spint<-cbind(phen_spint,sp)
head(phen_spint)
###############################################
# Now combining the phenology and the trait data

fake_spint<- data.frame(phen_spint,mm)
head(fake_spint)

#Note Stan doesn't like dots so changing X.Intercept name

colnames(fake_spint)[colnames(fake_spint)=="X.Intercept."] <- "Intercept"

####### Uncentered data #############
# phen <- rnorm(n = length(site), mean = mm %*% coeff, sd = 1) # This code works but the values are HUGE 
# phen

#right now I am not including site in the model, so the df needs to be trimmed
#what role does the intercept column play?
#fake_spint_cent<-fake_spint_cent[,c(1:3,5:9)]; head(fake_spint_cent) #includes intercept
# Below line of code not working ... 
fake_spint_cent<-fake_spint[,c(1:2,5:9)]; head(fake_spint_cent) #does NOT include intercept
