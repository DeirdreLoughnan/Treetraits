## Started Sept 8, 2019 ##

## DL writing code for fake data to test for relationships between functional traits and phenology under variable climate conditions

rm(list=ls()) 
options(stringsAsFactors = FALSE)

#For this model, I will use a single trait with C, F, and P
#partial pooling

setwd("~/Documents/github/Treetraits/Rcode/testdata")

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
# Set up: 4 sites, 28 species in eastern dataset, two levels each of warming and photoperiod, and three levels of chilling. 
# Data from western dataset will have 26 species, two levels each of warming and photoperiod, and two levels of chilling
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

nsite = 4 #number of sites with trait data
nsp = 28 #number of species
nind = 12

rep = 1 # making it greater for the test data

#number of treatment levels
nwarm = 2
nphoto = 2
nchill = 3

ntot<-nsite*nwarm*nphoto*nchill #48 
ntotsp<-nsite*nwarm*nphoto*nchill*nsp
# code below to loop through species and indiviudals 

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
#building the required dataframe 
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

site = gl(nsite, rep, length = ntot)
warm = gl(nwarm, rep*nsite, length = ntot)
photo = gl(nphoto, rep*nsite*nwarm, length = ntot)
chill = gl(nchill, rep*nsite*nwarm*nphoto, length = ntot)

# to start I am going to focus on one trait only
slav=rnorm(ntotsp, 5, 1)
# htv=rnorm(ntot, 11, 3) 
# cnv=rnorm(ntot, 10, 2)
# woodv=rnorm(ntot, 0.7, 0.05)

treatcombo = paste(warm, photo, chill, sep = "_")

d <- data.frame(treatcombo, site, warm, photo, slav)

###### Set up differences for each level
sitediff =2
sitediff2 = 2 
sitediff3 = 2 
sitediff4 = 2 
warmdiff = -20 # days earlier from 1 to 2
photodiff = -14
chilldiff = -20 # need to think more about this number and what would be reasonable
sladiff = -0.5

# htdiff=0.5
# cndiff=-0.5
# wooddiff=1

# interactions. 9 two-way interactions
sitewarm = 0
sitewarm2 = 0
sitewarm3 = 0
sitewarm4 = 0

sitephoto = 0
sitephoto2 = 0
sitephoto3 = 0
sitephoto4 = 0
warmphoto = 3.5
warmchill = -7 #negative, warmer days after longer chilling will advance budburst 
chillphoto = -3.5 #negative, but smaller, longer photoperiod will result in earlier budburst, but to a lesser extent


####### SD for each variable

sitediff.sd = 1.5 
warmdiff.sd = 1 
photodiff.sd = 1
chilldiff.sd =1
sitewarm.sd = 1
site2photo.sd = 1
site3photo.sd = 1
site4photo.sd = 1
warmphoto.sd = 1
warmchill.sd = 1
chillphoto.sd = 1

mm<-model.matrix(~(site+warm+photo)^2+slav, data.frame(site, warm, photo,slav))
head(mm)
colnames(mm)

fake<- vector()
for (i in 1:nsp){
  coeff <- c(1, sitediff2,sitediff3,sitediff4, warmdiff, photodiff, sladiff,
             sitewarm2,sitewarm3,sitewarm4, sitephoto2, sitephoto3,sitephoto4,
             warmphoto)
  bb <- rnorm(n = length(warm), mean = mm %*% coeff, sd = 1)
  faket<-merge(bb,mm)
  fakex<-data.frame(faket, sp=i)
  fake<-rbind(fake, fakex)
}

#I think it worked! I need to think more about incorporating the trait data, right now every site has a unique trait value. 