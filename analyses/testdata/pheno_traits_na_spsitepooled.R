## Started Sept 8, 2019 ##

## DL writing code for fake data to test for relationships between functional traits and phenology under variable climate conditions

## Oct 5 update: trying to partially pool across both species and sites; for purpose of fake data, increasing the number of sites to 6

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
#nsite=6 
#nsp = 28 #number of species
nsp = 5 #I reduced this temporarily for the purpose of building the fake data, just to be sure it is doing what I think it is
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


#### To start, keeping it simiple and ignoring interacitons between climate factors
#### Set up differences for each level

#GC climate differences
warmdiff = -20 # days earlier from 1 to 2
photodiff = -14
chill2diff = -20 # need to think more about this number and what would be reasonable
chill3diff = -10

#trait differences
sladiff = -0.5
# htdiff=0.5
# cndiff=-0.5
# wooddiff=1

####### SD for each variable

warmdiff.sd = 1 
photodiff.sd = 1
chill2diff.sd =1
chill3diff.sd =3

sladiff.sd = 1

warmphoto.sd = 1
warmchill2.sd = 1
warmchill3.sd = 1

chill2photo.sd = 1
chill3photo.sd = 1


mm<-model.matrix(~(warm+photo+chill)+slav, data.frame( warm, photo,chill,slav))
head(mm)
colnames(mm)

fake<- vector()
for (i in 1:nsp){
  coeff <- c(1, warmdiff, photodiff,chill2diff,chill3diff, sladiff)
  bb <- rnorm(n = length(warm), mean = mm %*% coeff, sd = 1)
  faket<-merge(bb,mm)
  #faket<-data.frame(bb,mm)
  fakesp<-data.frame(faket, sp=i)
  fake<-rbind(fake, faksp)
}
str(fake)
unique(fake$sp)
#I think it worked! I need to think more about incorporating the trait data, right now every site has a different trait value but all GC treatments have the same value. This makes sense to me, since we might see differnces across the 4 pops, but I don't have the data to know how the traits would change with warming, chilling & photo periods.
names(fake)
length(fake$x);length(fake$site2)
summary(lm(x ~ (warm2+photo2+chill2+chill3+slav), data = fake))

#########################################################################################################################
# Adding site and species differences:
baseint = 35
spint<-baseint +c(1:nsp)-mean(1:nsp) #chaning the intercept for each indiv sp
sitint<-baseint +c(1:nsite)-mean(1:nsite)

indivint<-sitint[1]+spint[1] +1:nind-mean(1:nind)
indivint

for(j in 1:nsp){
  for (i in sitint) {
    indivint<-sitint[1]+spint[1] # this might be for indiviudal tree id, +1:nind-mean(1:nind)
fake<-vector()
length(nsite)
  coeff <- c(spint[1], 
             rnorm(1, warmdiff, warmdiff.sd),
             rnorm(1, photodiff, photodiff.sd), 
             rnorm(1, chilldiff2, chill2diff.sd),
             rnorm(1, chilldiff3, chill3diff.sd),
             rnorm(1, sladiff, sladiff.sd)
  )
  bb <- rnorm(n = length(warm), mean = mm %*% coeff, sd = 0.1)
  
  fakex <- data.frame(bb, sp = 1:j, mm)
  fakex2<-data.frame(fakex, site =1:nsite)
  
  fake <- rbind(fake, fakex2)  
  }
}

#240
12*5*4
View(fake)

#########################################################################################################################
# Trying to make things a bit more complicated, accounting for differences in responses across sites...this might not be the direction to actually head given I only have data for 10 individuals for each population
###### Set up differences for each level
sitediff =2
sitediff2 = 2 
sitediff3 = 2 
sitediff4 = 2 
warmdiff = -20 # days earlier from 1 to 2
photodiff = -14
chilldiff2 = -20 # need to think more about this number and what would be reasonable
chilldiff3 = -10
sladiff = -0.5

# htdiff=0.5
# cndiff=-0.5
# wooddiff=1

# interactions. 9 two-way interactions
sitewarm = 0
sitewarm2 = 0
sitewarm3 = 0
sitewarm4 = 0

sitechill22 = 0
sitechill23 = 0
sitechill24 = 0
sitechill32 = 0
sitechill33 = 0
sitechill34 = 0

sitephoto = 0
sitephoto2 = 0
sitephoto3 = 0
sitephoto4 = 0
warmphoto = 3.5
warmchill2 = -7 #negative, warmer days after longer chilling will advance budburst 
warmchill3= -10 #additive effects of longer chilling and warm spring
chillphoto2 = -3.5
chillphoto3 = -1.5#negative, but smaller, longer photoperiod and greater chilling will result in earlier budburst, but to a lesser extent
sitesla2 = 0
sitesla3 = 0
sitesla4= 0
warmsla =0 
photosla=0

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


mm<-model.matrix(~(site+warm+photo+chill)^2+slav, data.frame(site, warm, photo,chill,slav))
head(mm)
colnames(mm)

fake<- vector()
for (i in 1:nsp){
  coeff <- c(1, sitediff2,sitediff3,sitediff4, warmdiff, photodiff,chilldiff2,chilldiff3, sladiff,
             sitewarm2,sitewarm3,sitewarm4, sitephoto2, sitephoto3,sitephoto4,sitechill22,sitechill23,sitechill24,
             sitechill32,sitechill33,sitechill34,
             warmphoto,warmchill2, warmchill3,chillphoto2,chillphoto3)
  bb <- rnorm(n = length(warm), mean = mm %*% coeff, sd = 1)
  faket<-merge(bb,mm)
  #faket<-data.frame(bb,mm)
  fakex<-data.frame(faket, sp=i)
  fake<-rbind(fake, fakex)
}
str(fake)
#I think it worked! I need to think more about incorporating the trait data, right now every site has a different trait value but all GC treatments have the same value. This makes sense to me, since we might see differnces across the 4 pops, but I don't have the data to know how the traits would change with warming, chilling & photo periods.
length(fake$slav); length(fake$warm2); length(fake$x)
fake$warm2
names(fake)
summary(lm(x~ site2+ site3 + warm2, data= fake))

names(fake)
