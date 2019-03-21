### Started January 29 2019 ###

## DL writing dummy data for 507 Project ##

#The aim of this code is to build test data for modeling the relationship between phenology and functional traits
#Update: Added partial pooling across species & also changed data to not exclude sites with no phenology data

rm(list=ls()) 
options(stringsAsFactors = FALSE)

#For this early model, I will have the following variables:
#Traits: SLA, height (ht), m.dbh, wood density, C:N, & stomatal density
#Note: After looking for collinearity in the x varibales, I decided to exclude DBH and just use ht

setwd("~/Documents/github/Treetraits")
source('Rcode/source/Cleaning_compiling_data.R') #Here I combined the budburst and trait data, averaging over individual trees for a given treatment. -->stopped giving me the final dataset, unclear why

str(comb)
length(unique(comb$sp))

cor<-comb[,10:15]
pairs(cor)
#Starting by looking for high correlations between traits

# Correlation panel
#there are many ways to do this, but this provides a busy but I think useful visual
# library("PerformanceAnalytics")
# chart.Correlation(comb[,c(10:15)], histogram=FALSE, pch=19)

###################################################
# Making test data
###################################################
#To start with a simple model, I am going to develop a linear model with just the five traits


int<-11 # days into the experiment of bb 
sigma<-2

nsite = 2 #number of sites
nsp = 28 #number of species
rep = 50 # I think the maximum number is 10, but to increase the size of the dataset I have made it 50
ntot<-nsite*nsp*rep #2800

int_sp<-rnorm(nsp, mean=int, sd=sigma)
#Building the required dataframe for all of the replication and traits
#Assigning a variable to site **This will be returned to later when partially pooling across sites
site=gl(nsite, rep, length=ntot)

#randomly generating numbers for each trait based on what I think a small tree would be
#sitev=rnorm(ntot,5,2) #Ultimatley I will be pooling across sites, so this is not necessary at this point in time

# Trait values (trait name + v for value)
slav= rnorm(ntot, 5, 1); 
htv=rnorm(ntot, 11, 3) #sigma should be really high (3m) bc it includes trees and woody shrubs
cnv=rnorm(ntot, 10, 2)
woodv=rnorm(ntot, 0.5, 0.05)
stomv=rnorm(ntot, 20, 5)

#effect sizes
site.diff=2 # I think site 2 will be delayed by 2 days due to the 5 degree diff in lat
sladiff= -0.5
htdiff=0.5
cndiff=-0.5
wooddiff=0.3
stomdiff=1

#making a matrix of the x dummy variables
 mm <- model.matrix(~(site+slav+htv+cnv+woodv+stomv), data.frame(site,slav,htv,cnv,woodv,stomv))
 mm

#making a simple data frame of the x dummy variables
df<-data.frame(site,slav,htv,cnv,woodv,stomv)

################################################
#To resolve initial issues of negative values and unusual predictions, the below code was used to test variables individually

#phen<-vector()
# #stomatal density on its own
# for (i in 1:ntot){
#   phen[i]<-int+stomdiff*stomv[i]+rnorm(1,0, sigma)
# }
# phen
# range(phen)
# #
# #cn on its own
# for (i in 1:ntot){
#   phen[i]<-int+wooddiff*woodv[i]+rnorm(1,0, sigma)
# }
# range(phen)

#cn on its own
# for (i in 1:ntot){
#   phen[i]<-int+cndiff*cnv[i]+rnorm(1,0, sigma)
#   
# }
# range(phen)
# 
# #SLA on its own
# for (i in 1:ntot){
#   phen[i]<-int+sladiff*slav[i]+rnorm(1,0, sigma)
# }
# range(phen)
# summary(lm(phen~slav))

#ht on its own
# for (i in 1:ntot){
#   phen[i]<-int+htdiff*htv[i]+rnorm(1,0, sigma)
# }
# # phen
# range(phen)
#
#Finally we can use all variables to produce y dummy data
phenfull<-vector()
for (i in 1:ntot){
  phenfull[i]<-int+sladiff*slav[i]+htdiff*htv[i]+cndiff*cnv[i]+wooddiff*woodv[i]+stomdiff*stomv[i]+rnorm(1, 0, sigma)
}
range(phenfull) #great, no negative values

summary(lm(phenfull~slav+htv+cnv+woodv+stomv))
#These values seem to be as good as I can get without choosing values I think are not ecologically relevant 

par(mfrow=c(1,1))
require(rethinking)
simplehist(phenfull, xlab="Day of budburst")


#<><><><><><><><><><><><><><><><><><>
# Now adding separate intercepts for each species
# Now adding separate intercepts for each species
# For centered data
baseinter=11
sp_int<-baseinter+c(1:nsp)-mean(1:nsp)
sp_int

fake_spint<-vector()
sp_mean <- rnorm(nsp, 11, 3)
for (j in c(1:length(sp_mean))){
  phen_spint<- int+sladiff*slav+htdiff*htv+cndiff*cnv+wooddiff*woodv+stomdiff*stomv+rnorm(rep*nsite, sp_mean[j], sigma)
  temp<-data.frame(phen_spint)
  phen_spint_cent<- rbind(fake_spint, temp)
}

#Add species 
sp=gl(nsp, rep*nsite, length= ntot); species

head(phen_spint)
range(phen_spint)

#Combine the phenology data with the species identity 
phen_spint<-cbind(phen_spint,sp)
head(phen_spint)
###############################################
# Creating centered trait data:
# centering the x dummy data
mm.cent<-scale(mm[,3:7]);head(mm.cent)
#pairs(mm.cent) # These plots are just clouds of points and not helpful
mm.cent<-cbind(mm[,1:2], mm.cent)
head(mm.cent)

##############################################
#Calculating the y dummy data using  matrix multiplication: I will use an alternative method for now. 

# coeff <- c(1, sladiff, htdiff, cndiff, wooddiff, stomdiff)
# 
# #centered
# phencent <- rnorm(n = length(site), mean = mm.cent %*% coeff, sd = 1) # This code works but the values are HUGE
# phencent
# 
# #uncentered
# phen <- rnorm(n = length(site), mean = mm %*% coeff, sd = 1) # This code works but the values are HUGE
# range(phen) #why are some of these values negative?
##############################################

#Creating centered data
slav<-scale(slav)
htv<-scale(htv)
cnv<-scale(cnv)
woodv<-scale(woodv)
stomv<-scale(stomv)

#Testing whether the issue is just stomatal density:
phencent<-vector()
for (i in 1:ntot){
  phencent[i]<-int+sladiff*slav[i]+htdiff*htv[i]+cndiff*cnv[i]+wooddiff*woodv[i]+stomdiff*stomv[i]+rnorm(1,0, sigma)
}
head(phencent)

#hist of sample distribution
simplehist(phencent, xlab="Budburst daty")

#<><><><><><><><><><><><><><><><><><>
# Now adding separate intercepts for each species
# For centered data
baseinter=11
sp_int<-baseinter+c(1:nsp)-mean(1:nsp)
sp_int

fake_spint<-vector()
sp_mean <- rnorm(nsp, 11, 3)
for (j in c(1:length(sp_mean))){
  phen_spint_cent<- int+sladiff*slav+htdiff*htv+cndiff*cnv+wooddiff*woodv+stomdiff*stomv+rnorm(rep*nsite, sp_mean[j], sigma)
  temp<-data.frame(phen_spint_cent
                   #, sp=c(1:length(sp_mean)
                   )
  phen_spint_cent<- rbind(fake_spint, temp)
}

#Add species 
sp=gl(nsp, rep*nsite, length= ntot); species

head(phen_spint_cent)
range(phen_spint_cent)

#Combine the phenology data with the species identity 
phen_spint_cent<-cbind(phen_spint_cent,sp)
head(phen_spint_cent)

####################################################################################
#Now creating the combined datasets of the fake data for use with map
fake<- data.frame(phenfull,mm)
head(fake)

fake_spint<- data.frame(phen_spint,mm)
head(fake_spint)

fakecent<- data.frame(phencent,mm.cent)
head(fakecent)

fake_spint_cent<- data.frame(phen_spint_cent,mm.cent)
head(fake_spint_cent)

#Note Stan doesn't like dots so changing X.Intercept name
colnames(fake)[colnames(fake)=="X.Intercept."] <- "Intercept"
colnames(fake_spint)[colnames(fake_spint)=="X.Intercept."] <- "Intercept"
colnames(fakecent)[colnames(fakecent)=="X.Intercept."] <- "Intercept"
colnames(fake_spint_cent)[colnames(fake_spint_cent)=="X.Intercept."] <- "Intercept"


# summary(lm(phenfull ~ (slav+htv+cnv+woodv+stomv), data = fake)) 
# 
# summary(lm(phencent ~ (slav+htv+cnv+woodv+stomv), data = fakecent))

####### Uncentered data #############
# phen <- rnorm(n = length(site), mean = mm %*% coeff, sd = 1) # This code works but the values are HUGE 
# phen

#right now I am not including site in the model, so the df needs to be trimmed
#what role does the intercept column play?
#fake_spint_cent<-fake_spint_cent[,c(1:3,5:9)]; head(fake_spint_cent) #includes intercept
fake_spint_cent<-fake_spint_cent[,c(1:2,5:9)]; head(fake_spint_cent) #does NOT include intercept

###### MAP2Stan MODEL ###########
#Simpiler model with just sla

#Starting from the very basics witht his new dataset 

int.m<- map2stan(
  alist(
    phencent~dnorm(mu, sigma),
    mu<-a[sp],
    a[sp]~dnorm(0, 100),
    sigma~dunif(0,1)
  ),
  data=fake_spint_cent, chains=4)


precis(int.m)

plot(int.m)
#It makes sense that this model has issues, the data was generated with individual intercepts in mind, but I did not expect it to be this bad
#has 56 divergent transitions among other issues

# <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <> #
# Model including individual species intercepts
int.m<- map2stan(
   alist(
    phencent~dnorm(mu, sigma),
    mu<-a[sp],
    a[sp]~dnorm(0, 100),
    sigma~dunif(0,100)
       ),
  data=fake_spint_cent)

int.m<- map2stan(
  alist(
    phencent~dnorm(mu, sigma),
    mu<-a+a_sp[sp],
    a_sp[sp]~dnorm(a,sigma),
    a~dnorm(0, 100),
    sigma~dunif(0,100)
  ),
  data=fake_spint_cent,
  control = list(max_treedepth = 15))

precis(int.m, depth=2)

plot(int.m)

head(fake_spint_cent)

slam <- map2stan(
  alist(
    phencent~ dnorm( mu , sigma ) ,
    mu<-a[sp]+bsla*slav,
    a[sp] ~ dnorm( 0,10) ,
    bsla~dnorm(0, 10),
    sigma ~ dunif(0,1)
  ), data=fake_spint_cent #, iter=4000 , chains=4 
  )

precis(slam)
plot(slam)
#plot the posterior
plot( phencent ~ slav , data=fakecent)
abline( a=coef(sla.m)["a"] , b=coef(sla.m)["bsla"] )

#To get the uncertainty plot many of these lines
post <- extract.samples( sla.m, n=1000)

# display raw data and sample size
plot( fake$phenfull , fake$slav ,
      xlim=range(fake$phenfull) , ylim=range(fake$slav) ,
      col=rangi2 , xlab="SLA" , ylab="Budburst day" )

# plot the lines, with transparency
for ( i in 1:20 )
  abline( a=post$a[i] , b=post$b[i] , col=col.alpha("black",0.3) )
#######################################################
#Now trying adding all variables
full_m<- map(
  alist(
    phencent~dnorm(mu, sigmahere),
    mu<-a+bsla*slav+bht*htv+bcn*cnv+bstom*stomv+bwood*woodv,
    a~dnorm(0, 10),
    bsla~dnorm(0, 10),
    bht~dnorm(0, 10),
    bcn~dnorm(0, 10),
    bstom~dnorm(0, 10),
    bwood~dnorm(0, 10),
    sigmahere~dunif(0,10)
  ),
  data=fakecent)
precis(full_m)


#Very small intervals when you plot the MAP values
plot(precis(full_m))

#Calculating the variance covariance matrix
vcmatrix<-vcov(full_m)
round(vcmatrix, digits=6)
#These values seem very very small

#Calculating the correlation matrix, again getting very very small values
cormatrix<-cov2cor(vcov(full_m))
round(cormatrix, digits=6)

precis( full_m , corr=TRUE )


head(fakecent)
#Now using map2stan
full_stan<- map2stan(
  alist(
    phencent~dnorm(mu, sigmahere),
    mu<-a+bsla*slav+bht*htv+bcn*cnv+bstom*stomv+bwood*woodv,
    a~dnorm(0, 10),
    bsla~dnorm(0, 10),
    bht~dnorm(0, 10),
    bcn~dnorm(0, 10),
    bstom~dnorm(0, 10),
    bwood~dnorm(0, 10),
    sigmahere~dunif(0,10)
  ),
  data=fakecent)
plot(full_stan)
precis(full_stan)

plot(full_stan)

post <- extract.samples( full_stan )
str(post)

pairs(post)
pairs(full_stan)

################################################################
#Testing whether the issue is just stomatal density:
nstom<-vector()
for (i in 1:ntot){
  nstom[i]<-int+sladiff*slav[i]+htdiff*htv[i]+cndiff*cnv[i]+wooddiff*woodv[i]+rnorm(1,0, sigma)
}

#Now creating the combined datasets of the fake data for use with map
fake_nstom<- data.frame(nstom,mm.cent)
head(fake_nstom)
#Note Stan doesn't like dots so changing X.Intercept name
colnames(fake_nstom)[colnames(fake_nstom)=="X.Intercept."] <- "Intercept"
head(fake_nstom)

full_nstom<- map(
  alist(
    phencent~dnorm(mu, sigmahere),
    mu<-a+bsla*slav+bht*htv+bcn*cnv+bwood*woodv,
    a~dnorm(0, 10),
    bsla~dnorm(0, 10),
    bht~dnorm(0, 10),
    bcn~dnorm(0, 10),
    bstom~dnorm(0, 10),
    sigmahere~dunif(0,10)
  ),
  data=fake_nstom)
precis(full_nstom)

#Testing whether the issue is just wood density: 
nwood<-vector()
for (i in 1:ntot){
  nwood[i]<-int+sladiff*slav[i]+htdiff*htv[i]+cndiff*cnv[i]+stomdiff*stomv[i]+rnorm(1,0, sigma)
}

#Now creating the combined datasets of the fake data for use with map
fake_nwood<- data.frame(nwood,mm.cent)
head(fake_nwood)
#Note Stan doesn't like dots so changing X.Intercept name
colnames(fake_nwood)[colnames(fake_nwood)=="X.Intercept."] <- "Intercept"
head(fake_nwood)

full_nwood<- map(
  alist(
    phencent~dnorm(mu, sigmahere),
    mu<-a+bsla*slav+bht*htv+bcn*cnv+bstom*stomv,
    a~dnorm(0, 10),
    bsla~dnorm(0, 10),
    bht~dnorm(0, 10),
    bcn~dnorm(0, 10),
    bstom~dnorm(0, 10),
    sigmahere~dunif(0,10)
  ),
  data=fake_nwood)
precis(full_nwood)
##########################################
# PLOTS
##########################################

par(mfrow=c(2,3))
dens(slav)
dens(htv)
dens(cnv)
dens(woodv)
dens(stomv)

par(mfrow=c(2,3))
plot(phencent~slav, data=fakecent)
plot(phencent~htv, data=fakecent)
plot(phencent~cnv, data=fakecent)
plot(phencent~woodv, data=fakecent)
plot(phencent~stomv, data=fakecent)

#looking at the data
#calculate the mean values for line calculations
msla<-mean(fakecent$slav)
mstom<-mean(fakecent$stomv)
mht<-mean(fakecent$htv)
mcn<-mean(fakecent$cnv)
mwood<-mean(fakecent$woodv)

#assign variable names to coefficients * probably not the most efficient way to do this
test<-precis(full_m)
a=coef(full_m)["a"] ; a
bsla=coef(full_m)["bsla"]; bsla
bht=coef(full_m)["bht"]; bht
bcn=coef(full_m)["bcn"]; bcn
bstom=coef(full_m)["bstom"]; bsla
bwood=coef(full_m)["bwood"]; bwood

range(fakecent$phencent)
y<-seq(-5,5,0.1);length(y)

par(mfrow=c(2,3))
plot(phencent~slav, data=fakecent, ylim=
       range(fakecent$phencent), xlab="SLA", ylab="Day of budburst")
lines(y,(a+bsla*y+bht*mht+bcn*mcn+bstom*mstom+bwood*mwood), col="red", lwd=3, lty=1)


plot(phencent~htv, data=fakecent, ylim=
       range(fakecent$phencent), xlab="Height", ylab="Day of budburst")
  lines(y,(a+bsla*msla+bht*y+bcn*mcn+bstom*mstom+bwood*mwood), col="red", lwd=3, lty=1)
  
plot(phencent~cnv, data=fakecent, ylim=
       range(fakecent$phencent), xlab="C:N", ylab="Day of budburst")
  lines(y,(a+bsla*msla+bht*mht+bcn*y+bstom*mstom+bwood*mwood), col="red", lwd=3, lty=1)

plot(phencent~stomv, data=fakecent, ylim=
       range(fakecent$phencent), xlab="Stomatal density", ylab="Day of budburst")
lines(y,(a+bsla*msla+bht*mht+bcn*mcn+bstom*y+bwood*mwood), col="red", lwd=3, lty=1)

plot(phencent~woodv, data=fakecent, ylim=
       range(fakecent$phencent), xlab="Wood Density", ylab="Day of budburst")
  lines(y,(a+bsla*msla+bht*mht+bcn*mcn+bstom*mstom*+bwood*y), col="red", lwd=3, lty=1)
  
###########################################################
# counterfactual plot for sla

  A.avg <- mean( fakecent$slav)
  R.seq <- seq( from=-5 , to=5 , length.out=2800 )
  pred.data <- data.frame(
    htv=R.seq,
    wood=R.seq,
    stomv=R.seq,
    cnv=R.seq,
    slav=A.avg
  )

mu <- link( full_m , data=fakecent )
mu.mean <- apply( mu , 2 , mean ); length(mu.mean)
mu.PI <- apply( mu , 2 , PI )
A.sim <- sim( full_m , data=fakecent , n=1e4 )
A.PI <- apply( A.sim , 2 , PI )
plot( phencent ~ slav , data=fakecent , type="n" )

length(R.seq); length(mu.mean)
lines( R.seq , mu.mean )
shade( mu.PI , R.seq )
shade( A.PI , R.seq )

#Counterfactual plot for ht
# 

A.avg <- mean( fakecent$htv)
R.seq <- seq( from=-5 , to=5 , length.out=2800 )
pred.data <- data.frame(
  htv=A.avg,
  wood=R.seq,
  stomv=R.seq,
  cnv=R.seq,
  slav=R.seq
)

mu <- link( full_m , data=fakecent )
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI )
A.sim <- sim( full_m , data=fakecent , n=1e4 )
A.PI <- apply( A.sim , 2 , PI )
plot( phencent ~ htv , data=fakecent , type="n" )

length(A.seq); length(mu.mean)
lines( R.seq , mu.mean )
shade( mu.PI , R.seq )
shade( A.PI , R.seq )

### ** These plots look the same and both look like fuzzy caterpillars **
#######################################################
#PPC as per r code 5.11 & 5.12
par(mfrow=c(1,1))

mu_ppc <- link( full_m)
mu.mean <- apply( mu_ppc , 2 , mean )
mu.PI <- apply( mu_ppc , 2 , PI )

#simulate data
phen.sim <- sim( full_m , n=1e4 )
phen.PI <- apply( phen.sim , 2 , PI )

plot( mu.mean ~ fakecent$phencent , col="darkgreen" , ylim=range(mu.PI) ,
      xlab="Observed phen" , ylab="Predicted phen" )
abline( a=0 , b=1 , lty=2, lwd=4)
for ( i in 1:nrow(fakecent) )
  lines( rep(fakecent$phencent[i],2) , c(mu.PI[1,i],mu.PI[2,i]) ,
         col="darkgreen" )

###############

# compute residuals
phen.resid <- fakecent$phencent - mu.mean
# get ordering by divorce rate
ord <- order(phen.resid)
# make the plot
dotchart( phen.resid[ord] , labels=d$Loc[ord]  , cex=0.6 )
abline( v=0 , col=col.alpha("black",0.2) )
for ( i in 1:nrow(fakecent) ) {
  j <- ord[i] # which State in order
  lines( fakecent$phencent[j]-c(mu.PI[1,j],mu.PI[2,j]) , rep(i,2) )
  points( fakecent$phencent[j]-c(phen.PI[1,j],phen.PI[2,j]) , rep(i,2),
          pch=3 , cex=0.6 , col="gray" )
}



