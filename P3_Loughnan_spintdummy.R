### Started January 29 2019 ###

## DL writing dummy data for 507 Project ##

#The aim of this code is to build test data for modeling the relationship between phenology and functional traits
#Update: Added partial pooling across species & also changed data to not exclude sites with no phenology data

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

int<-11 # days into the experiment of bb 
sigma<-2

nsite = 2 #number of sites
nsp = 28 #number of species
rep = 50 # I think the maximum number is 10, but to increase the size of the dataset I have made it 50
ntot<-nsite*nsp*rep #2800

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

#####################################################################
# Now adding separate intercepts for each species ** Starting very simple with just SLA

fake_spint<-vector()
sp_mean <- rnorm(nsp, mean=int, sd=sigma)
for (j in c(1:length(sp_mean))){
  # note that in the below you model the intercept twice ... given that you have int
  # at the start I would expect your species effects to be offsets (centered around zero)
  #I interpret this and your comment below to mean I should either remove the int at the beginning   # or that I should have another a parameter below to account for the offset
  #phen_spint<- int+sladiff*slav+rnorm(rep*nsite, sp_mean[j], sigma)
  phen_spint<- sladiff*slav+htdiff*htv+cndiff*cnv+wooddiff*woodv+stomdiff*stomv+rnorm(rep*nsite, sp_mean[j], sigma)
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

###### MAP2Stan MODEL ###########
#Simpiler model with just sla

#Starting from the very basics witht his new dataset 

int.m<- map2stan(
  alist(
    phen_spint~dnorm(mu, sigma), # changed your response variable ... did not see it in your dataframe
    mu<-a[sp],
    a[sp]~dnorm(0, 20),
    sigma~dnorm(0,10) # plot hist(runif(1000,0,1)) to see what dunif(0,1) prior looks like
  ),
  data=fake_spint, chains=4)

#hist(runif(1000,0,1))
precis(int.m, depth=2)

plot(int.m)
#Wow this model is bad

plot(precis(int.m, depth=2))
# <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <>  <> #
# Model including individual species intercepts

slam <- map2stan(
  alist(
    phen_spint~ dnorm( mu , sigma ) ,
    mu<-a[sp]+bsla*slav, 
    a[sp] ~ dnorm( 0,20) ,
    bsla~dnorm(0, 10),
    sigma ~ dnorm(0,10)),
  data=fake_spint, iter=4000 , chains=4 
  )

precis(slam, depth=2)

plot(slam)

plot(precis(slam, depth=2))

#plot the posterior
par(mfrow=c(1,1))
plot( phen_spint ~ slav , data=fake_spint)
#abline( a=coef(slam)["a"] , b=coef(slam)["bsla"] )

#To get the uncertainty plot many of these lines
post <- extract.samples( slam, n=1000)

# display raw data and sample size
plot( fake_spint$phen_spint , fake_spint$slav ,
      xlim=range(fake_spint$phen_spint) , ylim=range(fake_spint$slav) ,
      col=rangi2 , xlab="SLA" , ylab="Budburst day" )

# plot the lines, with transparency
for ( i in 1:20 )
  abline( a=post$a[i] , b=post$b[i] , col=col.alpha("black",0.3) )

#########################
# Plots
#########################
par(mfrow=c(1,1))

mu_ppc <- link( slam)
mu.mean <- apply( mu_ppc , 2 , mean )
mu.PI <- apply( mu_ppc , 2 , PI )

#simulate data
phen.sim <- sim( slam , n=1000 )
phen.PI <- apply( phen.sim , 2 , PI )

plot( mu.mean ~ fake_spint$phen_spint , col="darkgreen" , ylim=range(mu.PI) ,
      xlab="Observed phen" , ylab="Predicted phen" )
abline( a=0 , b=1 , lty=2, lwd=4)
for ( i in 1:nrow(fake_spint) )
  lines( rep(fake_spint$phen_spint[i],2) , c(mu.PI[1,i],mu.PI[2,i]) ,
         col="darkgreen" )

###############

# compute residuals
phen.resid <- fake_spint$phen_spint - mu.mean
# get ordering by divorce rate
ord <- order(phen.resid)
# make the plot
dotchart( phen.resid[ord] , labels=fake_spint$sp[ord]  , cex=0.6 )
abline( v=0 , col=col.alpha("black",0.2) )
for ( i in 1:nrow(fake_spint) ) {
  j <- ord[i] # which State in order
  lines( fake_spint$phen_spint[j]-c(mu.PI[1,j],mu.PI[2,j]) , rep(i,2) )
  points( fake_spint$phen_spint[j]-c(phen.PI[1,j],phen.PI[2,j]) , rep(i,2),
          pch=3 , cex=0.6 , col="gray" )
}


# Counterfactual plot
sla.avg <- mean( fake_spint$slav)
R.seq <- seq( from=-5 , to=5 , length.out=2800 )
pred.data <- data.frame(
  htv=R.seq,
  wood=R.seq,
  stomv=R.seq,
  cnv=R.seq,
  slav=A.avg
)

mu <- link( slam , data=fakecent )
mu.mean <- apply( mu , 2 , mean ); length(mu.mean)
mu.PI <- apply( mu , 2 , PI )
A.sim <- sim( full_m
              , data=fakecent , n=1e4 )
A.PI <- apply( A.sim , 2 , PI )
plot( phencent ~ slav , data=fakecent , type="n" )

length(R.seq); length(mu.mean)
lines( R.seq , mu.mean )
shade( mu.PI , R.seq )
shade( A.PI , R.seq )

##########################################
#FULL model
##########################################
full<- map2stan(
  alist(
    phen_spint~ dnorm( mu , sigma ) ,
    mu<-a[sp]+bsla*slav+bht*htv+bcn*cnv+bstom*stomv+bwood*woodv, 
    a[sp] ~ dnorm( 0,20) ,
    bsla~dnorm(0, 10),
    bht~dnorm(0, 50),
    bcn~dnorm(0, 15),
    bstom~dnorm(0, 50),
    bwood~dnorm(0, 10),
    sigma ~ dnorm(0,10)),
  data=fake_spint, iter=4000 , chains=4 
)

precis(full, depth=2)

plot(full)

par(mfrow=c(1,1))
plot(precis(full,depth=2))

#########################
# Plots: based on code from CH 5 box 5.11, 5.12
#########################

par(mfrow=c(1,1))

mu_ppc <- link( full)
mu.mean <- apply( mu_ppc , 2 , mean )
mu.PI <- apply( mu_ppc , 2 , PI )

#simulate data
phen.sim <- sim( full , n=1000 )
phen.PI <- apply( phen.sim , 2 , PI )

plot( mu.mean ~ fake_spint$phen_spint , col="darkgreen" , ylim=range(mu.PI) ,
      xlab="Observed phen" , ylab="Predicted phen" )
abline( a=0 , b=1 , lty=2, lwd=4)
for ( i in 1:nrow(fake_spint) )
  lines( rep(fake_spint$phen_spint[i],2) , c(mu.PI[1,i],mu.PI[2,i]) ,
         col="darkgreen" )
# This looks pretty good, the data is fairly aligned with the perfect fit line#

##########################

# compute residuals
phen.resid <- fake_spint$phen_spint - mu.mean
# get ordering by divorce rate
ord <- order(phen.resid)
# make the plot
dotchart( phen.resid[ord] , labels=fake_spint$sp[ord]  , cex=0.6 )
abline( v=0 , col=col.alpha("black",0.2) )
for ( i in 1:nrow(fake_spint) ) {
  j <- ord[i] # which State in order
  lines( fake_spint$phen_spint[j]-c(mu.PI[1,j],mu.PI[2,j]) , rep(i,2) )
  points( fake_spint$phen_spint[j]-c(phen.PI[1,j],phen.PI[2,j]) , rep(i,2),
          pch=3 , cex=0.6 , col="gray" )
}

##########################

# Counterfactual plot
sla.avg <- mean( fake_spint$slav)
R.seq <- seq( from=-10 , to=10 , length.out=2800 )
pred.data <- data.frame(
  htv=R.seq,
  wood=R.seq,
  stomv=R.seq,
  cnv=R.seq,
  slav=sla.avg
)

mu <- link(full , data=fake_spint )
mu.mean <- apply( mu , 2 , mean ); length(mu.mean)
mu.PI <- apply( mu , 2 , PI )
A.sim <- sim( full
              , data=fake_spint , n=1000 )
A.PI <- apply( A.sim , 2 , PI )
plot( phen_spint ~ slav , data=fake_spint , type="n" )

length(R.seq); length(mu.mean)
lines( R.seq , mu.mean )
shade( mu.PI , R.seq )
shade( A.PI , R.seq )
