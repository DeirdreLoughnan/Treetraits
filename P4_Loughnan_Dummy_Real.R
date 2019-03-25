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
##################################################
# Making test data
###################################################
#Now developing a model with all traits and partial pooling across species

int<-10 # days into the experiment of bb 
sigma<- 0.5
nsite = 2 #number of sites
nsp = 28 #number of species
rep = 50 #
ntot<-nsite*nsp*rep #2800

#Building the required dataframe for all of the replication and traits
#Assigning a variable to site **This will be returned to later when partially pooling across sites
#site=gl(nsite, rep, length=ntot)

#randomly generating numbers for each trait based on what I think a small tree would be
#sitev=rnorm(ntot,5,2) #Ultimatley I will be pooling across sites, so this is not necessary at this point in time

# Trait values (trait name + v for value)

slav=rnorm(ntot, 5, 1)
htv=rnorm(ntot, 11, 3) 
cnv=rnorm(ntot, 10, 2)
woodv=rnorm(ntot, 0.7, 0.05)
stomv=rnorm(ntot, 200, 50)


#effect sizes
site.diff=2 
sladiff= -0.5
htdiff=0.5
cndiff=-0.5
wooddiff=1
stomdiff=1


#####################################################################
# Building up to the full model
#####################################################################

#NOW working with CENTERED DATA
#Centering trait values
stomc=scale(rnorm(rep*nsite, 200, 50))
slac= scale(rnorm(rep*nsite, 5, 1)) 
htc=scale(rnorm(rep*nsite, 11, 3)) 
cnc=scale(rnorm(rep*nsite, 10, 2))
woodc=scale(rnorm(rep*nsite, 0.7, 0.05))

# Here's what I think the loop should be doing.... 
fake_cent<-vector()
sp_mean <- rnorm(nsp, mean=int, sd=5);
sp_mean<-sp_mean[order=T]
for (j in c(1:length(sp_mean))){
  phen_stom<- sp_mean[j]+ rnorm(rep*nsite, 0, sigma)+stomdiff*stomc+sladiff*slac+htdiff*htc+cndiff*cnc+wooddiff*woodc
  temp<-data.frame(phen_stom)
  fake_cent<- rbind(fake_cent, temp) # changed so not overwriting fake_stom every time ... 
}

range(fake_cent)

dim(fake_cent)

#Add species - there is probably a better way to do this
sp=gl(nsp, rep*nsite, length= ntot, order=TRUE)

#sp<-(c(1:28, rep*nsite)); sp

#Combine the phenology data with the species identity 
phen_cent<-cbind(fake_cent,sp)
head(phen_cent)
hist(phen_cent$phen_stom)
#Based on R.code 5.40: this seems to work and keep the order of the intercepts consistent
phen_cent$sp_id<-as.integer(phen_cent$sp)

#This line is taking the mean phen_stom for each species
#checkcode <- aggregate(phen_stom["phen_stom"], phen_stom["sp"], FUN=mean) ; checkcode

###############################################
# Now combining the phenology and the trait data
head(phen_cent)
fake_cent<- data.frame(phen_cent,stomc,slac,htc,cnc,woodc)
head(fake_cent)

#fake_stom$sp_id <- coerce_index(fake_stom$sp)
fake_cent<-fake_cent[,c(1,3:8)]
head(fake_cent)

full_m <- map2stan(
  alist(
    phen_cent ~ dnorm(mu, sigma) , # you have two sigmas in your fake data, should have two here, but coded only one (so I changed)
    mu <- a_sp[sp_id]+bstom*stomc+bsla*slac+bht*htc+bcn*cnc+bwood*woodc, 
    a_sp[sp_id] ~ dnorm(a, sigma_sp) , # line 65 gives sigma for sp as 5 
    a~dnorm(0, 100),
    bstom~dnorm(0, 100),
    bsla~dnorm(0, 100),
    bht~dnorm(0,100),
    bcn~dnorm(0,100),
    bwood~dnorm(0,100),
    sigma ~ dnorm(0,10),
    sigma_sp ~ dnorm(0,5)),
  data=fake_cent, iter=4000 , chains=1 , control = list(max_treedepth = 15)
)

#playing around with different distributions
# full_m <- map2stan(
#   alist(
#     phen_cent ~ dnorm(mu, sigma) , # you have two sigmas in your fake data, should have two here, but coded only one (so I changed)
#     mu <- a_sp[sp_id]+bstom*stomc+bsla*slac+bht*htc+bcn*cnc+bwood*woodc, 
#     a_sp[sp_id] ~ dgamma2(a, sigma_sp) , # line 65 gives sigma for sp as 5 
#     a~dgamma(0, 10),
#     bstom~dgamma(0, 10),
#     bsla~dgamma(0, 10),
#     bht~dgamma(0,10),
#     bcn~dgamma(0,10),
#     bwood~dgamma(0,10),
#     sigma ~ dgamma(0,1),
#     sigma_sp ~ dgamma(0,5)),
  
  a[sp] ~ dnorm( 0,20) ,
      bsla~dnorm(20,10),
      bht~dnorm(15,3),
      bcn~dnorm(20,10),
      bwood~dnorm(2.5,0.6),
      bstom~dnorm(400, 100),# the intercepts all become uniform when this variable is added
      sigma ~ dnorm(0,10))
  data=fake_cent, iter=4000 , chains=4 
)

plot(full_m)

sp_mean
stomdiff
sladiff
htdiff
cndiff
wooddiff
sigma
int

#a and sigma_sp tend to be a little off, but not terribly 
precis(full_m, depth=2)

####################################################################
# working with the real data
####################################################################

source('Rcode/source/Cleaning_compiling_data.R') #Here I combined the budburst and trait data, averaging over individual trees for a given treatment
head(comb)
length(unique(comb$sp))
d<-comb[,c(2,12:13, 15:18)]
d <- d[ complete.cases(d), ] #there are a lot of missing cases for height
length(unique(d$sp))

d$sp_in<-coerce_index(d$sp)

head(d)
d_in<-d[,c(2:8)]
head(d_in)
length(unique(d_in$sp_in))
# range(d$cn, na.rm=TRUE)
# hist(rnorm(1000, 20, 10)) # sla values in data range from 11 to to 54
# hist(rnorm(1000, 15, 3)) # hts range from 0.20 to 20.75
# hist(rnorm(1000, 400, 100)) # stom_d values range from 30 to 770
# hist(rnorm(1000, 2.5, 0.6)) # wood_den range from 0.69 to 4.2
# hist(rnorm(1000, 20, 10)) # cn values in data range from 10.59 to 41.47

unique(d_in$sp_in) 

# full_m <- map2stan(
#   alist(
#     bday ~ dnorm(mu, sigma) , # you have two sigmas in your fake data, should have two here, but coded only one (so I changed)
#     mu <- a_sp[sp]+bstom*stom_d+bsla*sla+bht*Height+bcn*cn+bwood*wood_den, 
#     a_sp[sp] ~ dnorm(a, sigma_sp) , # line 65 gives sigma for sp as 5 
#     a~dnorm(0, 10),
#     bstom~dnorm(0, 10),
#     bsla~dnorm(0, 10),
#     bht~dnorm(0,10),
#     bcn~dnorm(0,10),
#     bwood~dnorm(0,10),
#     sigma ~ dnorm(0,1),
#     sigma_sp ~ dnorm(0,5)),
#   data=d, iter=4000 , chains=4 
# )

#issue had something to do with making the index (w/ 28) and then removing the ones with no heights, left only 26st
full_in <- map2stan(
  alist(
    bday ~ dnorm(mu, sigma) , # you have two sigmas in your fake data, should have two here, but coded only one (so I changed)
    mu <- a[sp_in]+bstom*stom_d+bsla*sla+bht*Height+bcn*cn+bwood*wood_den,
    a[sp_in] ~ dnorm(mu_a, sigma_a) , # line 65 gives sigma for sp as 5
    mu_a~dnorm(0, 10),
    bstom~dnorm(0, 10),
    bsla~dnorm(0, 10),
    bht~dnorm(0,10),
    bcn~dnorm(0,10),
    bwood~dnorm(0,10),
    sigma ~ dnorm(0,1),
    sigma_a ~ dnorm(0,5)),
  data=d_in, iter=4000 , chains=4
)
precis(full_in, depth=2)

plot(full_in)

par(mfrow=c(1,1))
plot(precis(full_in,depth=2))

#<><><><><><><><><><><><><><><><><>
# # the estimates are mucg better EXCEPT for wood density and the intercepts are almost all the same! 
# infprior<- map2stan(
#   alist(
#     phen_spint~ dnorm( mu , sigma ) ,
#     mu<-a[sp]+bsla*slav+bht*htv+bcn*cnv+bwood*woodv
#     +bstom*stomv
#     , 
#     a[sp] ~ dnorm( 0,20) ,
#     bsla~dnorm(20,10),
#     bht~dnorm(15,3),
#     bcn~dnorm(20,10),
#     bwood~dnorm(2.5,0.6),
#     bstom~dnorm(400, 100),# the intercepts all become uniform when this variable is added
#     sigma ~ dnorm(0,10)),
#   data=fake_spint, iter=4000 , chains=4 
# )
# 
# precis(infprior, depth=2)
# plot(precis(infprior, depth=2))
# 

#########################
# Plots: based on code from CH 5 box 5.11, 5.12
#########################

par(mfrow=c(1,1))

mu_ppc <- link( full_m)
mu.mean <- apply( mu_ppc , 2 , mean )
mu.PI <- apply( mu_ppc , 2 , PI )

#simulate data
phen.sim <- sim( full_m , n=1000 )
phen.PI <- apply( phen.sim , 2 , PI )

plot( mu.mean ~ d$bday , col="darkgreen" , ylim=range(mu.PI) ,
      xlab="Observed phen" , ylab="Predicted phen" )
for ( i in 1:nrow(d) )
  lines( rep(d$bday[i],2) , c(mu.PI[1,i],mu.PI[2,i]) ,
         col="darkgreen" )
abline( a=0 , b=1 , lty=2, lwd=4)
# This looks pretty good, the data is fairly aligned with the perfect fit line#

##########################

# compute residuals
phen.resid <- d$bday - mu.mean
# get ordering by divorce rate
ord <- order(phen.resid)
# make the plot
dotchart( phen.resid[ord] , labels=d$sp[ord]  , cex=0.6 )
abline( v=0 , col=col.alpha("black",0.2) )
for ( i in 1:nrow(fake_spint) ) {
  j <- ord[i] # which State in order
  lines( fake_spint$phen_spint[j]-c(mu.PI[1,j],mu.PI[2,j]) , rep(i,2) )
  points( fake_spint$phen_spint[j]-c(phen.PI[1,j],phen.PI[2,j]) , rep(i,2),
          pch=3 , cex=0.6 , col="gray" )
}

####################################################
# Visualize the partial pooling ...
# See https://www.tjmahr.com/plotting-partial-pooling-in-mixed-effects-models/
##
coerce_index(unique(d$sp))
d_sort <- d_in[order(d_in$sp_in),]

post <- extract.samples(full_in) # grab the posterior

nsp=26
J <- nsp #once NA's are removed there are only 26

com_pool_mod <- lm(bday ~ 1, d_sort)
species <- unique(d_sort$sp_id)
no_pool <- data.frame(sp_in=unique(d_sort$sp_in),
                      intercept=rep(NA, J), 
                      slope=rep(NA, J))
head(d_sort)
with_pool <- data.frame(sp_in=unique(d_sort$sp_in),
                        intercept=rep(NA, J), 
                        slope=rep(NA, J))
df.gravity <- no_pool[1:2,2:3]
with_pool$model <- "partial pooling"
no_pool$model <- "no pooling"

d_sort
data<-subset(d_sort, sp_in==unique(d_sort$sp_in)[1])
data

#data=subset(vino.formodel.alt.sort, judge==unique(vino.formodel.alt.sort$judge)[j]))
for (j in 1:J){
  modhere <- lm(bday~1,data=subset(d_sort, sp_in==unique(d_sort$sp_in)[j]))
  no_pool$intercept[j] <- coef(modhere)[1]
  no_pool$slope[j] <- j # cheat for no-slope model
}

for (j in 1:J){
  with_pool$intercept[j] <- mean(post$a[1:1000, j])
  with_pool$slope[j] <- j # cheat for no-slope model
}

# head(no_pool)
# head(with_pool)

df.pulled <- rbind(no_pool, with_pool)

df.gravity$model <- NA
df.gravity$intercept[1] <- coef(com_pool_mod)["(Intercept)"]
df.gravity$slope[1] <- 4.5 # cheat for no-slope model: middle of plot
df.gravity$model[1] <- "complete pooling"
df.gravity$intercept[2] <- mean(post$mu_a[1:1000])
df.gravity$slope[2] <- 4.5 # cheat for no-slope model: middle of plot
df.gravity$model[2] <- "partial pooling (mu)"

library(ggplot2)
#pdf("modelscompare_pp.pdf", width = 9, height = 6)
ggplot(df.pulled) + 
  aes(x = intercept, y = slope, color = model) + 
  geom_point(size = 2) + 
  geom_point(data = df.gravity, size = 5) + 
  # Draw an arrow connecting the observations between models
  geom_path(aes(group = as.character(sp_in), color = NULL), 
            arrow = arrow(length = unit(.02, "npc"))) + 
  # Use ggrepel to jitter the labels away from the points
  ggrepel::geom_text_repel(
    aes(label = sp_in, color = NULL), 
    data = no_pool, size=2) + 
  theme(legend.position = "bottom") + 
  ggtitle("Pooling of regression parameters: Intercept-only") + 
  xlab("Intercept estimate") + 
  ylab("") + 
  scale_color_brewer(palette = "Dark2") 
#dev.off()

######################################################################
######################################################################
# Posterior Predictive Checks:

pppost <- extract.samples(full_in)
totaln <- nrow(d_sort)

# Make one new dataset from model parameters and means from posterior
newsp <- rnorm(nsp, mean(pppost$mu_a[1:1000]), mean(pppost$sigma_a[1:1000]))
newspdat <- rnorm(totaln, rep(newsp, each=nsite), mean(pppost$sigma[1:1000]))

par(mfrow=c(1,2))
hist(d_sort$bday, main="real data")
hist(newspdat, main="PPC")

# Now do it 10 times (still using means from posterior)
newspmat <- matrix(nrow=10, ncol=totaln)
for (i in 1:10){
  newsp <- rnorm(nsp, mean(pppost$mu_a[1:1000]), mean(pppost$sigma_a[1:1000]))
  newspdat <- rnorm(totaln, rep(newsp, each=nsite), mean(pppost$sigma[1:1000]))
  newspmat[i,] <- newspdat
}

par(mfrow=c(2,6))
hist(d_sort$bday, main="real data")
for (i in 1:10){
  hist(newspmat[i,], main="PPC")
}

# Similarly ...
library(bayesplot)
color_scheme_set("brightblue")
ppc_dens_overlay(d_sort$bday, newspmat[1:5, ])


