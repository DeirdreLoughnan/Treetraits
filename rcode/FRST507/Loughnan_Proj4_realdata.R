### Started March 30 2019 ###

## DL writing final code for analysis of real data for 507 Project IV ##

#The aim of this code is to apply my current model to the real trait and phenology data, testing the relationship between phenology and functional traits
#Update: Added partial pooling across species and multiple PPC

rm(list=ls()) 
options(stringsAsFactors = FALSE)

#For this early model, I will have the following variables:
#Traits: SLA, height (ht), m.dbh, wood density, C:N, & stomatal density

library(rethinking)

setwd("~/Documents/github/Treetraits")

#source('Rcode/Source/Cleaning_compiling_data.R') #Here I combined the budburst and trait data, averaging over individual trees for a given treatment

comb<-read.csv("input/combined_traitpheno.csv")
str(comb)
length(unique(comb$sp))
d<-comb[,c(3,13,14, 16:19)]
d <- d[ complete.cases(d), ] #there are a lot of missing cases for height & many indivudals that have phenology, but not trait data
str(d)
length(unique(d$sp))

#Note, data is not centered, attmepts to center (using both scalle funciton and mannually) resulted in the following error messages:
# 1: In map2stan(alist(bday ~ dnorm(mu, sigma), mu <- a[sp_in] + bstom *  :
# Stripping scale attributes from variable ht

#d_cent<-d
# d_cent$sla<-scale(d_cent$sla)
# d_cent$ht<-scale(d_cent$ht)
# d_cent$stom_d<-scale(d_cent$stom_d)
# d_cent$cn<-scale(d_cent$cn)
# d_cent$wood_den<-scale(d_cent$wood_den)

# d$sla<-(d$sla-mean(d$sla))/(sd(d$sla))
# d$ht<-(d$ht-mean(d$ht))/(sd(d$ht))
# d$stom_d<-(d$stom_d-mean(d$stom_d))/(sd(d$stom_d))
# d$cn<-(d$cn-mean(d$cn))/(sd(d$cn))
# d$wood_den<-(d$wood_den-mean(d$wood_den))/(sd(d$wood_den))

d$sp_in<-coerce_index(d$sp)

d_in<-d[,c(2:8)]
str(d_in)
length(unique(d_in$sp_in))

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#Final model: 5 traits with partial pooling across species
#Table 2: Model output for real data 

full_in <- map2stan(
  alist(
    bday ~ dnorm(mu, sigma) , 
    mu <- a[sp_in]+bstom*stom_d+bsla*sla+bht*ht+bcn*cn+bwood*wood_den,
    a[sp_in] ~ dnorm(mu_a, sigma_a) , 
    mu_a~dnorm(0, 50),
    bstom~dnorm(0, 50),
    bsla~dnorm(0, 50),
    bht~dnorm(0,50),
    bcn~dnorm(0,50),
    bwood~dnorm(0,50),
    sigma ~ dnorm(0,10),
    sigma_a ~ dnorm(0,10)),
  data=d_in, iter=4000 , chains=4
)

precis(full_in, depth=2)

plot(full_in)

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#Figure 1: Model estimates using real phenology and trait data 

par(mfrow=c(1,1))
plot(precis(full_in,depth=2))

####################################################
# Plots: based on code from CH 5 box 5.11, 5.12
####################################################

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#Figure 2: PPC of observed versus predicted days to budburst

par(mfrow=c(1,1))

mu_ppc <- link( full_in)
mu.mean <- apply( mu_ppc , 2 , mean )
mu.PI <- apply( mu_ppc , 2 , PI )

#simulate data
phen.sim <- sim( full_in , n=1000 )
phen.PI <- apply( phen.sim , 2 , PI )

plot( mu.mean ~ d_in$bday , col="darkgreen" , ylim=range(mu.PI) ,
      xlab="Observed phen" , ylab="Predicted phen" )
for ( i in 1:nrow(d) )
  lines( rep(d$bday[i],2) , c(mu.PI[1,i],mu.PI[2,i]) ,
         col="darkgreen" )
abline( a=0 , b=1 , lty=2, lwd=4)
# This only looks fair, most of the data is close to the perfect fit line, but it does seem that the model is underpredicting values, particularly the higher values


# #########################################################################
# #PPC based on code from:
# # https://cran.r-project.org/web/packages/bayesplot/vignettes/graphical-ppcs.html
# #########################################################################
# 
pppost <- extract.samples(full_in)
totaln <- nrow(d_in)
nsite=2
# Make one new dataset from model parameters and means from posterior
newsp <- rnorm(nsp, mean(pppost$mu_a[1:1000]), mean(pppost$sigma_a[1:1000]))
newspdat <- rnorm(totaln, rep(newsp, each=nsite), mean(pppost$sigma[1:1000]))

par(mfrow=c(1,2))
hist(d_in$bday, main="real data")
hist(newspdat, main="PPC")

# Now do it 10 times (still using means from posterior)
newspmat <- matrix(nrow=50, ncol=totaln)
for (i in 1:50){
  newsp <- rnorm(nsp, mean(pppost$mu_a[1:1000]), mean(pppost$sigma_a[1:1000]))
  newspdat <- rnorm(totaln, rep(newsp, each=nsite), mean(pppost$sigma[1:1000]))
  newspmat[i,] <- newspdat
}

par(mfrow=c(2,6))
hist(d_in$bday, main="real data")
for (i in 1:10){
  hist(newspmat[i,], main="PPC")
}

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#Figure 3: Distribution of observed and simulated days to budburst

library(bayesplot)
color_scheme_set("green")
ppc_dens_overlay(d_in$bday, newspmat[1:50, ])

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#Figure 4: Distribution of predicted mean day of budburst per species (light green) compared to observed mean (dark green).
ppc_stat_grouped(d_in$bday, newspmat, binwidth=2,group = d_in$sp_in, stat = "mean")



#########################################################################
# Visualize the partial pooling ...
# See https://www.tjmahr.com/plotting-partial-pooling-in-mixed-effects-models/
#########################################################################

coerce_index(unique(d$sp))
#d_sort <- fake[order(fake$sp_id),]
post <- extract.samples(full_in) # get the posterior
str(post)

unique(d$sp)
nsp=26
J <- nsp 

com.pool.mod <- lm(bday ~ stom_d+sla+ht+cn+wood_den, d_in)
spp <- sort(unique(d$sp_in)); spp


no.pool <- data.frame(sp_in=rep(NA, length(spp)),
                      intercept=rep(NA, length(spp)), 
                      stom_d=rep(NA, length(spp)),
                      ht=rep(NA, length(spp)), 
                      sla=rep(NA, length(spp)),
                      cn=rep(NA, length(spp)), 
                      wood_den=rep(NA, length(spp)),
                      model=rep(NA, length(spp))
)
head(no.pool)

with.pool <- no.pool
df.gravity <- no.pool[1:2,2:7]  
with.pool$model <- "partial pooling"
no.pool$model <- "no pooling"

for (sp in c(1:length(spp))){
  no.pool$sp_in[sp] <- spp[sp]
  subby <- subset(d_in,  sp_in==spp[sp])
  lmfit <- lm(bday~ stom_d+sla+ht+cn+wood_den, data=subby)
  no.pool$intercept[sp] <- coef(lmfit)["(Intercept)"]
  no.pool$stom_d[sp] <- coef(lmfit)["stom_d"]
  no.pool$sla[sp] <- coef(lmfit)["sla"]
  no.pool$ht[sp] <- coef(lmfit)["ht"]
  no.pool$cn[sp] <- coef(lmfit)["cn"]
  no.pool$wood_den[sp] <- coef(lmfit)["wood_den"]
}

# com.pool.mod <- lm(bday ~ stom_d+sla+ht+cn+wood_den, d_in)

spp <- sort(unique(d$sp_in)); spp
no.pool <- data.frame(sp_in=rep(NA, length(spp[3])),
                      intercept=rep(NA, length(spp[3])),
                      stom=rep(NA, length(spp[3])),
                      ht=rep(NA, length(spp[3])),
                      sla=rep(NA, length(spp[3])),
                      cn=rep(NA, length(spp[3])),
                      wood_den=rep(NA, length(spp[3])),
                      model=rep(NA, length(spp[3]))
)
  no.pool$sp_in[3] <- spp[3]
  subby <- subset(d_in,  sp_in==spp[3])
  lmfit <- lm(bday~ stom_d+sla+ht+cn+wood_den, data=subby)
  no.pool$intercept[3] <- coef(lmfit)["(Intercept)"]
  no.pool$stom_d[3] <- coef(lmfit)["stom_d"]
  no.pool$sla[3] <- coef(lmfit)["sla"]
  no.pool$ht[3] <- coef(lmfit)["ht"]
  no.pool$cn[3] <- coef(lmfit)["cn"]
  no.pool$wood_den[3] <- coef(lmfit)["wood_den"]
  64.88713


modhere<-full_in

for (sp in c(1:length(spp))){
  with.pool$sp_in[sp] <- spp[sp]
  with.pool$intercept[spp] <- coef(modhere)[1:26]
  with.pool$stom[spp] <- coef(modhere)["bstom"]
  with.pool$sla[spp] <- coef(modhere)["bsla"]
  with.pool$ht[spp] <-  coef(modhere)["bht"]
  with.pool$cn[spp] <- coef(modhere)["bcn"]
  with.pool$wood[spp] <- coef(modhere)["bwood"]
}

no.pool$sp_no <- sort(unique(d$sp_in))
with.pool$sp_no <- sort(unique(d$sp_in))
df.pulled <- rbind(no.pool, with.pool)
#df.pulled <- bind_rows(no.pool, with.pool)

df.gravity$model <- NA
df.gravity$intercept[1] <-coef(com.pool.mod)["(Intercept)"]
df.gravity$stom[1] <-coef(com.pool.mod)["stom_d"]
df.gravity$sla[1] <-coef(com.pool.mod)["sla"]
df.gravity$ht[1] <-coef(com.pool.mod)["ht"]
df.gravity$cn[1] <-coef(com.pool.mod)["cn"]
df.gravity$wood[1] <-coef(com.pool.mod)["wood_den"]
df.gravity$model[1] <- "complete pooling"

df.gravity$intercept[2] <- coef(modhere)["mu_a"]
df.gravity$stom[2] <- coef(modhere)["bstom"]
df.gravity$sla[2] <- coef(modhere)["bsla"]
df.gravity$htc[2] <- coef(modhere)["bht"]
df.gravity$cnc[2] <- coef(modhere)["bcn"]
df.gravity$woodc[2] <- coef(modhere)["bwood"]
df.gravity$model[2] <- "partial pooling (mu)"

#Plot stomatal density
ggplot(df.pulled) + 
  aes(x = intercept, y = stom, color = model)+ 
  geom_point(size = 3) + 
  geom_point(data = df.gravity, size = 2) + 
  #Draw an arrow connecting the observations between models
  # geom_path(aes(group = as.character(sp_in), color = NULL), 
  #           arrow = arrow(length = unit(.02, "npc"))) + 
  # Use ggrepel to jitter the labels away from the points
  ggrepel::geom_text_repel(
    aes(label = sp_in, color = NULL), 
    data = no.pool, size=2) + 
  theme(legend.position = "bottom") + 
  ggtitle("Pooling of regression parameters: Int and stomatal density") + 
  xlab("Intercept estimate") + 
  ylab("Slope estimate") + 
  scale_color_brewer(palette = "Dark2") 

range(df.pulled$intercept)
#Plot stomatal density
test<-subset(df.gravity, intercept<200 & intercept>-250)
test2<-subset(df.pulled, intercept<200& intercept>-250 )
ggplot(test2) + 
  aes(x = intercept, y = stom, color = model)+ 
  geom_point(size = 2) + 
  geom_point(data = test, size = 5) + 
  # Draw an arrow connecting the observations between models
  #geom_path(aes(group = as.character(sp_in), color = NULL), 
  #arrow = arrow(length = unit(.02, "npc"))) + 
  # Use ggrepel to jitter the labels away from the points
  #ggrepel::geom_text_repel(
  #  aes(label = sp_in, color = NULL), 
  #  data = no.pool, size=2) + 
  theme(legend.position = "bottom") + 
  ggtitle("Pooling of regression parameters: Int and stomatal density") + 
  xlab("Intercept estimate") + 
  ylab("Slope estimate") + 
  scale_color_brewer(palette = "Dark2") 

