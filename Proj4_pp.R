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

full_in <- map2stan(
  alist(
    bday ~ dnorm(mu, sigma) , 
    mu <- a[sp_in]+bstom*stom_d+bsla*sla+bht*Height+bcn*cn+bwood*wood_den,
    a[sp_in] ~ dnorm(mu_a, sigma_a) , # line 65 gives sigma for sp as 5
    mu_a~dnorm(0, 50),
    bstom~dnorm(0, 50),
    bsla~dnorm(0, 50),
    bht~dnorm(0,50),
    bcn~dnorm(0,50),
    bwood~dnorm(0,50),
    sigma ~ dnorm(0,0.1),
    sigma_a ~ dnorm(0,5)),
  data=d_in, iter=4000 , chains=4
)

coerce_index(unique(d$sp))
#d_sort <- fake[order(fake$sp_id),]
post <- extract.samples(full_in) # grab the posterior
str(post)

unique(d$sp)
nsp=26
J <- nsp 

com.pool.mod <- lm(bday ~ stom_d+sla+Height, d_in)
spp <- sort(unique(d$sp_in)); spp
no.pool <- data.frame(sp_in=rep(NA, length(spp)),
                      intercept=rep(NA, length(spp)), 
                      stom=rep(NA, length(spp)),
                      ht=rep(NA, length(spp)), 
                      sla=rep(NA, length(spp)),
                      cn=rep(NA, length(spp)), 
                      wood=rep(NA, length(spp))
)
no.pool

with.pool <- no.pool
df.gravity <- no.pool[1:2,2:4] #why is it 1:2?
with.pool$model <- "partial pooling"
no.pool$model <- "no pooling"

for (sp in c(1:length(spp))){
  no.pool$sp_in[sp] <- spp[sp]
  subby <- subset(d_in,  sp_in==spp[sp])
  lmfit <- lm(bday~ stom_d+sla+Height+cn+wood_den, data=subby)
  no.pool$intercept[sp] <- coef(lmfit)["(Intercept)"]
  no.pool$stom[sp] <- coef(lmfit)["stom_d"]
  no.pool$sla[sp] <- coef(lmfit)["sla"]
  no.pool$ht[sp] <- coef(lmfit)["Height"]
  no.pool$cn[sp] <- coef(lmfit)["cn"]
  no.pool$wood[sp] <- coef(lmfit)["wood_den"]
}

modhere<-full_in
coef(modhere)[1:26]
for (sp in c(1:length(spp))){
  with.pool$sp_in[sp] <- spp[sp]
  with.pool$intercept[spp] <- coef(modhere)[1:26]
  with.pool$stom[spp] <- coef(modhere)["bstom"]
  with.pool$sla[spp] <- coef(modhere)["bsla"]
  with.pool$ht[spp] <-  coef(modhere)["bht"]
  with.pool$cn[spp] <- coef(modhere)["bcn"]
  with.pool$wood[spp] <- coef(modhere)["bwood"]
}

dim(with.pool)
dim(no.pool)
head(no.pool)
head(with.pool)
no.pool$sp_no <- sort(unique(d$sp_in))
with.pool$sp_no <- sort(unique(d$sp_in))
df.pulled <- rbind(no.pool, with.pool)
#df.pulled <- bind_rows(no.pool, with.pool)

df.gravity$model <- NA
df.gravity$intercept[1] <-coef(com.pool.mod)["(Intercept)"]
df.gravity$stom[1] <-coef(com.pool.mod)["stom_d"]
df.gravity$sla[1] <-coef(com.pool.mod)["sla"]
df.gravity$ht[1] <-coef(com.pool.mod)["Height"]
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
  geom_point(size = 2) + 
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
test<-subset(df.gravity, intercept<200 )
test2<-subset(df.pulled, intercept<200 )
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

