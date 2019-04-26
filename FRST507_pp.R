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

d$sp_in<-coerce_index(d$sp)

d_in<-d[,c(2:8)]
str(d_in)
length(unique(d_in$sp_in))
#write.csv(d_in, "frst507_standata.csv")
##################################################
N<-as.numeric(nrow(d))
N_sp_in<-26
bday<-d_in$bday
ht<-d_in$ht
stom_d<-d_in$stom_d
sla<-d_in$sla
wood_den<-d_in$wood_den
cn<-d_in$cn
sp_in<-d_in$sp_in

#md = stan('507m_stan.stan', data = c("N","N_sp_in","bday","ht","stom_d", "sla","wood_den","cn","sp_in"),iter = 2500, warmup=1500,control = list(adapt_delta = 0.99))

save(md, file="md_output.Rda")


##################################################
load("md_output.Rda")

sumer.md <- summary(full_in)$summary
sumer.md[grep("a", rownames(sumer.md)),] #calls the mu values, but why?


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

com.pool.mod <- lm(bday ~ stom_d+sla+ht+cn+wood_den, d_in)

######################################################################
######################################################################
spp <- sort(unique(d$sp_in)); spp[3]
no.pool <- data.frame(sp_in=rep(NA, length(spp[3])),
                      intercept=rep(NA, length(spp[3])),
                      stom=rep(NA, length(spp[3])),
                      ht=rep(NA, length(spp[3])),
                      sla=rep(NA, length(spp[3])),
                      cn=rep(NA, length(spp[3])),
                      wood_den=rep(NA, length(spp[3])),
                      model=rep(NA, length(spp[3]))
)
 no.pool$sp_in <- spp[1]
 subby <- subset(d_in,  sp_in==spp[1])
 lmfit <- lm(bday~ ht+stom_d+sla+cn+wood_den, data=subby);lmfit
no.pool$intercept[14] <- coef(lmfit)["(Intercept)"]
no.pool$stom_d[14] <- coef(lmfit)["stom_d"]
no.pool$sla[14] <- coef(lmfit)["sla"]
no.pool$ht[14] <- coef(lmfit)["ht"]
no.pool$cn[14] <- coef(lmfit)["cn"]
no.pool$wood_den[14] <- coef(lmfit)["wood_den"]
64.88713

######################################################################
######################################################################
modhere<-sumer.md

for (sp in c(1:length(spp))){
  with.pool$sp_in[sp] <- spp[sp]
  with.pool$intercept[sp] <- coef(modhere)[1:26]
  with.pool$stom[sp] <- coef(modhere)["bstom"]
  with.pool$sla[sp] <- coef(modhere)["bsla"]
  with.pool$ht[sp] <-  coef(modhere)["bht"]
  with.pool$cn[sp] <- coef(modhere)["bcn"]
  with.pool$wood[sp] <- coef(modhere)["bwood"]
}

# for (sp in c(1:length(spp))){
#   with.pool$sp_in[sp] <- spp[sp]
#   with.pool$intercept[sp] <- modhere[grep("mu_a", rownames(modhere)),1][spp[sp]+2]
#   with.pool$stom[sp] <- modhere[grep("bstom", rownames(modhere)),1][spp[sp]+2]
#   with.pool$sla[sp] <- modhere[grep("bsla", rownames(modhere)),1][spp[sp]+2]
#   with.pool$ht[sp] <-  modhere[grep("bht", rownames(modhere)),1][spp[sp]+2]
#   with.pool$cn[sp] <- modhere[grep("bcn", rownames(modhere)),1][spp[sp]+2]
#   with.pool$wood[sp] <- modhere[grep("bwood", rownames(modhere)),1][spp[sp]+2]
# }

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
################################################################
################################################################
################################################################
#Adapting the code from the tutorial to use for my own model
rm(list=ls()) 
options(stringsAsFactors = FALSE)

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

d$sp_in<-coerce_index(d$sp)

d_in<-d[,c(2:8)]
str(d_in)
length(unique(d_in$sp_in))


library(lme4)
library(dplyr)
library(tibble)

#based on tutorial code
#no pooling 
head(d)

df_no_pooling<-lmList(bday ~ (sla+ht)|sp_in,data=d) %>% 
  coef()%>%
  add_column(Model = "No pooling") %>% 
  rownames_to_column("sp_in")


df_no_pooling

head(df_no_pooling)
colnames(df_no_pooling)[colnames(df_no_pooling)=="(Intercept)"] <- "Intercept"

#complete pooling
m_pooled <- lm(bday ~ sla, d) 

# Repeat the intercept and slope terms for each participant
df_pooled <- data_frame(
  Model = "Complete pooling",
  Subject = unique(d$sp_in),
  Intercept = coef(m_pooled)[1], 
  sla = coef(m_pooled)[2])

head(df_pooled)

m <- lmer(bday ~ 1+ sla + (1 +sla | sp_in), d)
arm::display(m)      

df_partial_pooling <- coef(m)[["sp_in"]] %>% 
  rownames_to_column("sp_in") %>% 
  as_tibble() %>% 
  add_column(Model = "Partial pooling")

head(df_partial_pooling)
colnames(df_partial_pooling)[colnames(df_partial_pooling)=="(Intercept)"] <- "Intercept"

# Also visualize the point for the fixed effects
df_fixef <- data_frame(
  Model = "Partial pooling (average)",
  Intercept = fixef(m)[1],
 sla = fixef(m)[2])

# Complete pooling / fixed effects are center of gravity in the plot
df_gravity <- df_pooled %>% 
  distinct(Model, Intercept, sla) %>% 
  bind_rows(df_fixef)
df_gravity
#> # A tibble: 2 x 3
#>                       Model Intercept Slope_Days
#>                       <chr>     <dbl>      <dbl>
#> 1          Complete pooling  252.3207   10.32766
#> 2 Partial pooling (average)  252.5426   10.45212

df_pulled <- bind_rows(df_no_pooling, df_partial_pooling)


pdf("test.pdf")
ggplot(df_pulled) + 
  aes(x = Intercept, y = sla, color = Model) + 
  geom_point(size = 2) + 
  geom_point(data = df_gravity, size = 5) +
  # Draw an arrow connecting the observations between models
  geom_path(aes(group = sp_in, color = NULL),
            arrow = arrow(length = unit(.02, "npc"))) +
#   # Use ggrepel to jitter the labels away from the points
  ggrepel::geom_text_repel(
    aes(label = sp_in, color = NULL),
    data = df_no_pooling) +
  theme(legend.position = "bottom") +
  ggtitle("Pooling of regression parameters") +
  xlab("Intercept estimate") +
  ylab("Slope estimate") +
  scale_color_brewer(palette = "Dark2")
dev.off()
################################################

head(d)
par(mfrow=c(1,1))
boxplot(bday~sp, data=d)
boxplot(sla~sp, data=d)
boxplot(ht~sp, data=d)
boxplot(cn~sp, data=d)
boxplot(stom_d~sp, data=d)
boxplot(wood_den~sp, data=d)
