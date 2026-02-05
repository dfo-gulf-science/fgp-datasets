##########################################
##########################################

############ Annotated Code ##############
###### ARTIFICIAL OYSTER REEFS CAN #######
####### FACILITATE THE RECOVERY OF ####### 
######## LOST ECOSYSTEM FUNCTION #########
##### IN FRAGMENTED SEAGRASS HABITAT #####

##########################################
##########################################

######## Get packages ########
library(dplyr)
library(multcomp)
library(stringi)
library(ggfortify)
library(car)
library(emmeans)
library(lmerTest)
library(patchwork)
library(data.table)
library(tidyverse)
library(Rmisc)
library(forcats)
library(plyr)
library(blmeco)
library(AER)
library(viridis)
library(ggtext)
library(tidymv)
library(devtools)
library(mgcv)
library(tidyverse)
library(broom)
library(plyr)
library(lme4)
library(ggpubr)
library(MASS)
library(cowplot)
library(blmeco)
library(MASS)
library(performance)
library(pscl)
library(optimx)
library(rstanarm)
library(vegan)
library(BiodiversityR)
library(readxl)
library(ggsci)
library(ggrepel)
library(ggforce)

-----------------------------

##############################
###### BIOLOGICAL DATA #######
##############################

#### STATISTICAL ANALYSES ####
## Species richness ####
# Upload dataset
sp.rich<-read.csv(file.choose()) #Filename: 'richness_biomass.csv'
attach(sp.rich)
summary(sp.rich)
sp.rich$year<-as.factor(sp.rich$year)
sp.rich$year.class<-as.factor(sp.rich$year.class)
sp.rich$pile.year<-as.factor(sp.rich$pile.year)
sp.rich$brick.all<-as.factor(sp.rich$brick.all)
sp.rich$brick.year<-as.factor(sp.rich$brick.year)

#Build model for species richness (use GLM with Poisson dist'n for count data)
m1<-glm(richness~year.class,family=poisson,data=sp.rich)

#Check for overdispersion
check_overdispersion(m1) #Not overdispersed; underdispersion due to 

#Get fixed effects
Anova(m1) #Significant effect of year

#Get pairwise results
rich.pw<-emmeans(m1,~year.class,data=sp.rich,adjustment="Holm")
pairs(rich.pw)
#0-year different from 1- and 2-year; 1- and 2-year similar

## Community biomass ####
#Use 'sp.rich' dataset from "Species richness" analysis above

#Build model for biomass (use LM for continuous data)
m2<-lm(biomass~year.class,data=sp.rich)

#Check assumptions
qqnorm(resid(m2))
qqline(resid(m2))
plot(m2)
#Assumptions violated

#Try log transform
log.biomass<-log(biomass)

#build model with log transformed data
m2a<-lm(log.biomass~year,data=sp.rich)

#Check assumptions
qqnorm(resid(m2a))
qqline(resid(m2a))
plot(m2a)

#still violated; use GLM with gamma dist'n

#build GLM model 
m2b<-glm(biomass~year.class,family=Gamma,data=sp.rich)

#Check for overdispersion
check_overdispersion(m2b) #Not overdispersed

#Get fixed effects
Anova(m2b) #Significant effect of year

#Get pairwise results
rich.pw<-emmeans(m2b,~year.class,data=sp.rich,adjustment="Holm")
pairs(rich.pw) #Everything significantly different

## Community composition dbRDA ####
#Upload data for whole-community biomass
biomass.dbrda<-read.csv(file.choose()) #Filename: 'community biomass_RDA.csv'
attach(biomass.dbrda)
summary(biomass.dbrda)
biomass.dbrda$year<-as.factor(biomass.dbrda$year)
biomass.dbrda$years.deployed<-as.factor(biomass.dbrda$years.deployed)
biomass.dbrda$brick<-as.factor(biomass.dbrda$brick)
biomass.dbrda$pile<-as.factor(biomass.dbrda$pile)
biomass.dbrda$obs<-as.factor(biomass.dbrda$obs)

#Upload data for sessile animal-only biomass
biomass.dbrda.sessile<-read.csv(file.choose()) #Filename: 'sessile biomass_RDA.csv'
attach(biomass.dbrda.sessile)
summary(biomass.dbrda.sessile)
biomass.dbrda.sessile$year<-as.factor(biomass.dbrda.sessile$year)
biomass.dbrda.sessile$years.deployed<-as.factor(biomass.dbrda.sessile$years.deployed)
biomass.dbrda.sessile$brick<-as.factor(biomass.dbrda.sessile$brick)
biomass.dbrda.sessile$pile<-as.factor(biomass.dbrda.sessile$pile)
biomass.dbrda.sessile$obs<-as.factor(biomass.dbrda.sessile$obs)

#Upload data for species counts
count.dbrda<-read.csv(file.choose()) #Filename: 'community abundance_RDA.csv'
attach(count.dbrda)
summary(count.dbrda)
count.dbrda$year<-as.factor(count.dbrda$year)
count.dbrda$years.deployed<-as.factor(count.dbrda$years.deployed)
count.dbrda$brick<-as.factor(count.dbrda$brick)
count.dbrda$pile<-as.factor(count.dbrda$pile)
count.dbrda$obs<-as.factor(count.dbrda$obs)

#Upload data for sessile only species counts
count.sessile.dbrda<-read.csv(file.choose()) #Filename: 'sessile abundance_RDA.csv'
attach(count.sessile.dbrda)
summary(count.sessile.dbrda)
count.sessile.dbrda$year<-as.factor(count.sessile.dbrda$year)
count.sessile.dbrda$years.deployed<-as.factor(count.sessile.dbrda$years.deployed)
count.sessile.dbrda$brick<-as.factor(count.sessile.dbrda$brick)
count.sessile.dbrda$pile<-as.factor(count.sessile.dbrda$pile)
count.sessile.dbrda$obs<-as.factor(count.sessile.dbrda$obs)

#Specify environmental + community data for each data set
fix(biomass.dbrda) #invokes edit on x an then assigns the new (edited) version of x in the user's workspace
row.names(biomass.dbrda)=biomass.dbrda$obs
biomass.species=biomass.dbrda[,6:42] # species start with Daphnia sp. at column 6
biomass.environment=biomass.dbrda[,2:5] # environment starts with pile at column 2
fix(biomass.species)
fix(biomass.environment)

fix(biomass.dbrda.sessile) #invokes edit on x an then assigns the new (edited) version of x in the user's workspace
row.names(biomass.dbrda.sessile)=biomass.dbrda.sessile$obs
biomass.species.sessile=biomass.dbrda.sessile[,6:19] # species start with Ciona at column 6
biomass.environment.sessile=biomass.dbrda.sessile[,2:5] # environment starts with pile at column 2
fix(biomass.species.sessile)
fix(biomass.environment.sessile)

fix(count.dbrda) #invokes edit on x an then assigns the new (edited) version of x in the user's workspace
row.names(count.dbrda)=count.dbrda$obs
count.species=count.dbrda[,6:33] # species start with Daphnia sp. at column 6
count.environment=count.dbrda[,2:5] # environment starts with pile at column 2
fix(count.species)
fix(count.environment)

fix(count.sessile.dbrda) #invokes edit on x an then assigns the new (edited) version of x in the user's workspace
row.names(count.sessile.dbrda)=count.sessile.dbrda$obs
count.sessile.species=count.sessile.dbrda[,6:19] # species start with Ciona at column 6
count.sessile.environment=count.sessile.dbrda[,2:5] # environment starts with pile at column 2
fix(count.sessile.species)
fix(count.sessile.environment)

#Add small number to zero values for each data set
biomass.species001= (biomass.species + 0.001)
fix(biomass.species001)

biomass.species.sessile001=(biomass.species.sessile + 0.001)
fix(biomass.species.sessile001)

count.species001= (count.species + 0.001)
fix(count.species001)

count.sessile.species001= (count.sessile.species + 0.001)
fix(count.sessile.species001)

#Transform species data for each data set
biomass.species.hell <- disttransform(biomass.species001, method='hellinger')
biomass.species.sessile.hell <- disttransform(biomass.species.sessile001, method='hellinger')
count.species.hell <- disttransform(count.species001, method='hellinger')
count.sessile.species.hell <- disttransform(count.sessile.species001, method='hellinger')

## Analyze and plot whole-community biomass dataset
#Build RDA model with years.deployed only
biomass.ord.m1 <- rda(biomass.species.hell ~ years.deployed, data=biomass.environment, scaling="biomass.species", method="braycurtis")

#View summary and get results and R2
summary(biomass.ord.m1)
anova(biomass.ord.m1) #significant
RsquareAdj(biomass.ord.m1) #high R2 (0.78) and Adj R2 (0.76)

#Build RDA model with years.deployed + pile (nested within years deployed)
biomass.ord.m2 <- rda(biomass.species.hell ~ years.deployed+pile/years.deployed, data=biomass.environment, scaling="biomass.species")

#View summary and get results and R2
summary(biomass.ord.m2)
anova(biomass.ord.m2) #significant
RsquareAdj(biomass.ord.m2) #high R2 (0.83) and Adj R2 (0.76)

#Compare models
anova(biomass.ord.m1,biomass.ord.m2) #not significantly different; use results from more parsimonious model (biomass.ord.m1)
anova(biomass.ord.m1)

#Get pairwise results for year classes
multiconstrained(method="rda", biomass.species.hell ~ years.deployed, data=biomass.environment,distance="euclidean",add=TRUE)

## Analyze and plot sessile biomass dataset
#Build RDA model with years.deployed only
biomass.sessile.ord.m1 <- rda(biomass.species.sessile.hell ~ years.deployed, data=biomass.environment.sessile, scaling="biomass.species.sessile")

#View summary and get results and R2
summary(biomass.sessile.ord.m1)
anova(biomass.sessile.ord.m1) #significant
RsquareAdj(biomass.sessile.ord.m1) #high R2 (0.85) and Adj R2 (0.83)

#Build RDA model with years.deployed + pile (nested within years deployed)
biomass.sessile.ord.m2 <- rda(biomass.species.sessile.hell ~ years.deployed+pile/years.deployed, data=biomass.environment.sessile, scaling="biomass.species")

#View summary and get results and R2
summary(biomass.sessile.ord.m2)
anova(biomass.sessile.ord.m2) #significant
RsquareAdj(biomass.sessile.ord.m2) #high R2 (0.88) and Adj R2 (0.82)

#Compare models
anova(biomass.sessile.ord.m1,biomass.sessile.ord.m2) #not significantly different; use results from more parsimonious model (biomass.ord.m1)
anova(biomass.sessile.ord.m1)

#Get pairwise results for year classes
multiconstrained(method="rda", biomass.species.sessile.hell ~ years.deployed, data=biomass.environment.sessile,distance="euclidean",add=TRUE)

## Analyze and plot count dataset
#Build RDA model with years.deployed only
count.ord.m1 <- rda(count.species.hell ~ years.deployed, data=count.environment, scaling="count.species")

#View summary and get results and R2
summary(count.ord.m1)
anova(count.ord.m1) #significant
RsquareAdj(count.ord.m1) #decent R2 (0.47) and Adj R2 (0.43)

#Build RDA model with years.deployed + pile (nested within years deployed)
count.ord.m2 <- rda(count.species.hell ~ years.deployed+pile/years.deployed, data=count.environment, scaling="count.species")

#View summary and get results and R2
summary(count.ord.m2)
anova(count.ord.m2) #significant
RsquareAdj(count.ord.m2) #decent R2 (0.55) and Adj R2 (0.35)

#Compare models
anova(count.ord.m1,count.ord.m2) #not significantly different; use results from more parsimonious model (count.ord.m1)
anova(count.ord.m1)

#Get pairwise results for year classes
multiconstrained(method="rda", count.species.hell ~ years.deployed, data=count.environment,distance="euclidean",add=TRUE)

## Analyze and plot sessile only count dataset
#Build RDA model with years.deployed only
count.sessile.ord.m1 <- rda(count.sessile.species.hell ~ years.deployed, data=count.sessile.environment, scaling="count.sessile.species")

#View summary and get results and R2
summary(count.sessile.ord.m1)
anova(count.sessile.ord.m1) #significant
RsquareAdj(count.sessile.ord.m1) #decent R2 (0.52) and Adj R2 (0.48)

#Build RDA model with years.deployed + pile (nested within years deployed)
count.sessile.ord.m2 <- rda(count.sessile.species.hell ~ years.deployed+pile/years.deployed, data=count.sessile.environment, scaling="count.species")

#View summary and get results and R2
summary(count.sessile.ord.m2)
anova(count.sessile.ord.m2) #significant
RsquareAdj(count.sessile.ord.m2) #decent R2 (0.63) and Adj R2 (0.46)

#Compare models
anova(count.sessile.ord.m1,count.sessile.ord.m2) #not significantly different; use results from more parsimonious model (count.ord.m1)
anova(count.sessile.ord.m1)

#Get pairwise results for year classes
multiconstrained(method="rda", count.sessile.species.hell ~ years.deployed, data=count.sessile.environment,distance="euclidean",add=TRUE)

## Oyster abundance & biomass ####
#Upload dataset
oysters<-read.csv(file.choose()) #Filename: 'oyster_abundance_biomass.csv'
attach(oysters)
summary(oysters)
oysters$year<-as.factor(oysters$year)
oysters$year.class<-as.factor(oysters$year.class)
oysters$brick<-as.factor(oysters$brick)
oysters$pile.year<-as.factor(oysters$pile.year)

#Split dataset by "adult.spat" factor levels
adults <- filter(oysters, adult.spat == "Adult")[,match(c("brick","pile.year","year","year.class","biomass","count"), colnames (oysters))]
spat <- filter(oysters, adult.spat =="Spat")[,match(c("brick","pile.year","year","year.class","biomass","count"), colnames (oysters))]
summary(adults)
summary(spat)

#Build model for adult counts (use GLM with Poisson dist'n for count data)
adult.m1<-glm(count~year.class,family=poisson,data=adults,link="log")

stan_glm(data = adults, count~year.class,family = neg_binomial_2)

zeroinfl(count~year.class,data = adults, dist = "negbin")

#Check for overdispersion
check_overdispersion(adult.m1) #Overdispersed

#Use negative binomial GLMER
adult.m2<-glm.nb(count~year.class,data=adults)

#Check for overdispersion
check_overdispersion(adult.m2) #Not overdispersed

#Get results
Anova(adult.m2) #Significant year effect

#Get pairwise results
adult.count.pw<-emmeans(adult.m2,~year.class,data=adults,adjustment="Holm")
pairs(adult.count.pw) #1-year signiciantly different from 2-year; non-significant differences between 0-year driven by all 0 values for 0-year adults, so we can assume significant differences

#Build model for adult biomass (use GLM with Poisson dist'n)
adult.m3<-glm(biomass~year.class,family=poisson,data=adults)

#Check for overdispersion
check_overdispersion(adult.m3) #Overdispersed

#Build negative binomial model for adult biomass
adult.m4<-glm.nb(biomass~year.class,data=adults)

#Check for overdispersion
check_overdispersion(adult.m4) #Not overdispersed

#Get results
Anova(adult.m4) #Significant year effect

#Get pairwise results
adult.biomass.pw<-emmeans(adult.m4,~year.class,data=adults,adjustment="Holm")
pairs(adult.biomass.pw) #1-year signiciantly different from 2-year; non-significant differences between 0-year driven by all 0 values for 0-year adults, so we can assume significant differences

#Build model for spat counts (use GLM with Poisson dist'n for count data)
spat.m1<-glm(count~year.class,family=poisson,data=spat)

#Check for overdispersion
check_overdispersion(spat.m1) #Overdispersed

#Build negative binomial model for spat counts
spat.m2<-glm.nb(count~year.class,data=spat)

#Check for overdispersion
check_overdispersion(spat.m2) #Not overdispersed

#Get results
Anova(spat.m2) #Significant year effect

#Get pairwise results
spat.count.pw<-emmeans(spat.m2,~year.class,data=adults,adjustment="Holm")
pairs(spat.count.pw) #1-year significantly different from 0-year and 2-year; marginally non-significant difference between 0-year and 2-year

#Build model for spat biomass (use GLM with Poisson dist'n)
spat.m3<-glm(biomass~year.class,family=poisson,data=spat)

#Check for overdispersion
check_overdispersion(spat.m3) #Not overdispersed

#Get results
Anova(spat.m3) #Significant year effect

#Get pairwise results
spat.biomass.pw<-emmeans(spat.m3,~year.class,data=adults,adjustment="Holm")
pairs(spat.biomass.pw) #0-year significantly different from 1-year and 2-year; no difference between 1-year and 2-year

#### FIGURES ####
## Figure 1 (map) ####
No code required; figure generated using existing maps and PowerPoint
## Figure 2 (species richness list) ####
No code required; figure generated using Excel and PowerPoint
## Figure 3 (Species richness & biomass boxplots) ####
#Generate boxplots
rich.plot<-ggplot(sp.rich,aes(x=year.class, y=richness))+stat_boxplot(geom ='errorbar',width=0.25)+ geom_boxplot(aes(fill=year.class),color="black",outlier.shape = "") + geom_jitter(aes(fill=year.class,shape=pile.year),size=3,alpha=0.4,width=0.36)+ scale_fill_manual(values = c("#5a8ac6","lightgray","#f8696b"))+ scale_color_manual(values = c("#5a8ac6","lightgray","#f8696b"))+scale_y_continuous(limits=c(0,24),breaks=seq(0, 24, 5)) + labs(x = "Reef year class", y="Taxonomic richness (count)")+ theme_bw(12)+ theme(legend.position="none")+ theme(panel.grid.minor = element_blank()) +theme(panel.grid.major = element_blank())+  annotate("text", x="2 year", y=23.2, label= "b")+  annotate("text", x="1 year", y=24, label= "b")+  annotate("text", x="0 year", y=10, label= "a")+  annotate("text", x="0 year", y=23, label= "A",hjust=2.2,size=6)
biomass.plot<-ggplot(sp.rich,aes(x=year.class, y=biomass))+stat_boxplot(geom ='errorbar',width=0.25)+ geom_boxplot(aes(fill=year.class),color="black",outlier.shape = "") + geom_jitter(aes(fill=year.class,shape=pile.year),size=3,alpha=0.4,width=0.36)+ scale_fill_manual(values = c("#5a8ac6","lightgray","#f8696b"))+ scale_color_manual(values = c("#5a8ac6","lightgray","#f8696b"))+ labs(x = "Reef year class", y="Total biomass (g)")+ theme_bw(12)+ theme(legend.position="none")+ theme(panel.grid.minor = element_blank()) +theme(panel.grid.major = element_blank())+  annotate("text", x="2 year", y=2210, label= "c")+  annotate("text", x="1 year", y=560, label= "b")+  annotate("text", x="0 year", y=180, label= "a") + annotate("text", x="0 year", y=2100, label= "B",hjust=2.2,size=6)

full.plot<-rich.plot+biomass.plot
full.plot

tiff("Oyster resto richness biomass FINAL2_lowercase.tiff", width = 7, height = 4, units = "in", res = 800)
full.plot
dev.off()

pdf("Oyster resto richness biomass FINAL2_lowercase.pdf", width = 7, height = 4)
full.plot
dev.off()

## Figure 4 (Ordination plots) ####
#Whole community biomass
biomass.p1<- ordiplot(biomass.ord.m1, choices=c(1,2))

sites.long1 <- sites.long(biomass.p1, env.data=biomass.environment)
head(sites.long1)

species.long1 <- species.long(biomass.p1)
species.long1

axis.long1 <- axis.long(biomass.ord.m1, choices=c(1, 2))
axis.long1

spec.envfit <- envfit(biomass.p1, env=biomass.species.hell)
spec.data.envfit <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)

species.long2 <- species.long(biomass.p1, spec.data=spec.data.envfit)
species.long2

species.long3 <- species.long2[species.long2$r >= 0.4, ]
species.long3

biplot.biomass <- ggplot() +  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) + geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  xlab(axis.long1[1, "label"]) +  ylab(axis.long1[2, "label"]) +  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +  geom_point(data=sites.long1,aes(x=axis1, y=axis2, colour=years.deployed), size=5,alpha=0.5) + geom_segment(data=species.long2, aes(x=0, y=0, xend=axis1*4, yend=axis2*4),colour="tan1", size=0.7, arrow=arrow()) +  geom_text_repel(data=species.long2, aes(x=axis1*4, y=axis2*4, label=labels),colour="gray30") +  coord_fixed(ratio=1)+theme_bw()+theme(panel.grid = element_blank(),legend.position = c(0.2,0.15))+theme(panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line("gray25"), text = element_text(size = 12), axis.text = element_text(size = 10, colour = "gray25"),axis.title = element_text(size = 14, colour = "gray25"),legend.title = element_text(size = 14),legend.text = element_text(size = 14),legend.key = element_blank())+ scale_color_manual(values = c("#5a8ac6","gray70","#f8696b"))+ scale_fill_manual(values = c("#5a8ac6","gray70","#f8696b"))+labs(colour="Reef year class")+scale_y_continuous(limits=c(-4.5,4.5),breaks=seq(-4.5, 4.5, 1))+scale_x_continuous(limits=c(-4.5,4.5),breaks=seq(-4.5, 4.5, 1))+annotate("text", x=-4.5, y=4.5, label= "A", size=8)+annotate("text", x=4.5, y=-4.5, label= "Whole-community biomass", size=4.5,hjust=1,fontface=3)
biplot.biomass

#Sessile biomass
biomass.sessile.p1<- ordiplot(biomass.sessile.ord.m1, choices=c(1,2))

sites.sessile.long1 <- sites.long(biomass.sessile.p1, env.data=biomass.environment.sessile)
head(sites.sessile.long1)

species.sessile.long1 <- species.long(biomass.sessile.p1)
species.sessile.long1

axis.sessile.long1 <- axis.long(biomass.sessile.ord.m1, choices=c(1, 2))
axis.sessile.long1

spec.envfit.sessile <- envfit(biomass.sessile.p1, env=biomass.species.sessile.hell)
spec.data.envfit.sessile <- data.frame(r=spec.envfit.sessile$vectors$r, p=spec.envfit.sessile$vectors$pvals)

species.sessile.long2 <- species.long(biomass.sessile.p1, spec.data=spec.data.envfit.sessile)
species.sessile.long2

species.sessile.long3 <- species.sessile.long2[species.sessile.long2$r >= 0.4, ]
species.sessile.long3

biplot.sessile.biomass <- ggplot() +  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) + geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  xlab(axis.sessile.long1[1, "label"]) +  ylab(axis.sessile.long1[2, "label"]) +  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +  geom_point(data=sites.sessile.long1,aes(x=axis1, y=axis2, colour=years.deployed), size=5,alpha=0.5) + geom_segment(data=species.sessile.long2, aes(x=0, y=0, xend=axis1*4, yend=axis2*4),colour="tan1", size=0.7, arrow=arrow()) +  geom_text_repel(data=species.sessile.long2, aes(x=axis1*4, y=axis2*4, label=labels),colour="gray30") +  coord_fixed(ratio=1)+theme_bw()+theme(panel.grid = element_blank(),legend.position = "none")+theme(panel.background = element_blank(),panel.grid = element_blank(), axis.line = element_line("gray25"), text = element_text(size = 12), axis.text = element_text(size = 10, colour = "gray25"),axis.title = element_text(size = 14, colour = "gray25"),legend.title = element_text(size = 14),legend.text = element_text(size = 14),legend.key = element_blank())+ scale_color_manual(values = c("#5a8ac6","gray70","#f8696b"))+ scale_fill_manual(values = c("#5a8ac6","gray70","#f8696b"))+labs(colour="Reef year class")+scale_y_continuous(limits=c(-4.5,4.5),breaks=seq(-4.5, 4.5, 1))+annotate("text", x=4.5, y=4.5, label= "B", size=8)+ scale_x_reverse()+annotate("text", x=-4.5, y=-4.5, label= "Sessile species biomass", size=4.5,hjust=1,fontface=3)+ scale_x_reverse(breaks=c(seq(-4.5,4.5,by=1)))
biplot.sessile.biomass

#Whole community abundance
count.p1<- ordiplot(count.ord.m1, choices=c(1,2))

sites.count.long1 <- sites.long(count.p1, env.data=count.environment)
head(sites.count.long1)

species.count.long1 <- species.long(count.p1)
species.count.long1

axis.count.long1 <- axis.long(count.ord.m1, choices=c(1, 2))
axis.count.long1

spec.envfit.count <- envfit(count.p1, env=count.species.hell)
spec.data.envfit.count <- data.frame(r=spec.envfit.count$vectors$r, p=spec.envfit.count$vectors$pvals)

species.count.long2 <- species.long(count.p1, spec.data=spec.data.envfit.count)
species.count.long2

species.count.long3 <- species.count.long2[species.count.long2$r >= 0.4, ]
species.count.long3

biplot.count <- ggplot() +  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) + geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  xlab(axis.count.long1[1, "label"]) +  ylab(axis.count.long1[2, "label"]) +  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +  geom_point(data=sites.count.long1,aes(x=axis1, y=axis2, colour=years.deployed), size=5,alpha=0.5) + geom_segment(data=species.count.long2, aes(x=0, y=0, xend=axis1*4, yend=axis2*4),colour="tan1", size=0.7, arrow=arrow()) +  geom_text_repel(data=species.count.long2, aes(x=axis1*4, y=axis2*4, label=labels),colour="gray30") +  coord_fixed(ratio=1)+theme_bw()+theme(panel.grid = element_blank(),legend.position = "none")+theme(panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line("gray25"), text = element_text(size = 12), axis.text = element_text(size = 10, colour = "gray25"),axis.title = element_text(size = 14, colour = "gray25"),legend.title = element_text(size = 14),legend.text = element_text(size = 14),legend.key = element_blank())+ scale_color_manual(values = c("#5a8ac6","gray70","#f8696b"))+ scale_fill_manual(values = c("#5a8ac6","gray70","#f8696b"))+labs(colour="Reef year class")+scale_y_continuous(limits=c(-2.5,2.5),breaks=seq(-2.5, 2.5, 0.5))+scale_x_continuous(limits=c(-2.5,2.5),breaks=seq(-2.5, 2.5, 0.5))+annotate("text", x=-2.5, y=2.5, label= "C", size=8)+annotate("text", x=2.5, y=-2.5, label= "Whole-community abundance", size=4.5,hjust=1,fontface=3)
biplot.count

#Sessile abundance
count.sessile.p1<- ordiplot(count.sessile.ord.m1, choices=c(1,2))

sites.count.sessile.long1 <- sites.long(count.sessile.p1, env.data=count.sessile.environment)
head(sites.count.sessile.long1)

species.count.sessile.long1 <- species.long(count.sessile.p1)
species.count.sessile.long1

axis.count.sessile.long1 <- axis.long(count.sessile.ord.m1, choices=c(1, 2))
axis.count.sessile.long1

spec.envfit.count.sessile <- envfit(count.sessile.p1, env=count.sessile.species.hell)
spec.data.envfit.count.sessile <- data.frame(r=spec.envfit.count.sessile$vectors$r, p=spec.envfit.count.sessile$vectors$pvals)

species.count.sessile.long2 <- species.long(count.sessile.p1, spec.data=spec.data.envfit.count.sessile)
species.count.sessile.long2

species.count.sessile.long3 <- species.count.sessile.long2[species.count.sessile.long2$r >= 0.4, ]
species.count.sessile.long3

biplot.count.sessile <- ggplot() +  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) + geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  xlab(axis.count.long1[1, "label"]) +  ylab(axis.count.long1[2, "label"]) +  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +  geom_point(data=sites.count.sessile.long1,aes(x=axis1, y=axis2, colour=years.deployed), size=5,alpha=0.5) + geom_segment(data=species.count.sessile.long2, aes(x=0, y=0, xend=axis1*4, yend=axis2*4),colour="tan1", size=0.7, arrow=arrow()) +  geom_text_repel(data=species.count.sessile.long2, aes(x=axis1*4, y=axis2*4, label=labels),colour="gray30") +  coord_fixed(ratio=1)+theme_bw()+theme(panel.grid = element_blank(),legend.position = "none")+theme(panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line("gray25"), text = element_text(size = 12), axis.text = element_text(size = 10, colour = "gray25"),axis.title = element_text(size = 14, colour = "gray25"),legend.title = element_text(size = 14),legend.text = element_text(size = 14),legend.key = element_blank())+ scale_color_manual(values = c("#5a8ac6","gray70","#f8696b"))+ scale_fill_manual(values = c("#5a8ac6","gray70","#f8696b"))+labs(colour="Reef year class")+scale_y_continuous(limits=c(-2.5,2.5),breaks=seq(-2.5, 2.5, 0.5))+annotate("text", x=2.5, y=2.5, label= "D", size=8)+annotate("text", x=-2.5, y=-2.5, label= "Sessile species abundance", size=4.5,hjust=1,fontface=3)+ scale_x_reverse(breaks=c(seq(-2.5,2.5,by=0.5)))
biplot.count.sessile

#Combine biplots
biplot.combo<-(biplot.biomass+biplot.sessile.biomass)/(biplot.count+biplot.count.sessile)
biplot.combo

#Export figure
tiff("Habitat resto biplots FINAL_JUNE 2025.tiff", width = 11, height = 11, units = "in", res = 800)
biplot.combo
dev.off()

pdf("Habitat resto biplots FINAL_JUNE 2025.pdf", width = 11, height = 11)
biplot.combo
dev.off()

## Figure 5 (Oyster abundance & biomass boxplots) ####
adult.count.plot<-ggplot(adults,aes(x=year.class, y=count))+stat_boxplot(geom ='errorbar',width=0.25)+ geom_boxplot(aes(fill=year.class),color="black",outlier.shape="") + geom_jitter(aes(fill=year.class,shape=pile.year),size=3,alpha=0.4,width=0.35)+ scale_fill_manual(values = c("#5a8ac6","lightgray","#f8696b"))+ scale_color_manual(values = c("#5a8ac6","lightgray","#f8696b"))+scale_y_continuous(limits=c(-1,320),breaks=seq(0, 300, 50)) + labs(x = "Reef year class", y="Oyster count")+ theme_bw(12)+ theme(legend.position="none") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_blank(),plot.title = element_text(size=18,face="bold"))+  annotate("text", x="0 year", y=20, label= "a")+  annotate("text", x="1 year", y=40, label= "b")+  annotate("text", x="2 year", y=290, label= "c")+  annotate("text", x="0 year", y=310, label= "B",size=6,hjust=2.2)+  annotate("text", x="2 year", y=320, label= "Adult count",size=3.5,hjust=0.48,fontface=3)

adult.biomass.plot<-ggplot(adults,aes(x=year.class, y=biomass))+stat_boxplot(geom ='errorbar',width=0.25)+ geom_boxplot(aes(fill=year.class),color="black",outlier.shape = "") + geom_jitter(aes(fill=year.class,shape=pile.year),size=3,alpha=0.4,width=0.35)+ scale_fill_manual(values = c("#5a8ac6","lightgray","#f8696b"))+ scale_color_manual(values = c("#5a8ac6","lightgray","#f8696b"))+ scale_y_continuous(limits=c(-1,1800),breaks=seq(0, 1750, 250))+labs(x = "Reef year class", y="Oyster biomass (g)")+ theme_bw(12)+ theme(legend.position="none")+ theme(panel.grid.minor = element_blank(), axis.title.y = element_blank()) +theme(panel.grid.major = element_blank())+  annotate("text", x="0 year", y=100, label= "a")+  annotate("text", x="1 year", y=120, label= "b")+  annotate("text", x="2 year", y=1600, label= "c")+  annotate("text", x="0 year", y=1750, label= "D",size=6,hjust=2.2)+  annotate("text", x="2 year", y=1800, label= "Adult biomass",size=3.5,hjust=0.6,fontface=3)

spat.count.plot<-ggplot(spat,aes(x=year.class, y=count))+stat_boxplot(geom ='errorbar',width=0.25)+ geom_boxplot(aes(fill=year.class),color="black",outlier.shape="") + geom_jitter(aes(fill=year.class,shape=pile.year),size=3,alpha=0.4,width=0.35)+ scale_fill_manual(values = c("#5a8ac6","lightgray","#f8696b"))+ scale_color_manual(values = c("#5a8ac6","lightgray","#f8696b"))+scale_y_continuous(limits=c(-1,2500),breaks=seq(0, 2500, 500)) + labs(x = "Reef year class", y="Oyster count")+ theme_bw(12)+ theme(legend.position="none") + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(), axis.title.x = element_blank(),axis.text.x = element_blank(),plot.title = element_text(size=18,face="bold"))+  annotate("text", x="0 year", y=1600, label= "a")+  annotate("text", x="1 year", y=350, label= "b")+  annotate("text", x="2 year", y=650, label= "a")+  annotate("text", x="0 year", y=2425, label= "A",size=6,hjust=2.2)+  annotate("text", x="2 year", y=2500, label= "Spat count",size=3.5,hjust=0.48,fontface=3)

spat.biomass.plot<-ggplot(spat,aes(x=year.class, y=biomass))+stat_boxplot(geom ='errorbar',width=0.25)+ geom_boxplot(aes(fill=year.class),color="black",outlier.shape = "") + geom_jitter(aes(fill=year.class,shape=pile.year),size=3,alpha=0.4,width=0.35)+ scale_fill_manual(values = c("#5a8ac6","lightgray","#f8696b"))+ scale_color_manual(values = c("#5a8ac6","lightgray","#f8696b"))+ scale_y_continuous(limits=c(-0.5,15),breaks=seq(0, 14, 2))+labs(x = "Reef year class", y="Oyster biomass (g)")+ theme_bw(12)+ theme(legend.position="none")+ theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+  annotate("text", x="0 year", y=12.2, label= "a")+  annotate("text", x="1 year", y=4, label= "b")+  annotate("text", x="2 year", y=5.5, label= "b")+  annotate("text", x="0 year", y=14.5, label= "C",size=6,hjust=2.2)+  annotate("text", x="2 year", y=15, label= "Spat biomass",size=3.5,hjust=0.6,fontface=3)

full.oyster.plot<- (spat.count.plot + adult.count.plot) / (spat.biomass.plot + adult.biomass.plot)
full.oyster.plot

tiff("Oyster counts biomass habitat restoration FINAL_lowercase.tiff", width = 7, height = 7, units = "in", res = 800)
full.oyster.plot
dev.off()

pdf("Oyster counts biomass habitat restoration FINAL_lowercase.pdf", width = 7, height = 7)
full.oyster.plot
dev.off()

## Figure 6 (Fish use) ####
No code required; figure generated using Excel and PowerPoint



#############################
####### ABIOTIC DATA ########
#############################
##### Upload datasets ####
hourly.abiotic<-read.csv(file.choose()) #Filename: abiotic_hourly.csv
attach(hourly.abiotic)
summary(hourly.abiotic)
hourly.abiotic$stand<-as.factor(hourly.abiotic$stand)
hourly.abiotic$date<-as.factor(hourly.abiotic$date)
hourly.abiotic$julian.date<-as.factor(hourly.abiotic$julian.date)
hourly.abiotic$hour<-as.factor(hourly.abiotic$hour)
hourly.abiotic$time.of.day<-as.factor(hourly.abiotic$time.of.day)
summary(hourly.abiotic)

#No 'Reef' data for DO; upload separate dataset
hourly.do<-read.csv(file.choose()) #Filename: do_hourly.csv
attach(hourly.do)
summary(hourly.do)
hourly.do$stand<-as.factor(hourly.do$stand)
hourly.do$date<-as.factor(hourly.do$date)
hourly.do$julian.date<-as.factor(hourly.do$julian.date)
hourly.do$hour<-as.factor(hourly.do$hour)
hourly.do$time.of.day<-as.factor(hourly.do$time.of.day)
summary(hourly.do)

#Missing turbidity data after Aug 13 for 'Reef'; upload separate dataset
partial.turb<-read.csv(file.choose()) #Filename: turbidity_partial.csv
attach(partial.turb)
summary(partial.turb)
partial.turb$stand<-as.factor(partial.turb$stand)
partial.turb$date<-as.factor(partial.turb$date)
partial.turb$julian.date<-as.factor(partial.turb$julian.date)
partial.turb$hour<-as.factor(partial.turb$hour)
partial.turb$time.of.day<-as.factor(partial.turb$time.of.day)
summary(partial.turb)

#### STATISTICAL ANALYSES ####
## Dissolved oxygen ####
#note: no data for 'Reef'; logger malfunctioned
do.mod<-lm(do~stand*time.of.day+(julian.date*stand),data=hourly.do)
anova(do.mod) #Significant day effect (unsurprising) and Stand X Time of Day effect; look for differences between stands, as well as differences between time of day within stands

#Get pairwise results
#Stand only
stand.do.pw<-emmeans(do.mod,~stand,data=hourly.do,adjustment="Holm")
pairs(stand.do.pw) #DO Inside reefs significantly higher than Outside of reefs.

#Stand x Time of Day interaction
ixn.do.pw<-emmeans(do.mod,~time.of.day|stand,data=hourly.do,adjustment="Holm")
pairs(ixn.do.pw) #Daytime DO higher than Night for Inside but not Outside; Reefs generate day-night DO cycle whereas bare sediment patches do not. Likely due to biological activity. Coincides with finding re. pH below

## pH ####
ph.mod<-lm(sw.ph~stand*time.of.day+(julian.date*stand),data=hourly.abiotic)
anova(ph.mod) #Significant day effect (unsurprising) and Stand X Time of Day effect; look for differences between stands, as well as differences between time of day within stands

#Get pairwise results
#Stand only
stand.ph.pw<-emmeans(ph.mod,~stand,data=hourly.abiotic,adjustment="Holm")
pairs(stand.ph.pw) #Everything significantly different; Inside lower than Outside, Reef lower than both Inside and Outside; pH declines as you move closer to the reefs

#Stand x Time of Day nteraction
ixn.ph.pw<-emmeans(ph.mod,~time.of.day|stand,data=hourly.abiotic,adjustment="Holm")
pairs(ixn.ph.pw) #Daytime pH higher than Night for Inside and Reef, but not Outside; Reefs generate day-night pH cycle whereas bare sediment patches do not. Likely due to biological activity.

## Temperature ####
temp.mod<-lm(temp~stand*time.of.day+(julian.date*stand),data=hourly.abiotic)
anova(temp.mod) #Significant day effect (unsurprising) and Stand X Time of Day effect; look for differences between stands, as well as differences between time of day within stands

#Get pairwise results
#Stand only
stand.temp.pw<-emmeans(temp.mod,~stand,data=hourly.abiotic,adjustment="Holm")
pairs(stand.temp.pw) #Inside and Reef temperatures significantly lower than Outside; cooling effect of reefs, perhaps due to shading from vegetation.

#Stand x Time of Day nteraction
ixn.temp.pw<-emmeans(temp.mod,~time.of.day|stand,data=hourly.abiotic,adjustment="Holm")
pairs(ixn.temp.pw) #Daytime temps significantly higher than Night for Inside and Reef, but not Outside; Reefs generate day-night temperature cycle whereas bare sediment patches do not. Similar to pH and DO trends.

## Turbidity ####
#NOTE: Missing data after Aug 13 for 'Outside' reference area; analyze separate dataset up to Aug 13 (19 days)
#NOTE: Removed one outlier (Outside, Aug 2, 18:00PM; value was >4X larger tan next highest value)
turbidity.mod<-glm(turbidity~stand*time.of.day+(julian.date*stand),family=Gamma,data=partial.turb)
Anova(turbidity.mod) #All significant

#Get pairwise results
#Stand only
stand.turbidity.pw<-emmeans(turbidity.mod,~stand,data=partial.turb,adjustment="Holm")
pairs(stand.turbidity.pw) #Reef significantly different from Inside and Outside; Inside-Outside not different

#Stand x Time of Day nteraction
ixn.turbidity.pw<-emmeans(turbidity.mod,~time.of.day|stand,data=partial.turb,adjustment="Holm")
pairs(ixn.turbidity.pw) #Daytime turbidity higher than Night at the outside reference site, but not near reefs.

#### FIGURES ####
## Figure 7 (Abiotic time series) ####
hourly.abiotic$stand<-factor(hourly.abiotic$stand,levels=c("Reef","Inside","Outside"))
hourly.do$stand<-factor(hourly.do$stand,levels=c("Reef","Inside","Outside"))
partial.turb$stand<-factor(partial.turb$stand,levels=c("Reef","Inside","Outside"))
do.timeseries.plot<-ggplot(hourly.do,aes(x=julian.date, y=do,group=stand))+ geom_point(aes(fill=stand,color=stand,group=stand,shape=stand),outlier.shape=NA,size=3,alpha=0.2)+geom_smooth(aes(fill=stand,color=stand), method="loess", se = TRUE,alpha=0.45)+scale_color_manual(values=c("#3D399F","#B94DAC", "#F9A669"))+ scale_fill_manual(values=c("#3D399F","#B94DAC", "#F9A669"))  + labs(x = "Julian date", y="DO (mg/L)")+ theme_bw(12)+ theme(legend.position="none",legend.title = element_blank())+ theme(panel.grid.minor = element_blank()) +theme(panel.grid.major = element_blank(),axis.text.x = element_blank(),axis.title.x=element_blank())+scale_y_continuous(limits=c(5.5, 15.5),breaks=seq(5.5, 15.5, 1))+annotate("text", x="207", y=15.5, label= "A", size=6,hjust=0,vjust=1)

ph.timeseries.plot<-ggplot(hourly.abiotic,aes(x=julian.date, y=sw.ph,group=stand))+ geom_point(aes(fill=stand,color=stand,group=stand,shape=stand),outlier.shape=NA,size=3,alpha=0.2)+geom_smooth(aes(fill=stand,color=stand), method="loess", se = TRUE,alpha=0.45)+scale_color_manual(values=c("#3D399F","#B94DAC", "#F9A669"))+ scale_fill_manual(values=c("#3D399F","#B94DAC", "#F9A669")) + labs(x = "Julian date", y="pH")+ theme_bw(12)+ theme(legend.position="none",legend.title = element_blank())+ theme(panel.grid.minor = element_blank()) +theme(panel.grid.major = element_blank(),axis.text.x = element_blank(),axis.title.x=element_blank())+scale_y_continuous(limits=c(7.7, 8.5),breaks=seq(7.7, 8.5, 0.1))+annotate("text", x="207", y=8.5, label= "B", size=6,hjust=0,vjust=1)

temp.timeseries.plot<-ggplot(hourly.abiotic,aes(x=julian.date, y=temp,group=stand))+ geom_point(aes(fill=stand,color=stand,group=stand,shape=stand),outlier.shape=NA,size=3,alpha=0.2)+geom_smooth(aes(fill=stand,color=stand), method="loess", se = TRUE,alpha=0.45)+scale_color_manual(values=c("#3D399F","#B94DAC", "#F9A669"))+ scale_fill_manual(values=c("#3D399F","#B94DAC", "#F9A669"))  + labs(x = "Julian date", y="Temperature (°C)")+ theme_bw(12)+ theme(legend.position="none",legend.title = element_blank())+ theme(panel.grid.minor = element_blank()) +theme(panel.grid.major = element_blank(),axis.text.x = element_blank(),axis.title.x=element_blank())+scale_y_continuous(limits=c(18, 28),breaks=seq(18, 28, 1))+annotate("text", x="207", y=28, label= "C", size=6,hjust=0,vjust=1)

turb.timeseries.plot<-ggplot(hourly.abiotic,aes(x=julian.date, y=turbidity,group=stand))+ geom_point(aes(fill=stand,color=stand,group=stand,shape=stand),outlier.shape=NA,size=3,alpha=0.2)+geom_smooth(aes(fill=stand,color=stand), method="loess", se = TRUE,alpha=0.45)+scale_color_manual(values=c("#3D399F","#B94DAC", "#F9A669"))+ scale_fill_manual(values=c("#3D399F","#B94DAC", "#F9A669"))  + labs(x = "Julian date", y="Turbidity (NTU)")+ theme_bw(12)+ theme(legend.position=c(0.9,0.78),legend.title = element_blank())+ theme(panel.grid.minor = element_blank()) +theme(panel.grid.major = element_blank())+scale_y_continuous(limits=c(0, 100),breaks=seq(0, 100, 10))+annotate("text", x="207", y=100, label= "D", size=6,hjust=0,vjust=1)

#Combine plots
combo.timeseries.plot<-do.timeseries.plot/ph.timeseries.plot/temp.timeseries.plot/turb.timeseries.plot
combo.timeseries.plot

#Export figure
tiff("Abiotic timeseries_JAN 2026.tiff", width = 11, height = 13, units = "in", res = 800)
combo.timeseries.plot
dev.off()

pdf("Abiotic timeseries_JAN 2026.pdf", width = 11, height = 13)
combo.timeseries.plot
dev.off()


## Figure 8 (Boxplots by stand; left panels) ####
do.boxplot<-ggplot(hourly.abiotic,aes(x=stand, y=do))+ geom_jitter(aes(color=stand,fill=stand),size=2.5,alpha=0.2,width=0.2)+stat_boxplot(geom ='errorbar',width=0.25)+ geom_boxplot(aes(fill=stand),color="black",outlier.shape=NA) +scale_color_manual(values=c("#3D399F","#B94DAC","#F9A669"))+ scale_fill_manual(values=c("#3D399F","#B94DAC","#F9A669"))  + labs(x = "Location", y="Dissolved oxygen (mg/L)")+ theme_bw(12)+ theme(legend.position="none")+ theme(panel.grid.minor = element_blank()) +theme(panel.grid.major = element_blank(),axis.text.x=element_blank(),axis.title.x = element_blank())+scale_y_continuous(limits=c(5.5, 15.5),breaks=seq(5.5, 15.5, 1))+annotate("text", x="Reef", y=5.52, label= "No data", size=5,fontface=3)+annotate("text", x="Inside", y=10.5, label= "a", size=5,hjust=-4)+annotate("text", x="Outside", y=10.05, label= "b", size=5,hjust=-4)

ph.boxplot<-ggplot(hourly.abiotic,aes(x=stand, y=sw.ph))+ geom_jitter(aes(color=stand,fill=stand),size=2.5,alpha=0.2,width=0.2)+stat_boxplot(geom ='errorbar',width=0.25)+ geom_boxplot(aes(fill=stand),color="black",outlier.shape=NA) +scale_color_manual(values=c("#3D399F","#B94DAC","#F9A669"))+ scale_fill_manual(values=c("#3D399F","#B94DAC","#F9A669"))  + labs(x = "Location", y="pH")+ theme_bw(12)+ theme(legend.position="none")+ theme(panel.grid.minor = element_blank()) +theme(panel.grid.major = element_blank(),axis.text.x=element_blank(),axis.title.x = element_blank())+scale_y_continuous(limits=c(7.7, 8.5),breaks=seq(7.7, 8.5, 0.1))+annotate("text", x="Reef", y=8.13, label= "a", size=5,hjust=-4)+annotate("text", x="Inside", y=8.19, label= "b", size=5,hjust=-4)+annotate("text", x="Outside", y=8.25, label= "c", size=5,hjust=-4)

temp.boxplot<-ggplot(hourly.abiotic,aes(x=stand, y=temp))+ geom_jitter(aes(color=stand,fill=stand),size=2.5,alpha=0.2,width=0.2)+stat_boxplot(geom ='errorbar',width=0.25)+ geom_boxplot(aes(fill=stand),color="black",outlier.shape=NA) +scale_color_manual(values=c("#3D399F","#B94DAC","#F9A669"))+ scale_fill_manual(values=c("#3D399F","#B94DAC","#F9A669"))  + labs(x = "Location", y="Temperature (°C)")+ theme_bw(12)+ theme(legend.position="none")+ theme(panel.grid.minor = element_blank()) +theme(panel.grid.major = element_blank(),axis.text.x=element_blank(),axis.title.x = element_blank())+scale_y_continuous(limits=c(18, 28),breaks=seq(18, 28, 1))+annotate("text", x="Reef", y=23.6, label= "a", size=5,hjust=-4)+annotate("text", x="Inside", y=23.6, label= "a", size=5,hjust=-4)+annotate("text", x="Outside", y=23.6, label= "b", size=5,hjust=-4)

turb.boxplot<-ggplot(partial.turb,aes(x=stand, y=turbidity))+ geom_jitter(aes(color=stand,fill=stand),size=2.5,alpha=0.2,width=0.2)+stat_boxplot(geom ='errorbar',width=0.25)+ geom_boxplot(aes(fill=stand),color="black",outlier.shape=NA) +scale_color_manual(values=c("#3D399F","#B94DAC","#F9A669"))+ scale_fill_manual(values=c("#3D399F","#B94DAC","#F9A669"))  + labs(x = "Location", y="Turbidity (NTU)")+ theme_bw(12)+ theme(legend.position="none")+ theme(panel.grid.minor = element_blank()) +theme(panel.grid.major = element_blank())+scale_y_continuous(limits=c(0, 100),breaks=seq(0, 100, 10))+annotate("text", x="Reef", y=6, label= "a", size=5,hjust=-4)+annotate("text", x="Inside", y=9, label= "b", size=5,hjust=-4)+annotate("text", x="Outside", y=8, label= "c", size=5,hjust=-4)

## Figure 8 (Violin plots day-night; right panels) ####
do.violin<-ggplot(hourly.abiotic,aes(x=time.of.day, y=do))+ geom_violin(aes(fill=time.of.day),color="black",outlier.shape=NA)+geom_boxplot(width=0.2,outlier.shape=NA)+ facet_grid(~stand)+scale_color_manual(values=c("#00AFBB","#FC4E07"))+ scale_fill_manual(values=c("#00AFBB","#FC4E07"))  + labs(x = "Location", y="Dissolved oxygen (mg/L)")+ theme_bw(12)+ theme(legend.position="none")+ theme(panel.grid.minor = element_blank(),axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()) +theme(panel.grid.major = element_blank())+scale_y_continuous(limits=c(5.5, 15.5),breaks=seq(5.5, 15.5, 1))

ph.violin<-ggplot(hourly.abiotic,aes(x=time.of.day, y=sw.ph))+ geom_violin(aes(fill=time.of.day),color="black",outlier.shape=NA)+geom_boxplot(width=0.2,outlier.shape=NA)+ facet_grid(~stand)+scale_color_manual(values=c("#00AFBB","#FC4E07"))+ scale_fill_manual(values=c("#00AFBB","#FC4E07"))  + labs(x = "Location", y="pH")+ theme_bw(12)+ theme(legend.position="none")+ theme(panel.grid.minor = element_blank()) +theme(panel.grid.major = element_blank(),axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank())+scale_y_continuous(limits=c(7.7, 8.5),breaks=seq(7.7, 8.5, 0.1))

temp.violin<-ggplot(hourly.abiotic,aes(x=time.of.day, y=temp))+ geom_violin(aes(fill=time.of.day),color="black",outlier.shape=NA)+geom_boxplot(width=0.2,outlier.shape=NA)+ facet_grid(~stand)+scale_color_manual(values=c("#00AFBB","#FC4E07"))+ scale_fill_manual(values=c("#00AFBB","#FC4E07"))  + labs(x = "Location", y="Temperature (°C)")+ theme_bw(12)+ theme(legend.position="none")+ theme(panel.grid.minor = element_blank()) +theme(panel.grid.major = element_blank(),axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank())+scale_y_continuous(limits=c(18, 28),breaks=seq(18, 28, 1))

turb.violin<-ggplot(partial.turb,aes(x=time.of.day, y=turbidity))+ geom_violin(aes(fill=time.of.day),color="black",outlier.shape=NA)+geom_boxplot(width=0.2,outlier.shape=NA)+ facet_grid(~stand)+scale_color_manual(values=c("#00AFBB","#FC4E07"))+ scale_fill_manual(values=c("#00AFBB","#FC4E07"))  + labs(x = "Time of day", y="Turbidity (NTU)")+ theme_bw(12)+ theme(legend.position="none")+ theme(panel.grid.minor = element_blank()) +theme(panel.grid.major = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank())+scale_y_continuous(limits=c(0, 100),breaks=seq(0, 100, 10))

#Combine plots
combo_boxplot_violin_abiotic<-(do.boxplot+do.violin)/(ph.boxplot+ph.violin)/(temp.boxplot+temp.violin)/(turb.boxplot+turb.violin)
combo_boxplot_violin_abiotic

#Export figure
tiff("Abiotic boxplot violin graphs_JUNE 2025_2.tiff", width = 11, height = 13, units = "in", res = 800)
combo_boxplot_violin_abiotic
dev.off()

pdf("Abiotic boxplot violin graphs_JUNE 2025_2.pdf", width = 11, height = 13)
combo_boxplot_violin_abiotic
dev.off()
