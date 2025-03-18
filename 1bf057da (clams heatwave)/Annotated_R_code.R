#########################################################################################################################################################
#########################################################################################################################################################
---------------------------------------------------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------- SUPPLEMENTARY CODE FOR CLEMENTS ET AL. ---------------------------------------------------------
------- FISHING DURING EXTREME HEATWAVES ALTERS ECOLOGICAL INTERACTIONS AND INCREASES INDIRECT FISHING MORTALITY IN A UBIQUITOUS NEARSHORE SYSTEM -------
---------------------------------------------------------------- COMMUNICATIONS BIOLOGY -----------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------------------------------------------------
#########################################################################################################################################################
#########################################################################################################################################################

#### PACKAGES ####
library(devtools)
library(ggplot2)
library(patchwork)
library(lmerTest)
library(car)
library(forcats)
library(tidyverse)
library(dplyr)
library(socviz)
library(viridis)
library(MASS)
library(emmeans)
library(betareg)
library(lmtest)
library(blme)
library(brglm2)
library(logistf)
library(jpeg)
library(RCurl)
library(grid)
library(Rmisc)
library(mgcv)

#### STATISTICAL ANALYSIS: EFFECTS OF EXPERIMENT (MONTH), TIME SINCE RELEASE, TIDE LEVEL, AND PREDATOR TREATMENT ON BURROWING AND MORTALITY####
#Upload datafile for burrowing and mortality after 24 and 48 hours
burrowing.mortality<-read.csv(file.choose()) #datafile = 'S-01 data.csv'
attach(burrowing.mortality)
summary(burrowing.mortality)
burrowing.mortality$plot.glob<-as.factor(burrowing.mortality$plot.glob)
burrowing.mortality$month<-as.factor(burrowing.mortality$month)
burrowing.mortality$time.since.deploy<-as.factor(burrowing.mortality$time.since.deploy)
burrowing.mortality$tide.level<-as.factor(burrowing.mortality$tide.level)
burrowing.mortality$pred.treat<-as.factor(burrowing.mortality$pred.treat)
burrowing.mortality$avg.temp<-as.numeric(burrowing.mortality$avg.temp)
burrowing.mortality$julian.date<-as.numeric(burrowing.mortality$julian.date)

#Build BGLMER models for burrowing and mortality at 24 and 48 hours at each tide level, for each predator treatment, in each experiment (month). Include plot to account for repeated measures over two time points. Note that there is complete separation in the data because some treatment levels have all 0 values. As such, use Bayesian GLMER approach [bglmer() as per Bolker (2018)), specifying zero-mean normal priors to account for separation. 
#We initially tried variance prior of 9 but this resulted in an error for the mortalitymod, most likely due to quasi-complete separation. We thus used a slightly larger variance prior, which resolved the issue with reasonable model estimates. Diagonal matrix value of 60 was chosen to match the number of terms in the models. 
bglmer_burrowmod<-bglmer(prop.burrowed~month*pred.treat*time.since.deploy*tide.level+(1|plot.glob),family=binomial,data=burrowing.mortality, fixef.prior = normal(cov = diag(9,60))) #Ignore convergence warnings; model estimates still OK
bglmer_mortalitymod<-bglmer(prop.dead~month*pred.treat*time.since.deploy*tide.level+(1|plot.glob),family=binomial,data=burrowing.mortality, fixef.prior = normal(cov = diag(9,60))) #Error; slightly increase variance prior value

bglmer_burrowmod2<-bglmer(prop.burrowed~month*pred.treat*time.since.deploy*tide.level+(1|plot.glob),family=binomial,data=burrowing.mortality, fixef.prior = normal(cov = diag(10,60))) #Ignore convergence warnings; model estimates still OK
bglmer_mortalitymod2<-bglmer(prop.dead~month*pred.treat*time.since.deploy*tide.level+(1|plot.glob),family=binomial,data=burrowing.mortality, fixef.prior = normal(cov = diag(10,60))) #Error; slightly increase variance prior value

#Get results
Anova(bglmer_burrowmod2,type=3)#significant effect of experiment (month) on proportion of clams burrowed
Anova(bglmer_mortalitymod2,type=3) #Marginally non-significant experiment x predator treatment x time effect. Interpret 3-way interactive effect

#Get pairwise results
#Burrowing model
burrow.pairwise<-emmeans(bglmer_burrowmod,~month,adjustment="Holm") #differences between experiments; pooled across pred treats, times, and tide levels
pairs(burrow.pairwise) 

#Mortality model
mortality.pairwise1<-emmeans(bglmer_mortalitymod,~month|pred.treat|time.since.deploy,data=burrowing.mortality,adjustment="Holm") #differences between experiments for each pred treat and time; pooled across tide levels
pairs(mortality.pairwise1) 

mortality.pairwise2<-emmeans(bglmer_mortalitymod,~pred.treat|month|time.since.deploy,data=burrowing.mortality,adjustment="Holm") #differences between pred treats for each experiment and time; pooled across tide levels
pairs(mortality.pairwise2)

mortality.pairwise3<-emmeans(bglmer_mortalitymod,~time.since.deploy|month|pred.treat,data=burrowing.mortality,adjustment="Holm") #differences between times for each experiment and pred treat ; pooled across tide levels
pairs(mortality.pairwise3)

## GAM analysis for robustness
#Build GAMMs and get output for reburrowing and mortality to look at effects of average temperature (i.e., heatwave effect), time since deployment for each predator treatment (to determine if burrowing/mortality within each predator treatment differed between 24 and 48h), predator treatment
burrow.gam <- gamm(prop.burrowed~s(avg.temp,k=3)+pred.treat*time.since.deploy*tide.level+s(julian.date,bs='re'),data=burrowing.mortality, family=nb, method="REML")
anova(burrow.gam$gam)

mort.gam <- gamm(prop.dead~s(avg.temp,k=3)+pred.treat*time.since.deploy*tide.level+s(julian.date,bs='re'),data=burrowing.mortality, family=nb, method="REML")
anova(mort.gam$gam)


#### MAIN FIGURES ####
#### Figure 1 ####
#### Air temp time series
#Upload datafile 
temp.ts<-read.csv(file.choose()) #datafile = 'S-02 data.csv'
attach(temp.ts)
summary(temp.ts)

#Max air temp
temp.ts1<-ggplot(data=temp.ts,aes(x=julian.date, y=max.temp)) + geom_line(aes(color=max.temp)) + geom_point(aes(fill=max.temp,color=max.temp), size=2)+  scale_color_gradient(low = "gray", high = "#ff4500",name="Max temp (°C)")+  scale_fill_gradient(low = "gray", high = "#ff4500",name="Max temp (°C)")+theme_bw(7) + theme(legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.title.x=element_blank(),axis.text.x = element_blank())+scale_y_continuous(limits=c(0,40),breaks=seq(0, 40, 5))+scale_x_continuous(limits=c(122,263),breaks=seq(122,263, 5))+ labs(x = "Julian date", y="Max air temp (°C)")+annotate(geom = "rect", xmin = 142, xmax = 146,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4)+annotate(geom = "rect", xmin = 171, xmax = 175,ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.2)+annotate(geom = "rect", xmin = 184, xmax = 188,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4)+annotate(geom = "rect", xmin = 231, xmax = 235,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4)+annotate(geom = "rect", xmin = 259, xmax = 263,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4)+annotate(geom = "text", x = 144, y = 0, label = "May", color = "black",angle = 90,hjust=0,size=2)+annotate(geom = "text", x = 173, y = 0, label = "June", color = "black",angle = 90,hjust=0,size=2)+annotate(geom = "text", x = 186, y = 0, label = "July", color = "black",angle = 90,hjust=0,size=2)+annotate(geom = "text", x = 233, y = 0, label = "Aug", color = "black",angle = 90,hjust=0,size=2)+annotate(geom = "text", x = 261, y = 0, label = "Sept", color = "black",angle = 90,hjust=0,size=2)

#Mean air temp
temp.ts2<-ggplot(data=temp.ts,aes(x=julian.date, y=mean.temp)) + geom_line(aes(color=mean.temp)) + geom_point(aes(fill=mean.temp,color=mean.temp), size=2)+  scale_color_gradient(low = "lightgray", high = "#ff4500",name="Mean temp (°C)")+  scale_fill_gradient(low = "lightgray", high = "#ff4500",name="Mean temp (°C)")+theme_bw(7) + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),legend.position = "none")+scale_y_continuous(limits=c(0,30),breaks=seq(0, 30, 5))+scale_x_continuous(limits=c(122,263),breaks=seq(122,263, 5))+ labs(x = "Julian date", y="Mean air temp (°C)")+annotate(geom = "rect", xmin = 142, xmax = 146,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4)+annotate(geom = "rect", xmin = 171, xmax = 175,ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.2)+annotate(geom = "rect", xmin = 184, xmax = 188,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4)+annotate(geom = "rect", xmin = 231, xmax = 235,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4)+annotate(geom = "rect", xmin = 259, xmax = 263,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4)+annotate(geom = "text", x = 144, y = 0, label = "May", color = "black",angle = 90,hjust=0,size=2)+annotate(geom = "text", x = 173, y = 0, label = "June", color = "black",angle = 90,hjust=0,size=2)+annotate(geom = "text", x = 186, y = 0, label = "July", color = "black",angle = 90,hjust=0,size=2)+annotate(geom = "text", x = 233, y = 0, label = "Aug", color = "black",angle = 90,hjust=0,size=2)+annotate(geom = "text", x = 261, y = 0, label = "Sept", color = "black",angle = 90,hjust=0,size=2)

#### Temp during clam digging
#Upload datafile 
temp.dig<-read.csv(file.choose()) #datafile = 'S-03 data.csv'
attach(temp.dig)
summary(temp.dig)

#First, reorder "month" variable so that months are ordered chrionologically
temp.dig$month<-factor(temp.dig$month,levels=c("May","June","July","Aug","Sept"))

#Next, generate summary statistic values to apply in plots. NOTE: need to run base code at the bottom of this .R file prior to running summarySE() commands below.
temp.month <- summarySE(temp.dig, measurevar="temp", groupvars=c("month"))
humidex.month <- summarySE(temp.dig, measurevar="humidex", groupvars=c("month"))

#Create avg. temp and humidex data columns for each month to colour geom_jitter and geom_points properly
temp.dig$avg.temp <- c(25.025,25.025,25.025,25.025,32.575,32.575,32.575,32.575,28.825,28.825,28.825,28.825,27.375,27.375,27.375,27.375,28.625,28.625,28.625,28.625)
temp.dig$avg.humidex <- c(28.5,28.5,28.5,28.5,39.75,39.75,39.75,39.75,30.25,30.25,30.25,30.25,35,35,35,35,33.75,33.75,33.75,33.75)

#Temp
temp.plot<-ggplot(temp.dig,aes(x=month, y=temp))+geom_jitter(aes(color=avg.temp,fill=avg.temp),position = position_jitter(width=0.2,height=0),alpha=0.3,size=1)+geom_errorbar(data=temp.month,aes(ymin=temp-sd, ymax=temp+sd,color=temp),width=.2)+ geom_point(data=temp.month,aes(fill=temp,color=temp),size=2)+  scale_color_gradient(low = "gray", high = "#ff4500",name="Temp (°C)")+  scale_fill_gradient(low = "gray", high = "#ff4500",name="Temp (°C)")+scale_y_continuous(limits=c(18,42),breaks=seq(20,40,5))+ labs(x = "Trial", y="Air temperature (°C)")+ theme_bw(7)+ theme(legend.position="none")+ theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())

#Humidex
humidex.plot<-ggplot(temp.dig,aes(x=month, y=humidex))+geom_jitter(aes(color=avg.humidex,fill=avg.humidex),position = position_jitter(width=0.3,height=0),alpha=0.2,size=1)+geom_errorbar(data=humidex.month,aes(ymin=humidex-sd, ymax=humidex+sd,color=humidex),width=.2)+ geom_point(data=humidex.month,aes(fill=humidex,color=humidex),size=2)+  scale_color_gradient(low = "gray", high = "#ff4500",name="Humidex (°C)")+  scale_fill_gradient(low = "gray", high = "#ff4500",name="Humidex (°C)")+scale_y_continuous(limits=c(18,42),breaks=seq(20,40,5))+ labs(x = "Trial", y="Humidex (°C)")+ theme_bw(7)+ theme(legend.position="none")+ theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())

#Combine time series and experiment plots
temp.ts.full<-temp.ts1/temp.ts2/(temp.plot+humidex.plot)
temp.ts.full

#Export full figure
tiff("R2_FSCP Fig 1.tiff", width = 180, height = 200, units = "mm", res = 800)
temp.ts.full
dev.off()

#### Figure 2 ####
#Use same dataset as used in the statistical analyses above ("burrowing.mortality"; from datafile 'S-01 data.csv')
#Create plots
#First, reorder "month" variable so that months are ordered chronologically
burrowing.mortality$month<-factor(burrowing.mortality$month,levels=c("May","June","July","Aug","Sept"))

#Next, generate summary statistic values to apply in plots. NOTE: need to run base code at the bottom of this .R file prior to running summarySE() commands below.
burrow.sum <- summarySE(burrowing.mortality, measurevar="prop.burrowed", groupvars=c("month","time.since.deploy","pred.treat","avg.temp"))
mort.sum <- summarySE(burrowing.mortality, measurevar="prop.dead", groupvars=c("month","time.since.deploy","pred.treat","avg.temp"))

#Generate and view base plot for proportion reburrowed
burrow.plot<-ggplot(burrowing.mortality,aes(x=month, y=prop.burrowed))+geom_jitter(aes(color=avg.temp,fill=avg.temp),position = position_jitter(width=0.2,height=0),alpha=0.2,size=0.5)+geom_errorbar(data=burrow.sum,aes(ymin=prop.burrowed-sd, ymax=prop.burrowed+sd,color=avg.temp),width=.2) + geom_point(data=burrow.sum,aes(fill=avg.temp,color=avg.temp),size=2) +  scale_color_gradient(low = "gray", high = "#ff4500",name="Mean air temp (°C)")+  scale_fill_gradient(low = "gray", high = "#ff4500",name="Mean air (°C)")+scale_y_continuous(limits=c(-0.15,1.15),breaks=seq(0, 1, 0.2))+ facet_grid(~pred.treat~time.since.deploy) + labs(x = "Trial", y="Proportion burrowed")+ theme_bw(7)+ theme(legend.position="none")+ theme(panel.grid = element_blank())
burrow.plot

#Create data frames for adding text for pairwise results
b.may.pe.24 <- data.frame(prop.burrowed = 1.03,month = "May",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
b.june.pe.24 <- data.frame(prop.burrowed = 0.48,month = "June",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
b.july.pe.24 <- data.frame(prop.burrowed = 1.1,month = "July",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
b.aug.pe.24 <- data.frame(prop.burrowed = 1.07,month = "Aug",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
b.sept.pe.24 <- data.frame(prop.burrowed = 1.1,month = "Sept",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
b.may.pi.24 <- data.frame(prop.burrowed = 1.05,month = "May",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
b.june.pi.24 <- data.frame(prop.burrowed = 0.12,month = "June",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
b.july.pi.24 <- data.frame(prop.burrowed = 1.11,month = "July",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
b.aug.pi.24 <- data.frame(prop.burrowed = 1.08,month = "Aug",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
b.sept.pi.24 <- data.frame(prop.burrowed = 1.13,month = "Sept",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
b.may.pe.48 <- data.frame(prop.burrowed = 1.06,month = "May",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
b.june.pe.48 <- data.frame(prop.burrowed = 0.31,month = "June",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
b.july.pe.48 <- data.frame(prop.burrowed = 1.1,month = "July",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
b.aug.pe.48 <- data.frame(prop.burrowed = 1.12,month = "Aug",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
b.sept.pe.48 <- data.frame(prop.burrowed = 1.12,month = "Sept",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
b.may.pi.48 <- data.frame(prop.burrowed = 1.08,month = "May",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
b.june.pi.48 <- data.frame(prop.burrowed = 0.09,month = "June",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
b.july.pi.48 <- data.frame(prop.burrowed = 1.11,month = "July",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
b.aug.pi.48 <- data.frame(prop.burrowed = 0.99,month = "Aug",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
b.sept.pi.48 <- data.frame(prop.burrowed = 1.11,month = "Sept",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))

#Generate and view final plot for reburrowing with pairwise labels
burrow.plot.pairwise<-burrow.plot+geom_text(data = b.may.pe.24,label = "A",size=2)+geom_text(data = b.june.pe.24,label = "B",size=2)+geom_text(data = b.july.pe.24,label = "A",size=2)+geom_text(data = b.aug.pe.24,label = "A",size=2)+geom_text(data = b.sept.pe.24,label = "A",size=2)+geom_text(data = b.may.pi.24,label = "A",size=2)+geom_text(data = b.june.pi.24,label = "B",size=2)+geom_text(data = b.july.pi.24,label = "A",size=2)+geom_text(data = b.aug.pi.24,label = "A",size=2)+geom_text(data = b.sept.pi.24,label = "A",size=2)+geom_text(data = b.may.pe.48,label = "A",size=2)+geom_text(data = b.june.pe.48,label = "B",size=2)+geom_text(data = b.july.pe.48,label = "A",size=2)+geom_text(data = b.aug.pe.48,label = "A",size=2)+geom_text(data = b.sept.pe.48,label = "A",size=2)+geom_text(data = b.may.pi.48,label = "A",size=2)+geom_text(data = b.june.pi.48,label = "B",size=2)+geom_text(data = b.july.pi.48,label = "A",size=2)+geom_text(data = b.aug.pi.48,label = "A",size=2)+geom_text(data = b.sept.pi.48,label = "A",size=2)
burrow.plot.pairwise

#Generate and view base plot for proportion dead
mort.plot<-ggplot(data=burrowing.mortality,aes(x=month, y=prop.dead))+geom_jitter(aes(color=avg.temp,fill=avg.temp),position = position_jitter(width=0.2,height=0),alpha=0.2,size=0.5)+geom_errorbar(data=mort.sum,aes(ymin=prop.dead-sd, ymax=prop.dead+sd,color=avg.temp),width=.2) + geom_point(data=mort.sum,aes(fill=avg.temp,color=avg.temp),size=2)  +  scale_color_gradient(low = "gray", high = "#ff4500",name="Mean air temp (°C)")+  scale_fill_gradient(low = "gray", high = "#ff4500",name="Mean air (°C)")+scale_y_continuous(limits=c(-0.15,1.15),breaks=seq(0, 1, 0.2))+ facet_grid(~pred.treat~time.since.deploy) + labs(x = "Trial", y="Proportion dead")+ theme_bw(7)+ theme(legend.position="top")+ theme(panel.grid = element_blank())+guides(fill="none")+  theme(legend.key.size = unit(4, 'mm'),legend.text=element_text(size=4),legend.title=element_text(size=5))
mort.plot

#Create data frames for adding text for pairwise results
m.may.pe.24 <- data.frame(prop.dead = 0.12,month = "May",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.june.pe.24 <- data.frame(prop.dead = 0.58,month = "June",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.july.pe.24 <- data.frame(prop.dead = 0.12,month = "July",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.aug.pe.24 <- data.frame(prop.dead = 0.09,month = "Aug",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.sept.pe.24 <- data.frame(prop.dead = 0.09,month = "Sept",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.may.pi.24 <- data.frame(prop.dead = 0.09,month = "May",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.june.pi.24 <- data.frame(prop.dead = 1.06,month = "June",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.july.pi.24 <- data.frame(prop.dead = 0.26,month = "July",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.aug.pi.24 <- data.frame(prop.dead = 0.24,month = "Aug",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.sept.pi.24 <- data.frame(prop.dead = 0.33,month = "Sept",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.may.pe.48 <- data.frame(prop.dead = 0.12,month = "May",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.june.pe.48 <- data.frame(prop.dead = 1.06,month = "June",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.july.pe.48 <- data.frame(prop.dead = 0.15,month = "July",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.aug.pe.48 <- data.frame(prop.dead = 0.15,month = "Aug",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.sept.pe.48 <- data.frame(prop.dead = 0.22,month = "Sept",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.may.pi.48 <- data.frame(prop.dead = 0.09,month = "May",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.june.pi.48 <- data.frame(prop.dead = 1.09,month = "June",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.july.pi.48 <- data.frame(prop.dead = 0.35,month = "July",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.aug.pi.48 <- data.frame(prop.dead = 0.52,month = "Aug",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.sept.pi.48 <- data.frame(prop.dead = 0.4,month = "Sept",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.june.24.pe.star <- data.frame(prop.dead = 0.29,month = "June",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.june.24.pi.star <- data.frame(prop.dead = 0.82,month = "June",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.aug.24.pe.star <- data.frame(prop.dead = 0.05,month = "Aug",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.aug.24.pi.star <- data.frame(prop.dead = 0.12,month = "Aug",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.sept.24.pe.star <- data.frame(prop.dead = 0.05,month = "Sept",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.sept.24.pi.star <- data.frame(prop.dead = 0.14,month = "Sept",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.july.48.pe.star <- data.frame(prop.dead = 0.07,month = "July",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.july.48.pi.star <- data.frame(prop.dead = 0.16,month = "July",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.aug.48.pe.star <- data.frame(prop.dead = 0.07,month = "Aug",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.aug.48.pi.star <- data.frame(prop.dead = 0.3,month = "Aug",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.sept.48.pe.star <- data.frame(prop.dead = 0.09,month = "Sept",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.sept.48.pi.star <- data.frame(prop.dead = 0.2,month = "Sept",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.june.24.pe.plus <- data.frame(prop.dead = 0.2,month = "June",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.june.48.pe.plus <- data.frame(prop.dead = 0.77,month = "June",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PE",levels = c("PE","PI")))
m.aug.24.pi.plus <- data.frame(prop.dead = 0.02,month = "Aug",lab = "Text", time.since.deploy = factor("24 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))
m.aug.48.pi.plus <- data.frame(prop.dead = 0.2,month = "Aug",lab = "Text", time.since.deploy = factor("48 h",levels = c("24 h","48 h")), pred.treat = factor("PI",levels = c("PE","PI")))

#Build plot with pairwise labels
mort.plot.pairwise<-mort.plot+geom_text(data = m.may.pe.24,label = "AB",size=2)+geom_text(data = m.june.pe.24,label = "B",size=2)+geom_text(data = m.july.pe.24,label = "A",size=2)+geom_text(data = m.aug.pe.24,label = "A",size=2)+geom_text(data = m.sept.pe.24,label = "A",size=2)+geom_text(data = m.may.pi.24,label = "A",size=2)+geom_text(data = m.june.pi.24,label = "B",size=2)+geom_text(data = m.july.pi.24,label = "A",size=2)+geom_text(data = m.aug.pi.24,label = "A",size=2)+geom_text(data = m.sept.pi.24,label = "A",size=2)+geom_text(data = m.may.pe.48,label = "A",size=2)+geom_text(data = m.june.pe.48,label = "B",size=2)+geom_text(data = m.july.pe.48,label = "A",size=2)+geom_text(data = m.aug.pe.48,label = "A",size=2)+geom_text(data = m.sept.pe.48,label = "C",size=2)+geom_text(data = m.may.pi.48,label = "A",size=2)+geom_text(data = m.june.pi.48,label = "B",size=2)+geom_text(data = m.july.pi.48,label = "AB",size=2)+geom_text(data = m.aug.pi.48,label = "AB",size=2)+geom_text(data = m.sept.pi.48,label = "AB",size=2)+geom_text(data = m.june.24.pe.plus,label = "+",hjust=-0.8,size=2)+geom_text(data = m.june.48.pe.plus,label = "+",hjust=-0.8,size=2)+geom_text(data = m.aug.24.pi.plus,label = "+",hjust=-0.8,size=2)+geom_text(data = m.aug.48.pi.plus,label = "+",hjust=-0.8,size=2)+geom_text(data = m.june.24.pe.star,label = "*",size=3,hjust=-0.8)+geom_text(data = m.june.24.pi.star,label = "*",size=3,hjust=-0.8)+geom_text(data = m.aug.24.pe.star,label = "*",size=3,hjust=-0.8)+geom_text(data = m.aug.24.pi.star,label = "*",size=3,hjust=-0.8)+geom_text(data = m.sept.24.pe.star,label = "*",size=3,hjust=-0.8)+geom_text(data = m.sept.24.pi.star,label = "*",size=3,hjust=-0.8)+geom_text(data = m.july.48.pe.star,label = "*",size=3,hjust=-0.8)+geom_text(data = m.july.48.pi.star,label = "*",size=3,hjust=-0.8)+geom_text(data = m.aug.48.pe.star,label = "*",size=3,hjust=-0.8)+geom_text(data = m.aug.48.pi.star,label = "*",size=3,hjust=-0.8)+geom_text(data = m.sept.48.pe.star,label = "*",size=3,hjust=-0.8)+geom_text(data = m.sept.48.pi.star,label = "*",size=3,hjust=-0.8)
mort.plot.pairwise

fscp.burrowing.mort.combo.horiz<-burrow.plot.pairwise+mort.plot.pairwise
fscp.burrowing.mort.combo.horiz

#Export final plot
tiff("R1_FSCP Fig 2_jitter.tiff", width = 180, height = 90, units = "mm", res = 800)
fscp.burrowing.mort.combo.horiz
dev.off()

#NOTE: Clam images added in Powerpoint and final figure exported as high res TIF file

#### Figure 3 ####
pred<-read.csv(file.choose()) #datafile = 'S-04 data.csv'
attach(pred)
summary(pred)

#First, reorder "month" variable so that months are ordered chrionologically
pred$month<-factor(pred$month,levels=c("May","June","July","Aug","Sept"))

#Create total predator index variable
pred$pred.index<-pred$crab.count+pred$mudsnail.buckets.count

#Create crab plot
crab.barplot<-ggplot(data=pred, aes(x=month, y=crab.count)) + geom_bar(stat="identity", aes(fill=avg.air.temp),color="black",size=0.2) +  scale_fill_gradient(low = "gray", high = "#ff4500",name="Mean air temp (°C)")+theme_bw(5)+scale_y_continuous(limits=c(0,16),breaks=seq(0, 16, 2))+ theme(panel.grid.minor = element_blank(),legend.position = "top",axis.title.x=element_blank(),axis.text.x=element_blank(),panel.grid=element_blank())+ labs(x = "Trial", y="Number of crabs in mesocosms")+  theme(legend.key.size = unit(4, 'mm'),legend.text=element_text(size=4),legend.title=element_text(size=5))
crab.barplot

#Create mudsnail plot
mudsnail.barplot<-ggplot(data=pred, aes(x=month, y=mudsnail.buckets.count)) + geom_bar(stat="identity", aes(fill=avg.air.temp),color="black",size=0.2) +  scale_fill_gradient(low = "gray", high = "#ff4500",name="Mean air temp (°C)")+theme_bw(5)+scale_y_continuous(limits=c(0,26),breaks=seq(0, 24, 4))+ theme(panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),axis.text.x=element_blank(),panel.grid = element_blank())+ labs(x = "Trial", y="Number of mesocosms with mudsnails")+  theme(legend.key.size = unit(4, 'mm'),legend.text=element_text(size=4),legend.title=element_text(size=5))
mudsnail.barplot

#Create predator activity index plot (# crabs + number mesocosms with mudsnails)
index.barplot<-ggplot(data=pred, aes(x=month, y=pred.index)) + geom_bar(stat="identity", aes(fill=avg.air.temp),color="black",size=0.2) +  scale_fill_gradient(low = "gray", high = "#ff4500",name="Mean air temp (°C)")+theme_bw(5)+scale_y_continuous(limits=c(0,42),breaks=seq(0, 40, 5))+ theme(panel.grid.minor = element_blank(),legend.position = "none",panel.grid = element_blank())+ labs(x = "Trial", y="Predator activity index")+  theme(legend.key.size = unit(4, 'mm'),legend.text=element_text(size=4),legend.title=element_text(size=5))
index.barplot

#Combine predatpr activity plots
pred.plot<-crab.barplot/mudsnail.barplot/index.barplot
pred.plot

#Export plot
tiff("R1_FSCP Fig 3.tiff", width = 58, height = 130, units = "mm", res = 800)
pred.plot
dev.off()

#NOTE: Crab and mudsnail images added in Powerpoint and final figure exported as high res TIF file

#### Figure 4 ####
#Upload datafile for burrowing and mortality
temp.curve<-read.csv(file.choose()) #datafile = 'S-05 data.csv'
attach(temp.curve)
summary(temp.curve)

#Datafile for predator activity is the same as for Fugre 3 (data frame "pred", from datafile'S4 data.csv')

#Generate and view plot
temp.curve.final<-ggplot()+geom_smooth(data=temp.curve,aes(x=exp.temp, y=proportion,group=metric,color=metric,fill=metric), method="loess",span=1, se = TRUE,alpha=0.2,size=0.6)+geom_smooth(data=pred,aes(x=avg.air.temp, y=rel.pred.activity,color="#BE562F", fill="#BE562F"), method="loess",span=1, se = FALSE,alpha=0.2,linetype = "dashed",size=0.6)+scale_x_continuous(limits=c(24,33),breaks=seq(24,33,1))+scale_y_continuous(limits=c(-0.05,1.05),breaks=seq(0, 1, 0.2))+ scale_color_manual(values = c("#BE562F","#5B90B4","#E59A56"))+ scale_fill_manual(values = c("#BE562F","#5B90B4","#E59A56"))+labs(y="Proportion",x="Air temperature during fishing (°C)")+theme_bw(7)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(legend.position = "none",legend.title = element_blank())+annotate(geom = "text", x = 26, y = 0.92, label = "Reburrowing", color = "#5B90B4",hjust=0,size=2)+annotate(geom = "text", x = 26, y = 0, label = "Mortality", color = "#E59A56",hjust=0,size=2)+annotate(geom = "text", x = 25, y = 0.3, label = "Relative predator activity", color = "#BE562F",hjust=0,size=2)
temp.curve.final
  
#Export plot
tiff("FSCP Fig 4.tiff", width = 160, height = 75, units = "mm", res = 1600)
temp.curve.final
dev.off()

#Add images to exported figure using Adobe Photoshop


#### SUPPLEMENTARY FIGURES ####
#### Figure S1 ####
#Upload datafile for 2021 and 2024 reburrowing proportions
crab.trapping<-read.csv(file.choose()) #datafile = 'S-06 data.csv'
attach(crab.trapping)
summary(crab.trapping)
crab.trapping$julian.date<-as.numeric(crab.trapping$julian.date)
crab.trapping$year<-as.factor(crab.trapping$year)

#Generate and view base plot for proportion reburrowed
crab.cpue.plot<-ggplot(crab.trapping,aes(x=julian.date, y=avg.cpue))+annotate(geom = "rect", xmin = 142, xmax = 146,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4)+annotate(geom = "rect", xmin = 171, xmax = 175,ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.2)+annotate(geom = "rect", xmin = 184, xmax = 188,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4)+annotate(geom = "rect", xmin = 231, xmax = 235,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4)+annotate(geom = "rect", xmin = 259, xmax = 263,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.4)+annotate(geom = "text", x = 144, y = 3, label = "May", color = "black",angle = 90,hjust=1,vjust=0.5,size=3.5)+annotate(geom = "text", x = 173, y = 3, label = "June", color = "black",angle = 90,hjust=1,vjust=0.5,size=3.5)+annotate(geom = "text", x = 186, y = 3, label = "July", color = "black",angle = 90,hjust=1,vjust=0.5,size=3.5)+annotate(geom = "text", x = 233, y = 3, label = "Aug", color = "black",angle = 90,hjust=1,vjust=0.5,size=3.5)+annotate(geom = "text", x = 261, y = 3, label = "Sept", color = "black",angle = 90,hjust=1,vjust=0.5,size=3.5)+ geom_line(aes(x=julian.date,y=avg.cpue,group=year,color=year),size=0.75)+ geom_point(aes(x=julian.date,y=avg.cpue,group=year,fill=year,shape=year),size=5,color="black")+ scale_fill_manual(values = c("lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","#FF4500"))+ scale_color_manual(values = c("lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","#FF4500"))+scale_shape_manual(values=c(22,23,24,25,21,21))+scale_y_continuous(limits=c(0,3),breaks=seq(0, 3, 0.5))+scale_x_continuous(limits=c(155,270),breaks=seq(155, 270, 5)) + labs(x = "Julian date", y="Crabs per day (all species)")+ theme_bw(12)+ theme(legend.position=c(0.08,0.75))+ theme(panel.grid = element_blank())+ guides(fill=guide_legend(title="Year"),color=guide_legend(title="Year"),shape=guide_legend(title="Year"))

#Export plot
tiff("R1_FSCP Fig S1.tiff", width = 10, height = 5, units = "in", res = 800)
crab.cpue.plot
dev.off()

#### Figure S2 ####
#Upload datafile for 2021 and 2024 reburrowing proportions
burrow.compare<-read.csv(file.choose()) #datafile = 'S-07 data.csv'
attach(burrow.compare)
summary(burrow.compare)
burrow.compare$day<-as.numeric(burrow.compare$day)
burrow.compare$year<-as.factor(burrow.compare$year)
burrow.compare$site<-as.factor(burrow.compare$site)
burrow.compare$bucket<-as.factor(burrow.compare$bucket)

#Next, generate summary statistic values to apply in plots. NOTE: need to run base code at the bottom of this .R file prior to running summarySE() commands below.
burrow.comp.sum <- summarySE(burrow.compare, measurevar="prop.burrowed", groupvars=c("day","year"))

#Generate and view base plot for proportion reburrowed
burrow.comp.plot<-ggplot(burrow.comp.sum,aes(x=day, y=prop.burrowed))+geom_errorbar(aes(ymin=prop.burrowed-se, ymax=prop.burrowed+se,color=year),width=2)+annotate(geom = "rect", xmin = 170, xmax = 175,ymin = -Inf, ymax = Inf, fill = "lightgrey", alpha = 0.4)+ geom_point(aes(fill=year,color=year),size=5)+ geom_smooth(aes(x=day,y=prop.burrowed,group=year,color=year),se=FALSE,method='loess',span=0.5,size=0.75,linetype = "dashed")+scale_fill_manual(values = c("darkgrey", "#ff4500"))+ scale_color_manual(values = c("darkgrey", "#ff4500"))+scale_y_continuous(limits=c(0,1),breaks=seq(0, 1, 0.2))+scale_x_continuous(limits=c(135,280),breaks=seq(135, 280, 5)) + labs(x = "Julian Date", y="Proportion burrowed")+ theme_bw(12)+ theme(legend.position=c(0.9,0.14))+ theme(panel.grid = element_blank())+ guides(fill=guide_legend(title="Year"),color=guide_legend(title="Year"))+annotate(geom = "text", x = 170, y = 0.089, label = "June 21", color = "#FF4500",hjust=1,size=4,fontface=2)+annotate(geom = "text", x = 171, y = 0.79, label = "June 22", color = "black",hjust=1,size=4,,fontface=2)+annotate(geom = "text", x = 174, y = 0.08, label = "Heatwave", color = "#FF4500",hjust=0,size=4,fontface=3)+annotate(geom = "text", x = 269, y = 0.45, label = "Low salinity", color = "black",hjust=1,size=4,fontface=3)
burrow.comp.plot

#Export plot
tiff("R1_FSCP Fig S2.tiff", width = 10, height = 6, units = "in", res = 800)
burrow.comp.plot
dev.off()

#### Figure S3 ####
#Upload datafile for June precipitation
precipitation<-read.csv(file.choose()) #datafile = 'S-08 data.csv'
attach(precipitation)
summary(precipitation)
precipitation$trial<-as.factor(precipitation$trial)

#Next, create separate datasets for each month
may <- filter(precipitation, month == "May")[,match(c("date","day.in.month","julian.date","precip.mm","trial"), colnames (precipitation))]
june <- filter(precipitation, month == "June")[,match(c("date","day.in.month","julian.date","precip.mm","trial"), colnames (precipitation))]
july <- filter(precipitation, month == "July")[,match(c("date","day.in.month","julian.date","precip.mm","trial"), colnames (precipitation))]
august <- filter(precipitation, month == "August")[,match(c("date","day.in.month","julian.date","precip.mm","trial"), colnames (precipitation))]
september <- filter(precipitation, month == "September")[,match(c("date","day.in.month","julian.date","precip.mm","trial"), colnames (precipitation))]

#Create precipitation plots for each month
precip.may.plot<-ggplot(data=may, aes(x=day.in.month, y=precip.mm))+geom_line()+ geom_point(aes(fill=trial,color=trial),size=4) + scale_fill_manual(values = c("#FF4500","grey"))+ scale_color_manual(values = c("#FF4500","grey"))+theme_bw(12)+scale_y_continuous(limits=c(0,60),breaks=seq(0, 60, 10))+scale_x_continuous(limits=c(1,31),breaks=seq(1, 31, 1))+ theme(panel.grid.minor = element_blank(),legend.position = "none",panel.grid=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank())+ labs(x = "Day of month", y="Total precipitation (mm)")+annotate(geom = "text", x = 1, y = 60, label = "May", color = "black",hjust=0,vjust=1.2,size=5)

precip.june.plot<-ggplot(data=june, aes(x=day.in.month, y=precip.mm))+geom_line()+ geom_point(aes(fill=trial,color=trial),size=4) + scale_fill_manual(values = c("#FF4500","grey"))+ scale_color_manual(values = c("#FF4500","grey"))+theme_bw(12)+scale_y_continuous(limits=c(0,60),breaks=seq(0, 60, 10))+scale_x_continuous(limits=c(1,31),breaks=seq(1, 31, 1))+ theme(panel.grid.minor = element_blank(),legend.position = "none",panel.grid=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank())+ labs(x = "Day of month", y="Total precipitation (mm)")+annotate(geom = "text", x = 1, y = 60, label = "June", color = "black",hjust=0,vjust=1.2,size=5)

precip.july.plot<-ggplot(data=july, aes(x=day.in.month, y=precip.mm))+geom_line()+ geom_point(aes(fill=trial,color=trial),size=4) + scale_fill_manual(values = c("#FF4500","grey"))+ scale_color_manual(values = c("#FF4500","grey"))+theme_bw(12)+scale_y_continuous(limits=c(0,60),breaks=seq(0, 60, 10))+scale_x_continuous(limits=c(1,31),breaks=seq(1, 31, 1))+ theme(panel.grid.minor = element_blank(),legend.position = "none",panel.grid=element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank())+ labs(x = "Day of month", y="Total precipitation (mm)")+annotate(geom = "text", x = 1, y = 60, label = "July", color = "black",hjust=0,vjust=1.2,size=5)

precip.aug.plot<-ggplot(data=august, aes(x=day.in.month, y=precip.mm))+geom_line()+ geom_point(aes(fill=trial,color=trial),size=4) + scale_fill_manual(values = c("#FF4500","grey"))+ scale_color_manual(values = c("#FF4500","grey"))+theme_bw(12)+scale_y_continuous(limits=c(0,60),breaks=seq(0, 60, 10))+scale_x_continuous(limits=c(1,31),breaks=seq(1, 31, 1))+ theme(panel.grid.minor = element_blank(),legend.position = "none",panel.grid=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank())+ labs(x = "Day of month", y="Total precipitation (mm)")+annotate(geom = "text", x = 1, y = 60, label = "August", color = "black",hjust=0,vjust=1.2,size=5)

precip.sept.plot<-ggplot(data=september, aes(x=day.in.month, y=precip.mm))+geom_line()+ geom_point(aes(fill=trial,color=trial),size=4) + scale_fill_manual(values = c("#FF4500","grey"))+ scale_color_manual(values = c("#FF4500","grey"))+theme_bw(12)+scale_y_continuous(limits=c(0,60),breaks=seq(0, 60, 10))+scale_x_continuous(limits=c(1,31),breaks=seq(1, 31, 1))+ theme(panel.grid.minor = element_blank(),legend.position = "none",panel.grid=element_blank(),axis.title.y=element_blank())+ labs(x = "Day of month", y="Total precipitation (mm)")+annotate(geom = "text", x = 1, y = 60, label = "September", color = "black",hjust=0,vjust=1.2,size=5)

#Combine plots
precip.plot<-precip.may.plot/precip.june.plot/precip.july.plot/precip.aug.plot/precip.sept.plot

#Export plot
tiff("R1_FSCP Fig S3.tiff", width = 10, height = 12, units = "in", res = 800)
precip.plot
dev.off()


#### Figure S4 ####
#Upload datafile for clam shell length
size<-read.csv(file.choose()) #datafile = 'S-09 data.csv'
attach(size)
summary(size)

#Reorder experiment levels
size$experiment<-factor(size$experiment,levels=c("May","June","July","August","September"))

#Generate and view boxplot for size (shell length) distributions for each experiment
size.plot<-ggplot(size,aes(x=experiment, y=shell.length))+stat_boxplot(geom ='errorbar',width=0.4)+ geom_boxplot(aes(fill=experiment),color="black",outlier.shape = NA) + geom_jitter(aes(fill=experiment, size=0.01,alpha=0.2),color="black",width=0.3,shape=21)+ scale_color_manual(values = c("#D3D3D3", "#FF4500", "#F79470", "#EEAD95","#F69775"))+scale_y_continuous(limits=c(25,55),breaks=seq(25, 55, 5))+ scale_fill_manual(values = c("#D3D3D3", "#FF4500", "#F79470", "#EEAD95","#F69775")) + labs(x = "Trial", y="Shell length (mm)")+ theme_bw(12)+ theme(legend.position="none")+ theme(panel.grid.minor = element_blank()) +theme(panel.grid.major = element_blank(),legend.title = element_blank())
size.plot

#Export plot
tiff("R1_FSCP Fig S4.tiff", width = 6, height = 6, units = "in", res = 800)
size.plot
dev.off()

#Add image to exported figure using Adobe Photoshop; also add colour gradient from previosu Figure (Fig 2,3)

#### Figure S5 ####
Not generated using code (image)

#### SUPPLEMENTARY ANALYSIS: DAY 1 REBURROWING ####
#This analysis pertains to the results provided in the supplementary file "Supplementary analysis"

#Upload datafile for burrowing and mortality after 24 and 48 hours
burrowing.day1<-read.csv(file.choose()) #datafile = 'S-10.csv'
attach(burrowing.day1)
summary(burrowing.day1)
burrowing.day1$plot.glob<-as.factor(burrowing.day1$plot.glob)
burrowing.day1$month<-as.factor(burrowing.day1$month)
burrowing.day1$time.since.deploy<-as.factor(burrowing.day1$time.since.deploy)
burrowing.day1$tide.level<-as.factor(burrowing.day1$tide.level)

#Build BGLMER models for reburrowing at each 15-min observation interval (2 hours total) at each tide level, in each experiment (month) on Day 1. Include plot to account for repeated measures over two time points. Note that there is complete separation in the data because some treatment levels have all 0 values. As such, use Bayesian GLMER approach [bglmer() as per Bolker (2018)), specifying zero-mean normal priors to account for separation. Did not include predator treatment in Day 1 burrowing because no predation or mortality was observed during this time and we would not expect predator treatment to affect burrowing.

#Construct model
#The lower variance prior was chosen due to the removal of predator treatment as a fixed factor and the increased number of replicates for each fixed factor combination level. Diagonal matrix value of 120 was chosen to match the number of terms in the models.
bglmer_burrowmod_full<-bglmer(prop.burrowed~month*time.since.deploy*tide.level+(1|plot.glob),family=binomial,data=burrowing.day1, fixef.prior = normal(cov = diag(6,120))) #Ignore convergence warnings; model estimates still OK

#Get results
Anova(bglmer_burrowmod_full,type=3) #significant effect of

#Get pairwise results
#Partial model
burrow.full.pairwise1<-emmeans(bglmer_burrowmod_full,~time.since.deploy|tide.level|month,adjustment="Holm")
pairs(burrow.full.pairwise1)


#### Figure SA1 ####
#Use data file from above statistical analysis (data frame "burrowing.day1", from datafile 'S-07 data')

#First, reorder "month" variable so that months are ordered chronologically
burrowing.day1$month<-factor(burrowing.day1$month,levels=c("May","June","July","Aug","Sept"))
burrowing.day1$time.since.deploy<-factor(burrowing.day1$time.since.deploy,levels=c("15 min","30 min","45 min","60 min","75 min","90 min","105 min","120 min"))

#Next, generate summary statistic values to apply in plots. NOTE: need to run base code at the bottom of this .R file prior to running summarySE() commands below.
burrow.sum.full <- summarySE(burrowing.day1, measurevar="prop.full", groupvars=c("month","time.since.deploy","tide.level"))

#Generate and view base plot for proportion reburrowed
burrow.full.plot<-ggplot(burrow.sum.full,aes(x=time.since.deploy, y=prop.burrowed))+geom_errorbar(aes(ymin=prop.full-sd, ymax=prop.full+sd,color=month),width=.2)+ geom_point(aes(fill=month,color=month),size=5)+ scale_fill_manual(values = c("#D3D3D3", "#FF4500", "#F79470", "#EEAD95","#F69775"))+ scale_color_manual(values = c("#D3D3D3", "#FF4500", "#F79470", "#EEAD95","#F69775"))+scale_y_continuous(limits=c(-0.15,1.15),breaks=seq(0, 1, 0.2))+ facet_grid(~tide.level~month) + labs(x = "Time since fishing", y="Proportion reburrowed")+ theme_bw(15)+ theme(legend.position="none")+ theme(panel.grid = element_blank())+ theme(axis.text.x = element_text(angle = 45,vjust=1.2,hjust=1.2))
burrow.full.plot

#Export plot
tiff("FSCP Fig SA1.tiff", width = 14, height = 10, units = "in", res = 800)
burrow.full.plot
dev.off()


#### BASE CODE FOR SUMMARY STATISTICS (FIG 2) ####
## Need to run prior to using summarySE() commands above ##
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}
## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}
