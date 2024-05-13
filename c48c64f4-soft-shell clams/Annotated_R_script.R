####Get packages####
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
library(ggplot2)  
library(grid)  

####Section 1: Analyses for site differences in temperature and sediment composition####
#upload and configure temperature dataset and sediment characteristics dataset
temp <-read.csv(file.choose()) #data file: "temperature.csv"
attach(temp)
temp$site<-as.factor(temp$site)
summary(temp)

sediment <-read.csv(file.choose()) #data file: "sediment.csv"
attach(sediment)
sediment$site<-as.factor(sediment$site)
summary(sediment)

#Build model for temperature
temp.mod<-lm(temp~site,data=temp)
anova(temp.mod)

#Get pairwise differences between sites
temppairwise<-emmeans(temp.mod,~site,adjustment="Bonferroni")
pairs(temppairwise)

#Build model for grain size (phi)
grainsize.mod<-lm(grainsize~site,data=sediment)
anova(grainsize.mod)

#Get pairwise differences between sites
grainsizepairwise<-emmeans(grainsize.mod,~site,adjustment="Bonferroni")
pairs(grainsizepairwise)

#Build model for % organic content
OC.mod<-lm(OC~site,data=sediment)
anova(OC.mod)

#Get pairwise differences between sites
OCpairwise<-emmeans(OC.mod,~site,adjustment="Bonferroni")
pairs(OCpairwise)

#Build model for relative moisture content
moisture.mod<-lm(moisture~site,data=sediment)
anova(moisture.mod)

#Get pairwise differences between sites
moisturepairwise<-emmeans(moisture.mod,~site,adjustment="Bonferroni")
pairs(moisturepairwise)

####Section 2: Analyses for predicting clam presence from siphon holes####
#upload and configure size calibration dataset
siphon.calib <-read.csv(file.choose()) #data file: "siphon hole size - clam size.csv"
attach(siphon.calib)
summary(siphon.calib)
siphon.calib$site<-as.factor(siphon.calib$site)
siphon.calib$correct.id<-as.factor(siphon.calib$correct.id)
summary(siphon.calib)

#Create separate datasets for each site
PowellsCove <- filter(siphon.calib, site == "Powell's Cove")[,match(c("correct.id","correct.id.binary","siphon.hole.length","clam.length","clam.weight"), colnames (siphon.calib))]

Maisonnette <- filter(siphon.calib, site == "Maisonnette")[,match(c("correct.id","correct.id.binary","siphon.hole.length","clam.length","clam.weight"), colnames (siphon.calib))]

Shemogue <- filter(siphon.calib, site == "Shemogue")[,match(c("correct.id","correct.id.binary","siphon.hole.length","clam.length","clam.weight"), colnames (siphon.calib))]

Kouchibouguac <- filter(siphon.calib, site == "Kouchibouguac")[,match(c("correct.id","correct.id.binary","siphon.hole.length","clam.length","clam.weight"), colnames (siphon.calib))]

#NOTE: Percentage of correct identifications computed in MS Excel. Maisonette and Powell's Cove had multiple incorrect IDs

#Chi square tests comparing observed percentages of correct:incorrect siphon hole identifications to expected 50:50% distribution for each site
#Manually enter data for each site
Chisq.maisonnette<-c(50,14)
Chisq.powells<-c(50,12)
Chisq.shemogue<-c(50,0)
Chisq.kouchi<-c(50,0)

#Build models and obtain results
maisonnette.x2.mod<- chisq.test(Chisq.maisonnette,p=c(1/2,1/2))
powells.x2.mod<- chisq.test(Chisq.powells,p=c(1/2,1/2))
shemogue.x2.mod<- chisq.test(Chisq.shemogue,p=c(1/2,1/2))
kouchi.x2.mod<- chisq.test(Chisq.kouchi,p=c(1/2,1/2))
maisonnette.x2.mod
powells.x2.mod
shemogue.x2.mod
kouchi.x2.mod

#Logistic regression for sites with multiple incorrect identifications  
#Build binary GLM models for Maisonette and Powell's Cove for effect of siphon hole length on correct ID binary
MAISONNETTE1<- glm(correct.id.binary ~ siphon.hole.length, Maisonnette, family = binomial)
POWELLS1 <- glm(correct.id.binary ~ siphon.hole.length, PowellsCove, family = binomial)

#Get significance for siphon hole size effect
Anova(MAISONNETTE1)
Anova(POWELLS1)

#obtain 50%, 75%, and 95% correct ID sizes for each site
dose.p(MAISONNETTE1, p = c(0.5,0.75,0.9))
dose.p(POWELLS1, p = c(0.5,0.75,0.9))


####Section 3: Analyses for siphon hole counts as a proxy for clam abundance#### 
#upload and configure count dataset
siphon.density <-read.csv(file.choose()) #data file: "siphon hole counts - clam abundance.csv"
attach(siphon.density)
summary(siphon.density)
siphon.density$site<-as.factor(siphon.density$site)
siphon.density$quadrat<-as.factor(siphon.density$quadrat)
siphon.density$observer<-as.factor(siphon.density$observer)
summary(siphon.density)

#build model for clam count ~ siphon hole density
density.mod<-lm(clams.above.20mm~hole.count*site, data=siphon.density)

#Check assumptions
plot(density.mod)
#Assumptions violated

#Try data transformation
log.clam.count<-log(1+clams.above.20mm)
log.hole.count<-log(1+hole.count)

#build model for siphon hole density with transformed data
log.density.mod<-lm(log.clam.count~hole.count*site, data=siphon.density)
plot(log.density.mod)

#Assumptions OK; slightly non-normal but OK

#Get results
anova(log.density.mod)

#build model for clam biomass ~ siphon hole density
biomass.mod<-lm(total.clam.weight~hole.count*site, data=siphon.density)

#Check assumptions
plot(biomass.mod)
#Assumptions OK; slightly non-normal but OK

#Get results
anova(biomass.mod)

#Create separate datasets for each site
PowellsCove.observer <- filter(siphon.density, site == "Powell's Cove")[,match(c("quadrat","observer","hole.count","clams.above.20mm","clams.below.20mm","total.clams","total.clam.weight","include.exclude"), colnames (siphon.density))]

Maisonnette.observer <- filter(siphon.density, site == "Maisonnette")[,match(c("quadrat","observer","hole.count","clams.above.20mm","clams.below.20mm","total.clams","total.clam.weight","include.exclude"), colnames (siphon.density))]

Shemogue.observer <- filter(siphon.density, site == "Shemogue")[,match(c("quadrat","observer","hole.count","clams.above.20mm","clams.below.20mm","total.clams","total.clam.weight","include.exclude"), colnames (siphon.density))]

Kouchibouguac.observer <- filter(siphon.density, site == "Kouchibouguac")[,match(c("quadrat","observer","hole.count","clams.above.20mm","clams.below.20mm","total.clams","total.clam.weight","include.exclude"), colnames (siphon.density))]

#Build lm for observer effect for each site and get oberver effect results
observer.maisonnette.mod<-lm(hole.count~observer, data=Maisonnette.observer)
plot(observer.maisonnette.mod)
anova(observer.maisonnette.mod)

observer.kouchi.mod<-lm(hole.count~observer, data=Kouchibouguac.observer)
plot(observer.kouchi.mod)
anova(observer.kouchi.mod)

observer.shemogue.mod<-lm(hole.count~observer, data=Shemogue.observer)
plot(observer.shemogue.mod)
anova(observer.shemogue.mod)

observer.powells.mod<-lm(hole.count~observer, data=PowellsCove.observer)
plot(observer.powells.mod)
anova(observer.powells.mod)

####Section 4: Analyses for siphon hole dimensions as a proxy for clam size####
#NOTE: use same data file that was uploaded for previous analyses for predicting clam presence from siphon holes" above (i.e., #data file: "siphon hole size - clam size.csv")
#build model for clam length and get statistical summaries
length.mod<-lm(clam.length~siphon.hole.length*site, data=siphon.calib)
plot(length.mod)
#Assumptions OK; proceed with ANCOVA
anova(length.mod)

#get pairwise differences between sites
lengthpairwise<-emmeans(length.mod,~site,adjustment="Bonferroni")
options(max.print=1000000)
pairs(lengthpairwise)

#Get relationships and R^2 for each site
powells.length.mod<-lm(clam.length~siphon.hole.length, data=PowellsCove)
summary(powells.length.mod)

maisonnette.length.mod<-lm(clam.length~siphon.hole.length, data=Maisonnette)
summary(maisonnette.length.mod)

shemogue.length.mod<-lm(clam.length~siphon.hole.length, data=Shemogue)
summary(shemogue.length.mod)

kouchi.length.mod<-lm(clam.length~siphon.hole.length, data=Kouchibouguac)
summary(kouchi.length.mod)

#build model for clam weight and get statistical summaries
weight.mod<-lm(clam.weight~siphon.hole.length*site, data=siphon.calib)
plot(weight.mod)
#Assumptions OK; slightly non-normal; proceed with ANCOVA
anova(weight.mod)

#get pairwise differences between sites
weightpairwise<-emmeans(weight.mod,~site|siphon.hole.length,adjustment="Bonferroni")
options(max.print=1000000)
pairs(weightpairwise)

#Get relationships and R^2 for each site
powells.weight.mod<-lm(clam.weight~siphon.hole.length, data=PowellsCove)
summary(powells.weight.mod)

maisonnette.weight.mod<-lm(clam.weight~siphon.hole.length, data=Maisonnette)
summary(maisonnette.weight.mod)

shemogue.weight.mod<-lm(clam.weight~siphon.hole.length, data=Shemogue)
summary(shemogue.weight.mod)

kouchi.weight.mod<-lm(clam.weight~siphon.hole.length, data=Kouchibouguac)
summary(kouchi.weight.mod)

####Section 5: Sensitivity analysis for siphon hole counts as a proxy for clam abundance for Kouchibouguac with the single large sample excluded#### 
#upload and configure Kouchibouguac count dataset
kouchi.count <-read.csv(file.choose()) #data file: "siphon hole counts clam abundance - kouchibouguac no large sample.csv"
attach(kouchi.count)
summary(kouchi.count)
kouchi.count$quadrat<-as.factor(kouchi.count$quadrat)
kouchi.count$observer<-as.factor(kouchi.count$observer)
summary(kouchi.count)

#build model for clam count ~ siphon hole density
kouchi.density.mod<-lm(clams.above.20mm~hole.count, data=kouchi.count)

#Get results
anova(kouchi.density.mod)
#Significant relationship

####Section 6: Create figures####
####Figure 1####
No code required

####Figure 2####
#**NOTE** NEeed to switch to R version 4.2.2 to run the code below and must download package 'rgdal' from CRAN Archive
#Get packages
library(geojsonio)
library(rmapshaper)
library(rgdal)
library(tidyverse)
library(sp)
library(socviz)
library(mapview)
library(ggplot2)
library(ggOceanMaps)
library(rnaturalearthdata)
library(raster)
library(proj4)
library(maptools)
library(mapproj)
library(ggspatial)
library(jpeg)
library(viridis)
library(raster)
library(lattice)
library(rgeos)

#Basecode for scalebar
library(maps) 
library(maptools)  

#Then, we need a function to get the scale bar coordinates:
#
# Result #
#--------#
# Return a list whose elements are :
#   - rectangle : a data.frame containing the coordinates to draw the first rectangle ;
#   - rectangle2 : a data.frame containing the coordinates to draw the second rectangle ;
#   - legend : a data.frame containing the coordinates of the legend texts, and the texts as well.
#
# Arguments : #
#-------------#
# lon, lat : longitude and latitude of the bottom left point of the first rectangle to draw ;
# distanceLon : length of each rectangle ;
# distanceLat : width of each rectangle ;
# distanceLegend : distance between rectangles and legend texts ;
# dist.units : units of distance "km" (kilometers) (default), "nm" (nautical miles), "mi" (statute miles). 
createScaleBar <-
  function(lon,lat,distanceLon,distanceLat,distanceLegend, dist.units =
             "km"){
    # First rectangle
    bottomRight <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distanceLon, dist.units = dist.units, model = "WGS84")
    
    topLeft <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distanceLat, dist.units = dist.units, model = "WGS84")
    rectangle <- cbind(lon=c(lon, lon, bottomRight[1,"long"], bottomRight[1,"long"], lon),
                       lat = c(lat, topLeft[1,"lat"], topLeft[1,"lat"],lat, lat))
    rectangle <- data.frame(rectangle, stringsAsFactors = FALSE)
    
    # Second rectangle t right of the first rectangle
    bottomRight2 <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distanceLon*2, dist.units = dist.units, model = "WGS84")
    rectangle2 <- cbind(lon = c(bottomRight[1,"long"], bottomRight[1,"long"], bottomRight2[1,"long"], bottomRight2[1,"long"],
                                bottomRight[1,"long"]),
                        lat=c(lat, topLeft[1,"lat"], topLeft[1,"lat"], lat, lat))
    rectangle2 <- data.frame(rectangle2, stringsAsFactors = FALSE)
    
    # Now let's deal with the text
    onTop <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distanceLegend, dist.units = dist.units, model = "WGS84")
    onTop2 <- onTop3 <- onTop
    onTop2[1,"long"] <- bottomRight[1,"long"]
    onTop3[1,"long"] <- bottomRight2[1,"long"]
    
    legend <- rbind(onTop, onTop2, onTop3)
    legend <- data.frame(cbind(legend, text = c(0, distanceLon, distanceLon*2)), stringsAsFactors = FALSE, row.names = NULL)
    return(list(rectangle = rectangle, rectangle2 = rectangle2, legend = legend)) } 

#We also need a function to obtain the coordinates of the North arrow:

#
# Result #
#--------#
# Returns a list containing :
#   - res : coordinates to draw an arrow ;
#   - coordinates of the middle of the arrow (where the "N" will be plotted).
#
# Arguments : #
#-------------#
# scaleBar : result of createScaleBar() ;
# length : desired length of the arrow ;
# distance : distance between legend rectangles and the bottom of the arrow ;
# dist.units : units of distance "km" (kilometers) (default), "nm" (nautical miles), "mi" (statute miles). 
createOrientationArrow <-
  function(scaleBar, length, distance = 1, dist.units = "km"){
    lon <- scaleBar$rectangle2[1,1]
    lat <- scaleBar$rectangle2[1,2]
    
    # Bottom point of the arrow
    begPoint <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance, dist.units = dist.units, model = "WGS84")
    lon <- begPoint[1,"long"]
    lat <- begPoint[1,"lat"]
    
    # Let us create the endpoint
    onTop <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = length, dist.units = dist.units, model = "WGS84")
    
    leftArrow <- gcDestination(lon = onTop[1,"long"], lat = onTop[1,"lat"], bearing = 225, dist = length/5, dist.units =
                                 dist.units, model = "WGS84")
    
    rightArrow <- gcDestination(lon = onTop[1,"long"], lat = onTop[1,"lat"], bearing = 135, dist = length/5, dist.units =
                                  dist.units, model = "WGS84")
    
    res <- rbind(
      cbind(x = lon, y = lat, xend = onTop[1,"long"], yend = onTop[1,"lat"]),
      cbind(x = leftArrow[1,"long"], y = leftArrow[1,"lat"], xend = onTop[1,"long"], yend = onTop[1,"lat"]),
      cbind(x = rightArrow[1,"long"], y = rightArrow[1,"lat"], xend = onTop[1,"long"], yend = onTop[1,"lat"]))
    
    res <- as.data.frame(res, stringsAsFactors = FALSE)
    
    # Coordinates from which "N" will be plotted
    coordsN <- cbind(x = lon, y = (lat + onTop[1,"lat"])/2)
    
    return(list(res = res, coordsN = coordsN)) } 
#The last function enables the user to draw the elements:

#
# Result #
#--------#
# This function enables to draw a scale bar on a ggplot object, and optionally an orientation arrow #
# Arguments : #
#-------------#
# lon, lat : longitude and latitude of the bottom left point of the first rectangle to draw ;
# distanceLon : length of each rectangle ;
# distanceLat : width of each rectangle ;
# distanceLegend : distance between rectangles and legend texts ;
# dist.units : units of distance "km" (kilometers) (by default), "nm" (nautical miles), "mi" (statute miles) ;
# rec.fill, rec2.fill : filling colour of the rectangles (default to white, and black, resp.);
# rec.colour, rec2.colour : colour of the rectangles (default to black for both);
# legend.colour : legend colour (default to black);
# legend.size : legend size (default to 3);
# orientation : (boolean) if TRUE (default), adds an orientation arrow to the plot ;
# arrow.length : length of the arrow (default to 500 km) ;
# arrow.distance : distance between the scale bar and the bottom of the arrow (default to 300 km) ;
# arrow.North.size : size of the "N" letter (default to 6). 
scaleBar <- function(lon, lat, distanceLon, distanceLat, distanceLegend,
                     dist.unit = "km", rec.fill = "white", rec.colour = "black", rec2.fill
                     = "black", rec2.colour = "black", legend.colour = "black", legend.size = 3, orientation = TRUE, arrow.length = 500, arrow.distance = 300, arrow.North.size = 6){
  laScaleBar <- createScaleBar(lon = lon, lat = lat, distanceLon = distanceLon, distanceLat = distanceLat, distanceLegend =
                                 distanceLegend, dist.unit = dist.unit)
  # First rectangle
  rectangle1 <- geom_polygon(data = laScaleBar$rectangle, aes(x = lon, y = lat), fill = rec.fill, colour = rec.colour)
  
  # Second rectangle
  rectangle2 <- geom_polygon(data = laScaleBar$rectangle2, aes(x = lon, y = lat), fill = rec2.fill, colour = rec2.colour)
  
  # Legend
  scaleBarLegend <- annotate("text", label = paste(laScaleBar$legend[,"text"], dist.unit, sep=""), x =
                               laScaleBar$legend[,"long"], y = laScaleBar$legend[,"lat"], size =
                               legend.size, colour = legend.colour)
  
  res <- list(rectangle1, rectangle2, scaleBarLegend)
  
  if(orientation){# Add an arrow pointing North
    coordsArrow <- createOrientationArrow(scaleBar = laScaleBar, length = arrow.length, distance = arrow.distance, dist.unit =
                                            dist.unit)
    arrow <- list(geom_segment(data = coordsArrow$res, aes(x = x, y = y, xend = xend, yend = yend)), annotate("text", label = "N", x =
                                                                                                                coordsArrow$coordsN[1,"x"], y = coordsArrow$coordsN[1,"y"], size =
                                                                                                                arrow.North.size, colour = "black"))
    res <- c(res, arrow)
  }
  return(res) }

#Create map of study sites
#Upload lat/long for sampling sites
#values are lat/long in decimal degrees
siphons<-read.csv(file.choose()) #data file: "siphon coords.csv"
attach(siphons)

#Upload shapefiles (shapefiles obtained from J. Barrell, DFO Gulf)
maritimes<- readOGR(dsn = "W:/Molluscan/People/ClementsJeffC/Maps - shapefiles/Maritimes_DecimalDegree", "Maritimes_NAD83_CSRS",stringsAsFactors = F)
summary(maritimes@data)
eastern<- readOGR(dsn = "W:/Molluscan/People/ClementsJeffC/Maps - shapefiles/Eastern Canada_DecimalDegree", "eastern_canada_wgs84",stringsAsFactors = F)
summary(eastern@data)

#Upload north arrow image
north<- readJPEG(file.choose())

#set axis label decimal places
scaleFUNx <- function(x) sprintf("%.2f", x)
scaleFUNy <- function(y) sprintf("%.2f", y)

#Create and export site map

siphons$Site<-factor(siphons$Site,levels=c("Maisonnette","Kouchibouguac","Shemogue","Powell's Cove"))
siphon_map <- ggplot() + geom_polygon(data = eastern, aes(x = long, y = lat,group=group), colour = "#A6A6A6",fill="lightgray")+coord_cartesian(xlim=c(-66.5,-59), ylim = c(43.5,48.5))+geom_point(data = siphons, aes(x = long, y = lat,fill=Site,shape=Site), size=5,color="black")+scale_shape_manual(values = c(21,22,24,23))+theme_bw(14)+scale_fill_viridis(discrete = TRUE ,alpha=0.8)+scale_color_viridis(discrete = TRUE)+ theme(axis.title.x = element_blank(), axis.title.y = element_blank(),axis.line = element_line(colour = "black",size=0.1), panel.border = element_rect(colour = "black", fill=NA, size=1),axis.ticks = element_line(colour="black"), axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+theme(legend.position = c(0.45,0.85))+ scaleBar(lon = -62.22, lat = 43.62, distanceLon = 100, distanceLat = 10, distanceLegend = -10, dist.unit = "km",arrow.length = 0,arrow.distance = 0, arrow.North.size = 0)+ scale_y_continuous(labels=scaleFUNy,breaks=seq(43,50,by=1))+ scale_x_continuous(labels=scaleFUNx,breaks=seq(-70,-56,by=2))+annotation_raster(north, xmin = -61.5, xmax = -60.5, ymin = 43.8, ymax = 44.6)
siphon_map

tiff("Siphon hole map.tiff", width = 8.5, height = 8, units = "in", res = 800)
siphon_map
dev.off()

#Figure aesthetic edited in Powerpoint, including resizing and placement of North American inset map (created previously for another project), placement of legend, placement of scale bar and compass, and addition of landmark names

####Figure 3####
#Create logistic regression plot for Maisonette
maisonnette.logreg<-ggplot(Maisonnette, aes(x=siphon.hole.length, y=correct.id.binary)) + 
  geom_point(aes(color="#440154",fill="#440154",size=2),alpha=0.5,shape=21) +stat_smooth(aes(color="#440154",fill="#440154"),method="glm", se=TRUE, fullrange=TRUE, method.args = list(family=binomial)) + ylab("Correct identification") + xlab("Siphon hole length (mm)")+ xlim(0, 20) +scale_color_manual(values = c("#440154"))+scale_fill_manual(values = c("#440154"))+theme_bw(14)+theme(legend.position = "none",panel.grid = element_blank())

#Create logistic regression plot for Powell's Cove
powells.logreg<-ggplot(PowellsCove, aes(x=siphon.hole.length, y=correct.id.binary)) + 
  geom_point(aes(color="#fde725",fill="#fde725",size=2),alpha=0.5,shape=23) +stat_smooth(aes(color="#fde725",fill="#fde725"),method="glm", se=TRUE, fullrange=TRUE, method.args = list(family=binomial)) + ylab("Correct identification") + xlab("Siphon hole length (mm)")+ xlim(0, 20) +scale_color_manual(values = c("#fde725"))+scale_fill_manual(values = c("#fde725"))+theme_bw(14)+theme(legend.position = "none",axis.title.y=element_blank(),axis.text.y=element_blank(),panel.grid = element_blank())

#Combine plots
logreg.combo.plot<- maisonnette.logreg+powells.logreg
logreg.combo.plot

#export plot
tiff("siphon hole logistic regressions.tiff", width = 8, height = 4, units = "in", res = 800)
logreg.combo.plot
dev.off()

#Aesthetic edited in Powerpoint, including addition of site names and 50%, 75%, and 90% values

####Figure 4####
#Subset hole-clam count dataset by site
PowellsCoveDensity <- filter(siphon.density, site == "Powell's Cove")[,match(c("quadrat","observer","hole.count","clams.above.20mm","clams.below.20mm","total.clams","total.clam.weight","include.exclude"), colnames (siphon.density))]

MaisonnetteDensity <- filter(siphon.density, site == "Maisonnette")[,match(c("quadrat","observer","hole.count","clams.above.20mm","clams.below.20mm","total.clams","total.clam.weight","include.exclude"), colnames (siphon.density))]

ShemogueDensity <- filter(siphon.density, site == "Shemogue")[,match(c("quadrat","observer","hole.count","clams.above.20mm","clams.below.20mm","total.clams","total.clam.weight","include.exclude"), colnames (siphon.density))]

KouchiDensity <- filter(siphon.density, site == "Kouchibouguac")[,match(c("quadrat","observer","hole.count","clams.above.20mm","clams.below.20mm","total.clams","total.clam.weight","include.exclude"), colnames (siphon.density))]

#Create plots of >20mm clam counts ~ hole counts for each site
powells.density.plot<-ggplot(PowellsCoveDensity, aes(x=hole.count, y=clams.above.20mm,fill=observer,shape=include.exclude)) +geom_point(aes(fill=observer,color=observer,shape=include.exclude), size = 4,alpha=0.6)+ guides(shape = FALSE)+scale_shape_manual(values = c(5,23))+scale_x_continuous(limits=c(0,16),breaks=seq(0, 16, 2))+geom_smooth(aes(x=hole.count, y=clams.above.20mm,fill=observer,color=observer),method="lm", se = TRUE,alpha=0.3)+scale_y_continuous(limits=c(0,16),breaks=seq(0, 16, 2))+scale_color_manual(values = c("#be0026", "#ffda75"))+scale_fill_manual(values = c("#be0026", "#ffda75"))+scale_shape_manual(values = c(5,23))+labs(y="Clam count",x="Siphon hole count")+theme_bw(15)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(legend.position = c(0.15,0.85),legend.title = element_blank(),axis.title.y=element_blank())+ geom_abline(intercept = 0, slope = 1,size=0.5,color="red",linetype="dashed")

maisonnette.density.plot<-ggplot(MaisonnetteDensity, aes(x=hole.count, y=clams.above.20mm,fill=observer,shape=include.exclude)) +geom_point(aes(fill=observer,color=observer,shape=include.exclude), size = 4,alpha=0.6)+ guides(shape = FALSE)+scale_x_continuous(limits=c(0,50),breaks=seq(0, 50, 5))+geom_smooth(aes(x=hole.count, y=clams.above.20mm,fill=observer,color=observer),method="lm", se = TRUE,alpha=0.3)+scale_y_continuous(limits=c(0,50),breaks=seq(0, 50, 5))+scale_color_manual(values = c("#071d58","#1d92c1", "#ffda75"))+scale_fill_manual(values = c("#071d58","#1d92c1", "#ffda75"))+scale_shape_manual(values = c(1,21))+labs(y="Clam count",x="Siphon hole count")+theme_bw(15)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(legend.position = c(0.15,0.85),legend.title = element_blank(),axis.title.x=element_blank())+ geom_abline(intercept = 0, slope = 1,size=0.5,color="red",linetype="dashed")

shem.density.plot<-ggplot(ShemogueDensity, aes(x=hole.count, y=clams.above.20mm,fill=observer)) +geom_point(aes(fill=observer,color=observer), size = 4,alpha=0.6,shape=24)+ guides(shape = FALSE)+scale_x_continuous(limits=c(0,16),breaks=seq(0, 16, 2))+scale_y_continuous(limits=c(0,16),breaks=seq(0, 16, 2))+geom_smooth(aes(x=hole.count, y=clams.above.20mm,fill=observer,color=observer),method="lm", se = TRUE,alpha=0.3)+scale_color_manual(values = c("#071d58","#1d92c1"))+scale_fill_manual(values = c("#071d58","#1d92c1"))+labs(y="Clam count",x="Siphon hole count")+theme_bw(15)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(legend.position = c(0.15,0.85),legend.title = element_blank())+ geom_abline(intercept = 0, slope = 1,size=0.5,color="red",linetype="dashed")+ guides(shape = FALSE)

kouchi.density.plot<-ggplot(KouchiDensity, aes(x=hole.count, y=clams.above.20mm,fill=observer)) +geom_point(aes(fill=observer,color=observer), size = 4,alpha=0.6,shape=22)+ guides(shape = FALSE)+scale_x_continuous(limits=c(0,120),breaks=seq(0, 120, 10))+geom_smooth(aes(x=hole.count, y=clams.above.20mm,fill=observer,color=observer),method="lm", se = TRUE,alpha=0.3)+scale_y_continuous(limits=c(0,120),breaks=seq(0, 120, 10))+scale_color_manual(values = c("#fd4e2a","#071d58","#1d92c1"))+scale_fill_manual(values = c("#fd4e2a","#071d58","#1d92c1"))+labs(y="Clam count",x="Siphon hole count")+theme_bw(15)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(legend.position = c(0.15,0.85),legend.title = element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())+ geom_abline(intercept = 0, slope = 1,size=0.5,color="red",linetype="dashed")

combo.density.plot<-maisonnette.density.plot+kouchi.density.plot+shem.density.plot+powells.density.plot
combo.density.plot

#export plot
tiff("siphon hole counts.tiff", width = 10, height = 10, units = "in", res = 800)
combo.density.plot
dev.off()

# Aesthetic edited in Powerpoint, including addition of site names and shapes in legend of Powell's Cove panel

####Figure 5####
#Create plots of total clam weights ~ hole counts for each site
powells.biomass.plot<-ggplot(PowellsCoveDensity, aes(x=hole.count, y=total.clam.weight,fill=observer,shape=include.exclude)) +geom_point(aes(fill=observer,color=observer,shape=include.exclude), size = 4,alpha=0.6)+ guides(shape = FALSE)+scale_x_continuous(limits=c(0,16),breaks=seq(0, 16, 2))+geom_smooth(aes(x=hole.count, y=total.clam.weight,fill=observer,color=observer),method="lm", se = TRUE,alpha=0.3)+scale_y_continuous(limits=c(0,400),breaks=seq(0, 400, 50))+scale_color_manual(values = c("#be0026", "#ffda75"))+scale_fill_manual(values = c("#be0026", "#ffda75"))+scale_shape_manual(values = c(5,23))+labs(y="Clam biomass (g)",x="Siphon hole count")+theme_bw(15)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(legend.position = c(0.15,0.85),legend.title = element_blank(),axis.title.y=element_blank())

maisonnette.biomass.plot<-ggplot(MaisonnetteDensity, aes(x=hole.count, y=total.clam.weight,fill=observer,shape=include.exclude)) +geom_point(aes(fill=observer,color=observer,shape=include.exclude), size = 4,alpha=0.6)+ guides(shape = FALSE)+scale_x_continuous(limits=c(0,50),breaks=seq(0, 50, 5))+geom_smooth(aes(x=hole.count, y=total.clam.weight,fill=observer,color=observer),method="lm", se = TRUE,alpha=0.3)+scale_y_continuous(limits=c(0,400),breaks=seq(0, 400, 50))+scale_color_manual(values = c("#071d58","#1d92c1", "#ffda75"))+scale_fill_manual(values = c("#071d58","#1d92c1", "#ffda75"))+scale_shape_manual(values = c(1,21))+labs(y="Clam biomass (g)",x="Siphon hole count")+theme_bw(15)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(legend.position = c(0.15,0.85),legend.title = element_blank(),axis.title.x=element_blank())

shem.biomass.plot<-ggplot(ShemogueDensity, aes(x=hole.count, y=total.clam.weight,fill=observer)) +geom_point(aes(fill=observer,color=observer), size = 4,alpha=0.6,shape=24)+ guides(shape = FALSE)+scale_x_continuous(limits=c(0,16),breaks=seq(0, 16, 2))+scale_y_continuous(limits=c(0,400),breaks=seq(0, 400, 50))+geom_smooth(aes(x=hole.count, y=total.clam.weight,fill=observer,color=observer),method="lm", se = TRUE,alpha=0.3)+scale_color_manual(values = c("#071d58","#1d92c1"))+scale_fill_manual(values = c("#071d58","#1d92c1"))+labs(y="Clam biomass (g)",x="Siphon hole count")+theme_bw(15)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(legend.position = c(0.15,0.85),legend.title = element_blank())+ guides(shape = FALSE)

kouchi.biomass.plot<-ggplot(KouchiDensity, aes(x=hole.count, y=total.clam.weight,fill=observer)) +geom_point(aes(fill=observer,color=observer), size = 4,alpha=0.6,shape=22)+ guides(shape = FALSE)+scale_x_continuous(limits=c(0,120),breaks=seq(0, 120, 10))+geom_smooth(aes(x=hole.count, y=total.clam.weight,fill=observer,color=observer),method="lm", se = TRUE,alpha=0.3)+scale_y_continuous(limits=c(0,400),breaks=seq(0, 400, 50))+scale_color_manual(values = c("#fd4e2a","#071d58","#1d92c1"))+scale_fill_manual(values = c("#fd4e2a","#071d58","#1d92c1"))+labs(y="Clam biomass (g)t",x="Siphon hole count")+theme_bw(15)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(legend.position = c(0.15,0.85),legend.title = element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())

combo.biomass.plot<-maisonnette.biomass.plot+kouchi.biomass.plot+shem.biomass.plot+powells.biomass.plot
combo.biomass.plot

#export plot
tiff("siphon hole counts biomass.tiff", width = 10, height = 10, units = "in", res = 800)
combo.biomass.plot
dev.off()

# Aesthetic edited in Powerpoint, including addition of site names and shapes in legend of Powell's Cove panel

####Figure 6####
#create plot for clam length as a function of siphon hole length for each site, and density plot for site effect
siphon.calib$site<-factor(siphon.calib$site,levels=c("Maisonnette","Kouchibouguac","Shemogue","Powell's Cove"))
length.plot<-ggplot(siphon.calib, aes(x=siphon.hole.length, y=clam.length,color=site,fill=site)) +geom_point(aes(fill=site,color = site,shape=site), size = 3,alpha=0.4)+geom_smooth(aes(fill=site,color=site,fill=site), method="lm", se = TRUE,alpha=0.2)+scale_x_continuous(limits=c(0,20),breaks=seq(0, 20, 2))+scale_y_continuous(limits=c(10,90),breaks=seq(10, 90, 10))+ scale_color_viridis(discrete = TRUE)+ scale_fill_viridis(discrete = TRUE,alpha=0.8)+labs(y="Clam length (mm)",x="Siphon hole length (mm)")+scale_shape_manual(values = c(21,22,24,23))+theme_bw(15)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.title.x = element_blank())+theme(legend.position = c(0.2,0.88),legend.title = element_blank())

#Create density plot for clam length
length.density <- ggplot(siphon.calib, aes(clam.length,fill=site))+geom_density(color="black",alpha=0.4)+scale_x_continuous(limits=c(10,90),breaks=seq(10, 90, 10))+ theme_void()+theme(legend.position = "none")+labs(x="Clam length (mm)",y="")+ scale_fill_viridis(discrete = TRUE,alpha=0.5)+coord_flip()

#Combine scatter and density plots for clam length
length.plot.final<-plot_grid(length.plot, length.density, ncol = 2,align="hv",rel_widths = c(2, 0.75))
length.plot.final

#create plot for clam weight as a function of siphon hole length for each site, and density plot for site effect
weight.plot<-ggplot(siphon.calib, aes(x=siphon.hole.length, y=clam.weight,color=site,fill=site)) +geom_point(aes(fill=site,color = site,shape=site), size = 3,alpha=0.4)+geom_smooth(aes(fill=site,color=site,fill=site), method="lm", se = TRUE,alpha=0.2)+scale_x_continuous(limits=c(0,20),breaks=seq(0, 20, 2))+scale_y_continuous(limits=c(0,65),breaks=seq(0, 65, 5))+ scale_color_viridis(discrete = TRUE)+scale_shape_manual(values = c(21,22,24,23))+ scale_fill_viridis(discrete = TRUE ,alpha=0.8)+labs(y="Clam weight (g)",x="Siphon hole length (mm)")+theme_bw(15)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(legend.position = "none",legend.title = element_blank())

#Create density plot for clam weight
weight.density <- ggplot(siphon.calib, aes(clam.weight,fill=site))+geom_density(color="black",alpha=0.4)+scale_x_continuous(limits=c(0,65),breaks=seq(0, 65, 5))+ theme_void()+theme(legend.position = "none")+labs(x="Clam weight (g)",y="")+ scale_fill_viridis(discrete = TRUE,alpha=0.5)+coord_flip()

#Combine scatter and density plots for clam weight
weight.plot.final<-plot_grid(weight.plot, weight.density, ncol = 2,align="hv",rel_widths = c(2, 0.75))

#Combine length and weight plots for final figure
length.weight.combo.plot<-length.plot.final/weight.plot.final
length.weight.combo.plot

#export plot
tiff("Length weight scatter density siphon hole.tiff", width = 9, height = 12, units = "in", res = 800)
length.weight.combo.plot 
dev.off()

# Aesthetic edited in Powerpoint, including box around legend, black lines around the bottom and sides of density plots, and the addition of equations and R^2 values

####Figure S1####
No code required

####Figure S2####
#Upload daily mean temperature dataset
dailytemp <-read.csv(file.choose()) #data file: "daily mean temperature.csv"
attach(dailytemp)
dailytemp$site<-as.factor(dailytemp$site)
summary(dailytemp)

#Create plot for daily mean temperature over time
dailytemp$site<-factor(dailytemp$site,levels=c("Maisonnette","Kouchibouguac","Shemogue","Powell's Cove"))
meantempplot<-ggplot(dailytemp, aes(x=julian.date, y=temp, color=site)) + geom_line(aes(color=site), size=1) + geom_point(aes(color=site,fill=site,shape=site), size=2)+scale_shape_manual(values = c(21,22,24,23))+ scale_color_viridis(discrete = TRUE,alpha=0.6)+ scale_fill_viridis(discrete = TRUE)+labs(y="Daily mean temperature (?C)",x="Julian date")+theme_bw()+theme(panel.grid = element_blank(), legend.position = c(0.9,0.8),legend.title = element_blank())+scale_x_continuous(limits=c(150,300),breaks=seq(150, 300, 10))
meantempplot

#export plot
tiff("temp plot 2.tiff", width = 10, height = 5, units = "in", res = 800)
meantempplot
dev.off()

####Figure S3####
#Create plots of total clam counts ~ hole counts for each site
powells.density.plot.total<-ggplot(PowellsCoveDensity, aes(x=hole.count, y=total.clams,fill=observer,shape=include.exclude)) +geom_point(aes(fill=observer,color=observer,shape=include.exclude), size = 4,alpha=0.6)+ guides(shape = FALSE)+scale_shape_manual(values = c(5,23))+scale_x_continuous(limits=c(0,16),breaks=seq(0, 16, 2))+geom_smooth(aes(x=hole.count, y=total.clams,fill=observer,color=observer),method="lm", se = TRUE,alpha=0.3)+scale_y_continuous(limits=c(0,16),breaks=seq(0, 16, 2))+scale_color_manual(values = c("#be0026", "#ffda75"))+scale_fill_manual(values = c("#be0026", "#ffda75"))+scale_shape_manual(values = c(5,23))+labs(y="Clam count",x="Siphon hole count")+theme_bw(15)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(legend.position = c(0.15,0.85),legend.title = element_blank(),axis.title.y=element_blank())+ geom_abline(intercept = 0, slope = 1,size=0.5,color="red",linetype="dashed")

maisonnette.density.plot.total<-ggplot(MaisonnetteDensity, aes(x=hole.count, y=total.clams,fill=observer,shape=include.exclude)) +geom_point(aes(fill=observer,color=observer,shape=include.exclude), size = 4,alpha=0.6)+ guides(shape = FALSE)+scale_x_continuous(limits=c(0,50),breaks=seq(0, 50, 5))+geom_smooth(aes(x=hole.count, y=total.clams,fill=observer,color=observer),method="lm", se = TRUE,alpha=0.3)+scale_y_continuous(limits=c(0,55),breaks=seq(0, 55, 5))+scale_color_manual(values = c("#071d58","#1d92c1", "#ffda75"))+scale_fill_manual(values = c("#071d58","#1d92c1", "#ffda75"))+scale_shape_manual(values = c(1,21))+labs(y="Clam count",x="Siphon hole count")+theme_bw(15)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(legend.position = c(0.15,0.85),legend.title = element_blank(),axis.title.x=element_blank())+ geom_abline(intercept = 0, slope = 1,size=0.5,color="red",linetype="dashed")

shem.density.plot.total<-ggplot(ShemogueDensity, aes(x=hole.count, y=total.clams,fill=observer)) +geom_point(aes(fill=observer,color=observer), size = 4,alpha=0.6,shape=24)+ guides(shape = FALSE)+scale_x_continuous(limits=c(0,16),breaks=seq(0, 16, 2))+scale_y_continuous(limits=c(0,16),breaks=seq(0, 16, 2))+geom_smooth(aes(x=hole.count, y=total.clams,fill=observer,color=observer),method="lm", se = TRUE,alpha=0.3)+scale_color_manual(values = c("#071d58","#1d92c1"))+scale_fill_manual(values = c("#071d58","#1d92c1"))+labs(y="Clam count",x="Siphon hole count")+theme_bw(15)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(legend.position = c(0.15,0.85),legend.title = element_blank())+ geom_abline(intercept = 0, slope = 1,size=0.5,color="red",linetype="dashed")+ guides(shape = FALSE)

kouchi.density.plot.total<-ggplot(KouchiDensity, aes(x=hole.count, y=total.clams,fill=observer)) +geom_point(aes(fill=observer,color=observer), size = 4,alpha=0.6,shape=22)+ guides(shape = FALSE)+scale_x_continuous(limits=c(0,120),breaks=seq(0, 120, 10))+geom_smooth(aes(x=hole.count, y=total.clams,fill=observer,color=observer),method="lm", se = TRUE,alpha=0.3)+scale_y_continuous(limits=c(0,120),breaks=seq(0, 120, 10))+scale_color_manual(values = c("#fd4e2a","#071d58","#1d92c1"))+scale_fill_manual(values = c("#fd4e2a","#071d58","#1d92c1"))+labs(y="Clam count",x="Siphon hole count")+theme_bw(15)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(legend.position = c(0.15,0.85),legend.title = element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())+ geom_abline(intercept = 0, slope = 1,size=0.5,color="red",linetype="dashed")

combo.density.plot.total<-maisonnette.density.plot.total+kouchi.density.plot.total+shem.density.plot.total+powells.density.plot.total
combo.density.plot.total

#export plot
tiff("siphon hole total counts.tiff", width = 10, height = 10, units = "in", res = 800)
combo.density.plot.total
dev.off()

# Aesthetic edited in Powerpoint, including addition of site names

####Figure S4####
No code required
####Figure S5####
kouchi.density.plot.largeremoved<-ggplot(kouchi.count, aes(x=hole.count, y=clams.above.20mm,fill=observer)) +geom_point(aes(fill=observer,color=observer), size = 4,alpha=0.6,shape=22)+ guides(shape = FALSE)+scale_x_continuous(limits=c(0,50),breaks=seq(0, 50, 5))+geom_smooth(aes(x=hole.count, y=clams.above.20mm,fill=observer,color=observer),method="lm", se = TRUE,alpha=0.3)+scale_y_continuous(limits=c(0,50),breaks=seq(0, 50, 5))+scale_color_manual(values = c("#fd4e2a","#071d58","#1d92c1"))+scale_fill_manual(values = c("#fd4e2a","#071d58","#1d92c1"))+labs(y="Clam count",x="Siphon hole count")+theme_bw(15)+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+theme(legend.position = c(0.15,0.85),legend.title = element_blank())+ geom_abline(intercept = 0, slope = 1,size=0.5,color="red",linetype="dashed")
kouchi.density.plot.largeremoved

#export plot
tiff("kouchi count plot.tiff", width = 6, height = 6, units = "in", res = 800)
kouchi.density.plot.largeremoved
dev.off()
