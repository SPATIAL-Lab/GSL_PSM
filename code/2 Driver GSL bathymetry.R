library(coda)
library(rjags)
library(R2jags)
library(mcmcplots)
library(bayestestR)
library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)
library(isoWater)
library(tidyr)

setwd("C:/Users/ydmag/Google Drive/U of U/GSL proxy/GSL_PSM")

source("code/1 Helper functions.R")

#Get bathymetry chart for GSL
GSL.bathy <- read.csv("data/GSL bathymetry.csv")

GSL.level <- GSL.bathy$Depth*0.3048 #convert to meters

GSL.area <- (GSL.bathy$South.Area..acres+GSL.bathy$North.Area..acres)/247.1 #convert to km2

GSL.volume <- (GSL.bathy$South.Volume..acre.feet+GSL.bathy$North.Volume..acre.feet)/810700 #convert to km3

###add data from levels above 1285m, data interpolated from Jewell, 2021
#Lake level is bounded by the threshold to the West Desert at ~ 1286 m (Jewell, 2021)

#check data trend
GSL.level.1286 <- c(GSL.level,1286)
GSL.area.1286 <- c(GSL.area,10000)
GSL.volume.1286 <- c(GSL.volume,54)

plot(GSL.level.1286,GSL.area.1286,type="l")
plot(GSL.level.1286,GSL.volume.1286, type="l")

#need to predict lake area from volume?

#estimated salt load of the GSL (Loving 2000, Mohammed 2012)
SL.GSL <- 4.05e3 

#saturation of salt at 355 g/L
sub.GSL <- subset(GSL.volume.1286,abs(SL.GSL/GSL.volume.1286)<355)

#salinity curve
GSL.sali <- c(rep(355,length(GSL.volume.1286)-length(sub.GSL)),SL.GSL/sub.GSL)

# GSL.sali.sg <- 0.99907 + GSL.sali/1000 #calculate specific gravity at 15.56 celcius

plot(GSL.level.1286,GSL.sali,type = "l", ylim=c(0,350),ylab="Salinity, g/L")
#check historic salinity measurement: in 1949, Lake level at 1279.1 m. Salinity at ~270 g/L
points(1279.1,270, col = "red") #plotting right on the line

#surface area of the water body will have an effect on the evaporation rates for areas of similar climates
1-0.029*log(100*GSL.area.1286)

#use empirical data to derive salt water evaporation coefficient 
#study #1 Arons and Kientzler, 20 degrees C
Es = 6.1078 * exp(17.269 * 20/(237.3 + 20)) #millibar

A.K.redu.hgmm <- c(0.08, 0.15, 0.32, 0.49, 0.67,0.90, 1.14, 1.36, 1.62,1.92, 2.2, 
                    2.55, 2.91, 3.28,3.63,4.01, 4.49)
A.K.redu <- A.K.redu.hgmm*1.33322/(6.1078 * exp(17.269 * 20/(237.3 + 20))) #seems ok
A.K.sg <- 0.9982 + c(5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160)*41.55/23/1000
A.K.slpm <- c(5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160)*41.55/23

plot(A.K.slpm,1-A.K.redu,ylab="Vapor pressure coefficient", xlim=c(0,355),ylim = c(0.2,1))

plot(A.K.sg,1-A.K.redu,ylab="Vapor pressure coefficient", xlim=c(1,1.35),ylim = c(0.2,1))
  
#study #2 Dickson 1965, 20 degrees C
Dickson.redu <- c(0.20, 0.23, 0.26, 0.30, 0.34, 0.39, 0.42)
Dickson.sg <- c(1.188, 1.198, 1.208, 1.218, 1.229, 1.239, 1.244)
Dickson.slpm <- c(247, 256, 264, 273, 281, 289, 292)
points(Dickson.slpm,1-Dickson.redu,col="blue")

#study #3 Salhotra 1985
Saihotra.redu <- 1-c(0.97, 0.87, 0.82, 0.83, 0.78,0.73)
Saihotra.slpm <- c((48+60)/2, (185+215)/2, (221+245)/2, (188+209)/2, (220+242)/2, (232+272)/2)
points(Saihotra.slpm,1-Saihotra.redu,col="green")

beta.sali <- c(1-A.K.redu,1-Dickson.redu, 1-Saihotra.redu,1)
exp.sali <- c(A.K.slpm, Dickson.slpm, Saihotra.slpm,0)
exp.data <- data.frame(exp.sali, beta.sali)

nls(beta.sali~ 1-a*exp.sali^b, data = exp.data, start = list(a=0.001, b = 0.1))
# a = 9.098e-07
# b = 2.267

#test plots
plot.x <- seq(from=0,to = 355,by = 1)
plot(plot.x, 1-0.01749*exp(0.01046*plot.x),type = "l", xlab = "Salinity g/L", 
     ylab = "Vapor pressure coefficient")
points(exp.sali,beta.sali)

plot(plot.x, 1-9.098e-07*plot.x^2.267,type = "l", xlab = "Salinity g/L", 
     ylab = "Vapor pressure coefficient")
points(exp.sali,beta.sali)

GSL.salicoeff <- 1-9.098e-07*GSL.sali^2.267

# nls(beta.sali~ 1-a*exp(b*exp.sali), data = exp.data, start = list(a=0.1, b = 0.01))
# # a = 0.01749
# # b = 0.01046
