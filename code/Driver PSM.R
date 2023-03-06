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

GSL.sali.sg <- 0.99907 + GSL.sali/1000 #calculate specific gravity at 15.56 celcius

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

nls(beta.sali~ 1-a*exp(b*exp.sali), data = exp.data, start = list(a=0.1, b = 0.01))
# a = 0.01749
# b = 0.01046

nls(beta.sali~ 1-a*exp.sali^b, data = exp.data, start = list(a=0.001, b = 0.1))
# a = 9.098e-07
# b = 2.267

plot.x <- seq(from=0,to = 355,by = 1)
plot(plot.x, 1-0.01749*exp(0.01046*plot.x),type = "l", xlab = "Salinity g/L", 
     ylab = "Vapor pressure coefficient")
points(exp.sali,beta.sali)

plot(plot.x, 1-9.098e-07*plot.x^2.267,type = "l", xlab = "Salinity g/L", 
     ylab = "Vapor pressure coefficient")
points(exp.sali,beta.sali)

GSL.salicoeff <- 1-9.098e-07*GSL.sali^2.267

#########################water and vapor isotopes###########################
#isoWater to get atm vapor: Not so useful
GSL.Vapor <- wiDB_data(minLat = 39.5, maxLat = 42.5, minLong = -113.5,
                      maxLong = -110.5, minElev = NULL, maxElev = NULL, minDate = NULL,
                      maxDate = NULL, countries = "US", states = "UT", types = "Vapor",
                      projects = NULL, fields = NULL, tmpdir = tempdir(), clean = TRUE)
GSL.Vapor.iso <- data.frame(GSL.Vapor$data$d2H, GSL.Vapor$data$d18O)
GSL.Vapor.iso <- na.omit(GSL.Vapor.iso)
GSL.Vapor.mwl <- mwl(GSL.Vapor.iso, plot = TRUE)
GSL.Vapor.mwl# slope = 7.09
GSL.Vapor.mwl[1]
GSL.Vapor.mwl[2]

GSL.Vapor.mean.d18O <- mean(GSL.Vapor.iso$GSL.Vapor.data.d18O)
GSL.Vapor.mean.d2H <- mean(GSL.Vapor.iso$GSL.Vapor.data.d2H)
points(GSL.Vapor.mean.d18O,GSL.Vapor.mean.d2H,col="red")
d18O.V.mean <- GSL.Vapor.mean.d18O
d2H.V.mean <- GSL.Vapor.mean.d2H

#try to limit vapor data to the warm season
GSL.Vapor.data <- GSL.Vapor$data
GSL.Vapor.dates <- separate(GSL.Vapor.data, "Collection_Date", c("Date", "Time"), " ") 

GSL.Vapor.dates.m <- separate(GSL.Vapor.dates, "Date", c("Year", "Month", "Day"), "-")

GSL.Vapor.dates.m$Month <- as.integer(GSL.Vapor.dates.m$Month)

#select warm season data 
GSL.Vapor.warm <- subset(GSL.Vapor.dates.m, GSL.Vapor.dates.m$Month >4 & GSL.Vapor.dates.m$Month <11)
GSL.Vapor.warm.iso <- data.frame(GSL.Vapor.warm$d2H, GSL.Vapor.warm$d18O)
GSL.Vapor.warm.iso <- na.omit(GSL.Vapor.warm.iso)
GSL.Vapor.warm.mwl <- mwl(GSL.Vapor.warm.iso, plot = TRUE)
GSL.Vapor.warm.mwl[1]
GSL.Vapor.warm.mwl[2]
GSL.Vapor.warm.d18O <- mean(GSL.Vapor.warm.iso$GSL.Vapor.warm.d18O)
GSL.Vapor.warm.d2H <- mean(GSL.Vapor.warm.iso$GSL.Vapor.warm.d2H)
GSL.Vapor.warm.d18O
GSL.Vapor.warm.d2H
points(GSL.Vapor.warm.d18O, GSL.Vapor.warm.d2H, col = "orange", pch = 15)
#GSL is highly seasonal, so evaporation has to go up

#isoWater to get LMWL: extract intercept and slope for the prior distribution of runoff
GSL.precip <- wiDB_data(minLat = 39.5, maxLat = 42.5, minLong = -113.5,
                        maxLong = -110.5, minElev = NULL, maxElev = NULL, minDate = NULL,
                        maxDate = NULL, countries = "US", states = "UT", types = "Precipitation",
                        projects = NULL, fields = NULL, tmpdir = tempdir(), clean = TRUE)
#extract H and O isotopes from data
GSL.precip.iso <- data.frame(GSL.precip$data$d2H, GSL.precip$data$d18O)
GSL.precip.iso <- na.omit(GSL.precip.iso)

GSL.precip.mean.d18O <- mean(GSL.precip.iso$GSL.precip.data.d18O)
GSL.precip.mean.d2H <- mean(GSL.precip.iso$GSL.precip.data.d2H)
points(GSL.precip.mean.d18O,GSL.precip.mean.d2H,col="blue",pch = 15)
d18O.prec.mean <- GSL.precip.mean.d18O
d2H.prec.mean <- GSL.precip.mean.d2H

GSL.precip.mwl <- mwl(GSL.precip.iso, plot = TRUE)
GSL.precip.mwl

#try to limit vapor data to the warm season
GSL.precip.data <- GSL.precip$data
GSL.precip.dates <- separate(GSL.precip.data, "Collection_Date", c("Date", "Time"), " ") 

GSL.precip.dates.m <- separate(GSL.precip.dates, "Date", c("Year", "Month", "Day"), "-")

GSL.precip.dates.m$Month <- as.integer(GSL.precip.dates.m$Month)

#select warm season data 
GSL.precip.warm <- subset(GSL.precip.dates.m, GSL.precip.dates.m$Month >4 & GSL.precip.dates.m$Month <11)
GSL.precip.warm.iso <- data.frame(GSL.precip.warm$d2H, GSL.precip.warm$d18O)
GSL.precip.warm.iso <- na.omit(GSL.precip.warm.iso)
GSL.precip.warm.mwl <- mwl(GSL.precip.warm.iso, plot = TRUE)
GSL.precip.warm.mwl[1]
GSL.precip.warm.mwl[2]
GSL.precip.warm.d18O <- mean(GSL.precip.warm.iso$GSL.precip.warm.d18O)
GSL.precip.warm.d2H <- mean(GSL.precip.warm.iso$GSL.precip.warm.d2H)
GSL.precip.warm.d18O
GSL.precip.warm.d2H
points(GSL.precip.warm.d18O, GSL.precip.warm.d2H, col = "blue", pch = 16)
#GSL is highly seasonal, so evaporation has to go up

#Get GSL lake water d18O and d2H to extract intercept and slope for the prior distribution of lake water
GSL.lake <- wiDB_data(minLat = 40.6, maxLat = 41.7, minLong = -113,
                        maxLong = -111.9, minElev = NULL, maxElev = NULL, minDate = NULL,
                        maxDate = NULL, countries = "US", states = "UT", types = "Lake",
                        projects = NULL, fields = NULL, tmpdir = tempdir(), clean = TRUE)
GSL.lake.iso <- data.frame(GSL.lake$data$d2H, GSL.lake$data$d18O)
GSL.lake.iso <- na.omit(GSL.lake.iso)

GSL.lake.mwl <- mwl(GSL.lake.iso, plot = TRUE)
GSL.lake.mwl #slope, intercept, meanO, sso, rmse, ns

GSL.mean.d18O <- mean(GSL.lake.iso$GSL.lake.data.d18O)
GSL.mean.d2H <- mean(GSL.lake.iso$GSL.lake.data.d2H)
GSL.sd.d18O <- sd(GSL.lake.iso$GSL.lake.data.d18O)
GSL.sd.d2H <- sd(GSL.lake.iso$GSL.lake.data.d2H)
GSL.OHcov <- cov(GSL.lake.iso$GSL.lake.data.d18O, GSL.lake.iso$GSL.lake.data.d2H)

#use the mean, sd, and covariance to summarize GLS lake water data
GSL.obs <- iso(GSL.mean.d2H, GSL.mean.d18O, GSL.sd.d2H, GSL.sd.d18O, GSL.OHcov)

#possible sources of GSL, this can be built into the model
GSL.source.2 <- mwlSource(obs=GSL.obs, MWL = GSL.precip.mwl, slope = c(GSL.lake.mwl[1],0.5), stype=2,
                        ngens=1000, ncores = 3)

points(GSL.source.2$results$source_d18O,GSL.source.2$results$source_d2H,pch=16,col="red")




#########################Get geochronology of GSL core#####################################
#use bacon
library(rbacon)


#use the posterior of GSL source as the prior for runoff, also prescribe cov

#check salinity factor
# plot(GSL.level,1 -(0.788*GSL.sali/(1000+0.63*GSL.sali)))
# plot(GSL.level,(1 -(0.788*GSL.sali/(1000+0.63*GSL.sali)))*GSL.area/1000*1.3)
# 
# #try activity coefficient curve
# f.residuw <- seq(0,1,by=0.01)
# act.coef <- 0.99931 - 0.018521*f.residuw - 0.000543*f.residuw^2
# plot(f.residuw, act.coef)
# #almost linear one