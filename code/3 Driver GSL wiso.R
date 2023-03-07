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

