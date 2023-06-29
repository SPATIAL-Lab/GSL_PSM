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

#source rbacon output
Bacon(coredir = mydir, core = "GSL_1B",run = F)

#########################Get GSL core data#####################################
#####1##### Get bulk carbonate d18O
GSL.carb <- read.csv("data/GSL carbonates.csv")

#get posterior ages for depths from rbacon
GSL.carb$Depth

GSL.carb.age <- sapply(GSL.carb$Depth,Bacon.Age.d)

dim(GSL.carb.age) #row = iteration, col = depths

#####2##### Get Artimia cyst geochem
GSL.cyst <- read.csv("data/GSL BS cysts.csv")
GSL.cyst <- GSL.cyst %>% arrange(GSL.cyst$Depth)

#get posterior ages for depths from rbacon
GSL.cyst$Depth

GSL.cyst.age <- sapply(GSL.cyst$Depth,Bacon.Age.d)

dim(GSL.cyst.age) #row = iteration, col = depths

#calculate accumulation rate at cyst depths


acc.rate.cyst <- rep(0, length(GSL.cyst$Depth))
for(i in 1:length(GSL.cyst$Depth)){
  acc.rate.cyst[i] <- median(accrate.depth(GSL.cyst$Depth[i]))
}

plot(GSL.cyst$Depth, acc.rate.cyst)

#yield = mg/g of sediment
#acc.rate = depth/year or g/year
#so scaled cyst yield should be yield * acc.rate = mg(cyst)/year

cyst.yield.scl <- GSL.cyst$yield * acc.rate.cyst

plot(GSL.cyst$Depth, log(cyst.yield.scl)) #cyst yield generally decreased over time, but there are skipes
