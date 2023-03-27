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

#get posterior ages for depths from rbacon
GSL.cyst$Depth

GSL.cyst.age <- sapply(GSL.cyst$Depth,Bacon.Age.d)

dim(GSL.cyst.age) #row = iteration, col = depths
