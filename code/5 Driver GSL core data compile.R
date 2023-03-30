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

#preliminary plots
par(mfrow=c(1,1))
plot(GSL.carb$Depth, GSL.carb$Ave.d18O,type="l",col = "blue", ylim=c(-9,19))
lines(GSL.cyst$Depth, GSL.cyst$Ave.d18O,col="green")

plot(GSL.cyst$Depth, GSL.cyst$Ave.d2H,type="l",col = "blue", ylim=c(-200,-100))
lines(GSL.alk.C17.summ$Depth_cm, GSL.alk.C17.summ$avg,col="green")

plot(GSL.cyst$Depth, GSL.cyst$Ave.d2H,type="l",col = "blue", ylim=c(-200,-100))
lines(GSL.alk.C19.summ$Depth_cm, GSL.alk.C19.summ$avg,col="green")

plot(GSL.cyst$Depth, GSL.cyst$Ave.d2H,type="l",col = "blue", ylim=c(-220,-100))
lines(GSL.alk.C25.summ$Depth_cm, GSL.alk.C25.summ$avg,col="green")

plot(GSL.carb$Depth, GSL.carb$Ave.d13C,type="l",col = "blue", ylim=c(-60,10)) #carbonate d13C: source of DIC terrestrial vs lacustrine
#carbon used as runoff indicator? Erossion of old carbon
#but oxygen record doesn't make sense!
lines(GSL.alk.C25.summ$Depth_cm, GSL.alk.C25.summ$avg-GSL.alk.C17.summ$avg,col="green")
