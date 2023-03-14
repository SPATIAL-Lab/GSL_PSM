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

#########################Get geochronology of GSL core#####################################
#use bacon
library(rbacon)

GSL.1B.14C <- read.csv("data/GSL_core1B_14C.csv")

#use the posterior of GSL source as the prior for runoff, also prescribe cov
