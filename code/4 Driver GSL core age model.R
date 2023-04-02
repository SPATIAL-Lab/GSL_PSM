library(rbacon)

setwd("C:/Users/ydmag/Google Drive/U of U/GSL proxy/GSL_PSM")

source("code/1 Helper functions.R")

#########################Get geochronology of GSL core#####################################
#use rbacon

GSL.1B.14C <- read.csv("data/GSL_core1B_14C.csv")

Bacon.cleanup()

mydir <- "C:/Users/ydmag/Google Drive/U of U/GSL proxy/GSL_PSM/data/"

max.depth = 1300 #cm

# Bacon(coredir = mydir, core = "GSL_1B",d.min = 33, d.max = max.depth, d.by=0.5, cc = 1)#age per 0.5 cm

Bacon(coredir = mydir, core = "GSL_1B",d.min = 33, d.max = max.depth, d.by=0.5, cc = 1,run = FALSE)#age per 0.5 cm
#simple solution: get posterior ages for each depth for all data types

#The posterior ages will be used as JAGS input, then for each data type, sample from the ages

#use the posterior of GSL source as the prior for runoff, also prescribe cov
