library(rbacon)

setwd("C:/Users/ydmag/Google Drive/U of U/GSL proxy/GSL_PSM")

source("code/1 Helper functions.R")

#########################Get geochronology of GSL core#####################################
#use rbacon

GSL.1B.14C <- read.csv("data/GSL_1B/GSL_1B.csv")

Bacon.cleanup()

mydir <- "C:/Users/ydmag/Google Drive/U of U/GSL proxy/GSL_PSM/data/"

max.depth = 1300 #cm

Bacon(coredir = mydir, core = "GSL_1B",d.min = 33, d.max = max.depth, d.by=0.5, cc = 1)#age per 0.5 cm
#select y for the popups

# Bacon(coredir = mydir, core = "GSL_1B",d.min = 33, d.max = max.depth, d.by=0.5, cc = 1,run = FALSE)#age per 0.5 cm

agedepth(rotate.axes=TRUE,rev.age=TRUE)

###compare to previous age model with intcal13 calibration
library(rintcal)

ccurve(cc = "IntCal13", postbomb = FALSE)

intcal <- intcal.read.data()

intcal.data(cal1,
            cal2,
            cc1 = "IntCal13")
myintcal <- tempfile()
intcal.write.data(intcal, myintcal)

#####use cc = 5 to recreate IntCal 13 results!

#simple solution: get posterior ages for each depth for all data types

#The posterior ages will be used as JAGS input, then for each data type, sample from the ages

#use the posterior of GSL source as the prior for runoff, also prescribe cov
