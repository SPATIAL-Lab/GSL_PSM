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
library(dplyr)

setwd("C:/Users/ydmag/Google Drive/U of U/GSL proxy/GSL_PSM")

########## Get organic geochem lipid d2H
GSL.alk.raw <- read.csv("data/GSL_alkanes.csv")

#short chains represent lake water
GSL.alk.C17 <- GSL.alk.raw %>% filter(compound=="C17")
#calculate mean and sd of the 3 repeats 
GSL.alk.C17.summ <- GSL.alk.C17%>% group_by(Depth_cm) %>% summarize(avg = mean(d2H),sd = sd(d2H))

GSL.alk.C19 <- GSL.alk.raw %>% filter(compound=="C19")
#calculate mean and sd of the 3 repeats 
GSL.alk.C19.summ <- GSL.alk.C19%>% group_by(Depth_cm) %>% summarize(avg = mean(d2H),sd = sd(d2H))

#mid chain C25 represent a mixture of lake water and runoff (mostly runoff)
GSL.alk.C25 <- GSL.alk.raw %>% filter(compound=="C25")
#calculate mean and sd of the 3 repeats 
GSL.alk.C25.summ <- GSL.alk.C25%>% group_by(Depth_cm) %>% summarize(avg = mean(d2H),sd = sd(d2H))

#long chain represent terrestrial plants from runoff
GSL.alk.C29 <- GSL.alk.raw %>% filter(compound=="C29")
#calculate mean and sd of the 3 repeats 
GSL.alk.C29.summ <- GSL.alk.C29%>% group_by(Depth_cm) %>% summarize(avg = mean(d2H),sd = sd(d2H))

#get posterior ages for depths from rbacon
GSL.alk.C17.summ$Depth_cm

GSL.alk.age <- sapply(GSL.alk.C17.summ$Depth_cm,Bacon.Age.d)

dim(GSL.cyst.age) #row = iteration, col = depths


#preliminary plot
plot(GSL.alk.C17.summ$Depth_cm, GSL.alk.C17.summ$avg,type="l",col = "blue", ylim=c(-250,-100))
lines(GSL.alk.C19.summ$Depth_cm, GSL.alk.C19.summ$avg,col="green")
lines(GSL.alk.C25.summ$Depth_cm, GSL.alk.C25.summ$avg,col="orange")
lines(GSL.alk.C29.summ$Depth_cm, GSL.alk.C29.summ$avg,col="red")