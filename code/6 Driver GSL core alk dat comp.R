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
GSL.alk.raw <- read.csv("data/GSL_alkanes_alt.csv")

# GSL.alk.no <- na.omit(GSL.alk.raw)

#short chains represent lake water
GSL.alk.C17 <- GSL.alk.raw %>% filter(compound=="C17")
#calculate mean and sd of the 3 repeats 
GSL.alk.C17.summ <- GSL.alk.C17%>% group_by(Depth_cm) %>% summarize(C17.avg = mean(d2H),C17.sd = sd(d2H))

GSL.alk.C17.naom <- na.omit(GSL.alk.C17.summ)

GSL.alk.C19 <- GSL.alk.raw %>% filter(compound=="C19")
#calculate mean and sd of the 3 repeats 
GSL.alk.C19.summ <- GSL.alk.C19%>% group_by(Depth_cm) %>% summarize(C19.avg = mean(d2H),C19.sd = sd(d2H))

GSL.alk.C19.naom <- na.omit(GSL.alk.C19.summ)

#mid chain C25 represent a mixture of lake water and runoff (mostly runoff)
GSL.alk.C25 <- GSL.alk.raw %>% filter(compound=="C25")
#calculate mean and sd of the 3 repeats 
GSL.alk.C25.summ <- GSL.alk.C25%>% group_by(Depth_cm) %>% summarize(C25.avg = mean(d2H),C25.sd = sd(d2H))

GSL.alk.C25.naom <- na.omit(GSL.alk.C25.summ)

#long chain represent terrestrial plants from runoff
GSL.alk.C29 <- GSL.alk.raw %>% filter(compound=="C29")
#calculate mean and sd of the 3 repeats 
GSL.alk.C29.summ <- GSL.alk.C29%>% group_by(Depth_cm) %>% summarize(C29.avg = mean(d2H),C29.sd = sd(d2H))

GSL.alk.C29.naom <- na.omit(GSL.alk.C29.summ)

#get posterior ages for depths from rbacon
GSL.scwax.age <- sapply(GSL.alk.C19.naom$Depth_cm,Bacon.Age.d)

GSL.lcwax.age <- sapply(GSL.alk.C29.naom$Depth_cm,Bacon.Age.d)

dim(GSL.scwax.age) #row = iteration, col = depths
dim(GSL.lcwax.age)

#preliminary plot
plot(GSL.alk.C17.summ$Depth_cm, GSL.alk.C17.summ$avg,type="l",col = "blue", ylim=c(-250,-100))
lines(GSL.alk.C19.summ$Depth_cm, GSL.alk.C19.summ$avg,col="green")
lines(GSL.alk.C25.summ$Depth_cm, GSL.alk.C25.summ$avg,col="orange")
lines(GSL.alk.C29.summ$Depth_cm, GSL.alk.C29.summ$avg,col="red")