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

########## Get organic geochem lipid d2H################
GSL.FAME.raw <- read.csv("data/GSL_FAME_alt.csv")

d2H.meOH <- -123.7 #d2H of added methyl group, Huang et al 2002
#short chains represent lake water

######C18
GSL.FAME.C18 <- GSL.FAME.raw %>% filter(compound=="C18")
C18 <- 18

#mass balance correction for FAME values
GSL.FAME.C18.corr <- GSL.FAME.C18 %>% 
  mutate(d2H.corr = (d2H * (2 * C18 + 2) - 3 * d2H.meOH)/(2 * C18 - 1))


GSL.FAME.C18.naom <- na.omit(GSL.FAME.C18.summ)

#long chain represent terrestrial plants from runoff
GSL.FAME.C28 <- GSL.FAME.raw %>% filter(compound=="C28")
C28 <- 28

#mass balance correction for FAME values
GSL.FAME.C28.corr <- GSL.FAME.C28 %>% 
  mutate(d2H.corr = (d2H * (2 * C28 + 2) - 3 * d2H.meOH)/(2 * C28 - 1))

#aggregate data by depth 
GSL.FAME.C18.summ <- GSL.FAME.C18.corr%>% group_by(Depth_cm) %>%
  summarize(C18.avg = mean(d2H.corr),C18.sd = sd(d2H.corr),PA.18 = mean(area2+area3))

#aggregate data by depth
GSL.FAME.C28.summ <- GSL.FAME.C28.corr%>% group_by(Depth_cm) %>%
  summarize(C28.avg = mean(d2H.corr),C28.sd = sd(d2H.corr),PA.28 = mean(area2+area3))

#check sd for each depth
print(cbind(GSL.FAME.C18.summ$C18.sd,GSL.FAME.C18.summ$Depth_cm))

GSL.FAME.C18.naom <- na.omit(GSL.FAME.C18.summ)

GSL.FAME.C28.naom <- na.omit(GSL.FAME.C28.summ)

#get posterior ages for depths from rbacon
GSL.scwax.age <- sapply(GSL.FAME.C16.naom$Depth_cm,Bacon.Age.d)

GSL.lcwax.age <- sapply(GSL.FAME.C28.naom$Depth_cm,Bacon.Age.d)

dim(GSL.scwax.age) #row = iteration, col = depths
dim(GSL.lcwax.age)


# GSL.FAME.C16 <- GSL.FAME.raw %>% filter(compound=="C16")
# C16 <- 16
# 
# #mass balance correction for FAME values
# GSL.FAME.C16.corr <- GSL.FAME.C16 %>% 
#   mutate(d2H.corr = (d2H * (2 * C16 + 2) - 3 * d2H.meOH)/(2 * C16 - 1))
# 
# GSL.FAME.C16.naom <- na.omit(GSL.FAME.C16.summ)

#aggregate data by depth
# GSL.FAME.C16.summ <- GSL.FAME.C16.corr%>% group_by(Depth_cm) %>%
#   summarize(C16.avg = mean(d2H.corr),C16.sd = sd(d2H.corr),PA.16 = mean(area2+area3))