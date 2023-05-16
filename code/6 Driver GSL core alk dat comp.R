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

d2H.meOH <- -123
#short chains represent lake water
GSL.FAME.C16 <- GSL.FAME.raw %>% filter(compound=="C16")
C16 <- 16

#mass balance correction for FAME values
GSL.FAME.C16.corr <- GSL.FAME.C16 %>% 
  mutate(d2H.corr = (d2H * (2 * C16 + 2) - 3 * d2H.meOH)/(2 * C16 - 1))

#calculate mean and sd of the 3 repeats 
GSL.FAME.C16.summ <- GSL.FAME.C16.corr%>% group_by(Depth_cm) %>%
  summarize(C16.avg = mean(d2H.corr),C16.sd = sd(d2H.corr),PA.16 = mean(area2+area3))

GSL.FAME.C16.naom <- na.omit(GSL.FAME.C16.summ)

######C18
GSL.FAME.C18 <- GSL.FAME.raw %>% filter(compound=="C18")
C18 <- 18

#mass balance correction for FAME values
GSL.FAME.C18.corr <- GSL.FAME.C18 %>% 
  mutate(d2H.corr = (d2H * (2 * C18 + 2) - 3 * d2H.meOH)/(2 * C18 - 1))

#calculate mean and sd of the 3 repeats 
GSL.FAME.C18.summ <- GSL.FAME.C18.corr%>% group_by(Depth_cm) %>%
  summarize(C18.avg = mean(d2H.corr),C18.sd = sd(d2H.corr),PA.18 = mean(area2+area3))

GSL.FAME.C18.naom <- na.omit(GSL.FAME.C18.summ)

#long chain represent terrestrial plants from runoff
GSL.FAME.C28 <- GSL.FAME.raw %>% filter(compound=="C28")
C28 <- 28

#mass balance correction for FAME values
GSL.FAME.C28.corr <- GSL.FAME.C28 %>% 
  mutate(d2H.corr = (d2H * (2 * C28 + 2) - 3 * d2H.meOH)/(2 * C28 - 1))

#calculate mean and sd of the 3 repeats 
GSL.FAME.C28.summ <- GSL.FAME.C28.corr%>% group_by(Depth_cm) %>%
  summarize(C28.avg = mean(d2H.corr),C28.sd = sd(d2H.corr),PA.28 = mean(area2+area3))

GSL.FAME.C28.naom <- na.omit(GSL.FAME.C28.summ)

#get posterior ages for depths from rbacon
GSL.scwax.age <- sapply(GSL.FAME.C16.naom$Depth_cm,Bacon.Age.d)

GSL.lcwax.age <- sapply(GSL.FAME.C28.naom$Depth_cm,Bacon.Age.d)

dim(GSL.scwax.age) #row = iteration, col = depths
dim(GSL.lcwax.age)

#preliminary plot
plot(GSL.FAME.C16.summ$Depth_cm, GSL.FAME.C16.summ$C16.avg,type="l",col = "blue", ylim=c(-220,-100))
lines(GSL.FAME.C28.summ$Depth_cm, GSL.FAME.C28.summ$C28.avg,col="red")
lines(GSL.FAME.C18.summ$Depth_cm, GSL.FAME.C18.summ$C18.avg,col="green")

plot(GSL.FAME.C16.naom$Depth_cm,ratio.17.n, type = "l")

# GSL.alk.raw <- read.csv("data/GSL_alkanes_alt.csv")
# 
# # GSL.alk.no <- na.omit(GSL.alk.raw)
# 
# #short chains represent lake water
# GSL.alk.C17 <- GSL.alk.raw %>% filter(compound=="C17")
# #calculate mean and sd of the 3 repeats 
# GSL.alk.C17.summ <- GSL.alk.C17%>% group_by(Depth_cm) %>%
#   summarize(C17.avg = mean(d2H),C17.sd = sd(d2H),PA.17 = mean(area2+area3))
# 
# GSL.alk.C17.naom <- na.omit(GSL.alk.C17.summ)
# 
# # GSL.alk.C19 <- GSL.alk.raw %>% filter(compound=="C19")
# # #calculate mean and sd of the 3 repeats 
# # GSL.alk.C19.summ <- GSL.alk.C19%>% group_by(Depth_cm) %>% 
# #   summarize(C19.avg = mean(d2H),C19.sd = sd(d2H),PA.19 = mean(area2+area3))
# 
# #long chain represent terrestrial plants from runoff
# GSL.alk.C29 <- GSL.alk.raw %>% filter(compound=="C29")
# #calculate mean and sd of the 3 repeats 
# GSL.alk.C29.summ <- GSL.alk.C29%>% group_by(Depth_cm) %>% 
#   summarize(C29.avg = mean(d2H),C29.sd = sd(d2H),PA.29 = mean(area2+area3))
# 
# # sc.lc.wax <- cbind(GSL.alk.C17.summ,GSL.alk.C19.summ,GSL.alk.C29.summ)
# 
# GSL.alk.C29.naom <- na.omit(GSL.alk.C29.summ)

# plot(GSL.alk.C17.summ$Depth_cm, GSL.alk.C17.summ$PA.17/GSL.alk.C19.summ$PA.19)#C17 to C19 ratio
# 
# ratio.17.19 <- GSL.alk.C17.naom$PA.17/GSL.alk.C19.naom$PA.19
# 
# ratio.17.n <- GSL.alk.C17.naom$PA.17/(GSL.alk.C19.naom$PA.19+GSL.alk.C17.naom$PA.17)
# 
# plot(ratio.17.n,GSL.alk.C17.naom$C17.avg-GSL.alk.C19.naom$C19.avg)
# 
# 
# plot(GSL.alk.C17.naom$C17.avg-GSL.alk.C19.naom$C19.avg,ratio.17.19)#/??
# 
# plot(GSL.alk.C19.naom$C19.avg,ratio.17.19)#/??
# 
# GSL.alk.C19.naom <- na.omit(GSL.alk.C19.summ)
# 
# sc.pw.d2H <- (GSL.alk.C19.naom$C19.avg*GSL.alk.C19.naom$PA.19 + GSL.alk.C17.naom$C17.avg*GSL.alk.C17.naom$PA.17)/
#   (GSL.alk.C17.naom$PA.17+GSL.alk.C19.naom$PA.19)

# #mid chain C23 represent a mixture of lake water and runoff (mostly runoff)
# GSL.alk.C23 <- GSL.alk.raw %>% filter(compound=="C23")
# #calculate mean and sd of the 3 repeats 
# GSL.alk.C23.summ <- GSL.alk.C23%>% group_by(Depth_cm) %>% 
#   summarize(C23.avg = mean(d2H),C23.sd = sd(d2H),PA.23 = mean(area2+area3))
# 
# GSL.alk.C23.naom <- na.omit(GSL.alk.C23.summ)
# 
# #mid chain C25 represent a mixture of lake water and runoff (mostly runoff)
# GSL.alk.C25 <- GSL.alk.raw %>% filter(compound=="C25")
# #calculate mean and sd of the 3 repeats 
# GSL.alk.C25.summ <- GSL.alk.C25%>% group_by(Depth_cm) %>% 
#   summarize(C25.avg = mean(d2H),C25.sd = sd(d2H),PA.25 = mean(area2+area3))
# 
# GSL.alk.C25.naom <- na.omit(GSL.alk.C25.summ)
# 
# GSL.alk.raw %>% group_by(Depth_cm) %>% filter(compound=="C29")
# 
# GSL.alk.raw %>% group_by(Depth_cm) %>% filter(compound=="C17")

# #long chain represent terrestrial plants from runoff
# GSL.alk.C31 <- GSL.alk.raw %>% filter(compound=="C31")
# #calculate mean and sd of the 3 repeats 
# GSL.alk.C31.summ <- GSL.alk.C31%>% group_by(Depth_cm) %>% 
#   summarize(C31.avg = mean(d2H),C31.sd = sd(d2H),PA.31 = mean(area2+area3))
# 
# GSL.alk.C31.naom <- na.omit(GSL.alk.C31.summ)

#try to assess wetland vs terrestrial input using Paq?

#get posterior ages for depths from rbacon
# GSL.scwax.age <- sapply(GSL.alk.C17.naom$Depth_cm,Bacon.Age.d)
# 
# GSL.lcwax.age <- sapply(GSL.alk.C29.naom$Depth_cm,Bacon.Age.d)
# 
# dim(GSL.scwax.age) #row = iteration, col = depths
# dim(GSL.lcwax.age)
# 
# #preliminary plot
# plot(GSL.alk.C17.summ$Depth_cm, GSL.alk.C17.summ$C17.avg,type="l",col = "blue", ylim=c(-250,-100))
# # lines(GSL.alk.C19.summ$Depth_cm, GSL.alk.C19.summ$C19.avg,col="green")
# # lines(GSL.alk.C25.summ$Depth_cm, GSL.alk.C25.summ$C25.avg,col="orange")
# lines(GSL.alk.C29.summ$Depth_cm, GSL.alk.C29.summ$C29.avg,col="red")
# 
# plot(GSL.alk.C23.naom$Depth_cm, GSL.alk.C23.naom$C23.avg,col="blue", type = "l")
# 
# plot(GSL.alk.C17.naom$Depth_cm,ratio.17.n, type = "l")
