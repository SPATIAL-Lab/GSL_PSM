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


################use synthetic data to test lake mass balance#####
#age resolution
age.res = 10 #this determins the number of time steps in the simulation

Bacon.Age.d(max.depth)#cm, this is consistent with rbacon input

max.age = max(Bacon.Age.d(max.depth)) #use the maximum age to constrain simulation 

#number of time steps in the simulation
t <- ceiling(max.age/age.res) + 200 #1006 time steps + 200 for buffer

#######use calibrations for short chain wax fractionation
# post.leng.scwax <- length(post)
# 
# post.scwax.alpha.sl <- post
# post.scwax.alpha.inc <- post

#######compile priors#######

#supply slopes and intercepts
MWL.slope <- GSL.lake.mwl[1]
MWL.intc <- GSL.lake.mwl[2]

V.intc <- GSL.Vapor.warm.mwl[1]
V.slope <- GSL.Vapor.warm.mwl[2]

Lw.slope <- GSL.lake.mwl[1]
Lw.intc <- GSL.lake.mwl[2]

#warm season vapor d18O
d18O.vap.warm <- GSL.Vapor.warm.d18O

#rbacon MCMC iternations
n.bac.it <- nrow(GSL.carb.age)

####### compile measured data#####

d18O.car.dat <- GSL.carb$Ave.d18O

d18O.car.sd <- GSL.carb$SD.d18O

lcwax.d2H.dat <- GSL.alk.C29.naom$C29.avg

lcwax.d2H.sd <- GSL.alk.C29.naom$C29.sd

scwax.d2H.dat <- GSL.alk.C19.naom$C19.avg

scwax.d2H.sd <- GSL.alk.C19.naom$C19.sd

BScyst.d18O.dat <- GSL.cyst$Ave.d18O

BScyst.d18O.sd <- rep(0.75, n.cyst)#fixed precision

BScyst.d2H.dat <- GSL.cyst$Ave.d2H

BScyst.d2H.sd <- rep(2, n.cyst) #fixed precision

#carbonate ages
post.age.carb <- GSL.carb.age

n.carb <- ncol(GSL.carb.age)

#cyst ages
post.age.cyst <- GSL.cyst.age

n.cyst <- ncol(GSL.cyst.age)

#lcwax ages
post.age.lcwax <- GSL.lcwax.age

n.lc.wax <- ncol(GSL.lcwax.age)

#lcwax ages
post.age.scwax <- GSL.scwax.age

n.sc.wax <- ncol(GSL.scwax.age)

#parameters to save
parameters <- c("L.level","rh", "nsws", "LST", "T.gap","AT","Runoff", "sal", "S.coeff",
                "E.rate", "Evap", "LV.A", "LST.cps.ac","Lw.d2H","Lw.d18O","f",
                "Ro.d2H", "Ro.d18O", "air.d2H","air.d18O","evap.d2H","evap.d18O",
                "BScyst.d2H", "BScyst.d18O", "scwax.d2H", "lcwax.d2H", "d18O.car",
                "BScyst.slope.lw.2H", "BScyst.slope.diet.2H", "BScyst.inc.2H",
                "BScyst.slope.lw.18O", "BScyst.slope.diet.18O", "BScyst.inc.18O",
                "epsilon.2H.carbohy", "epsilon.18O.carbohy", "epsilon.2H.AkAcid",
                "scwax.alpha.sl", "scwax.alpha.inc", "lcwax.d2H.inc", "lcwax.d2H.slope",
                "d2H.gap.MAP_Ro", "alpha2H.sc.alkane", "Ro.intc.m", "Ro.slope.m",
                "carb.age.gap")


dat = list(GSL.level = GSL.level.1286, GSL.area = GSL.area.1286, 
           GSL.volume = GSL.volume.1286, GSL.sali = GSL.sali, t = t, age.res = age.res,
           V.intc = V.intc, V.slope = V.slope, MWL.intc = MWL.intc, MWL.slope = MWL.slope,
           Lw.intc = Lw.intc, Lw.slope = Lw.slope,
           d18O.vap.warm = d18O.vap.warm,
           n.bac.it = n.bac.it, n.carb = n.carb, n.cyst = n.cyst, n.lc.wax = n.lc.wax, n.sc.wax = n.sc.wax,
           post.age.carb = post.age.carb, post.age.cyst = post.age.cyst, 
           d18O.car.dat = d18O.car.dat, lcwax.d2H.dat = lcwax.d2H.dat, 
           post.age.scwax = post.age.scwax, post.age.lcwax = post.age.lcwax,
           scwax.d2H.dat = scwax.d2H.dat, BScyst.d18O.dat = BScyst.d18O.dat, 
           BScyst.d2H.dat = BScyst.d2H.dat,
           d18O.car.sd = d18O.car.sd, lcwax.d2H.sd = lcwax.d2H.sd, 
           scwax.d2H.sd = scwax.d2H.sd, BScyst.d18O.sd = BScyst.d18O.sd, 
           BScyst.d2H.sd = BScyst.d2H.sd)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e3
n.burnin = 1e3
n.thin = 1

#Run it
SI.Arc = do.call(jags.parallel,list(model.file = "code/JAGS PSM Arc.R", 
                                     parameters.to.save = parameters, 
                                     data = dat, n.chains=5, n.iter = n.iter, 
                                     n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #150 seconds for t = 1206, n.iter = 2e3
#estimated time taken: t = 1206 (decadal variation), 1.4 hours for n.iter = 2e4

SI.Arc$BUGSoutput$summary

#check relative humidity
plot(density(SI.Arc$BUGSoutput$sims.list$rh)) #getting low

#check advection coefficient
plot(density(SI.Arc$BUGSoutput$sims.list$f)) #

#check wind speed
plot(density(SI.Arc$BUGSoutput$sims.list$nsws)) #normal

#check LST ac
plot(density(SI.Arc$BUGSoutput$sims.list$LST.cps.ac)) #all over

#check atm vapor d18O
plot(density(SI.Arc$BUGSoutput$sims.list$air.d18O)) 

#check LST
SI.Arc.LST <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$LST,0.89)
plot(1:t,SI.Arc.LST[[1]],type="l",ylim=c(15,30), main="LST")
lines(1:t,SI.Arc.LST[[2]],lty=2)
lines(1:t,SI.Arc.LST[[3]],lty=2)
#by far, temperature has the least constraint

#check AT
SI.Arc.AT <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$AT,0.89)
plot(1:t,test.lwmb.AT[[1]],type="l",ylim=c(15,30), main="LST")
lines(1:t,test.lwmb.AT[[2]],lty=2)
lines(1:t,test.lwmb.AT[[3]],lty=2)

#check Runoff and evaporation
SI.lwmb.Runoff <- MCMC.CI.bound(SI.lwmb$BUGSoutput$sims.list$Runoff,0.89)
plot(1:t,SI.lwmb.Runoff[[1]],type="l",ylim=c(0,12), main="Runoff and Evap")
lines(1:t,SI.lwmb.Runoff[[2]],lty=2)
lines(1:t,SI.lwmb.Runoff[[3]],lty=2)

#check evaporation
SI.lwmb.Evap <- MCMC.CI.bound(SI.lwmb$BUGSoutput$sims.list$Evap,0.89)
lines(1:t,SI.lwmb.Evap[[1]],col="red")
lines(1:t,SI.lwmb.Evap[[2]],lty=2,col="red")
lines(1:t,SI.lwmb.Evap[[3]],lty=2,col="red")

#check evaporation rate
SI.lwmb.Evaprate <- MCMC.CI.bound(SI.lwmb$BUGSoutput$sims.list$E.rate,0.89)
plot(1:t,SI.lwmb.Evaprate[[1]],type="l",main="Evap rate", ylim=c(0,1))
lines(1:t,SI.lwmb.Evaprate[[2]],lty=2)
lines(1:t,SI.lwmb.Evaprate[[3]],lty=2)

#check lake level
SI.lwmb.L.level <- MCMC.CI.bound(SI.lwmb$BUGSoutput$sims.list$L.level,0.89)
plot(1:t,SI.lwmb.L.level[[1]],type="l",ylim=c(1275,1285), main="Lake level")
lines(1:t,SI.lwmb.L.level[[2]],lty=2)
lines(1:t,SI.lwmb.L.level[[3]],lty=2)

#check salinity
SI.lwmb.sal <- MCMC.CI.bound(SI.lwmb$BUGSoutput$sims.list$sal,0.89)
plot(1:t,SI.lwmb.sal[[1]],type="l",ylim=c(50,355), main="Salinity")
lines(1:t,SI.lwmb.sal[[2]],lty=2)
lines(1:t,SI.lwmb.sal[[3]],lty=2)

#check runoff isotopes
SI.lwmb.Rod18O <- MCMC.CI.bound(SI.lwmb$BUGSoutput$sims.list$Ro.d18O,0.89)
plot(1:t,SI.lwmb.Rod18O[[1]],type="l",ylim=c(-30,-10), main="Runoff d18O")
lines(1:t,SI.lwmb.Rod18O[[2]],lty=2)
lines(1:t,SI.lwmb.Rod18O[[3]],lty=2)

#check evaporation isotopes
SI.lwmb.evapd18O <- MCMC.CI.bound(SI.lwmb$BUGSoutput$sims.list$evap.d18O,0.89)
plot(1:t,SI.lwmb.evapd18O[[1]],type="l",ylim=c(-30,-10), main="Evap d18O")
lines(1:t,SI.lwmb.evapd18O[[2]],lty=2)
lines(1:t,SI.lwmb.evapd18O[[3]],lty=2)

#check lake water isotopes
SI.lwmb.Lw.d18O <- MCMC.CI.bound(SI.lwmb$BUGSoutput$sims.list$Lw.d18O,0.89)
plot(1:t,SI.lwmb.Lw.d18O[[1]],type="l",ylim=c(-8,-2), main="Lake water d18O")
lines(1:t,SI.lwmb.Lw.d18O[[2]],lty=2)
lines(1:t,SI.lwmb.Lw.d18O[[3]],lty=2)