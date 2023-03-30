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
age.res = 100 #this determins the number of time steps in the simulation

Bacon.Age.d(1300)#cm, this is consistent with rbacon input

max.age = max(Bacon.Age.d(max.depth)) #use the maximum age to constrain simulation 

#number of time steps in the simulation
t <- ceiling(max.age/age.res) #+ 20 #101 time steps + 20 for buffer

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

####### compile measured data#####

d18O.car.dat <- GSL.carb$Ave.d18O

d18O.car.sd <- GSL.carb$SD.d18O

lcwax.d2H.dat <- GSL.alk.C29.naom$C29.avg

lcwax.d2H.sd <- rep(2.5, n.lc.wax)

scwax.d2H.dat <- GSL.alk.C19.naom$C19.avg

scwax.d2H.sd <- rep(2.5, n.sc.wax)

BScyst.d18O.dat <- GSL.cyst$Ave.d18O

BScyst.d18O.sd <- rep(0.75, n.cyst)#fixed precision

BScyst.d2H.dat <- GSL.cyst$Ave.d2H

BScyst.d2H.sd <- rep(2, n.cyst) #fixed precision



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
n.iter = 1e4
n.burnin = 5e3
n.thin = 2

#Run it
SI.Arc = do.call(jags.parallel,list(model.file = "code/JAGS PSM Arc.R", 
                                     parameters.to.save = parameters, 
                                     data = dat, n.chains=5, n.iter = n.iter, 
                                     n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 
#estimated time taken: t = 131 (centenial variation), 2.3 hours for n.iter = 2e4

SI.Arc$BUGSoutput$summary

par(mfrow=c(3,3))
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

#check brine shrimp cyst intercept
plot(density(SI.Arc$BUGSoutput$sims.list$BScyst.inc.2H)) 
abline(v=-92)

#check brine shrimp cyst slopes
plot(density(SI.Arc$BUGSoutput$sims.list$BScyst.slope.lw.2H)) 
abline(v=0.34)

#check brine shrimp cyst intercept
plot(density(SI.Arc$BUGSoutput$sims.list$BScyst.slope.diet.2H)) 
abline(v=0.26)

#check brine shrimp cyst intercept
plot(density(SI.Arc$BUGSoutput$sims.list$BScyst.inc.18O)) 
abline(v=15.9)
#check brine shrimp cyst slopes
plot(density(SI.Arc$BUGSoutput$sims.list$BScyst.slope.lw.18O)) 
abline(v=0.692)
#check brine shrimp cyst intercept
plot(density(SI.Arc$BUGSoutput$sims.list$BScyst.slope.diet.18O)) 
abline(v=0.101)

plot(density(SI.Arc$BUGSoutput$sims.list$d2H.gap.MAP_Ro)) 
abline(v=45)

#epsilon BS diet (smaller than pure algea, means that there is trophic enrichment?)
plot(density(SI.Arc$BUGSoutput$sims.list$epsilon.2H.carbohy)) 
abline(v=-100)

plot(density(SI.Arc$BUGSoutput$sims.list$epsilon.18O.carbohy)) 
abline(v=27)

plot(density(SI.Arc$BUGSoutput$sims.list$lcwax.d2H.inc)) 
abline(v=-129) #much smaller than -129, at -103!, mixing with a -170 per mil source?

plot(density(SI.Arc$BUGSoutput$sims.list$lcwax.d2H.slope))
abline(v=0.78)

plot(density(SI.Arc$BUGSoutput$sims.list$carb.age.gap))
abline(v=0.78)

plot(density(SI.Arc$BUGSoutput$sims.list$scwax.alpha.sl))#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
abline(v=0.0008)#this is very wrong! Use a different prior! Or use something else to constrain salinity

plot(density(SI.Arc$BUGSoutput$sims.list$scwax.alpha.inc))#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Or remove this relationship!

SI.Arc.scalpha <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$alpha2H.sc.alkane,0.89)
plot(1:t*age.res, SI.Arc.scalpha[[1]],type="l",ylim=c(0.6,1), main="sc alpha")
lines(1:t*age.res,SI.Arc.scalpha[[2]],lty=2)
lines(1:t*age.res,SI.Arc.scalpha[[3]],lty=2)

#check LST
SI.Arc.LST <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$LST,0.89)
plot(1:t*age.res,SI.Arc.LST[[1]],type="l",ylim=c(0,30), main="LST")
lines(1:t*age.res,SI.Arc.LST[[2]],lty=2)
lines(1:t*age.res,SI.Arc.LST[[3]],lty=2)
#by far, temperature has the least constraint

#check AT
SI.Arc.AT <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$AT,0.89)
plot(1:t*age.res,test.Arc.AT[[1]],type="l",ylim=c(15,30), main="LST")
lines(1:t*age.res,test.Arc.AT[[2]],lty=2)
lines(1:t*age.res,test.Arc.AT[[3]],lty=2)

#check Runoff and evaporation
SI.Arc.Runoff <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$Runoff,0.89)
plot(1:t*age.res,SI.Arc.Runoff[[1]],type="l",ylim=c(0,12), main="Runoff and Evap")
lines(1:t*age.res,SI.Arc.Runoff[[2]],lty=2)
lines(1:t*age.res,SI.Arc.Runoff[[3]],lty=2)

#check evaporation
SI.Arc.Evap <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$Evap,0.89)
lines(1:t*age.res,SI.Arc.Evap[[1]],col="red")
lines(1:t*age.res,SI.Arc.Evap[[2]],lty=2,col="red")
lines(1:t*age.res,SI.Arc.Evap[[3]],lty=2,col="red")

#check evaporation rate
SI.Arc.Evaprate <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$E.rate,0.89)
plot(1:t*age.res,SI.Arc.Evaprate[[1]],type="l",main="Evap rate", ylim=c(0,2))
lines(1:t*age.res,SI.Arc.Evaprate[[2]],lty=2)
lines(1:t*age.res,SI.Arc.Evaprate[[3]],lty=2)

#check lake level
SI.Arc.L.level <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$L.level,0.89)
plot(1:t*age.res,SI.Arc.L.level[[1]],type="l",ylim=c(1276,1286), main="Lake level")
lines(1:t*age.res,SI.Arc.L.level[[2]],lty=2)
lines(1:t*age.res,SI.Arc.L.level[[3]],lty=2)

#check runoff isotopes
SI.Arc.Rod18O <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$Ro.d18O,0.89)
plot(1:t*age.res,SI.Arc.Rod18O[[1]],type="l",ylim=c(-30,-10), main="Runoff d18O")
lines(1:t*age.res,SI.Arc.Rod18O[[2]],lty=2)
lines(1:t*age.res,SI.Arc.Rod18O[[3]],lty=2)

SI.Arc.Rod2H <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$Ro.d2H,0.89)
plot(1:t*age.res,SI.Arc.Rod2H[[1]],type="l",ylim=c(-200,-140), main="Runoff d2H")
lines(1:t*age.res,SI.Arc.Rod2H[[2]],lty=2)
lines(1:t*age.res,SI.Arc.Rod2H[[3]],lty=2)

#check evaporation isotopes
SI.Arc.evapd18O <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$evap.d18O,0.89)
plot(1:t*age.res,SI.Arc.evapd18O[[1]],type="l",ylim=c(-30,-10), main="Evap d18O")
lines(1:t*age.res,SI.Arc.evapd18O[[2]],lty=2)
lines(1:t*age.res,SI.Arc.evapd18O[[3]],lty=2)

#check lake water isotopes
SI.Arc.Lw.d18O <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$Lw.d18O,0.89)
plot(1:t*age.res,SI.Arc.Lw.d18O[[1]],type="l",ylim=c(-10,0), main="Lake water d18O")
lines(1:t*age.res,SI.Arc.Lw.d18O[[2]],lty=2)
lines(1:t*age.res,SI.Arc.Lw.d18O[[3]],lty=2)

SI.Arc.Lw.d2H <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$Lw.d2H,0.89)
plot(1:t*age.res,SI.Arc.Lw.d2H[[1]],type="l",ylim=c(-100,0), main="Lake water d2H")
lines(1:t*age.res,SI.Arc.Lw.d2H[[2]],lty=2)
lines(1:t*age.res,SI.Arc.Lw.d2H[[3]],lty=2)

#check the gap between runoff d2H and lake water d2H
SI.Arc.Rod2H <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$Ro.d2H,0.89)
SI.Arc.Lw.d2H<- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$Lw.d2H,0.89)
plot(1:t*age.res,SI.Arc.Rod2H[[1]],type="l",ylim=c(-200,0), main="Runoff vs lake d2H")
lines(1:t*age.res,SI.Arc.Lw.d2H[[1]],col="green",lwd =2)

#check lake water evaporation line
plot(SI.Arc.Lw.d18O[[1]],SI.Arc.Lw.d2H[[1]],xlim=c(-10,3),ylim=c(-100,0))
abline(a=-41.189633, b=5.031288) 

#check meteoric line
plot(SI.Arc.Rod18O[[1]],SI.Arc.Rod2H[[1]],xlim=c(-25,-15),ylim=c(-200,-120))
abline(a=0.7314362, b=7.5242923)