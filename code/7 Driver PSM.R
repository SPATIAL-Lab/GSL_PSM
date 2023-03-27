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
#use lake level history to derive precipitation amout
t <- 20
set.seed(1234)
LL.inc <- rnorm(t,0,0.3)
LL.mod <- rep(0,t)
LL.mod[1] <- rnorm(1,1280,0.5)
for(i in 2:t){
  LL.mod[i] <- LL.mod[i - 1] + LL.inc[i - 1]
}

Lw.d18O.inc <- rnorm(t,0,0.3)
Lw.d18O.mod <- rep(0,t)
Lw.d18O.mod[1] <- rnorm(1,-4,0.1)
for(i in 2:t){
  Lw.d18O.mod[i] <- Lw.d18O.mod[i - 1] + Lw.d18O.inc[i - 1]
}

Lw.d2H.mod <- GSL.lake.mwl[1] * Lw.d18O.mod + GSL.lake.mwl[2]

points(Lw.d18O.mod,Lw.d2H.mod, col = "green",pch = 16)


#######use calibrations for short chain wax fractionation
post.leng.scwax <- length(post)

post.scwax.alpha.sl <- post
post.scwax.alpha.inc <- post

#supply slopes and intercepts
MWL.slope <- GSL.lake.mwl[1]
MWL.intc <- GSL.lake.mwl[2]

V.intc <- GSL.Vapor.warm.mwl[1]
V.slope <- GSL.Vapor.warm.mwl[2]

Lw.slope <- GSL.lake.mwl[1]
Lw.intc <- GSL.lake.mwl[2]

#warm season vapor d18O
d18O.vap.warm <- GSL.Vapor.warm.d18O

#carbonate ages
GSL.carb.age

#cyst ages
GSL.carb.age

#parameters to save
parameters <- c("L.level","rh", "nsws", "LST", "T.gap","AT","Runoff", "sal", "S.coeff",
                "E.rate", "Evap", "LV.A", "LST.cps.ac","Lw.d2H","Lw.d18O","f",
                "Ro.d2H", "Ro.d18O", "air.d2H","air.d18O","evap.d2H","evap.d18O")

dat = list(L.level.rec = LL.mod, GSL.level = GSL.level.1286, GSL.area = GSL.area.1286, 
           GSL.volume = GSL.volume.1286, GSL.sali = GSL.sali, t = t, 
           Lw.d2H.dat = Lw.d2H.mod, Lw.d18O.dat = Lw.d18O.mod, 
           V.intc = V.intc, V.slope = V.slope, MWL.intc = MWL.intc, MWL.slope = MWL.slope,
           Lw.intc = Lw.intc, Lw.slope = Lw.slope,
           d18O.vap.warm = d18O.vap.warm)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 2e4
n.burnin = 1e4
n.thin = 1

#Run it
SI.lwmb = do.call(jags.parallel,list(model.file = "code/JAGS PSM Lwmb SI.R", 
                                     parameters.to.save = parameters, 
                                     data = dat, n.chains=5, n.iter = n.iter, 
                                     n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #150 seconds for t = 20, n.iter = 2e4
#estimated time taken: t = 800 (decadal variation), 1.4 hours for n.iter = 2e4

SI.lwmb$BUGSoutput$summary

#check relative humidity
plot(density(SI.lwmb$BUGSoutput$sims.list$rh)) #getting low

#check advection coefficient
plot(density(SI.lwmb$BUGSoutput$sims.list$f)) #

#check wind speed
plot(density(SI.lwmb$BUGSoutput$sims.list$nsws)) #normal

#check LST ac
plot(density(SI.lwmb$BUGSoutput$sims.list$LST.cps.ac)) #all over

#check atm vapor d18O
plot(density(SI.lwmb$BUGSoutput$sims.list$air.d18O)) 

#check LST
SI.lwmb.LST <- MCMC.CI.bound(SI.lwmb$BUGSoutput$sims.list$LST,0.89)
plot(1:t,test.lwmb.LST[[1]],type="l",ylim=c(15,30), main="LST")
lines(1:t,test.lwmb.LST[[2]],lty=2)
lines(1:t,test.lwmb.LST[[3]],lty=2)
#by far, temperature has the least constraint

#check AT
SI.lwmb.AT <- MCMC.CI.bound(SI.lwmb$BUGSoutput$sims.list$AT,0.89)
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