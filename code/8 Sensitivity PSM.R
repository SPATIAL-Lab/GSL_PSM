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

################Sensitivity of PSM with the following elements:#####
#error term using 1/2 pairwise difference 
#no covariance between LST and Ro.d18O
#sc.wax record with n-C18
#with proximal lake mixing
#with d18O and d2H organic matter exchange
#all flexible empirical relationships and epsilons
#evaporation number: 150

#age resolution
age.res = 100 #this determins the number of time steps in the simulation

Bacon.Age.d(1300)#cm, this is consistent with rbacon input

max.age = max(Bacon.Age.d(max.depth)) #use the maximum age to constrain simulation 

t.avg = 28 #number of time steps

#900 years for mean lc alkane offset, 
lcwax.offset = 900

# ~3000 years for carbonate offset
carb.offset = 3000

#number of time steps in the simulation
t <- ceiling(max.age/age.res) + t.avg #101 time steps + 28 for averaging

#######compile priors#######

#supply slopes and intercepts
MWL.slope <- GSL.precip.mwl[1]
MWL.intc <- GSL.precip.mwl[2]

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
#use half pairwise difference as uncertainty around the mean#

d18O.car.dat <- GSL.carb$Ave.d18O

d18O.car.sd <- rep(mean(abs(diff(GSL.carb$Ave.d18O))/2), n.carb)

lcwax.d2H.dat <- GSL.FAME.C28.naom$C28.avg

lcwax.d2H.sd <- rep(mean(abs(diff(GSL.FAME.C28.naom$C28.avg))/2), n.lc.wax)

scwax.d2H.dat <- GSL.FAME.C18.naom$C18.avg

scwax.d2H.sd <- rep(mean(abs(diff(GSL.FAME.C18.naom$C18.avg))/2), n.sc.wax) 

BScyst.d18O.dat <- GSL.cyst$Ave.d18O

BScyst.d18O.sd <- rep(mean(abs(diff(GSL.cyst$Ave.d18O))/2), n.cyst)

BScyst.d2H.dat <- GSL.cyst$Ave.d2H

BScyst.d2H.sd <- rep(mean(abs(diff(GSL.cyst$Ave.d2H))/2), n.cyst)

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
                "carb.age.gap", "Ro.d18O.cps","age.carb.d", "age.lc.d","age.cyst.d","age.sc.d",
                "ind.dep", "carb.bins", "lcwax.bins","scwax.bins", "cyst.bins",
                "prx.sal", "prx.L.d2H", "prx.L.d18O", "f.m.ro","w.avg.carb", "w.avg.lcwax",
                "d2H.C14d.lcwax","d18O.C14d.carb","f.C14d.carb", "f.C14d.lcwax",
                "T.cps.slope", "sl.cpsAT.18O", "r.exO", "r.exH", "alpha.exH", "alpha.exO",
                "epsilon.alk.acid","LST.pre","tl","ind.dep")


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
           BScyst.d2H.sd = BScyst.d2H.sd, t.avg=t.avg, carb.offset= carb.offset, lcwax.offset=lcwax.offset,
           n.d.evap = 150)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 5e3
n.thin = 5

#Run it
PSM.s.evap150 = do.call(jags.parallel,list(model.file = "code/JAGS PSM full.R", 
                                   parameters.to.save = parameters, 
                                   data = dat, n.chains=5, n.iter = n.iter, 
                                   n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #8.5 hours for 1e4
#estimated time taken: t = 129 (centenial variation), 2.3 hours for n.iter = 2e4

write.csv(PSM.s.evap150$BUGSoutput$summary, "out/PSM_s_evap150_sum.csv")

PSM.f$BUGSoutput$summary

par(mfrow=c(3,3))
#check relative humidity
plot(density(PSM.s.evap150$BUGSoutput$sims.list$rh),main="RH") #getting low

#check advection coefficient
plot(density(PSM.s.evap150$BUGSoutput$sims.list$f),main="F convection") #

#check wind speed
plot(density(PSM.s.evap150$BUGSoutput$sims.list$nsws),main="Wind") #normal

#check LST ac
plot(density(PSM.s.evap150$BUGSoutput$sims.list$LST.pre),main="LST pre") #very high!

#check Ro ac
plot(density(PSM.s.evap150$BUGSoutput$sims.list$Ro.d18O.cps),main="Ro ac")

#check atm vapor d18O
plot(density(PSM.s.evap150$BUGSoutput$sims.list$air.d18O),main="Atm d18O") 
abline(v=d18O.vap.warm)

# check brine shrimp cyst intercept
plot(density(PSM.s.evap150$BUGSoutput$sims.list$BScyst.inc.2H),main="Cyst d2H interc")
abline(v=-92) #ok

#check brine shrimp cyst slopes
plot(density(PSM.s.evap150$BUGSoutput$sims.list$BScyst.slope.lw.2H),main="Cyst lw d2H slope")
abline(v=0.34) #ok

#check brine shrimp cyst intercept
plot(density(PSM.s.evap150$BUGSoutput$sims.list$BScyst.slope.diet.2H),main="Cyst diet d2H slope")
abline(v=0.26) #0.2

#check brine shrimp cyst intercept
plot(density(PSM.s.evap150$BUGSoutput$sims.list$BScyst.inc.18O),main="Cyst d18O interc")
abline(v=15.9)#ok
#check brine shrimp cyst slopes
plot(density(PSM.s.evap150$BUGSoutput$sims.list$BScyst.slope.lw.18O),main="Cyst lw d18O slope")
abline(v=0.692)#ok
#check brine shrimp cyst intercept
plot(density(PSM.s.evap150$BUGSoutput$sims.list$BScyst.slope.diet.18O),main="Cyst diet d18O slope")
abline(v=0.101)#ok

plot(density(PSM.s.evap150$BUGSoutput$sims.list$d2H.gap.MAP_Ro),main="MAP Ro gap") 
abline(v=45)

#epsilon BS diet (smaller than pure algea, means that there is trophic enrichment?)
plot(density(PSM.s.evap150$BUGSoutput$sims.list$epsilon.2H.carbohy),main="2H carbohy ep") 
abline(v=-100)#likely a mixture of both algae and terrestrial/wetland POM

plot(density(PSM.s.evap150$BUGSoutput$sims.list$epsilon.18O.carbohy),main="18O carbohy ep") 
abline(v=27) #much lower for d18O!

plot(density(PSM.s.evap150$BUGSoutput$sims.list$lcwax.d2H.inc),main="lc MAP inc") 
abline(v=-125) #much smaller than -129, at -115!, mixing with a -170 per mil source?

plot(density(PSM.s.evap150$BUGSoutput$sims.list$lcwax.d2H.slope),main="lc MAP slope")
abline(v=0.62) #0.76

plot(density(PSM.s.evap150$BUGSoutput$sims.list$scwax.alpha.sl),main="sc alpha slope")
abline(v=0.0008)#0.0009

plot(density(PSM.s.evap150$BUGSoutput$sims.list$scwax.alpha.inc),main="sc alpha interc")
abline(v=0.80745) # Consider make this fixed
#Or remove this relationship!

plot(density(PSM.s.evap150$BUGSoutput$sims.list$f.m.ro),main="f.m.ro") #more like 0.4
abline(v=0.2) 

plot(density(PSM.s.evap150$BUGSoutput$sims.list$f.C14d.carb),main="f.C14d.carb") #ok
abline(v=0.05) 

plot(density(PSM.s.evap150$BUGSoutput$sims.list$f.C14d.lcwax),main="f.C14d.lcwax") #ok
abline(v=0.01) 

plot(density(PSM.s.evap150$BUGSoutput$sims.list$epsilon.alk.acid),main="epsilon.alk.acid") #ok
abline(v=25)

plot(density(PSM.s.evap150$BUGSoutput$sims.list$T.cps.slope), main="T.cps.slope")

plot(density(PSM.s.evap150$BUGSoutput$sims.list$tl), main="Date aligned")

PSM.s.evap150.scalpha <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$alpha2H.sc.alkane,0.89)
plot(t:1*age.res, PSM.s.evap150.scalpha[[1]],type="l",ylim=c(0.8,1), main="sc alpha")
lines(t:1*age.res,PSM.s.evap150.scalpha[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap150.scalpha[[3]],lty=2)

PSM.s.evap150.prx.sal <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$prx.sal,0.89)
plot(t:1*age.res, PSM.s.evap150.prx.sal[[1]],type="l",ylim=c(0,200), main="pre sal")
lines(t:1*age.res,PSM.s.evap150.prx.sal[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap150.prx.sal[[3]],lty=2)

PSM.s.evap150.prx.L.d2H <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$prx.L.d2H,0.89)
plot(t:1*age.res, PSM.s.evap150.prx.L.d2H[[1]],type="l",ylim=c(-150,0), main="prx L d2H")
lines(t:1*age.res,PSM.s.evap150.prx.L.d2H[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap150.prx.L.d2H[[3]],lty=2)

#check LST
PSM.s.evap150.LST <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$LST,0.89)
plot(t:1*age.res,PSM.s.evap150.LST[[1]],type="l",ylim=c(0,25), main="LST")
lines(t:1*age.res,PSM.s.evap150.LST[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap150.LST[[3]],lty=2)
#by far, temperature has the least constraint, still too high!

#check AT
PSM.s.evap150.AT <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$AT,0.89)
plot(t:1*age.res,PSM.s.evap150.AT[[1]],type="l",ylim=c(0,25), main="AT")
lines(t:1*age.res,PSM.s.evap150.AT[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap150.AT[[3]],lty=2)

#check Runoff and evaporation
PSM.s.evap150.Runoff <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$Runoff,0.89)
plot(t:1*age.res,PSM.s.evap150.Runoff[[1]],type="l",ylim=c(0,12), main="Runoff and Evap")
lines(t:1*age.res,PSM.s.evap150.Runoff[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap150.Runoff[[3]],lty=2)

#check evaporation
PSM.s.evap150.Evap <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$Evap,0.89)
lines(t:1*age.res,PSM.s.evap150.Evap[[1]],col="red")
lines(t:1*age.res,PSM.s.evap150.Evap[[2]],lty=2,col="red")
lines(t:1*age.res,PSM.s.evap150.Evap[[3]],lty=2,col="red")

#check evaporation rate
PSM.s.evap150.Evaprate <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$E.rate,0.89)
plot(t:1*age.res,PSM.s.evap150.Evaprate[[1]],type="l",main="Evap rate", ylim=c(0,2))
lines(t:1*age.res,PSM.s.evap150.Evaprate[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap150.Evaprate[[3]],lty=2)

#check lake level
PSM.s.evap150.L.level <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$L.level,0.89)
plot(t:1*age.res,PSM.s.evap150.L.level[[1]],type="l",ylim=c(1276,1286), main="Lake level")
lines(t:1*age.res,PSM.s.evap150.L.level[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap150.L.level[[3]],lty=2)

#check salinity
PSM.s.evap150.sal <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$sal,0.89)
plot(t:1*age.res,PSM.s.evap150.sal[[1]],type="l",ylim=c(50,360), main="salinity")
lines(t:1*age.res,PSM.s.evap150.sal[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap150.sal[[3]],lty=2)

#check runoff isotopes
PSM.s.evap150.Rod18O <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$Ro.d18O,0.89)
plot(t:1*age.res,PSM.s.evap150.Rod18O[[1]],type="l",ylim=c(-25,-10), main="Runoff d18O")
lines(t:1*age.res,PSM.s.evap150.Rod18O[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap150.Rod18O[[3]],lty=2)

PSM.s.evap150.Rod2H <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$Ro.d2H,0.89)
plot(t:1*age.res,PSM.s.evap150.Rod2H[[1]],type="l",ylim=c(-180,-80), main="Runoff d2H")
lines(t:1*age.res,PSM.s.evap150.Rod2H[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap150.Rod2H[[3]],lty=2)

#check evaporation isotopes
PSM.s.evap150.evap.d18O <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$evap.d18O,0.89)
plot(t:1*age.res,PSM.s.evap150.evap.d18O[[1]],type="l",ylim=c(-25,-10), main="Evap d18O")
lines(t:1*age.res,PSM.s.evap150.evap.d18O[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap150.evap.d18O[[3]],lty=2)

PSM.s.evap150.evap.d2H <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$evap.d2H,0.89)
plot(t:1*age.res,PSM.s.evap150.evap.d2H[[1]],type="l",ylim=c(-180,-80), main="Evap d2H")
lines(t:1*age.res,PSM.s.evap150.evap.d2H[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap150.evap.d2H[[3]],lty=2)

#check lake water isotopes
PSM.s.evap150.Lw.d18O <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$Lw.d18O,0.89)
plot(t:1*age.res,PSM.s.evap150.Lw.d18O[[1]],type="l",ylim=c(-10,0), main="Lake water d18O")
lines(t:1*age.res,PSM.s.evap150.Lw.d18O[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap150.Lw.d18O[[3]],lty=2)

PSM.s.evap150.Lw.d2H <- MCMC.CI.bound(PSM.s.evap150$BUGSoutput$sims.list$Lw.d2H,0.89)
plot(t:1*age.res,PSM.s.evap150.Lw.d2H[[1]],type="l",ylim=c(-100,0), main="Lake water d2H")
lines(t:1*age.res,PSM.s.evap150.Lw.d2H[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap150.Lw.d2H[[3]],lty=2)

plot(PSM.s.evap150.Lw.d18O[[1]],PSM.s.evap150.Lw.d2H[[1]],xlim=c(-25,5),ylim=c(-200,20))
abline(a=Lw.intc, b=Lw.slope) #ok
#check meteoric line
points(PSM.s.evap150.Rod18O[[1]],PSM.s.evap150.Rod2H[[1]])
abline(a=MWL.intc, b=MWL.slope)
points(PSM.s.evap150.evap.d18O[[1]],PSM.s.evap150.evap.d2H[[1]],col="red")

################Sensitivity of PSM with the following elements:#####
#error term using 1/2 pairwise difference 
#no covariance between LST and Ro.d18O
#sc.wax record with n-C18
#with proximal lake mixing
#with d18O and d2H organic matter exchange
#all flexible empirical relationships and epsilons
#evaporation number: 150

#age resolution
age.res = 100 #this determins the number of time steps in the simulation

Bacon.Age.d(1300)#cm, this is consistent with rbacon input

max.age = max(Bacon.Age.d(max.depth)) #use the maximum age to constrain simulation 

t.avg = 28 #number of time steps

#900 years for mean lc alkane offset, 
lcwax.offset = 900

# ~3000 years for carbonate offset
carb.offset = 3000

#number of time steps in the simulation
t <- ceiling(max.age/age.res) + t.avg #101 time steps + 28 for averaging

#######compile priors#######

#supply slopes and intercepts
MWL.slope <- GSL.precip.mwl[1]
MWL.intc <- GSL.precip.mwl[2]

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
#use half pairwise difference as uncertainty around the mean#

d18O.car.dat <- GSL.carb$Ave.d18O

d18O.car.sd <- rep(mean(abs(diff(GSL.carb$Ave.d18O))/2), n.carb)

lcwax.d2H.dat <- GSL.FAME.C28.naom$C28.avg

lcwax.d2H.sd <- rep(mean(abs(diff(GSL.FAME.C28.naom$C28.avg))/2), n.lc.wax)

scwax.d2H.dat <- GSL.FAME.C18.naom$C18.avg

scwax.d2H.sd <- rep(mean(abs(diff(GSL.FAME.C18.naom$C18.avg))/2), n.sc.wax) 

BScyst.d18O.dat <- GSL.cyst$Ave.d18O

BScyst.d18O.sd <- rep(mean(abs(diff(GSL.cyst$Ave.d18O))/2), n.cyst)

BScyst.d2H.dat <- GSL.cyst$Ave.d2H

BScyst.d2H.sd <- rep(mean(abs(diff(GSL.cyst$Ave.d2H))/2), n.cyst)

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
                "carb.age.gap", "Ro.d18O.cps","age.carb.d", "age.lc.d","age.cyst.d","age.sc.d",
                "ind.dep", "carb.bins", "lcwax.bins","scwax.bins", "cyst.bins",
                "prx.sal", "prx.L.d2H", "prx.L.d18O", "f.m.ro","w.avg.carb", "w.avg.lcwax",
                "d2H.C14d.lcwax","d18O.C14d.carb","f.C14d.carb", "f.C14d.lcwax",
                "T.cps.slope", "sl.cpsAT.18O", "r.exO", "r.exH", "alpha.exH", "alpha.exO",
                "epsilon.alk.acid","LST.pre","tl","ind.dep")


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
           BScyst.d2H.sd = BScyst.d2H.sd, t.avg=t.avg, carb.offset= carb.offset, lcwax.offset=lcwax.offset,
           n.d.evap = 120)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 5e3
n.thin = 5

#Run it
PSM.s.evap120 = do.call(jags.parallel,list(model.file = "code/JAGS PSM full.R", 
                                           parameters.to.save = parameters, 
                                           data = dat, n.chains=5, n.iter = n.iter, 
                                           n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #8.5 hours for 1e4
#estimated time taken: t = 129 (centenial variation), 2.3 hours for n.iter = 2e4

write.csv(PSM.s.evap120$BUGSoutput$summary, "out/PSM_s_evap120_sum.csv")

PSM.f$BUGSoutput$summary

par(mfrow=c(3,3))
#check relative humidity
plot(density(PSM.s.evap120$BUGSoutput$sims.list$rh),main="RH") #getting low

#check advection coefficient
plot(density(PSM.s.evap120$BUGSoutput$sims.list$f),main="F convection") #

#check wind speed
plot(density(PSM.s.evap120$BUGSoutput$sims.list$nsws),main="Wind") #normal

#check LST ac
plot(density(PSM.s.evap120$BUGSoutput$sims.list$LST.pre),main="LST pre") #very high!

#check Ro ac
plot(density(PSM.s.evap120$BUGSoutput$sims.list$Ro.d18O.cps),main="Ro ac")

#check atm vapor d18O
plot(density(PSM.s.evap120$BUGSoutput$sims.list$air.d18O),main="Atm d18O") 
abline(v=d18O.vap.warm)

# check brine shrimp cyst intercept
plot(density(PSM.s.evap120$BUGSoutput$sims.list$BScyst.inc.2H),main="Cyst d2H interc")
abline(v=-92) #ok

#check brine shrimp cyst slopes
plot(density(PSM.s.evap120$BUGSoutput$sims.list$BScyst.slope.lw.2H),main="Cyst lw d2H slope")
abline(v=0.34) #ok

#check brine shrimp cyst intercept
plot(density(PSM.s.evap120$BUGSoutput$sims.list$BScyst.slope.diet.2H),main="Cyst diet d2H slope")
abline(v=0.26) #0.2

#check brine shrimp cyst intercept
plot(density(PSM.s.evap120$BUGSoutput$sims.list$BScyst.inc.18O),main="Cyst d18O interc")
abline(v=15.9)#ok
#check brine shrimp cyst slopes
plot(density(PSM.s.evap120$BUGSoutput$sims.list$BScyst.slope.lw.18O),main="Cyst lw d18O slope")
abline(v=0.692)#ok
#check brine shrimp cyst intercept
plot(density(PSM.s.evap120$BUGSoutput$sims.list$BScyst.slope.diet.18O),main="Cyst diet d18O slope")
abline(v=0.101)#ok

plot(density(PSM.s.evap120$BUGSoutput$sims.list$d2H.gap.MAP_Ro),main="MAP Ro gap") 
abline(v=45)

#epsilon BS diet (smaller than pure algea, means that there is trophic enrichment?)
plot(density(PSM.s.evap120$BUGSoutput$sims.list$epsilon.2H.carbohy),main="2H carbohy ep") 
abline(v=-100)#likely a mixture of both algae and terrestrial/wetland POM

plot(density(PSM.s.evap120$BUGSoutput$sims.list$epsilon.18O.carbohy),main="18O carbohy ep") 
abline(v=27) #much lower for d18O!

plot(density(PSM.s.evap120$BUGSoutput$sims.list$lcwax.d2H.inc),main="lc MAP inc") 
abline(v=-125) #much smaller than -129, at -115!, mixing with a -170 per mil source?

plot(density(PSM.s.evap120$BUGSoutput$sims.list$lcwax.d2H.slope),main="lc MAP slope")
abline(v=0.62) #0.76

plot(density(PSM.s.evap120$BUGSoutput$sims.list$scwax.alpha.sl),main="sc alpha slope")
abline(v=0.0008)#0.0009

plot(density(PSM.s.evap120$BUGSoutput$sims.list$scwax.alpha.inc),main="sc alpha interc")
abline(v=0.80745) # Consider make this fixed
#Or remove this relationship!

plot(density(PSM.s.evap120$BUGSoutput$sims.list$f.m.ro),main="f.m.ro") #more like 0.4
abline(v=0.2) 

plot(density(PSM.s.evap120$BUGSoutput$sims.list$f.C14d.carb),main="f.C14d.carb") #ok
abline(v=0.05) 

plot(density(PSM.s.evap120$BUGSoutput$sims.list$f.C14d.lcwax),main="f.C14d.lcwax") #ok
abline(v=0.01) 

plot(density(PSM.s.evap120$BUGSoutput$sims.list$epsilon.alk.acid),main="epsilon.alk.acid") #ok
abline(v=25)

plot(density(PSM.s.evap120$BUGSoutput$sims.list$T.cps.slope), main="T.cps.slope")

plot(density(PSM.s.evap120$BUGSoutput$sims.list$tl), main="Date aligned")

PSM.s.evap120.scalpha <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$alpha2H.sc.alkane,0.89)
plot(t:1*age.res, PSM.s.evap120.scalpha[[1]],type="l",ylim=c(0.8,1), main="sc alpha")
lines(t:1*age.res,PSM.s.evap120.scalpha[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.scalpha[[3]],lty=2)

PSM.s.evap120.prx.sal <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$prx.sal,0.89)
plot(t:1*age.res, PSM.s.evap120.prx.sal[[1]],type="l",ylim=c(0,200), main="pre sal")
lines(t:1*age.res,PSM.s.evap120.prx.sal[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.prx.sal[[3]],lty=2)

PSM.s.evap120.prx.L.d2H <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$prx.L.d2H,0.89)
plot(t:1*age.res, PSM.s.evap120.prx.L.d2H[[1]],type="l",ylim=c(-150,0), main="prx L d2H")
lines(t:1*age.res,PSM.s.evap120.prx.L.d2H[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.prx.L.d2H[[3]],lty=2)

#check LST
PSM.s.evap120.LST <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$LST,0.89)
plot(t:1*age.res,PSM.s.evap120.LST[[1]],type="l",ylim=c(0,25), main="LST")
lines(t:1*age.res,PSM.s.evap120.LST[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.LST[[3]],lty=2)
#by far, temperature has the least constraint, still too high!

#check AT
PSM.s.evap120.AT <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$AT,0.89)
plot(t:1*age.res,PSM.s.evap120.AT[[1]],type="l",ylim=c(0,25), main="AT")
lines(t:1*age.res,PSM.s.evap120.AT[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.AT[[3]],lty=2)

#check Runoff and evaporation
PSM.s.evap120.Runoff <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$Runoff,0.89)
plot(t:1*age.res,PSM.s.evap120.Runoff[[1]],type="l",ylim=c(0,12), main="Runoff and Evap")
lines(t:1*age.res,PSM.s.evap120.Runoff[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.Runoff[[3]],lty=2)

#check evaporation
PSM.s.evap120.Evap <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$Evap,0.89)
lines(t:1*age.res,PSM.s.evap120.Evap[[1]],col="red")
lines(t:1*age.res,PSM.s.evap120.Evap[[2]],lty=2,col="red")
lines(t:1*age.res,PSM.s.evap120.Evap[[3]],lty=2,col="red")

#check evaporation rate
PSM.s.evap120.Evaprate <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$E.rate,0.89)
plot(t:1*age.res,PSM.s.evap120.Evaprate[[1]],type="l",main="Evap rate", ylim=c(0,2))
lines(t:1*age.res,PSM.s.evap120.Evaprate[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.Evaprate[[3]],lty=2)

#check lake level
PSM.s.evap120.L.level <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$L.level,0.89)
plot(t:1*age.res,PSM.s.evap120.L.level[[1]],type="l",ylim=c(1276,1286), main="Lake level")
lines(t:1*age.res,PSM.s.evap120.L.level[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.L.level[[3]],lty=2)

#check salinity
PSM.s.evap120.sal <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$sal,0.89)
plot(t:1*age.res,PSM.s.evap120.sal[[1]],type="l",ylim=c(50,360), main="salinity")
lines(t:1*age.res,PSM.s.evap120.sal[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.sal[[3]],lty=2)

#check runoff isotopes
PSM.s.evap120.Rod18O <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$Ro.d18O,0.89)
plot(t:1*age.res,PSM.s.evap120.Rod18O[[1]],type="l",ylim=c(-25,-10), main="Runoff d18O")
lines(t:1*age.res,PSM.s.evap120.Rod18O[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.Rod18O[[3]],lty=2)

PSM.s.evap120.Rod2H <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$Ro.d2H,0.89)
plot(t:1*age.res,PSM.s.evap120.Rod2H[[1]],type="l",ylim=c(-180,-80), main="Runoff d2H")
lines(t:1*age.res,PSM.s.evap120.Rod2H[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.Rod2H[[3]],lty=2)

#check evaporation isotopes
PSM.s.evap120.evap.d18O <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$evap.d18O,0.89)
plot(t:1*age.res,PSM.s.evap120.evap.d18O[[1]],type="l",ylim=c(-25,-10), main="Evap d18O")
lines(t:1*age.res,PSM.s.evap120.evap.d18O[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.evap.d18O[[3]],lty=2)

PSM.s.evap120.evap.d2H <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$evap.d2H,0.89)
plot(t:1*age.res,PSM.s.evap120.evap.d2H[[1]],type="l",ylim=c(-180,-80), main="Evap d2H")
lines(t:1*age.res,PSM.s.evap120.evap.d2H[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.evap.d2H[[3]],lty=2)

#check lake water isotopes
PSM.s.evap120.Lw.d18O <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$Lw.d18O,0.89)
plot(t:1*age.res,PSM.s.evap120.Lw.d18O[[1]],type="l",ylim=c(-10,0), main="Lake water d18O")
lines(t:1*age.res,PSM.s.evap120.Lw.d18O[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.Lw.d18O[[3]],lty=2)

PSM.s.evap120.Lw.d2H <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$Lw.d2H,0.89)
plot(t:1*age.res,PSM.s.evap120.Lw.d2H[[1]],type="l",ylim=c(-100,0), main="Lake water d2H")
lines(t:1*age.res,PSM.s.evap120.Lw.d2H[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.Lw.d2H[[3]],lty=2)

plot(PSM.s.evap120.Lw.d18O[[1]],PSM.s.evap120.Lw.d2H[[1]],xlim=c(-25,5),ylim=c(-200,20))
abline(a=Lw.intc, b=Lw.slope) #ok
#check meteoric line
points(PSM.s.evap120.Rod18O[[1]],PSM.s.evap120.Rod2H[[1]])
abline(a=MWL.intc, b=MWL.slope)
points(PSM.s.evap120.evap.d18O[[1]],PSM.s.evap120.evap.d2H[[1]],col="red")

################Sensitivity of PSM with the following elements:#####
#error term using 1/2 pairwise difference 
#no covariance between LST and Ro.d18O
#sc.wax record with n-C18
#with proximal lake mixing
# * without * d18O and d2H organic matter exchange
#all flexible empirical relationships and epsilons
#evaporation number: 150

#age resolution
age.res = 100 #this determins the number of time steps in the simulation

Bacon.Age.d(1300)#cm, this is consistent with rbacon input

max.age = max(Bacon.Age.d(max.depth)) #use the maximum age to constrain simulation 

t.avg = 28 #number of time steps

#900 years for mean lc alkane offset, 
lcwax.offset = 900

# ~3000 years for carbonate offset
carb.offset = 3000

#number of time steps in the simulation
t <- ceiling(max.age/age.res) + t.avg #101 time steps + 28 for averaging

#######compile priors#######

#supply slopes and intercepts
MWL.slope <- GSL.precip.mwl[1]
MWL.intc <- GSL.precip.mwl[2]

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
#use half pairwise difference as uncertainty around the mean#

d18O.car.dat <- GSL.carb$Ave.d18O

d18O.car.sd <- rep(mean(abs(diff(GSL.carb$Ave.d18O))/2), n.carb)

lcwax.d2H.dat <- GSL.FAME.C28.naom$C28.avg

lcwax.d2H.sd <- rep(mean(abs(diff(GSL.FAME.C28.naom$C28.avg))/2), n.lc.wax)

scwax.d2H.dat <- GSL.FAME.C18.naom$C18.avg

scwax.d2H.sd <- rep(mean(abs(diff(GSL.FAME.C18.naom$C18.avg))/2), n.sc.wax) 

BScyst.d18O.dat <- GSL.cyst$Ave.d18O

BScyst.d18O.sd <- rep(mean(abs(diff(GSL.cyst$Ave.d18O))/2), n.cyst)

BScyst.d2H.dat <- GSL.cyst$Ave.d2H

BScyst.d2H.sd <- rep(mean(abs(diff(GSL.cyst$Ave.d2H))/2), n.cyst)

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
                "carb.age.gap", "Ro.d18O.cps","age.carb.d", "age.lc.d","age.cyst.d","age.sc.d",
                "ind.dep", "carb.bins", "lcwax.bins","scwax.bins", "cyst.bins",
                "prx.sal", "prx.L.d2H", "prx.L.d18O", "f.m.ro","w.avg.carb", "w.avg.lcwax",
                "d2H.C14d.lcwax","d18O.C14d.carb","f.C14d.carb", "f.C14d.lcwax",
                "T.cps.slope", "sl.cpsAT.18O", "r.exO", "r.exH", "alpha.exH", "alpha.exO",
                "epsilon.alk.acid","LST.pre","tl","ind.dep")


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
           BScyst.d2H.sd = BScyst.d2H.sd, t.avg=t.avg, carb.offset= carb.offset, lcwax.offset=lcwax.offset)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 5e3
n.thin = 5

#Run it
PSM.s.wo.ex = do.call(jags.parallel,list(model.file = "code/JAGS PSM full wo om ex.R", 
                                           parameters.to.save = parameters, 
                                           data = dat, n.chains=5, n.iter = n.iter, 
                                           n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #8.5 hours for 1e4
#estimated time taken: t = 129 (centenial variation), 2.3 hours for n.iter = 2e4

write.csv(PSM.s.wo.ex$BUGSoutput$summary, "out/PSM.s.wo.ex_sum.csv")

PSM.s.wo.ex$BUGSoutput$summary

par(mfrow=c(3,3))
#check relative humidity
plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$rh),main="RH") #getting low

#check advection coefficient
plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$f),main="F convection") #

#check wind speed
plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$nsws),main="Wind") #normal

#check LST ac
plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$LST.pre),main="LST pre") #very high!

#check Ro ac
plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$Ro.d18O.cps),main="Ro ac")

#check atm vapor d18O
plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$air.d18O),main="Atm d18O") 
abline(v=d18O.vap.warm)

# check brine shrimp cyst intercept
plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$BScyst.inc.2H),main="Cyst d2H interc")
abline(v=-92) #ok

#check brine shrimp cyst slopes
plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$BScyst.slope.lw.2H),main="Cyst lw d2H slope")
abline(v=0.34) #ok

#check brine shrimp cyst intercept
plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$BScyst.slope.diet.2H),main="Cyst diet d2H slope")
abline(v=0.26) #0.2

#check brine shrimp cyst intercept
plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$BScyst.inc.18O),main="Cyst d18O interc")
abline(v=15.9)#ok
#check brine shrimp cyst slopes
plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$BScyst.slope.lw.18O),main="Cyst lw d18O slope")
abline(v=0.692)#ok
#check brine shrimp cyst intercept
plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$BScyst.slope.diet.18O),main="Cyst diet d18O slope")
abline(v=0.101)#ok

plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$d2H.gap.MAP_Ro),main="MAP Ro gap") 
abline(v=45)

#epsilon BS diet (smaller than pure algea, means that there is trophic enrichment?)
plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$epsilon.2H.carbohy),main="2H carbohy ep") 
abline(v=-100)#likely a mixture of both algae and terrestrial/wetland POM

plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$epsilon.18O.carbohy),main="18O carbohy ep") 
abline(v=27) #much lower for d18O!

plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$lcwax.d2H.inc),main="lc MAP inc") 
abline(v=-125) #much smaller than -129, at -115!, mixing with a -170 per mil source?

plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$lcwax.d2H.slope),main="lc MAP slope")
abline(v=0.62) #0.76

plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$scwax.alpha.sl),main="sc alpha slope")
abline(v=0.0008)#0.0009

plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$scwax.alpha.inc),main="sc alpha interc")
abline(v=0.80745) # Consider make this fixed
#Or remove this relationship!

plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$f.m.ro),main="f.m.ro") #more like 0.4
abline(v=0.2) 

plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$f.C14d.carb),main="f.C14d.carb") #ok
abline(v=0.05) 

plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$f.C14d.lcwax),main="f.C14d.lcwax") #ok
abline(v=0.01) 

plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$epsilon.alk.acid),main="epsilon.alk.acid") #ok
abline(v=25)

plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$T.cps.slope), main="T.cps.slope")

plot(density(PSM.s.wo.ex$BUGSoutput$sims.list$tl), main="Date aligned")

PSM.s.wo.ex.scalpha <- MCMC.CI.bound(PSM.s.wo.ex$BUGSoutput$sims.list$alpha2H.sc.alkane,0.89)
plot(t:1*age.res, PSM.s.wo.ex.scalpha[[1]],type="l",ylim=c(0.8,1), main="sc alpha")
lines(t:1*age.res,PSM.s.wo.ex.scalpha[[2]],lty=2)
lines(t:1*age.res,PSM.s.wo.ex.scalpha[[3]],lty=2)

PSM.s.wo.ex.prx.sal <- MCMC.CI.bound(PSM.s.wo.ex$BUGSoutput$sims.list$prx.sal,0.89)
plot(t:1*age.res, PSM.s.wo.ex.prx.sal[[1]],type="l",ylim=c(0,200), main="pre sal")
lines(t:1*age.res,PSM.s.wo.ex.prx.sal[[2]],lty=2)
lines(t:1*age.res,PSM.s.wo.ex.prx.sal[[3]],lty=2)

PSM.s.wo.ex.prx.L.d2H <- MCMC.CI.bound(PSM.s.wo.ex$BUGSoutput$sims.list$prx.L.d2H,0.89)
plot(t:1*age.res, PSM.s.wo.ex.prx.L.d2H[[1]],type="l",ylim=c(-150,0), main="prx L d2H")
lines(t:1*age.res,PSM.s.wo.ex.prx.L.d2H[[2]],lty=2)
lines(t:1*age.res,PSM.s.wo.ex.prx.L.d2H[[3]],lty=2)

#check LST
PSM.s.wo.ex.LST <- MCMC.CI.bound(PSM.s.wo.ex$BUGSoutput$sims.list$LST,0.89)
plot(t:1*age.res,PSM.s.wo.ex.LST[[1]],type="l",ylim=c(0,25), main="LST")
lines(t:1*age.res,PSM.s.wo.ex.LST[[2]],lty=2)
lines(t:1*age.res,PSM.s.wo.ex.LST[[3]],lty=2)
#by far, temperature has the least constraint, still too high!

#check AT
PSM.s.wo.ex.AT <- MCMC.CI.bound(PSM.s.wo.ex$BUGSoutput$sims.list$AT,0.89)
plot(t:1*age.res,PSM.s.wo.ex.AT[[1]],type="l",ylim=c(0,25), main="AT")
lines(t:1*age.res,PSM.s.wo.ex.AT[[2]],lty=2)
lines(t:1*age.res,PSM.s.wo.ex.AT[[3]],lty=2)

#check Runoff and evaporation
PSM.s.wo.ex.Runoff <- MCMC.CI.bound(PSM.s.wo.ex$BUGSoutput$sims.list$Runoff,0.89)
plot(t:1*age.res,PSM.s.wo.ex.Runoff[[1]],type="l",ylim=c(0,12), main="Runoff and Evap")
lines(t:1*age.res,PSM.s.wo.ex.Runoff[[2]],lty=2)
lines(t:1*age.res,PSM.s.wo.ex.Runoff[[3]],lty=2)

#check evaporation
PSM.s.wo.ex.Evap <- MCMC.CI.bound(PSM.s.wo.ex$BUGSoutput$sims.list$Evap,0.89)
lines(t:1*age.res,PSM.s.wo.ex.Evap[[1]],col="red")
lines(t:1*age.res,PSM.s.wo.ex.Evap[[2]],lty=2,col="red")
lines(t:1*age.res,PSM.s.wo.ex.Evap[[3]],lty=2,col="red")

#check evaporation rate
PSM.s.wo.ex.Evaprate <- MCMC.CI.bound(PSM.s.wo.ex$BUGSoutput$sims.list$E.rate,0.89)
plot(t:1*age.res,PSM.s.wo.ex.Evaprate[[1]],type="l",main="Evap rate", ylim=c(0,2))
lines(t:1*age.res,PSM.s.wo.ex.Evaprate[[2]],lty=2)
lines(t:1*age.res,PSM.s.wo.ex.Evaprate[[3]],lty=2)

#check lake level
PSM.s.wo.ex.L.level <- MCMC.CI.bound(PSM.s.wo.ex$BUGSoutput$sims.list$L.level,0.89)
plot(t:1*age.res,PSM.s.wo.ex.L.level[[1]],type="l",ylim=c(1276,1286), main="Lake level")
lines(t:1*age.res,PSM.s.wo.ex.L.level[[2]],lty=2)
lines(t:1*age.res,PSM.s.wo.ex.L.level[[3]],lty=2)

#check salinity
PSM.s.wo.ex.sal <- MCMC.CI.bound(PSM.s.wo.ex$BUGSoutput$sims.list$sal,0.89)
plot(t:1*age.res,PSM.s.wo.ex.sal[[1]],type="l",ylim=c(50,360), main="salinity")
lines(t:1*age.res,PSM.s.wo.ex.sal[[2]],lty=2)
lines(t:1*age.res,PSM.s.wo.ex.sal[[3]],lty=2)

#check runoff isotopes
PSM.s.wo.ex.Rod18O <- MCMC.CI.bound(PSM.s.wo.ex$BUGSoutput$sims.list$Ro.d18O,0.89)
plot(t:1*age.res,PSM.s.wo.ex.Rod18O[[1]],type="l",ylim=c(-25,-10), main="Runoff d18O")
lines(t:1*age.res,PSM.s.wo.ex.Rod18O[[2]],lty=2)
lines(t:1*age.res,PSM.s.wo.ex.Rod18O[[3]],lty=2)

PSM.s.wo.ex.Rod2H <- MCMC.CI.bound(PSM.s.wo.ex$BUGSoutput$sims.list$Ro.d2H,0.89)
plot(t:1*age.res,PSM.s.wo.ex.Rod2H[[1]],type="l",ylim=c(-180,-80), main="Runoff d2H")
lines(t:1*age.res,PSM.s.wo.ex.Rod2H[[2]],lty=2)
lines(t:1*age.res,PSM.s.wo.ex.Rod2H[[3]],lty=2)

#check evaporation isotopes
PSM.s.wo.ex.evap.d18O <- MCMC.CI.bound(PSM.s.wo.ex$BUGSoutput$sims.list$evap.d18O,0.89)
plot(t:1*age.res,PSM.s.evap120.evap.d18O[[1]],type="l",ylim=c(-25,-10), main="Evap d18O")
lines(t:1*age.res,PSM.s.evap120.evap.d18O[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.evap.d18O[[3]],lty=2)

PSM.s.evap120.evap.d2H <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$evap.d2H,0.89)
plot(t:1*age.res,PSM.s.evap120.evap.d2H[[1]],type="l",ylim=c(-180,-80), main="Evap d2H")
lines(t:1*age.res,PSM.s.evap120.evap.d2H[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.evap.d2H[[3]],lty=2)

#check lake water isotopes
PSM.s.evap120.Lw.d18O <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$Lw.d18O,0.89)
plot(t:1*age.res,PSM.s.evap120.Lw.d18O[[1]],type="l",ylim=c(-10,0), main="Lake water d18O")
lines(t:1*age.res,PSM.s.evap120.Lw.d18O[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.Lw.d18O[[3]],lty=2)

PSM.s.evap120.Lw.d2H <- MCMC.CI.bound(PSM.s.evap120$BUGSoutput$sims.list$Lw.d2H,0.89)
plot(t:1*age.res,PSM.s.evap120.Lw.d2H[[1]],type="l",ylim=c(-100,0), main="Lake water d2H")
lines(t:1*age.res,PSM.s.evap120.Lw.d2H[[2]],lty=2)
lines(t:1*age.res,PSM.s.evap120.Lw.d2H[[3]],lty=2)

plot(PSM.s.evap120.Lw.d18O[[1]],PSM.s.evap120.Lw.d2H[[1]],xlim=c(-25,5),ylim=c(-200,20))
abline(a=Lw.intc, b=Lw.slope) #ok
#check meteoric line
points(PSM.s.evap120.Rod18O[[1]],PSM.s.evap120.Rod2H[[1]])
abline(a=MWL.intc, b=MWL.slope)
points(PSM.s.evap120.evap.d18O[[1]],PSM.s.evap120.evap.d2H[[1]],col="red")