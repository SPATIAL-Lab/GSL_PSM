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

t.avg = 28 #number of time steps

#900 years for mean lc alkane offset, 
lcwax.offset = 900

# ~3000 years for carbonate offset
carb.offset = 3000

#number of time steps in the simulation
t <- ceiling(max.age/age.res) + t.avg #101 time steps + 28 for averaging

#######use calibrations for short chain wax fractionation
# post.leng.scwax <- length(post)
# 
# post.scwax.alpha.sl <- post
# post.scwax.alpha.inc <- post

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

d18O.car.dat <- GSL.carb$Ave.d18O

d18O.car.sd <- rep(0.2, n.carb)#fixed precision

lcwax.d2H.dat <- GSL.alk.C29.naom$C29.avg

lcwax.d2H.sd <- rep(5, n.lc.wax)#fixed precision

# scwax.d2H.dat <- sc.pw.d2H

# scwax.d2H.dat <- GSL.alk.C17.naom$C17.avg
# 
# scwax.d2H.sd <- rep(5, n.sc.wax)#fixed precision

C17.d2H.dat <- GSL.alk.C17.naom$C17.avg

C17.d2H.sd <- rep(5, n.sc.wax)#fixed precision

C19.d2H.dat<- GSL.alk.C19.naom$C19.avg

C19.d2H.sd <- rep(5, n.sc.wax)#fixed precision

r.C17.n <- ratio.17.n

BScyst.d18O.dat <- GSL.cyst$Ave.d18O

BScyst.d18O.sd <- rep(1, n.cyst)#fixed precision

BScyst.d2H.dat <- GSL.cyst$Ave.d2H

BScyst.d2H.sd <- rep(3, n.cyst) #fixed precision

#parameters to save
parameters <- c("L.level","rh", "nsws", "LST", "T.gap","AT","Runoff", "sal.A", "S.coeff",
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
                "alpha2H.gr.alkane", "alpha2H.cyan.alkane","E.I.ratio")


# dat = list(GSL.level = GSL.level.1286, GSL.area = GSL.area.1286, 
#            GSL.volume = GSL.volume.1286, GSL.sali = GSL.sali, t = t, age.res = age.res,
#            V.intc = V.intc, V.slope = V.slope, MWL.intc = MWL.intc, MWL.slope = MWL.slope,
#            Lw.intc = Lw.intc, Lw.slope = Lw.slope,
#            d18O.vap.warm = d18O.vap.warm,
#            n.bac.it = n.bac.it, n.carb = n.carb, n.cyst = n.cyst, n.lc.wax = n.lc.wax, n.sc.wax = n.sc.wax,
#            post.age.carb = post.age.carb, post.age.cyst = post.age.cyst, 
#            d18O.car.dat = d18O.car.dat, lcwax.d2H.dat = lcwax.d2H.dat, 
#            post.age.scwax = post.age.scwax, post.age.lcwax = post.age.lcwax,
#            scwax.d2H.dat = scwax.d2H.dat, BScyst.d18O.dat = BScyst.d18O.dat, 
#            BScyst.d2H.dat = BScyst.d2H.dat, r.C17.n=r.C17.n,
#            d18O.car.sd = d18O.car.sd, lcwax.d2H.sd = lcwax.d2H.sd, 
#            scwax.d2H.sd = scwax.d2H.sd, BScyst.d18O.sd = BScyst.d18O.sd, 
#            BScyst.d2H.sd = BScyst.d2H.sd, t.avg=t.avg, carb.offset= carb.offset, lcwax.offset=lcwax.offset)


dat = list(GSL.level = GSL.level.1286, GSL.area = GSL.area.1286, 
           GSL.volume = GSL.volume.1286, GSL.sali = GSL.sali, t = t, age.res = age.res,
           V.intc = V.intc, V.slope = V.slope, MWL.intc = MWL.intc, MWL.slope = MWL.slope,
           Lw.intc = Lw.intc, Lw.slope = Lw.slope,
           d18O.vap.warm = d18O.vap.warm,
           n.bac.it = n.bac.it, n.carb = n.carb, n.cyst = n.cyst, n.lc.wax = n.lc.wax, n.sc.wax = n.sc.wax,
           post.age.carb = post.age.carb, post.age.cyst = post.age.cyst, 
           d18O.car.dat = d18O.car.dat, lcwax.d2H.dat = lcwax.d2H.dat, 
           post.age.scwax = post.age.scwax, post.age.lcwax = post.age.lcwax,
           C17.d2H.dat = C17.d2H.dat, C19.d2H.dat = C19.d2H.dat, BScyst.d18O.dat = BScyst.d18O.dat, 
           BScyst.d2H.dat = BScyst.d2H.dat, r.C17.n=r.C17.n,
           d18O.car.sd = d18O.car.sd, lcwax.d2H.sd = lcwax.d2H.sd, 
           C17.d2H.sd = C17.d2H.sd, C19.d2H.sd = C19.d2H.sd, BScyst.d18O.sd = BScyst.d18O.sd, 
           BScyst.d2H.sd = BScyst.d2H.sd, t.avg=t.avg, carb.offset= carb.offset, lcwax.offset=lcwax.offset)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 5e3
n.thin = 5

#Run it
SI.f2 = do.call(jags.parallel,list(model.file = "code/JAGS PSM full2.R", 
                                     parameters.to.save = parameters, 
                                     data = dat, n.chains=5, n.iter = n.iter, 
                                     n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #8.5 hours for 1e4
#estimated time taken: t = 129 (centenial variation), 2.3 hours for n.iter = 2e4

SI.Arc$BUGSoutput$summary

par(mfrow=c(3,3))
#check relative humidity
plot(density(SI.f2$BUGSoutput$sims.list$rh),main="RH") #getting low

#check advection coefficient
# plot(density(SI.Arc$BUGSoutput$sims.list$f),main="F convection") #

#check wind speed
plot(density(SI.f2$BUGSoutput$sims.list$nsws),main="Wind") #normal

#check LST ac
plot(density(SI.f2$BUGSoutput$sims.list$LST.cps.ac),main="LST ac") #very high!

#check Ro ac
plot(density(SI.f2$BUGSoutput$sims.list$Ro.d18O.cps),main="Ro ac")

#check atm vapor d18O
plot(density(SI.f2$BUGSoutput$sims.list$air.d18O),main="Atm d18O") 
abline(v=d18O.vap.warm)

# check brine shrimp cyst intercept
plot(density(SI.f2$BUGSoutput$sims.list$BScyst.inc.2H),main="Cyst d2H interc")
abline(v=-92) #ok

#check brine shrimp cyst slopes
plot(density(SI.f2$BUGSoutput$sims.list$BScyst.slope.lw.2H),main="Cyst lw d2H slope")
abline(v=0.34) #ok

#check brine shrimp cyst intercept
plot(density(SI.f2$BUGSoutput$sims.list$BScyst.slope.diet.2H),main="Cyst diet d2H slope")
abline(v=0.26) #0.2

#check brine shrimp cyst intercept
plot(density(SI.f2$BUGSoutput$sims.list$BScyst.inc.18O),main="Cyst d18O interc")
abline(v=15.9)#ok
#check brine shrimp cyst slopes
plot(density(SI.f2$BUGSoutput$sims.list$BScyst.slope.lw.18O),main="Cyst lw d18O slope")
abline(v=0.692)#ok
#check brine shrimp cyst intercept
plot(density(SI.f2$BUGSoutput$sims.list$BScyst.slope.diet.18O),main="Cyst diet d18O slope")
abline(v=0.101)#ok

plot(density(SI.f2$BUGSoutput$sims.list$d2H.gap.MAP_Ro),main="MAP Ro gap") 
abline(v=45)

#epsilon BS diet (smaller than pure algea, means that there is trophic enrichment?)
plot(density(SI.f2$BUGSoutput$sims.list$epsilon.2H.carbohy),main="2H carbohy ep") 
abline(v=-100)#likely a mixture of both algae and terrestrial/wetland POM

plot(density(SI.f2$BUGSoutput$sims.list$epsilon.18O.carbohy),main="18O carbohy ep") 
abline(v=30) #much lower for d18O!

plot(density(SI.f2$BUGSoutput$sims.list$lcwax.d2H.inc),main="lc MAP inc") 
abline(v=-129) #much smaller than -129, at -115!, mixing with a -170 per mil source?

plot(density(SI.f2$BUGSoutput$sims.list$lcwax.d2H.slope),main="lc MAP slope")
abline(v=0.78) #0.76

plot(density(SI.f2$BUGSoutput$sims.list$scwax.alpha.sl),main="sc alpha slope")
abline(v=0.0008)#much lower fractionation for short chain alkane!, could also be inflacted salinity

plot(density(SI.f2$BUGSoutput$sims.list$scwax.alpha.inc),main="sc alpha interc")
abline(v=0.80745) # Consider make this fixed
#Or remove this relationship!

plot(density(SI.f2$BUGSoutput$sims.list$f.m.ro),main="f.m.ro") #more like 1, no premix?
abline(v=0.2) 

plot(density(SI.f2$BUGSoutput$sims.list$f.C14d.carb),main="f.C14d.carb") #ok
abline(v=0.05) 

plot(density(SI.f2$BUGSoutput$sims.list$f.C14d.lcwax),main="f.C14d.lcwax") #ok
abline(v=0.01) 

SI.f2.scalpha <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$alpha2H.sc.alkane,0.89)
plot(t:1*age.res, SI.f2.scalpha[[1]],type="l",ylim=c(0.6,1), main="sc alpha")
lines(t:1*age.res,SI.f2.scalpha[[2]],lty=2)
lines(t:1*age.res,SI.f2.scalpha[[3]],lty=2)

SI.f2.prx.sal <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$prx.sal,0.89)
plot(t:1*age.res, SI.f2.prx.sal[[1]],type="l",ylim=c(0,300), main="pre sal")
lines(t:1*age.res,SI.f2.prx.sal[[2]],lty=2)
lines(t:1*age.res,SI.f2.prx.sal[[3]],lty=2)

SI.f2.prx.L.d2H <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$prx.L.d2H,0.89)
plot(t:1*age.res, SI.f2.prx.L.d2H[[1]],type="l",ylim=c(-150,0), main="prx L d2H")
lines(t:1*age.res,SI.f2.prx.L.d2H[[2]],lty=2)
lines(t:1*age.res,SI.f2.prx.L.d2H[[3]],lty=2)

#check LST
SI.f2.LST <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$LST,0.89) 
plot(t:1*age.res,SI.f2.LST[[1]],type="l",ylim=c(10,40), main="LST")
lines(t:1*age.res,SI.f2.LST[[2]],lty=2)
lines(t:1*age.res,SI.f2.LST[[3]],lty=2)
#by far, Now this is more reasonable!

#check AT
SI.f2.AT <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$AT,0.89)
plot(t:1*age.res,SI.f2.AT[[1]],type="l",ylim=c(20,50), main="AT")
lines(t:1*age.res,SI.f2.AT[[2]],lty=2)
lines(t:1*age.res,SI.f2.AT[[3]],lty=2)

#check Runoff and evaporation
SI.f2.Runoff <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$Runoff,0.89)
plot(t:1*age.res,SI.f2.Runoff[[1]],type="l",ylim=c(0,12), main="Runoff and Evap")
lines(t:1*age.res,SI.f2.Runoff[[2]],lty=2)
lines(t:1*age.res,SI.f2.Runoff[[3]],lty=2)

#check evaporation
SI.f2.Evap <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$Evap,0.89)
lines(t:1*age.res,SI.f2.Evap[[1]],col="red")
lines(t:1*age.res,SI.f2.Evap[[2]],lty=2,col="red")
lines(t:1*age.res,SI.f2.Evap[[3]],lty=2,col="red")

#check evaporation rate
SI.f2.Evaprate <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$E.rate,0.89)
plot(t:1*age.res,SI.f2.Evaprate[[1]],type="l",main="Evap rate", ylim=c(0,2))
lines(t:1*age.res,SI.f2.Evaprate[[2]],lty=2)
lines(t:1*age.res,SI.f2.Evaprate[[3]],lty=2)

#check lake level
SI.f2.L.level <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$L.level,0.89)
plot(t:1*age.res,SI.f2.L.level[[1]],type="l",ylim=c(1276,1286), main="Lake level")
lines(t:1*age.res,SI.f2.L.level[[2]],lty=2)
lines(t:1*age.res,SI.f2.L.level[[3]],lty=2)

#check salinity
SI.f2.sal.A <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$sal.A,0.89)
plot(t:1*age.res,SI.f2.sal.A[[1]],type="l",ylim=c(50,360), main="salinity")
lines(t:1*age.res,SI.f2.sal.A[[2]],lty=2)
lines(t:1*age.res,SI.f2.sal.A[[3]],lty=2)

#check runoff isotopes
SI.f2.Rod18O <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$Ro.d18O,0.89)
plot(t:1*age.res,SI.f2.Rod18O[[1]],type="l",ylim=c(-25,-5), main="Runoff d18O")
lines(t:1*age.res,SI.f2.Rod18O[[2]],lty=2)
lines(t:1*age.res,SI.f2.Rod18O[[3]],lty=2)

SI.f2.Rod2H <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$Ro.d2H,0.89)
plot(t:1*age.res,SI.f2.Rod2H[[1]],type="l",ylim=c(-200,-50), main="Runoff d2H")
lines(t:1*age.res,SI.f2.Rod2H[[2]],lty=2)
lines(t:1*age.res,SI.f2.Rod2H[[3]],lty=2)

#check evaporation isotopes
SI.f2.evapd18O <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$evap.d18O,0.89)
plot(t:1*age.res,SI.f2.evapd18O[[1]],type="l",ylim=c(-30,-10), main="Evap d18O")
lines(t:1*age.res,SI.f2.evapd18O[[2]],lty=2)
lines(t:1*age.res,SI.f2.evapd18O[[3]],lty=2)

#check lake water isotopes
SI.f2.Lw.d18O <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$Lw.d18O,0.89)
plot(t:1*age.res,SI.f2.Lw.d18O[[1]],type="l",ylim=c(-10,0), main="Lake water d18O")
lines(t:1*age.res,SI.f2.Lw.d18O[[2]],lty=2)
lines(t:1*age.res,SI.f2.Lw.d18O[[3]],lty=2)

SI.f2.Lw.d2H <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$Lw.d2H,0.89)
plot(t:1*age.res,SI.f2.Lw.d2H[[1]],type="l",ylim=c(-80,10), main="Lake water d2H")
lines(t:1*age.res,SI.f2.Lw.d2H[[2]],lty=2)
lines(t:1*age.res,SI.f2.Lw.d2H[[3]],lty=2)

#check the gap between runoff d2H and lake water d2H
SI.f2.Rod2H <- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$Ro.d2H,0.89)
SI.f2.Lw.d2H<- MCMC.CI.bound(SI.f2$BUGSoutput$sims.list$Lw.d2H,0.89)
plot(t:1*age.res,SI.f2.Rod2H[[1]],type="l",ylim=c(-200,0), main="Runoff vs lake d2H")
lines(t:1*age.res,SI.f2.Lw.d2H[[1]],col="green",lwd =2)

#check lake water evaporation line
plot(SI.Arc.Lw.d18O[[1]],SI.Arc.Lw.d2H[[1]],xlim=c(-10,3),ylim=c(-100,30))
abline(a=Lw.intc, b=Lw.slope) #no good...

#check meteoric line
plot(SI.Arc.Rod18O[[1]],SI.Arc.Rod2H[[1]],xlim=c(-25,-10),ylim=c(-180,-80))
abline(a=MWL.intc, b=MWL.slope)

SI.Arc.w.carb <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$w.avg.carb,0.89)
plot(t.avg:1 -1,SI.Arc.w.carb[[1]],type="l", main="carbonate weights")
lines(t.avg:1 -1,SI.Arc.w.carb[[2]],lty=2)
lines(t.avg:1 -1,SI.Arc.w.carb[[3]],lty=2)

SI.Arc.w.lcwax <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$w.avg.lcwax,0.89)
plot(t.avg:1 -1,SI.Arc.w.lcwax[[1]],type="l", main="lc wax weights")
lines(t.avg:1 -1,SI.Arc.w.lcwax[[2]],lty=2)
lines(t.avg:1 -1,SI.Arc.w.lcwax[[3]],lty=2)


#check age-depth model with different substrates
#cysts
SI.Arc.age.cyst.d <- MCMC.CI.bound(SI.Arc$BUGSoutput$sims.list$age.cyst.d,0.89)
plot(GSL.cyst$Depth,SI.Arc.age.cyst.d[[1]],type="l",ylim=c(0,1e4), main="cyst age")
lines(GSL.cyst$Depth,SI.Arc.age.cyst.d[[2]],lty=2)
lines(GSL.cyst$Depth,SI.Arc.age.cyst.d[[3]],lty=2)

