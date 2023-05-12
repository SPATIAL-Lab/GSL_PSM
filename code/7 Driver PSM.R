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
#use half pairwise difference as uncertainty around the mean#

d18O.car.dat <- GSL.carb$Ave.d18O

d18O.car.sd <- rep(mean(abs(diff(GSL.carb$Ave.d18O))/2), n.carb)

lcwax.d2H.dat <- GSL.FAME.C28.naom$C28.avg

lcwax.d2H.sd <- rep(mean(abs(diff(GSL.FAME.C28.naom$C28.avg))/2), n.lc.wax)

scwax.d2H.dat <- GSL.FAME.C16.naom$C16.avg

scwax.d2H.sd <- rep(mean(abs(diff(GSL.FAME.C16.naom$C16.avg))/2), n.sc.wax) 

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
                "epsilon.alk.acid","LST.pre","tl")


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
PSM.f = do.call(jags.parallel,list(model.file = "code/JAGS PSM full.R", 
                                     parameters.to.save = parameters, 
                                     data = dat, n.chains=5, n.iter = n.iter, 
                                     n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #8.5 hours for 1e4
#estimated time taken: t = 129 (centenial variation), 2.3 hours for n.iter = 2e4

save(PSM.f, file = "out/PSM.f.RData")