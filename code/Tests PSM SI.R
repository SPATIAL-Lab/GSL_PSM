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

setwd("C:/Users/ydmag/Google Drive/U of U/GSL proxy/GSL_PSM")

##############test evaporation rate calculation####

#if using simplified penman equation, may need to use annual average?
rh <- 0.4
LST <- 20
AT <- LST + 10
Vpfs.20 = 0.61078 * exp(17.269 * AT/(237.3 + AT))

E.rate.20 =  2.909* 5.8 * (GSL.area.1286*1e6)^-0.05 * (GSL.salicoeff - rh) * Vpfs.20/(2.501-0.002361* LST)/1e6 * 183*1000

E.vol.20 <- E.rate.20 * GSL.area.1286/1000
#reaching mass balance between 1279.703m and 1279.855m
LST <- 25
AT <- LST + 10
Vpfs.25 = 0.61078 * exp(17.269 * AT/(237.3 + AT))

E.rate.25 =  2.909* 5.8 * (GSL.area.1286*1e6)^-0.05 * (GSL.salicoeff - rh) * Vpfs.25/(2.501-0.002361* LST)/1e6 * 183*1000

E.vol.25 <- E.rate.25 * GSL.area.1286/1000

LST <- 15
AT <- LST + 10
Vpfs.15 = 0.61078 * exp(17.269 * AT/(237.3 + AT))

E.rate.15 =  2.909* 5.8 * (GSL.area.1286*1e6)^-0.05 * (GSL.salicoeff - rh) * Vpfs.15/(2.501-0.002361* LST)/1e6 * 183*1000

E.vol.15 <- E.rate.15 * GSL.area.1286/1000

plot(GSL.level.1286,E.vol.25,type="l",col="red",ylim=c(0,20),
     main = "Evaporation at rh = 0.4",xlab="Lake level",ylab="Evaporation")
lines(GSL.level.1286,E.vol.20,col="green")
lines(GSL.level.1286,E.vol.15,col="blue")
abline(h=3.3,lty=2,lwd = 2)
legend(1272,20,c("LST = 25 degrees","LST = 20 degrees","LST = 15 degrees","mean inflow"),
       col=c("red","green","blue","black"), lwd=c(1,1,1,2), lty=c(1,1,1,2))

#annual mean lake surface temperature is 16 degrees, max = 28 degrees warm season mean is probably ~22 degrees

rh <- 0.35
LST <- 20
AT <- LST + 10
Vpfs.20 = 0.61078 * exp(17.269 * AT/(237.3 + AT))

E.rate.20 =  2.909* 5.8 * (GSL.area.1286*1e6)^-0.05 * (GSL.salicoeff - rh) * Vpfs.20/(2.501-0.002361* LST)/1e6 * 183*1000 

E.vol.20 <- E.rate.20 * GSL.area.1286/1000
#reaching mass balance between 1279.703m and 1279.855m 
LST <- 25
AT <- LST + 10
Vpfs.25 = 0.61078 * exp(17.269 * AT/(237.3 + AT))

E.rate.25 =  2.909* 5.8 * (GSL.area.1286*1e6)^-0.05 * (GSL.salicoeff - rh) * Vpfs.25/(2.501-0.002361* LST)/1e6 * 183*1000 

E.vol.25 <- E.rate.25 * GSL.area.1286/1000

LST <- 15
AT <- LST + 10
Vpfs.15 = 0.61078 * exp(17.269 * AT/(237.3 + AT))

E.rate.15 =  2.909* 5.8 * (GSL.area.1286*1e6)^-0.05 * (GSL.salicoeff - rh) * Vpfs.15/(2.501-0.002361* LST)/1e6 * 183*1000 

E.vol.15 <- E.rate.15 * GSL.area.1286/1000

plot(GSL.level.1286,E.vol.25,type="l",col="red",ylim=c(0,20),
     main = "Evaporation at rh = 0.35",xlab="Lake level",ylab="Evaporation")
lines(GSL.level.1286,E.vol.20,col="green")
lines(GSL.level.1286,E.vol.15,col="blue")
abline(h=3.3,lty=2,lwd = 2)
legend(1272,20,c("LST = 25 degrees","LST = 20 degrees","LST = 15 degrees","mean inflow"),
       col=c("red","green","blue","black"), lwd=c(1,1,1,2), lty=c(1,1,1,2))

#evaporation occurs primarily in the warmer half of the year 
################use synthetic data to test lake mass balance#####
#use lake level history to derive precipitation amout
t <- 20
set.seed(1234)
LL.inc <- rnorm(20,0,0.3)
LL.mod <- rep(0,20)
LL.mod[1] <- rnorm(1,1280,0.5)
for(i in 2:t){
  LL.mod[i] <- LL.mod[i - 1] + LL.inc[i - 1]
}

plot(1:20,LL.mod, type = "l")

#parameters to save
parameters <- c("L.level","rh", "nsws", "LST", "Runoff", "sal", "S.coeff",
                "E.rate", "Evap", "LV", "LA", "LST.cps.ac")

dat = list(L.level.rec = LL.mod, GSL.level = GSL.level.1286, GSL.area = GSL.area.1286, GSL.volume = GSL.volume.1286,
           GSL.sali = GSL.sali, t = t)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 5e3
n.thin = 1

#Run it
test.lwmb = do.call(jags.parallel,list(model.file = "code/JAGS PSM Lwmb test.R", 
                                       parameters.to.save = parameters, 
                                       data = dat, n.chains=5, n.iter = n.iter, 
                                       n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #22 second

test.lwmb$BUGSoutput$summary

#check relative humidity
plot(density(test.lwmb$BUGSoutput$sims.list$rh)) #getting low

#check wind speed
plot(density(test.lwmb$BUGSoutput$sims.list$nsws)) #normal

#check LST ac
plot(density(test.lwmb$BUGSoutput$sims.list$LST.cps.ac)) #all over

#check LST
test.lwmb.LST <- MCMC.CI.bound(test.lwmb$BUGSoutput$sims.list$LST,0.89)
plot(1:20,test.lwmb.LST[[1]],type="l",ylim=c(15,30))
lines(1:20,test.lwmb.LST[[2]],lty=2)
lines(1:20,test.lwmb.LST[[3]],lty=2)
#by far, temperature has the least constraint

#check Runoff
test.lwmb.Runoff <- MCMC.CI.bound(test.lwmb$BUGSoutput$sims.list$Runoff,0.89)
plot(1:20,test.lwmb.Runoff[[1]],type="l",ylim=c(0,12))
lines(1:20,test.lwmb.Runoff[[2]],lty=2)
lines(1:20,test.lwmb.Runoff[[3]],lty=2)

#check evaporation
test.lwmb.Evap <- MCMC.CI.bound(test.lwmb$BUGSoutput$sims.list$Evap,0.89)
plot(1:20,test.lwmb.Evap[[1]],type="l")
lines(1:20,test.lwmb.Evap[[2]],lty=2)
lines(1:20,test.lwmb.Evap[[3]],lty=2)

###########################################test 2, with water isotopes####################
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

#supply slopes and intercepts
MWL.slope <- GSL.lake.mwl[1]
MWL.intc <- GSL.lake.mwl[2]

V.intc <- GSL.Vapor.warm.mwl[1]
V.slope <- GSL.Vapor.warm.mwl[2]

Lw.slope <- GSL.lake.mwl[1]
Lw.intc <- GSL.lake.mwl[2]

#warm season vapor d18O
d18O.vap.warm <- GSL.Vapor.warm.d18O

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

#check the gap between runoff d2H and lake water d2H
SI.lwmb.Rod2H <- MCMC.CI.bound(SI.lwmb$BUGSoutput$sims.list$Ro.d2H,0.89)
SI.lwmb.Lw.d2H<- MCMC.CI.bound(SI.lwmb$BUGSoutput$sims.list$Lw.d2H,0.89)
plot(1:t,SI.lwmb.Rod2H[[1]],type="l",ylim=c(-200,-40), main="Runoff vs lake d2H")
lines(1:t,SI.lwmb.Lw.d2H[[1]],col="green",lwd =2)
#the gap is around 80 permil, measured gap can be as low as 0, or as high as 40

#if epsilon higher plants are indeed -120/-100 for shrubs, then expected wax should be around -270/250!
#so there are additional processes involved!

#consider the mixture of three sources: autotropic and heterotrophic bacteria, 
#and higher plants (only small fraction of the mixture?), or mixing with an oil source (fixed marine source would lead to !

###############test isotope fractionation equations######

#first identify the source of lake water!
#slope of lake water, used to determine a combination of parameters such as rh, t and delta P?
5.031288

#equilibrium vapor from runoff
#-23.06999
#-167.39

#equilibrium fractionation, try 25 degrees C
LST.k <- 273.15+25

sal.corr.d18O = (NaCl.mc.O * f.NaCl + MgCl2.mc.O * f.MgCl2)
sal.corr.d2H = (NaCl.mc.H * f.NaCl + MgCl2.mc.H * f.MgCl2)

H.frac.coef <- 1- sal.corr.d2H/1000*(300/sal.mol.ms)
O.frac.coef <- 1- sal.corr.d18O/1000*(300/sal.mol.ms)

Alpha.2H = exp(24844/(LST.k^2)- 76.248/LST.k + 0.05261) * H.frac.coef

Alpha.18O = exp(1137/(LST.k^2)- 0.4156/LST.k - 0.00207) * O.frac.coef 

epsilon.2H = 1000*(Alpha.2H-1)

epsilon.18O = 1000*(Alpha.18O-1)


#0.5 equilibrium vapor from annual precip 
annp.eq.vap.d18O <- (d18O.prec.mean - 0.3*epsilon.18O)/(1+1e-3*epsilon.18O)

annp.eq.vap.d2H <- (d2H.prec.mean - 0.3*epsilon.2H)/(1+1e-3*epsilon.2H)

points(annp.eq.vap.d18O, annp.eq.vap.d2H, col = "red",pch=15) # yes, this is fine

#equilibrium from warm season precipitation
LST.k <- 273.15+25 #temperature is same as LST

Alpha.2H = exp(24844/(LST.k^2)- 76.248/LST.k + 0.05261)

Alpha.18O = exp(1137/(LST.k^2)- 0.4156/LST.k - 0.00207)

epsilon.2H = 1000*(Alpha.2H-1)

epsilon.18O = 1000*(Alpha.18O-1)

warm.eq.vap.d18O <- (GSL.precip.warm.d18O - epsilon.18O)/(1+1e-3*epsilon.18O)

warm.eq.vap.d2H <- (GSL.precip.warm.d2H - epsilon.2H)/(1+1e-3*epsilon.2H)

points(warm.eq.vap.d18O, warm.eq.vap.d2H, col = "blue",pch=19) # yes, this is fine

#mixing with vapor from lake water
Alpha.2H.co = exp(24844/(LST.k^2)- 76.248/LST.k + 0.05261) * H.frac.coef

Alpha.18O.co = exp(1137/(LST.k^2)- 0.4156/LST.k - 0.00207) * O.frac.coef 

epsilon.2H.co = 1000*(Alpha.2H.co-1)

epsilon.18O.co = 1000*(Alpha.18O.co-1)

lw.eq.vap.d18O <- (lw.d18O - epsilon.18O.co)/(1+1e-3*epsilon.18O.co)

lw.eq.vap.d2H <- (lw.d2H - epsilon.2H.co)/(1+1e-3*epsilon.2H.co)

f.v <- 0.1 #fraction of vapor that is from the lake, same as turbulence?

air.d2H.int =  f.v* lw.eq.vap.d2H +  (1-f.v)*warm.eq.vap.d2H
air.d18O.int =  f.v* lw.eq.vap.d18O +  (1-f.v)*warm.eq.vap.d18O

#this can be use as the initial value of air vapor 

#runoff 
Ro.d18O.int <- -17

Ro.d2H.int<- GSL.precip.mwl[1]*Ro.d18O.int +GSL.precip.mwl[2]

rh.n <- seq(0.5,1,by=0.05)
#normalized relative humidity

#equilibrium fractionation, try 25 degrees C
LST.k <- 273.15+25

sal.corr.d18O = (NaCl.mc.O * f.NaCl + MgCl2.mc.O * f.MgCl2)
sal.corr.d2H = (NaCl.mc.H * f.NaCl + MgCl2.mc.H * f.MgCl2)

H.frac.coef <- 1- sal.corr.d2H/1000*(300/sal.mol.ms)
O.frac.coef <- 1- sal.corr.d18O/1000*(300/sal.mol.ms)

Alpha.2H = exp(24844/(LST.k^2)- 76.248/LST.k + 0.05261) * H.frac.coef

Alpha.18O = exp(1137/(LST.k^2)- 0.4156/LST.k - 0.00207) * O.frac.coef 

epsilon.2H = 1000*(Alpha.2H-1)

epsilon.18O = 1000*(Alpha.18O-1)

f.advc <- 0.5 #fraction of advaction, lower advc, lower slope

#here do not use normalized rh! because it describes the kinetic fractionation due to diffusion
rh.10 <- rh.n #ambient air relative humidity, plus advection 
epsilon.kH <- 12.5*(1-rh.10)*(1-f.advc)
epsilon.kO <- 14.2*(1-rh.10)*(1-f.advc)

#test evaporation slope of GSL
SLEL.num <- (rh.n * (air.d2H.int-Ro.d2H.int)+(1+Ro.d2H.int*1e-3)*(epsilon.kH+epsilon.2H/Alpha.2H))/(rh.n * 1e3-epsilon.kH-epsilon.2H/Alpha.2H)

SLEL.denm <- (rh.n * (air.d18O.int-Ro.d18O.int)+(1+Ro.d18O.int*1e-3)*(epsilon.kO+epsilon.18O/Alpha.18O))/(rh.n * 1e3-epsilon.kO-epsilon.18O/Alpha.18O)

#alternative: annual precipitation, more negative

plot(rh.n,SLEL.num/SLEL.denm,ylim=c(0,10))
abline(h=5)
#ambient vapor isotope ratios
ra2H <- (air.d2H.int*1e-03) + 1 #-134.4 per mil
ra18O <- (air.d18O.int*1e-03) + 1 #-18.4 per mil

#alternative air: half equilibrium with annual precip

ra2H <- (GSL.Vapor.warm.d2H*1e-03) + 1 #-134.4 per mil
ra18O <- (GSL.Vapor.warm.d18O*1e-03) + 1 #-18.4 per mil

#alternative using warm season eq vapor
ra2H <- (warm.eq.vap.d2H*1e-03) + 1 #-134.4 per mil
ra18O <- (warm.eq.vap.d18O*1e-03) + 1 #-18.4 per mil

f <- 0.25 #turbulence coefficient (range: 0.3 - 0.2)

rh.10 <- 0.3 #ambient air relative humidity (range: 0.4-0.3)

#lake water

lw.d18O <- -5
lw.d2H <- GSL.lake.mwl[1]*lw.d18O + GSL.lake.mwl[2]

points(lw.d18O, lw.d2H, col = "blue",pch =16) # yes, this is fine
abline(a = GSL.precip.mwl[2],b= GSL.precip.mwl[1],col="blue")
abline(a = GSL.lake.mwl[2],b= GSL.lake.mwl[1],col="green")

rlw2H <- (lw.d2H*1e-03) + 1 #-60 per mil

rlw18O <- (lw.d18O*1e-03) + 1 #-4 per mil

epsilon.kH <- 12.5*(1-rh.10)*(1-f)
epsilon.kO <- 14.2*(1-rh.10)*(1-f)

alphak.2H = 1-epsilon.kH*1e-3
alphak.18O = 1-epsilon.kO*1e-3

LV <- 30 #km3
Ro <- 3.3 #km3
Evap <- c(5,10,15) #km3

#runoff 
Ro.d18O.int <- -17

Ro.d2H.int<- GSL.precip.mwl[1]*Ro.d18O.int +GSL.precip.mwl[2]

rRo.2H <- (Ro.d2H.int*1e-03) + 1 #-60 per mil
rRo.18O <- (Ro.d18O.int*1e-03) + 1 #-60 per mil

rlw2H.m = (rlw2H*LV + rRo.2H*Ro)/(LV+Ro)
rlw18O.m = (rlw18O*LV + rRo.18O*Ro)/(LV+Ro)

#back to delta values
rl.d2H.m = (rlw2H.m - 1) * 1e3
rl.d18O.m = (rlw18O.m - 1) * 1e3

re2H = (rlw2H.m/Alpha.2H-rh.10*f*ra2H) / (((1-rh.10)/alphak.2H) + (rh.10 * (1 - f))) #Dee 2015

re18O = (rlw18O.m/Alpha.18O-rh.10*f*ra18O) / (((1-rh.10)/alphak.18O) + (rh.10 * (1 - f))) #Dee 2015
#back to delta values
re.d2H = (re2H - 1) * 1e3
re.d18O = (re18O - 1) * 1e3

####alternative equation
# re.d2H = ((lw.d2H-epsilon.2H)/Alpha.2H-rh.n*air.d2H.int - epsilon.kH)/(1-rh.n+epsilon.kH*1e-3)
# 
# re.d18O = ((lw.d18O-epsilon.18O)/Alpha.18O-rh.n*air.d18O.int - epsilon.kO)/(1-rh.n+epsilon.kO*1e-3)

#evaporated water
re.d2H
# -127.1229, when rh.n = 0.8, f =0.5
re.d18O
# -11.78631, when rh.n = 0.8, f =0.5

points(re.d18O, re.d2H, col = "black",pch=20)

#runoff
points(Ro.d18O.int, Ro.d2H.int, col = "blue",pch=17)



points(lw.d18O, lw.d2H, col = "green",pch=15)

lw.d18O.n = (Ro * Ro.d18O.int - Evap * re.d18O + LV * lw.d18O)/(LV+Ro-Evap)
lw.d2H.n = (Ro * Ro.d2H.int - Evap * re.d2H + LV * lw.d2H)/(LV+Ro-Evap)

lw.d18O.n
lw.d2H.n

points(lw.d18O.n, lw.d2H.n, col = "black",pch=15)

