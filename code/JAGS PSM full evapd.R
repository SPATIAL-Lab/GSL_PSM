model {
  
  #evaluation, use carb.bins
  for (i in 1:n.carb){
    #carbonate d18O
    d18O.car.dat[i] ~ dnorm(carb.d18O.dc[carb.bins[i]],1/d18O.car.sd[i]^2)
    
  }
  
  #evaluation, use cyst.bins
  for (i in 1:n.cyst){
    #Brine Shrimp cysts
    BScyst.d2H.dat[i] ~ dnorm(BScyst.d2H[cyst.bins[i]], 1/BScyst.d2H.sd[i]^2)
    
    BScyst.d18O.dat[i] ~ dnorm(BScyst.d18O[cyst.bins[i]], 1/BScyst.d18O.sd[i]^2)
    
  }
  
  #evaluation, use scwax and lcwax bins
  for (i in 1:n.sc.wax){
    #short chain wax
    scwax.d2H.dat[i] ~ dnorm(scwax.d2H[scwax.bins[i]], 1/scwax.d2H.sd[i]^2)
    
  }
  
  for (i in 1:n.lc.wax){
    #long chain wax
    lcwax.d2H.dat[i] ~ dnorm(lcwax.d2H.dc[lcwax.bins[i]], 1/lcwax.d2H.sd[i]^2)
    
  }
  
  ###SAMPLE MODEL###
  
  #################mapping depth onto t################# 
  
  #assign age bins
  carb.bins = t - trunc(age.carb.d/age.res) #pure age offset
  
  age.carb.d = post.age.carb[ind.dep,]
  
  lcwax.bins = t - trunc(age.lc.d/age.res)
  
  scwax.bins = t - trunc(age.sc.d/age.res)
  
  age.lc.d = post.age.lcwax[ind.dep,]
  
  age.sc.d = post.age.scwax[ind.dep,]
  
  cyst.bins = t - trunc(age.cyst.d/age.res)
  
  age.cyst.d = post.age.cyst[ind.dep,]
  
  ind.dep ~ dcat(rep(1, n.bac.it))
  
  # ind.dep = trunc(ind.dep.doub)
  # 
  # ind.dep.doub ~ dunif(1, n.bac.it+1)
  
  ###ARCH MODEL###
  #################signal attenuation with certain types of data#######################
  
  #modeling storage effect of long chain wax, creating a simulated downcore record
  #with averaging window t.avg.car, and averaging weights w.avg.lcwax for lc wax and
  
  #get weighted averaged lc wax
  # for (i in 1:(t - t.avg)){
  #   lcwax.d2H.dc[i] = sum(lcwax.d2H[i:(i + t.avg - 1)] * w.avg.lcwax)*(1 - f.C14d.lcwax) + f.C14d.lcwax* d2H.C14d.lcwax
  # }
  
  #get weighted averaged lc wax
  for (i in t.avg:t ){
    lcwax.d2H.dc[i] = sum(lcwax.d2H[(i - t.avg+1):i] * w.avg.lcwax)*(1 - f.C14d.lcwax) + f.C14d.lcwax* d2H.C14d.lcwax
  }
  
  d2H.C14d.lcwax ~ dnorm(-150,1/15^2)
  
  #evaluate lc wax 14C offset
  lcwax.offset ~ dnorm(lcwax.offset.m, 1/200^2) #uncertainty of 200 years for long chain wax 
  
  lcwax.offset.m = sum(abs((1 - t.avg):0) * age.res * w.avg.lcwax)*(1 - f.C14d.lcwax) + 50000 * f.C14d.lcwax
  
  #normalized weights
  w.avg.lcwax <- (w.trnpt.lcwax+w.stor.lcwax)/sum(w.trnpt.lcwax+w.stor.lcwax)
  
  w.trnpt.lcwax = dnorm(abs((1 - t.avg):0), trnpt.lcwax.mean, 1/10^2)
  
  trnpt.lcwax.mean ~ dnorm(15,1/5^2) T(1+5, t.avg-5) #a reasonable prior /dunif(1,t.avg-5) 
  
  # w.stor.lcwax = dexp(0:(t.avg - 1), stor.par.lcwax)
  
  w.stor.lcwax = dexp(abs((1 - t.avg):0), stor.par.lcwax)
  
  stor.par.lcwax ~ dbeta(5,40) #1/stor.par ~0.1, longer storage for lc wax
  
  #very little carbon dead lc alkanes
  f.C14d.lcwax ~ dbeta(2,100) #~0.5% dead material
  
  for (i in t.avg:t ){
    carb.d18O.dc[i] = sum(carb.d18O[(i - t.avg+1):i] * w.avg.carb)*(1 - f.C14d.carb) + d18O.C14d.carb * f.C14d.carb
  }
  d18O.C14d.carb ~ dnorm(-5,1/1^2)#lake carbonate has d18O near -5 per mil, sd = 1!!!!!!!!!!!!!!!!!!!!!!!!!
  
  #evaluate carbonate 14C offset
  carb.offset ~ dnorm(carb.offset.m, 1/200^2) #uncertainty of 200 years for carbonates 
  
  carb.offset.m = sum(abs((1 - t.avg):0) * age.res * w.avg.carb)*(1 - f.C14d.carb) + 50000 * f.C14d.carb
  
  #normalized weights
  w.avg.carb = (w.trnpt.carb+w.stor.carb)/sum(w.trnpt.carb+w.stor.carb)
  
  w.trnpt.carb = dnorm(abs((1 - t.avg):0), trnpt.carb.mean, 1/5^2)
  
  trnpt.carb.mean ~ dnorm(15,1/5^2) T(1+5, t.avg-5) #a reasonable prior
  
  w.stor.carb = dexp(abs((1 - t.avg):0), stor.par.carb)
  
  stor.par.carb ~ dbeta(5,20) #1/stor.par ~0.3, shorter storage for carbonates
  
  f.C14d.carb ~ dbeta(10,200) #~5% carbon dead material
  
  ###SEN MODEL###
  #convert vsmow to vpdb
  carb.d18O = 0.97001 * carb.d18O.vsmow - 29.99
  
  #convert to delta
  scwax.d2H = (rscwax.2H - 1) * 1e3
  
  for (i in 1:t){
    #BS cysts
    BScyst.d2H[i] = BScyst.slope.lw.2H * Lw.d2H[i] + BScyst.slope.diet.2H * BSdiet.d2H[i] + BScyst.inc.2H
    
    BScyst.d18O[i] = BScyst.slope.lw.18O * Lw.d18O[i] + BScyst.slope.diet.18O * BSdiet.d18O[i] + BScyst.inc.18O
    
    # BScyst.d2H[i] = 0.34 * Lw.d2H[i] + 0.26 * BSdiet.d2H[i] -92
    # 
    # BScyst.d18O[i] = 0.692 * Lw.d18O[i] + 0.101 * BSdiet.d18O[i] + 15.9
    
    #convert back to delta values
    BSdiet.d2H[i] = (rBSdiet.2H[i]-1) * 1e3
    
    BSdiet.d18O[i] = (rBSdiet.18O[i]-1) * 1e3
    
    #model algal cellulose d2H and d18O 
    #consider exchangable and non-exchangable H and O in cellulose
    #Filot 2006 for H isotope exchange of cellulose
    rBSdiet.2H[i] = rprx.L.2H[i] * (epsilon.2H.carbohy/1000 + 1) #* (1 - r.exH) + r.exH * rlw2H[i] * alpha.exH
    #Cernusak 2005 for exchangeable Oxygen in cellulose
    rBSdiet.18O[i] = rprx.L.18O[i] * (epsilon.18O.carbohy/1000 + 1) #*(1 - r.exO) + r.exO * rlw18O[i] * alpha.exO
    
    #short chain wax from proximal lake water and n-C18 acid
    rscwax.2H[i] = rprx.L.2H[i] * alpha2H.sc.alkane[i] * alpha.alk.acid
    
    
    #short chain wax n-C17 epsilon scale with salinity, Sachse and Sachs
    #Here try to use alpha for better consistency, and alpha cannot be larger than 1
    alpha2H.sc.alkane[i] = ifelse(prx.sal[i] * scwax.alpha.sl + scwax.alpha.inc > 1,1, prx.sal[i] * scwax.alpha.sl + scwax.alpha.inc)
    
    #long chain wax 
    #use the equation in McFarlin et al 2019, with slope and intercept
    lcwax.d2H[i] = lcwax.d2H.slope * d2H.MAP[i] + lcwax.d2H.inc
    
    d2H.MAP[i] = Ro.d2H[i] + d2H.gap.MAP_Ro
    
    #lake carbonate (aragonite) d18O (VPDB) from lake water d18O (VSMOW) #Kim et al 2007 
    carb.d18O.vsmow[i] = Lw.d18O[i] + 17.88 * (1000/LST.k[i]) - 31.14 + d18O.Mg.ef[i]#this is vsmow!
    
  }
  #Mg effect on d18O araganite precipitation (passive degassing)
  d18O.Mg.ef = -1 * sal*f.MgCl2/95.21*1000/902.5 #Kim et al 2007
  
  #n-alkane to n-alkanoic acid fractionation: Chikaraishi & Naraoka 2007
  
  alpha.alk.acid = epsilon.alk.acid * 1e-03 + 1 
  
  epsilon.alk.acid ~ dnorm (25, 1/16^2) T(0,45) #sd is large, but n-acid should be more 2H enriched
  
  #cellulose exchange ratio and alpha
  
  # r.exH ~ dnorm(0.37, 1/ 0.01^2) #Nielson and Bowen 2010: slope of 0.63 for algae grown in water
  # 
  # alpha.exH ~ dnorm(1.082, 1/0.014^2) #Filot et al 2006
  # 
  # r.exO ~ dnorm(0.31, 1/ 0.01^2) #Nielson and Bowen 2010: slope of 0.69 for algae grown in water
  # 
  # alpha.exO ~ dnorm(1.027, 1/0.001^2) #Cernusak et al 2005
  
  # Brine shrimp cyst slopes and intercepts, Nielson and Bowen 2010
  BScyst.slope.lw.2H ~ dnorm(0.34, 1/0.019^2)
  
  BScyst.slope.diet.2H ~ dnorm(0.26, 1/0.025^2)
  
  BScyst.inc.2H ~ dnorm(-92, 1/3.6^2)
  
  BScyst.slope.lw.18O ~ dnorm(0.692, 1/0.013^2)
  
  BScyst.slope.diet.18O ~ dnorm(0.101, 1/0.017^2)
  
  BScyst.inc.18O ~ dnorm(15.9, 1/0.31^2)
  
  #carbohydrate 2H fractionation for Brine Shrimp diet (algae), Estep and Hoering 1980 
  epsilon.2H.carbohy ~ dnorm(-100, 1/15^2)
  
  #carbohydrate 18O fractionation for Brine Shrimp diet (algae), DeNiro and Epstine 1981 
  epsilon.18O.carbohy ~ dnorm(27, 1/5^2)
  
  #simple version: use linear relationship by sachse and sachs 2008, but it does not apply to high salinity
  scwax.alpha.sl ~ dnorm(0.00080, 1/0.00005^2) T (0.0007,0.0009)
  
  scwax.alpha.inc ~ dnorm(0.80745, 1/0.005^2) T (0.79,0.83)
  
  # #use the equation in McFarlin et al 2019, with slope and intercept for n-acid
  lcwax.d2H.inc ~ dnorm(-125, 1/20^2) T(-150,-100)
  lcwax.d2H.slope ~ dnorm(0.62, 1/0.01^2)
  
  #Assuming a gap between MAP d2H and runoff, as a prescribed covariance
  d2H.gap.MAP_Ro ~ dnorm(d2H.gap.MAP_Ro.mean, 1/1^2) #allow some variation
  d2H.gap.MAP_Ro.mean ~ dnorm(45,1/5^2) #gap mean = 50 +-5, informed by modern value
  
  #get mid-chain alkanes epsilon and water mixture.
  #proximal lake water mixing
  #mixing with runoff
  
  #proximal lake salinity
  prx.sal = sal * prx.L.v/(prx.L.v + Runoff) #diluted by runoff
  
  #convert back to delta values
  prx.L.d2H = (rprx.L.2H - 1) * 1e3
  
  prx.L.d18O = (rprx.L.18O- 1) * 1e3
  
  rprx.L.2H = (prx.L.v * rlw2H + Runoff * rRo.2H)/(prx.L.v + Runoff)
  
  rprx.L.18O = (prx.L.v * rlw18O + Runoff * rRo.18O)/(prx.L.v + Runoff)
  
  prx.L.v = LV * f.m.ro
  #fraction of lake water that is mixed with runoff, let the model explore
  f.m.ro ~ dbeta(10,40) #~0.2
  
  ###EVN MODEL###
  
  #results: lake water d18O, d2H, salinity
  
  #evaporation choices? T, wind speed, albeno, and Rs: radiation?????????????????????????????????????????????
  #penman's equation (energy-balance equation), simplified by Valiantzas 2006 eq 32
  
  #convert back to delta values
  Lw.d2H = (rlw2H - 1) * 1e3 
  Lw.d18O = (rlw18O - 1) * 1e3 
  
  evap.d2H = 1e3 *(re2H - 1)
  evap.d18O = 1e3 *(re18O - 1)
  
  for (i in 2:t){
    #lake evaporation amount, km3
    Evap[i] = E.rate[i] * LA[i]/1000
    
    # E.rate[i] =  2.909* nsws * (LA[i]*1e6)^-0.05 * (S.coeff[i] - rh) * Vpfs[i]/(2.501-0.002361*LST[i])/1e6 * 210*1000
    
    #McMillan 1973 for wind function
    E.rate[i] =  (3.6 + 2.5 * nsws) * (5/LA[i])^0.05 * (S.coeff[i] - rh) * Vpfs[i]/(2.501-0.002361*LST[i])/1e6 * n.d.evap*1000
    #if S.coeff[i] - rh > 0, then evaporation happens, if not, condensation happens
    
    #fresh water saturation vapor pressure Teten's eq in kPa
    Vpfs[i] = 0.61078 * exp(17.269 * AT[i]/(237.3 + AT[i])) #saturated vapor pressure
    
    #evaporated water fractionation, ratio for water vapor in boundary layer air
    re2H[i] = (rlw2H[i]/Alpha.2H[i]-rh.n[i]*f*ra2H) / (((1-rh.n[i])/alphak.2H[i]) + (rh.n[i]*(1-f))) #Dee 2015

    re18O[i] = (rlw18O[i]/Alpha.18O[i]-rh.n[i]*f*ra18O) / (((1-rh.n[i])/alphak.18O[i]) + (rh.n[i]*(1-f))) #Dee 2015
    # 
    # re2H[i] = (rlw2H[i]/Alpha.2H[i]-rh.n[i]*ra2H) / (1-rh.n[i])*alphak.2H[i]
    # 
    # #re18O is the 18O ratio for water vapor in boundary layer air
    # re18O[i] = (rlw18O[i]/Alpha.18O[i]-rh.n[i]*ra18O) / (1-rh.n[i])*alphak.18O[i] 
    
    #equilibrium fractionations >1, temperature dependent, range : 0 to 100 degrees C
    Alpha.2H[i] = exp(24844/(LST.k[i]^2)- 76.248/LST.k[i] + 0.05261) * H.frac.coef[i] #Majoube 1971 + salt correction
    
    Alpha.18O[i] = exp(1137/(LST.k[i]^2)- 0.4156/LST.k[i] - 0.00207) * O.frac.coef[i] #Majoube 1971 + salt correction
    
    H.frac.coef[i] = 1- sal.corr.d2H/1000*(sal[i]/sal.mol.ms)
    
    O.frac.coef[i] = 1- sal.corr.d18O/1000*(sal[i]/sal.mol.ms)
    
    alphak.2H[i] = 1 - 12.5 * (1 - rh.n[i]) * (1 - f) * 1e-3 #wind speed lower than 7 m/s
    alphak.18O[i] = 1 - 14.2 * (1 - rh.n[i]) * (1 - f) * 1e-3
    
    rh.n[i] = rh/S.coeff[i] #salt water normalized vapor pressure
    
    #saline water vapor pressure coefficient using a combination of data
    S.coeff[i] = 1-9.098e-07*sal[i]^2.267
    
    #use salinity table to interpolate salinity, *GSL.volume and GSL.sali are inputs*
    sal[i] = interp.lin(LV[i], GSL.volume, GSL.sali)
    
    #lake area, use bathymetric table
    LA[i] = interp.lin(LV[i], GSL.volume, GSL.area)
    
    L.level[i] = interp.lin(LV[i], GSL.volume, GSL.level)
    
    # #Lake water isotope mass balance
    rlw2H[i] = (rlw2H[i - 1] * LV[i - 1]  + rRo.2H[i - 1] * Runoff[i - 1] - Evap[i - 1] * re2H[i - 1])/LV[i]
    
    rlw18O[i] = (rlw18O[i - 1] * LV[i - 1]  + rRo.18O[i - 1] * Runoff[i - 1] - Evap[i - 1] * re18O[i - 1])/LV[i]
    
    # #lake volume, use bathymetric table
    LV[i] = LV[i - 1] + Runoff[i - 1] - Evap[i - 1]
    
  }
  
  #define t = 1 for priming of the parameters
  #lake water isotope ratios after evaporation
  
  Evap[1] = E.rate[1]*LA[1]/1000
  
  #Finch and Calver 2008 #acceptible rates with a warm season bias
  #lake evaporation is happening primarily between April and October, 7 months
  # E.rate[1] =  2.909* nsws * (LA[1]*1e6)^-0.05 * (S.coeff[1] - rh) * Vpfs[1]/(2.501-0.002361*LST[1])/1e6 * 210*1000
  
  E.rate[1] =  (3.6 + 2.5 * nsws) * (5/LA[1])^0.05 * (S.coeff[1] - rh) * Vpfs[1]/(2.501-0.002361*LST[1])/1e6 * n.d.evap*1000
  
  #if S.coeff[i] - rh > 0, then evaporation happens, if not, condensation happens
  
  #fresh water saturation vapor pressure Teten's eq in kPa
  Vpfs[1] = 0.61078 * exp(17.269 * AT[1]/(237.3 + AT[1])) #saturated vapor pressure
  
  #####LST-dependent evaporation isotope calculations####
  
  
  re2H[1] = (rlw2H[1]/Alpha.2H[1]-rh.n[1]*f*ra2H) / (((1-rh.n[1])/alphak.2H[1]) + (rh.n[1] * (1 - f))) #Dee 2015

  #re18O is the 18O ratio for water vapor in boundary layer air
  re18O[1] = (rlw18O[1]/Alpha.18O[1]-rh.n[1]*f*ra18O) / (((1-rh.n[1])/alphak.18O[1]) + (rh.n[1] * (1 - f))) #Dee 2015
  
  # re2H[1] = (rlw2H[1]/Alpha.2H[1]-rh.n[1]*ra2H) / (1-rh.n[1])*alphak.2H[1]
  # 
  # #re18O is the 18O ratio for water vapor in boundary layer air
  # re18O[1] = (rlw18O[1]/Alpha.18O[1]-rh.n[1]*ra18O) / (1-rh.n[1])*alphak.18O[1] 
  
  #equilibrium fractionations >1, temperature dependent; Majoube, 1971, correction for salinity (Gibson 2016, Koehler 2013)
  Alpha.2H[1] = exp(24844/(LST.k[1]^2)- 76.248/LST.k[1] + 0.05261) * H.frac.coef[1]
  
  Alpha.18O[1] = exp(1137/(LST.k[1]^2)- 0.4156/LST.k[1] - 0.00207) * O.frac.coef[1]
  
  H.frac.coef[1] = 1- sal.corr.d2H/1000*(sal[1]/sal.mol.ms)
  
  O.frac.coef[1] = 1- sal.corr.d18O/1000*(sal[1]/sal.mol.ms)
  
  #kinetic fractionation (see Horita 2008)
  alphak.2H[1] = 1 - 12.5 * (1 - rh.n[1]) * (1 - f) * 1e-3 #wind speed lower than 7 m/s
  alphak.18O[1] = 1 - 14.2 * (1 - rh.n[1]) * (1 - f) * 1e-3
  
  rh.n[1] = rh/S.coeff[1] #salt water normalized vapor pressure
  
  S.coeff[1] = 1-9.098e-07*sal[1]^2.267 #salinity vapor pressure coefficent for evaporation
  
  sal[1] = interp.lin(L.level[1], GSL.level, GSL.sali)
  
  LA[1] = interp.lin(L.level[1], GSL.level, GSL.area)
  
  LV[1] = interp.lin(L.level[1], GSL.level, GSL.volume) #LV.A is the real lake volume at the end of the seasonal cycle
  
  # L.level[1] ~ dnorm(1280,1/3^2) T(1276,1284)#m
  
  L.level[1] ~ dunif(1276,1284)#m
  
  n.d.evap ~ dnorm(n.d.evap.int,1/10^2) 
  
  n.d.evap.int = 120  #120 days +- 10 days
  
  #####Water vapor isotopes initial values#####
  
  ##convert to water vapor ratios in air 
  ra2H = (air.d2H*1e-03) + 1 #convert to water vapor 2H ratio in air
  
  ra18O = (air.d18O*1e-03) + 1 #convert to water vapor 18O ratio in air
  
  air.d2H = air.d18O*V.slope.m + V.intc.m #calculate vapor isotopes
  
  #here V.intc and V.slope parameters are from *model input*
  V.intc.m ~ dnorm(V.intc, 1 / 1 ^ 2) #with some uncertainty
  
  V.slope.m ~ dnorm(V.slope, 1 / 0.1 ^ 2) #with some uncertainty
  
  air.d18O ~ dnorm(air.d18O.int, 1/0.5^2) #allowed some variation
  
  #an initial value that centers around estimates of modern warm season water vapor 
  air.d18O.int ~ dnorm(d18O.vap.warm , 1/2^2)
  
  #f: fraction of advected air over lake, Tanganyika is set at 0.3, GSL is set at a slightly higher value
  # f = 0.3
  f ~ dbeta(20, 50)
  
  #####Lake water isotopes initial values#####
  ##convert to lake water ratios
  rlw2H[1] = (Lw.d2H.int*1e-03) + 1 #convert to water vapor 2H ratio in air
  
  rlw18O[1] = (Lw.d18O.int*1e-03) + 1 #convert to water vapor 18O ratio in air
  
  #lake water d18O and d2H using a line
  #calculate Lw.d2H using d18O
  Lw.d2H.int = Lw.d18O.int * Lw.slope.m + Lw.intc.m
  
  #here parameters are from *model input*
  Lw.intc.m ~ dnorm(Lw.intc, 1 / 1 ^ 2) 
  
  Lw.slope.m ~ dnorm(5, 1 / 0.5 ^ 2) 
  
  #lake water initial value
  #an initial value that centers around modern estimates, uninformative prior
  Lw.d18O.int ~ dnorm(Lw.d18O.int.mean, 1/Lw.d18O.int.sd^2) 
  
  Lw.d18O.int.mean ~ dunif(-10,0)
  
  Lw.d18O.int.sd = 1
  
  #####Runoff isotopes time series####
  #runoff ratio from delta values
  rRo.2H = (Ro.d2H * 1e-3) + 1 #runoff 2H ratio from delta value
  
  rRo.18O = (Ro.d18O * 1e-3) + 1 #runoff 18O ratio from delta value
  
  for (i in 1:t){
    #runoff d18O and d2H are correlated and evolve along MWL, but not a time series
    Ro.d2H[i] = Ro.d18O[i]*Ro.slope.m + Ro.intc.m
    
    Ro.d18O[i] ~ dnorm(Ro.d18O.int.mean, 1/Ro.d18O.int.sd^2) T(-25,-9)
    
  }
  
  # for (i in 2:t){
  #   #runoff d18O and d2H are correlated and evolve along MWL, but not a time series
  #   Ro.d2H[i] = Ro.d18O[i]*Ro.slope.m + Ro.intc.m
  # 
  #   Ro.d18O[i] = Ro.d18O[i - 1] + Ro.d18O.cps[i]
  # 
  #   Ro.d18O.cps[i] ~ dnorm(LST.cps[i] * T.cps.slope * sl.cpsAT.18O, 1/1^2)
  # 
  # }
  # # covariance between runoff d18O and LST, use a normal distribution for the slope, applied to cps
  # Ro.d2H[1] = Ro.d18O[1]*Ro.slope.m + Ro.intc.m
  # 
  # Ro.d18O[1] = Ro.d18O.int.mean #initial value
  # 
  # Ro.d18O.cps[1] ~ dnorm(LST.cps[1] * T.cps.slope * sl.cpsAT.18O, 1/1^2)
  # 
  # #slope is from Sturm et al 2010, but consider LST being less variable than MAT
  # sl.cpsAT.18O ~ dnorm(0.5, 1/0.2^2) T(0,0.911) #lowest is 0 (not correlated)
  
  #parameters are from *model input*
  Ro.d18O.int.mean ~ dunif(-20,-15)
  
  Ro.d18O.int.sd = 2
  
  Ro.intc.m ~ dnorm(MWL.intc, 1 / 1 ^ 2)
  
  Ro.slope.m ~ dnorm(MWL.slope, 1 / 0.1 ^ 2)
  
  #####Runoff amount#####
  # not an autocorrelated time series
  for (i in 1:t){
    
    Runoff[i] ~ dlnorm(log(Runoff.int.mean), log(Runoff.pre)) #allow some variation
    
  }
  
  Runoff.int.mean ~ dnorm(3.5,1/0.5^2) T(2,)#3.5 +- 0.5, from Mohammed 2011
  
  Runoff.pre ~ dgamma(Runoff.pre.shp, Runoff.pre.rate)
  Runoff.pre.shp = 200
  Runoff.pre.rate = 2
  
  #####LST time series####
  #normal distribution, in degrees C
  for (i in 2:t){
    LST.k[i] = 273.15 + LST[i]
    
    AT[i] = LST[i - 1] + LST.cps[i] * T.cps.slope + T.gap #air temperature of the warm season is set at a constant offset
    
    LST[i] = LST[i - 1] + LST.cps[i]
    
    LST.cps[i] ~ dnorm(LST.cps[i - 1] * LST.cps.ac, LST.pre) T(-1,1)
    
  }
  
  AT[1] = LST[1] + LST.cps[1] * T.cps.slope + T.gap
  
  LST.cps[1] ~ dnorm(0, LST.pre) #allowed some variation
  
  LST.cps.ac ~ dunif(0.001, 0.8)
  
  #Air temperature covaries with LST, but at a slightly higher magnitude
  T.cps.slope ~ dnorm(1.5, 1/0.2^2) T(1,2)
  #temperature gap is modeled as stochastic
  T.gap ~ dnorm(5, 1/0.5^2) # allow some variation
  
  # initiate the series with an reasonable prior
  LST.k[1] = 273.15 + LST[1]
  LST[1] ~ dnorm(LST.int, LST.pre) #allowed some variation
  
  #an uninformative initial value: 20+-5 degrees C with a warm season bias, Steenburgh et al 2000
  LST.int ~ dunif(10, 30)  
  
  LST.pre ~ dgamma(LST.pre.shp, LST.pre.rate) # ~0.75 degrees error/100 years
  LST.pre.shp = 50
  LST.pre.rate = 2
  
  #####starting values####
  
  nsws ~ dnorm(5.8, 1/0.5^2) T(4,7)#wind speed data from Steenburgh, 2000
  
  #relative humidity ~0.35
  rh ~ dbeta(40, 72) T(0.1,0.5)
  
  #salt correction to equilibrium fractionation factors
  sal.mol.ms = (58.44*f.NaCl + 95.211*f.MgCl2)
  
  sal.corr.d18O = (NaCl.mc.O * f.NaCl + MgCl2.mc.O * f.MgCl2)
  sal.corr.d2H = (NaCl.mc.H * f.NaCl + MgCl2.mc.H * f.MgCl2)
  
  f.NaCl = 1 - f.MgCl2
  
  f.MgCl2 = 0.1092 #fraction from Dickson 1965
  
  #salt isotope correction, see Koehler et al 2013, mainly with NaCl and MgCl2
  NaCl.mc.O = 0.3
  
  MgCl2.mc.O = -1.03
  
  NaCl.mc.H = 1.8
  
  MgCl2.mc.H = 6.2
  
}