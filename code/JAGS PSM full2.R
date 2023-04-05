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
    C17.d2H.dat[i] ~ dnorm(C17.d2H[scwax.bins[i]], 1/C17.d2H.sd[i]^2)
    
    C19.d2H.dat[i] ~ dnorm(C19.d2H[scwax.bins[i]], 1/C19.d2H.sd[i]^2)
    
    r.C17.n[i] ~ dt(1,f.cyan[scwax.bins[i]], 1/0.2^2) #weak constraint
    
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
  
  #evaluate lc wax 14C offset
  # lcwax.offset ~ dnorm(lcwax.offset.m, 1/200^2) #uncertainty of 200 years for long chain wax 
  # 
  # lcwax.offset.m = sum(0:(t.avg - 1) * age.res * w.avg.lcwax)*(1 - f.C14d.lcwax) + 50000 * f.C14d.lcwax
  
  lcwax.offset ~ dnorm(lcwax.offset.m, 1/200^2) #uncertainty of 200 years for long chain wax 
  
  lcwax.offset.m = sum(abs((1 - t.avg):0) * age.res * w.avg.lcwax)*(1 - f.C14d.lcwax) + 50000 * f.C14d.lcwax
  
  #normalized weights
  # w.avg.lcwax <- (w.trnpt.lcwax+w.stor.lcwax)/sum(w.trnpt.lcwax+w.stor.lcwax)
  
  w.avg.lcwax <- (w.trnpt.lcwax+w.stor.lcwax)/sum(w.trnpt.lcwax+w.stor.lcwax)
  
  # w.trnpt.lcwax = dnorm(0:(t.avg - 1), trnpt.lcwax.mean, 1/10^2)
  
  w.trnpt.lcwax = dnorm(abs((1 - t.avg):0), trnpt.lcwax.mean, 1/10^2)
  
  trnpt.lcwax.mean ~ dunif(1,t.avg-5) #uninformative prior
  
  # w.stor.lcwax = dexp(0:(t.avg - 1), stor.par.lcwax)
  
  w.stor.lcwax = dexp(abs((1 - t.avg):0), stor.par.lcwax)
  
  stor.par.lcwax ~ dbeta(10,80) #1/stor.par ~0.1, longer storage for lc wax
  
  #very little carbon dead lc alkanes
  f.C14d.lcwax ~ dbeta(2,100) #~2% dead material
  
  d2H.C14d.lcwax ~ dnorm(-150,1/30^2)#mesozoic lc wax is !-150 per mil, sd = 30!!!!!!!!!!!!!!!!!!!!!!!!!
  
  #get weighted averaged carbonate d18O
  # for (i in 1:(t - t.avg)){
  #   carb.d18O.dc[i] = sum(carb.d18O[i:(i + t.avg - 1)] * w.avg.carb)*(1 - f.C14d.carb) + d18O.C14d.carb * f.C14d.carb 
  # }
  
  for (i in t.avg:t ){
    carb.d18O.dc[i] = sum(carb.d18O[(i - t.avg+1):i] * w.avg.carb)*(1 - f.C14d.carb) + d18O.C14d.carb * f.C14d.carb
  }
  
  #evaluate carbonate 14C offset
  # carb.offset ~ dnorm(carb.offset.m, 1/200^2) #uncertainty of 200 years for carbonates 
  # 
  # carb.offset.m = sum(0:(t.avg - 1) * age.res * w.avg.carb)*(1 - f.C14d.carb) + 50000 * f.C14d.carb
  
  #evaluate carbonate 14C offset
  carb.offset ~ dnorm(carb.offset.m, 1/200^2) #uncertainty of 200 years for carbonates 
  
  carb.offset.m = sum(abs((1 - t.avg):0) * age.res * w.avg.carb)*(1 - f.C14d.carb) + 50000 * f.C14d.carb
  
  #normalized weights
  # w.avg.carb = (w.trnpt.carb+w.stor.carb)/sum(w.trnpt.carb+w.stor.carb)
  
  w.avg.carb = (w.trnpt.carb+w.stor.carb)/sum(w.trnpt.carb+w.stor.carb)
  
  # w.trnpt.carb = dnorm(0:(t.avg - 1), trnpt.carb.mean, 1/5^2)
  
  w.trnpt.carb = dnorm(abs((1 - t.avg):0), trnpt.carb.mean, 1/5^2)
  
  trnpt.carb.mean ~ dunif(1,t.avg-5) #uninformative prior
  
  # w.stor.carb = dexp(0:(t.avg - 1), stor.par.carb)
  
  w.stor.carb = dexp(abs((1 - t.avg):0), stor.par.carb)
  
  stor.par.carb ~ dbeta(5,10) #1/stor.par ~0.3, shorter storage for carbonates
  
  f.C14d.carb ~ dbeta(10,200) #~5% carbon dead material
  
  d18O.C14d.carb ~ dnorm(-5,1/2^2)#lake carbonate has d18O near -5 per mil, sd = 3!!!!!!!!!!!!!!!!!!!!!!!!!

  ###SEN MODEL###
  #convert vsmow to vpdb
  carb.d18O = 0.97001 * carb.d18O.vsmow - 29.99
  
  C17.d2H =(rwax.cyan.2H-1)* 1e3
  
  C19.d2H =(rwax.gr.2H-1)* 1e3
  
  #convert back to delta values
  BSdiet.d2H = (rBSdiet.2H-1) * 1e3
  
  BSdiet.d18O = (rBSdiet.18O-1) * 1e3
  
  for (i in 1:t){
    #BS cysts
    BScyst.d2H[i] = BScyst.slope.lw.2H * Lw.d2H[i] + BScyst.slope.diet.2H * BSdiet.d2H[i] + BScyst.inc.2H

    BScyst.d18O[i] = BScyst.slope.lw.18O * Lw.d18O[i] + BScyst.slope.diet.18O * BSdiet.d18O[i] + BScyst.inc.18O

    # BScyst.d2H[i] = 0.34 * Lw.d2H[i] + 0.26 * BSdiet.d2H[i] -92
    # 
    # BScyst.d18O[i] = 0.692 * Lw.d18O[i] + 0.101 * BSdiet.d18O[i] + 15.9
    
    #model algal cellulose d2H and d18O 
    rBSdiet.2H[i] = rprx.L.2H[i] * (epsilon.2H.carbohy/1000 + 1)* f.cyan[i] + rlw2H.A[i]* (epsilon.2H.carbohy/1000 + 1)* (1 - f.cyan[i])
    
    rBSdiet.18O[i] = rprx.L.18O[i] * (epsilon.18O.carbohy/1000 + 1)* f.cyan[i]+ rlw18O.A[i]* (epsilon.18O.carbohy/1000 + 1)*(1 - f.cyan[i])
    
    
    f.cyan[i] ~ dbeta(1, 1) #algae mixing ratio for BS diet
    
    #short chain wax as a mixture between proximal lake water and lake center
    # rscwax.2H[i] = rprx.L.2H[i] * alpha2H.sc.alkane[i]
    
    # rscwax.2H[i] = rwax.cyan.2H[i] * f.cyan[i] + rwax.gr.2H[i] * (1 - f.cyan[i])
    
    #assume all C17 is produced by cyanobacteria
    #cyanobacteria from proximal lake water, use salinity correction
    rwax.cyan.2H[i] = rprx.L.2H[i] * alpha2H.cyan.alkane[i]
    
    #assume all C19 is produced by green algae
    #green algae from lake water, use salinity correction
    rwax.gr.2H[i] = rlw2H.A[i] * alpha2H.gr.alkane[i]
    
    #short chain wax epsilon scale with salinity, Sachse and Sachs
    #Here try to use alpha for better consistency, and alpha cannot be larger than 1
    alpha2H.gr.alkane[i] = ifelse(prx.sal[i] * scwax.alpha.sl + scwax.alpha.inc > 1,1, sal.A[i] * scwax.alpha.sl + scwax.alpha.inc)
    
    alpha2H.cyan.alkane[i] = ifelse(prx.sal[i] * scwax.alpha.sl + scwax.alpha.inc > 1,1, prx.sal[i] * scwax.alpha.sl + scwax.alpha.inc)

    #long chain wax 
    #use the equation in McFarlin et al 2019, with slope and intercept
    lcwax.d2H[i] = lcwax.d2H.slope * d2H.MAP[i] + lcwax.d2H.inc

    #assuming local wax sources using runoff ratio
    d2H.MAP[i] = Ro.d2H[i] + d2H.gap.MAP_Ro
    
    # lcwax.d2H[i] = (rlcwax2H[i]- 1) * 1e3 
    # rlcwax2H[i] = rRo.2H[i]*(1 + epsilon2H.lcwax/1e3) #this is a fixed epsilon relationship

    #lake carbonate d18O #have to use aragonite equation (Mckinzie 1987)
    carb.d18O.vsmow[i] = Lw.d18O[i] + 17.88 * (1000/LST.k[i]) - 31.14 + d18O.Mg.ef[i]#this is vsmow!

  }
  #Mg effect on d18O araganite precipitation (passive degassing)
  d18O.Mg.ef = -1 * sal.A*f.MgCl2/95.211*1000/902.5 #Kim et al 2007
  
  #Brine shrimp cyst slopes and intercepts, Nielson and Bowen 2010
  BScyst.slope.lw.2H ~ dnorm(0.34, 1/0.019^2)

  BScyst.slope.diet.2H ~ dnorm(0.26, 1/0.025^2)

  BScyst.inc.2H ~ dnorm(-92, 1/3.6^2)

  BScyst.slope.lw.18O ~ dnorm(0.692, 1/0.013^2)

  BScyst.slope.diet.18O ~ dnorm(0.101, 1/0.017^2)

  BScyst.inc.18O ~ dnorm(15.9, 1/0.31^2)
  
  #carbohydrate 2H fractionation for Brine Shrimp diet (algae), Estep and Hoering 1980 
  epsilon.2H.carbohy ~ dnorm(-100, 1/15^2)
  
  #carbohydrate 18O fractionation for Brine Shrimp diet (algae), Savage 2021 
  epsilon.18O.carbohy ~ dnorm(30, 1/5^2)
  
  #simple version: use linear relationship by sachse and sachs 2008, but it does not apply to high salinity
  scwax.alpha.sl ~ dnorm(0.00080, 1/0.0001^2) T (0.0002,0.0015)
  
  scwax.alpha.inc ~ dnorm(0.80745, 1/0.05^2) T (0.79,0.83)
  
  #more complicated version: use regression 
  #use posterior of the short chain wax calibration in the calculation 
  # scwax.alpha.sl = post.scwax.alpha.sl[indx.scwax]
  # 
  # scwax.alpha.inc = post.scwax.alpha.inc[indx.scwax]
  # 
  # scwax.slope[indx.scwax]
  # scwax.inc[indx.scwax]
  # indx.scwax ~ dcat(rep(1, post.leng.scwax))
  
  #also explore alternatives, such as all plants from mid latitudes 30 - 60 degrees Sensitivity test
  #set up epsilon.lcwax, using , Cn-29, taken from Konecky
  # epsilon2H.lcwax ~ dnorm(epsilon2H.lcwax.mean, 1/epsilon2H.lcwax.sd^2)
  # epsilon2H.lcwax.mean = -121 #boreal trees, Kocecky 2019
  # epsilon2H.lcwax.sd = 22
  
  #use the equation in McFarlin et al 2019, with slope and intercept, alkanes
  lcwax.d2H.inc ~ dnorm(-129, 1/15^2)
  lcwax.d2H.slope ~ dnorm(0.78, 1/0.01^2)
  
  # #use the equation in McFarlin et al 2019, with slope and intercept for n-acid
  # lcwax.d2H.inc ~ dnorm(-125, 1/20^2)
  # lcwax.d2H.slope ~ dnorm(0.62, 1/0.01^2)
  
  #Assuming a gap between MAP d2H and runoff, as a prescribed covariance
  d2H.gap.MAP_Ro ~ dnorm(d2H.gap.MAP_Ro.mean, 1/1^2) #allow some variation
  d2H.gap.MAP_Ro.mean ~ dnorm(45,1/5^2) #gap mean = 50 +-5, informed by modern value
  
  #get mid-chain alkanes epsilon and water mixture.
  #proximal lake water mixing
  #mixing with runoff
  
  #proximal lake salinity
  prx.sal = sal.A * prx.L.v/(prx.L.v + Runoff) #diluted by runoff
  
  #convert back to delta values
  prx.L.d2H = (rprx.L.2H - 1) * 1e3
  
  prx.L.d18O = (rprx.L.18O- 1) * 1e3

  rprx.L.2H = (prx.L.v * rlw2H.A + Runoff * rRo.2H)/(prx.L.v + Runoff)
  
  rprx.L.18O = (prx.L.v * rlw18O.A + Runoff * rRo.18O)/(prx.L.v + Runoff)
  
  prx.L.v = LV.A * f.m.ro
  #fraction of lake water that is mixed with runoff, let the model explore
  f.m.ro ~ dbeta(1,1)
  
  ###EVN MODEL###
  #results: lake water d18O, d2H, salinity
  
  #evaporation choices? T, wind speed, albeno, and Rs: radiation?????????????????????????????????????????????
  #penman's equation (energy-balance equation), simplified by Valiantzas 2006 eq 32

  for (i in 2:t){
    #convert back to delta values
    Lw.d2H[i] = (rlw2H.A[i] - 1) * 1e3 
    Lw.d18O[i] = (rlw18O.A[i] - 1) * 1e3 
    
    #Lake water isotope mass balance for after evaporation
    rlw18O.A[i] = ( LV.P[i] * rlw18O.P[i] - Evap[i] * re18O[i] )/LV.A[i]
    
    rlw2H.A[i] = ( LV.P[i] * rlw2H.P[i] - Evap[i] * re2H[i])/LV.A[i]
    
    sal.A[i] = interp.lin(L.level[i], GSL.level, GSL.sali)
    
    L.level[i] = interp.lin(LV.A[i], GSL.volume, GSL.level)
    
    LV.A[i] = LV.P[i] - Evap[i] #LV.A is the real lake volume at the end of the seasonal cycle
    
    #lake evaporation amount, km3
    Evap[i] = E.rate[i] * LA.P[i]/1000
    
    E.rate[i] =  2.909* nsws * (LA.P[i]*1e6)^-0.05 * (S.coeff[i] - rh) * Vpfs[i]/(2.501-0.002361*LST[i])/1e6 * 183*1000
    #if S.coeff[i] - rh > 0, then evaporation happens, if not, condensation happens
    
    #fresh water saturation vapor pressure Teten's eq in kPa
    Vpfs[i] = 0.61078 * exp(17.269 * AT[i]/(237.3 + AT[i])) #saturated vapor pressure
    
    evap.d2H[i] = 1e3 *(re2H[i] - 1)
    evap.d18O[i] = 1e3 *(re18O[i] - 1)
    
    #evaporated water fractionation, ratio for water vapor in boundary layer air
    re2H[i] = (rlw2H.P[i]/Alpha.2H[i]-rh*f*ra2H) / (((1-rh)/alphak.2H) + (rh*(1-f))) #Dee 2015
    
    re18O[i] = (rlw18O.P[i]/Alpha.18O[i]-rh*f*ra18O) / (((1-rh)/alphak.18O) + (rh*(1-f))) #Dee 2015
    
    #equilibrium fractionations >1, temperature dependent, range : 0 to 100 degrees C
    Alpha.2H[i] = exp(24844/(LST.k[i]^2)- 76.248/LST.k[i] + 0.05261) * H.frac.coef[i] #Majoube 1971 + salt correction
    
    Alpha.18O[i] = exp(1137/(LST.k[i]^2)- 0.4156/LST.k[i] - 0.00207) * O.frac.coef[i] #Majoube 1971 + salt correction
    
    H.frac.coef[i] = 1- sal.corr.d2H/1000*(sal.P[i]/sal.mol.ms)
    
    O.frac.coef[i] = 1- sal.corr.d18O/1000*(sal.P[i]/sal.mol.ms)
    
    #saline water vapor pressure coefficient using a combination of data
    S.coeff[i] = 1-9.098e-07*sal.P[i]^2.267
    
    #use salinity table to interpolate salinity, *GSL.volume and GSL.sali are inputs*
    sal.P[i] = interp.lin(LV.P[i], GSL.volume, GSL.sali)
    
    #calculate mixed lake and runoff ratio before evaporation
    rlw2H.P[i] = (rlw2H.A[i - 1] * LV.A[i - 1]  + rRo.2H[i] * Runoff[i])/LV.P[i]
    
    rlw18O.P[i] = (rlw18O.A[i - 1] * LV.A[i - 1]  + rRo.18O[i] * Runoff[i])/LV.P[i]
    
    #lake area, use bathymetric table
    LA.P[i] = interp.lin(LV.P[i], GSL.volume, GSL.area)
    #lake volume, use bathymetric table
    LV.P[i] = LV.A[i - 1] + Runoff[i]
    
  }
  
  #define t = 1 for priming of the parameters
  #lake water isotope ratios after evaporation
  rlw18O.A[1] = ( LV.P[1] * rlw18O.P[1] - Evap[1] * re18O[1] )/LV.A[1]
  
  rlw2H.A[1] = ( LV.P[1] * rlw2H.P[1] - Evap[1] * re2H[1])/LV.A[1]
  
  L.level[1] = interp.lin(LV.A[1], GSL.volume, GSL.level)
  
  LV.A[1] = LV.P[1] - Evap[1] #LV.A is the real lake volume at the end of the seasonal cycle
  
  Evap[1] = E.rate[1]*LA.P[1]/1000
  
  #Finch and Calver 2008 #acceptible rates with a warm season bias
  E.rate[1] =  2.909* nsws * (LA.P[1]*1e6)^-0.05 * (S.coeff[1] - rh) * Vpfs[1]/(2.501-0.002361*LST[1])/1e6 * 183*1000
  #if S.coeff[i] - rh > 0, then evaporation happens, if not, condensation happens
  
  #fresh water saturation vapor pressure Teten's eq in kPa
  Vpfs[1] = 0.61078 * exp(17.269 * AT[1]/(237.3 + AT[1])) #saturated vapor pressure
  
  #####LST-dependent evaporation isotope calculations####
  evap.d2H[1] = 1e3 *(re2H[1] - 1)
  evap.d18O[1] = 1e3 *(re18O[1] - 1)
  
  #alphas need further correction for salinity (Gibson 2016, Koehler 2013)
  re2H[1] = (rlw2H.P[1]/Alpha.2H[1]-rh*f*ra2H) / (((1-rh)/alphak.2H) + (rh * (1 - f))) #Dee 2015
  
  #re18O is the 18O ratio for water vapor in boundary layer air
  re18O[1] = (rlw18O.P[1]/Alpha.18O[1]-rh*f*ra18O) / (((1-rh)/alphak.18O) + (rh * (1 - f))) #Dee 2015
  
  #equilibrium fractionations >1, temperature dependent; Majoube, 1971
  Alpha.2H[1] = exp(24844/(LST.k[1]^2)- 76.248/LST.k[1] + 0.05261) * H.frac.coef[1]
  
  Alpha.18O[1] = exp(1137/(LST.k[1]^2)- 0.4156/LST.k[1] - 0.00207) * O.frac.coef[1]
  
  
  H.frac.coef[1] = 1- sal.corr.d2H/1000*(sal.P[1]/sal.mol.ms)
  
  O.frac.coef[1] = 1- sal.corr.d18O/1000*(sal.P[1]/sal.mol.ms)
  
  S.coeff[1] = 1-9.098e-07*sal.P[1]^2.267 #salinity vapor pressure coefficent for evaporation
  
  sal.P[1] = interp.lin(L.level.int, GSL.level, GSL.sali)
  
  sal.A[1] = interp.lin(L.level.int, GSL.level, GSL.sali)
  
  #####Water vapor isotopes initial values#####
  
  ##convert to water vapor ratios in air 
  ra2H = (air.d2H*1e-03) + 1 #convert to water vapor 2H ratio in air
  
  ra18O = (air.d18O*1e-03) + 1 #convert to water vapor 18O ratio in air
  
  air.d2H = air.d18O*V.slope.m + V.intc.m #calculate vapor isotopes
  
  #here V.intc and V.slope parameters are from *model input*
  V.intc.m ~ dnorm(V.intc, 1 / 1 ^ 2) #with some uncertainty
  
  V.slope.m ~ dnorm(V.slope, 1 / 0.1 ^ 2) #with some uncertainty
  
  air.d18O ~ dnorm(air.d18O.int, 1/0.2^2) #allowed some variation
  
  #an initial value that centers around estimates of modern warm season water vapor 
  air.d18O.int ~ dnorm(d18O.vap.warm , 1/0.5^2)
  
  #kinetic fractionation (see Horita 2008)
  alphak.2H = 1 - 12.5 * (1 - rh) * (1 - f) * 1e-3 #wind speed lower than 7 m/s
  alphak.18O = 1 - 14.2 * (1 - rh) * (1 - f) * 1e-3
  
  #f: fraction of advected air over lake, stochastic, Tanganyika is set at 0.3, GSL is set at a similar value
  # f = 0.3
  f ~ dbeta(20, 50)
  
  #####Lake water isotopes initial values#####
  ##convert to lake water ratios
  rlw2H.P[1] = (Lw.d2H[1]*1e-03) + 1 #convert to water vapor 2H ratio in air
  
  rlw18O.P[1] = (Lw.d18O[1]*1e-03) + 1 #convert to water vapor 18O ratio in air
  
  #lake water d18O and d2H using a line
  #calculate Lw.d2H using d18O
  Lw.d2H[1] = Lw.d18O[1]*Lw.slope.m + Lw.intc.m
  
  #here parameters are from *model input*
  Lw.intc.m ~ dnorm(Lw.intc, 1 / 1 ^ 2) 
  
  Lw.slope.m ~ dnorm(Lw.slope, 1 / 0.1 ^ 2) 
  #lake water initial value
  Lw.d18O[1] ~ dnorm(Lw.d18O.int, 1/0.2^2) #allowed some variation
  
  #an initial value that centers around modern estimates, uninformative prior
  Lw.d18O.int ~ dnorm(Lw.d18O.int.mean, 1/Lw.d18O.int.sd^2) #T(-10,5)
  
  Lw.d18O.int.mean = -5
  
  Lw.d18O.int.sd = 3
  
  #####Runoff isotopes time series####
  #runoff ratio from delta values
  rRo.2H[1:t] = (Ro.d2H[1:t] * 1e-3) + 1 #runoff 2H ratio from delta value
  
  rRo.18O[1:t] = (Ro.d18O[1:t] * 1e-3) + 1 #runoff 18O ratio from delta value
  
  for (i in 1:t){
    #runoff d18O and d2H are correlated and evolve along MWL, but not a time series
    Ro.d2H[i] = Ro.d18O[i]*Ro.slope.m + Ro.intc.m

    Ro.d18O[i] ~ dnorm(Ro.d18O.int.mean, 1/Ro.d18O.int.sd^2) T(-25,-10)
      
  }
  
  #parameters are from *model input*
  Ro.d18O.int.mean = -17

  Ro.d18O.int.sd = 3
  
  Ro.intc.m ~ dnorm(MWL.intc, 1 / 1 ^ 2)
  
  Ro.slope.m ~ dnorm(MWL.slope, 1 / 0.1 ^ 2)
  
  #lake area, use bathymetric table
  LA.P[1] = interp.lin(LV.P[1], GSL.volume, GSL.area)
  #lake volume, use bathymetric table
  LV.P[1] = interp.lin(L.level.int, GSL.level, GSL.volume) + Runoff[1]
  
  #####Runoff amount#####
  # not an autocorrelated time series
  for (i in 1:t){
    
    Runoff[i] ~ dlnorm(log(Runoff.mean), 2*log(Runoff.pre)) #allow some variation
    
  }
  
  Runoff.mean ~ dnorm(3.2,1/0.5^2)
  
  Runoff.pre ~ dgamma(Runoff.pre.shp, Runoff.pre.rate)
  Runoff.pre.shp = 200
  Runoff.pre.rate = 2
  
  #####LST time series####
  #normal distribution, in degrees C
  for (i in 2:t){
    LST.k[i] = 273.15 + LST[i]
    
    AT[i] = LST[i] + T.gap #air temperature of the warm season is set at a constant offset
    
    LST[i] = LST[i - 1] + LST.cps[i]
    
    LST.cps[i] ~ dnorm(LST.cps[i - 1] * LST.cps.ac, LST.pre) T(-1,1)
    
  }
  LST.cps[1] ~ dnorm(0, LST.pre) #allowed some variation
  
  LST.cps.ac ~ dunif(0.001, 0.8)
  
  AT[1] = LST[1] + T.gap
  #temperature gap is modeled as stochastic
  T.gap ~ dnorm(T.gap.mean, 1/0.5^2) # allow some variation
  T.gap.mean ~ dnorm(5, 1/1^2)# try 5 degrees higher temp
  # initiate the series with an reasonable prior
  LST.k[1] = 273.15 + LST[1]
  LST[1] ~ dnorm(LST.int, LST.pre) #allowed some variation
  
  #an uninformative initial value: 20+-2 degrees C with a warm season bias, Steenburgh et al 2000
  # LST.int ~ dnorm(20, 1/5^2) T(10, 30)  
  LST.int ~ dunif(10, 30) #uninformative prior
  
  LST.pre ~ dgamma(LST.pre.shp, LST.pre.rate) # ~0.25 degrees error/100 years
  LST.pre.shp = 50
  LST.pre.rate = 2
  
  #####starting values####
  nsws ~ dnorm(nsws.mean, 1/0.5^2) #allow some variation
  nsws.mean ~ dnorm(5.8, 1/0.2^2) #wind speed data from Steenburgh, 2000
  
  #relative humidity ~0.35 +- 0.05
  rh ~ dbeta(40, 72) T(0.2,0.5)

  #lake starting level: 1280 +-1 meters
  L.level.int ~ dnorm(1280,1/2^2) T(1275,1286)#m
  
  sal.mol.ms = (58.44*f.NaCl + 95.211*f.MgCl2)
  
  sal.corr.d18O = (NaCl.mc.O * f.NaCl + MgCl2.mc.O * f.MgCl2)
  sal.corr.d2H = (NaCl.mc.H * f.NaCl + MgCl2.mc.H * f.MgCl2)
  
  f.NaCl = 1 - f.MgCl2
  
  f.MgCl2 = 0.1092 #fraction from Dickson 1965
  
  #salt isotope correction, Koehler et al 2013, mainly with NaCl and MgCl2
  NaCl.mc.O = 0.3
  
  MgCl2.mc.O = -1.03
  
  NaCl.mc.H = 1.8
  
  MgCl2.mc.H = 6.2
  
}