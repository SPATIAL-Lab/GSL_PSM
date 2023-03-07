model {
  #evaluation
  
  ###SAMPLE MODEL###
  
  ###ARCH MODEL###
  
  ###SEN MODEL###
  # for (i in 1:t){
  #   
  #   
  #   #short chain wax from lake water
  #   Lw.d2H[i]
  #   
  #   #short chain wax epsilon scale with salinity 
  #   sal[i]
  #   
  #   #long chain wax
  #   #assuming local wax sources using runoff d2H
  #   Ro.d2H[i]
  #   
  #   #lake carbonate d18O (VPDB) from lake water d18O (VSMOW)
  #   d18O.car[i] = Lw.d18O[i] - 0.27 + 25.8 - sqrt(25.8^2 - 11.1(16.1 - LST[i])) #Dee et al. 2018
  # 
  # }
  
  ###EVN MODEL###
  #results: lake water d18O, d2H, salinity
  

  #evaporation choices? T, wind speed, albeno, and Rs: radiation?????????????????????????????????????????????
  #penman's equation (energy-balance equation), simplified by Valiantzas 2006 eq 32
  #evaluation
  for (i in 1:t){
    # L.level.rec[i] ~ dnorm(L.level[i],1/0.1^2)
    
    Lw.d2H.dat[i] ~ dnorm(Lw.d2H[i],1/1^2)
    
    Lw.d18O.dat[i] ~ dnorm(Lw.d18O[i],1/0.1^2)
  
  }

  for (i in 2:t){
    #convert back to delta values
    Lw.d2H[i] = (rlw2H.A[i] - 1) * 1e3 
    Lw.d18O[i] = (rlw18O.A[i] - 1) * 1e3 
    
    #Lake water isotope mass balance for after evaporation
    rlw18O.A[i] = ( LV.P[i] * rlw18O.P[i] - Evap[i] * re18O[i] )/LV.A[i]
    
    rlw2H.A[i] = ( LV.P[i] * rlw2H.P[i] - Evap[i] * re2H[i])/LV.A[i]
    
    #lake evaporation amount, km3
    Evap[i] = E.rate[i] * LA.P[i]/1000
    
    L.level[i] = interp.lin(LV.A[i], GSL.volume, GSL.level)
    
    LV.A[i] = LV.P[i] - Evap[i] #LV.A is the real lake volume at the end of the seasonal cycle
    
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
    
    H.frac.coef[i] = 1- sal.corr.d2H/1000*(sal[i]/sal.mol.ms)
    
    O.frac.coef[i] = 1- sal.corr.d18O/1000*(sal[i]/sal.mol.ms)
    
    #saline water vapor pressure coefficient using a combination of data
    S.coeff[i] = 1-9.098e-07*sal[i]^2.267
    
    #use salinity table to interpolate salinity, *GSL.level and GSL.sali are inputs*
    sal[i] = interp.lin(L.level[i-1], GSL.level, GSL.sali)
    
    #calculate mixed lake and runoff ratio before evaporation
    rlw2H.P[i] = (rlw2H.A[i - 1] * LV.A[i - 1]  + rRo.2H[i] * Runoff[i])/LV.P[i]
    
    rlw18O.P[i] = (rlw18O.A[i - 1] * LV.A[i - 1]  + rRo.18O[i] * Runoff[i])/LV.P[i]
    
    #lake area, use bathymetric table
    LA.P[i] = interp.lin(LV.P[i], GSL.volume, GSL.area)
    #lake volume, use bathymetric table
    LV.P[i] = interp.lin(L.level[i-1], GSL.level, GSL.volume) + Runoff[i]
    
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
  
  
  H.frac.coef[1] = 1- sal.corr.d2H/1000*(sal[1]/sal.mol.ms)
  
  O.frac.coef[1] = 1- sal.corr.d18O/1000*(sal[1]/sal.mol.ms)
  
  S.coeff[1] = 1-9.098e-07*sal[1]^2.267 #salinity vapor pressure coefficent for evaporation
  
  sal[1] = interp.lin(L.level.int, GSL.level, GSL.sali)
  
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
  f ~ dnorm(f.mean, 1/0.05^2) #allow some variation
  f.mean ~ dnorm(0.3, 1/0.02^2) #~0.3 +- 0.02 rh
  
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
  Lw.d18O.int ~ dnorm(Lw.d18O.int.mean, 1/Lw.d18O.int.sd^2) T(-15,5)
  
  Lw.d18O.int.mean = -5
  
  Lw.d18O.int.sd = 4
  
  #####Runoff isotopes time series####
  #runoff ratio from delta values
  rRo.2H[1:t] = (Ro.d2H[1:t] * 1e-3) + 1 #runoff 2H ratio from delta value
  
  rRo.18O[1:t] = (Ro.d18O[1:t] * 1e-3) + 1 #runoff 18O ratio from delta value

  for (i in 2:t){
    #runoff d18O and d2H are correlated and evolve along MWL
    Ro.d2H[i] = Ro.d18O[i]*Ro.slope.m + Ro.intc.m
    
    Ro.d18O[i] = Ro.d18O[i - 1] + Ro.d18O.cps[i]
    
    Ro.d18O.cps[i] ~ dnorm(Ro.d18O.cps[i - 1] * Ro.cps.ac, Ro.d18O.pre)
    
  }
  Ro.d18O.cps[1] ~ dnorm(0, Ro.d18O.pre) #centered around 0, allowed some variation
  
  Ro.cps.ac ~ dunif(0, 0.8)
  
  Ro.d18O.pre ~ dgamma(Runoff.pre.shp, Runoff.pre.rate) # ~0.5 per mil error/10 years
  Ro.d18O.pre.shp = 16
  Ro.d18O.pre.rate = 4
  
  #calculate Ro.d2H using d18O
  Ro.d2H[1] = Ro.d18O[1]*Ro.slope.m + Ro.intc.m
  
  #parameters are from *model input*
  Ro.intc.m ~ dnorm(MWL.intc, 1 / 1 ^ 2) 
  
  Ro.slope.m ~ dnorm(MWL.slope, 1 / 0.1 ^ 2) 
  
  #Runoff initial value, uninformative prior
  Ro.d18O[1] ~ dnorm(Ro.d18O.int, 1/0.2^2) #allowed some variation
  
  #an initial value that centers around modern estimates, but uninformative
  Ro.d18O.int ~ dnorm(Ro.d18O.int.mean, 1/Ro.d18O.int.sd^2) T(-25,-10)
  
  Ro.d18O.int.mean = -17
  
  Ro.d18O.int.sd = 2
  
  #lake area, use bathymetric table
  LA.P[1] = interp.lin(LV.P[1], GSL.volume, GSL.area)
  #lake volume, use bathymetric table
  LV.P[1] = interp.lin(L.level.int, GSL.level, GSL.volume) + Runoff[1]
  
  #####Runoff amount#####
  # not an autocorrelated time series
  for (i in 1:t){
    
    Runoff[i] ~ dlnorm(log(Runoff.int.mean), 2*log(Runoff.pre)) #allow some variation
    
  }
  
  Runoff.int.mean ~ dnorm(3.5,1/0.5^2)
  
  Runoff.pre ~ dgamma(Runoff.pre.shp, Runoff.pre.rate)
  Runoff.pre.shp = 100
  Runoff.pre.rate = 2
  
  #####LST time series####
  #normal distribution, in degrees C
  for (i in 2:t){
    LST.k[i] = 273.15 + LST[i]
    
    AT[i] = LST[i] + T.gap #air temperature of the warm season is set at a constant offset
    
    LST[i] = LST[i - 1] + LST.cps[i]
    
    LST.cps[i] ~ dnorm(LST.cps[i - 1] * LST.cps.ac, LST.pre)
    
  }
  LST.cps[1] ~ dnorm(0, LST.pre) #allowed some variation
  
  LST.cps.ac ~ dunif(0.001, 0.8)
  
  AT[1] = LST[1] + T.gap
  #temperature gap is modeled as stochastic
  T.gap ~ dnorm(T.gap.mean, 1/0.5^2) # allow some variation
  T.gap.mean ~ dnorm(10, 1/1^2)
  # initiate the series with an reasonable prior
  LST.k[1] = 273.15 + LST[1]
  LST[1] ~ dnorm(LST.int, LST.pre) #allowed some variation
  
  #an uninformative initial value: 20+-5 degrees C with a warm season bias, Steenburgh et al 2000
  LST.int ~ dnorm(20, 1/5^2) T(10, 30)  
  
  LST.pre ~ dgamma(LST.pre.shp, LST.pre.rate) # ~0.15 degrees error/10 years
  LST.pre.shp = 100
  LST.pre.rate = 2
  
  #####starting values####
  
  nsws ~ dnorm(nsws.mean, 1/0.5^2) #allow some variation
  nsws.mean ~ dnorm(5.8, 1/0.2^2) #wind speed data from Steenburgh, 2000
  
  rh ~ dnorm(rh.mean, 1/0.02^2) #allow some variation
  rh.mean ~ dnorm(0.35, 1/0.05^2) #~0.35 +- 0.05 rh
  
  #lake starting level: 1280 +-1 meters
  L.level.int ~ dnorm(1280,1/1^2) T(1275,1286)#m
  
  # #NaCl, MgCl2, 
  # salt.molar.mass = (58.44*0.7591 + 95.211*0.1092)
  # #NaCl, MgCl2, Na2SO4, KCl
  
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