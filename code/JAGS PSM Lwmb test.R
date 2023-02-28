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
  
  #input values: T, wind speed, albeno, and Rs: radiation

  #evaporation choices
  #penman's equation (energy-balance equation), simplified by Valiantzas 2006 eq 32
  
  #evaluation
  for (i in 1:t){
    L.level.rec[i] ~ dnorm(L.level[i],1/0.1^2)
  }

  

  for (i in 2:t){
    #convert back to delta values
    # Lw.d2H[i] = (rlw2H[i] - 1) * 1e3 #d18O of boundary layer air
    # Lw.d18O[i] = (rlw18O[i] - 1) * 1e3 #d18O of boundary layer air
    
    #lake evaporation amount, km3
    Evap[i] = E.rate[i] * LA[i]/1000
    
    #water evaporation rate equation, convert to meters
    #Finch and Calver 2008 # produces acceptable rates with a warm season bias
    E.rate[i] =  2.909* nsws * (LA[i]*1e6)^-0.05 * (S.coeff[i] - rh) * Vpfs[i]/(2.501-0.002361*LST[i])/1e6 * 183*1000
    #if S.coeff[i] - rh > 0, then evaporation happens, if not, condensation happens
    
    #fresh water saturation vapor pressure Teten's eq in kPa
    Vpfs[i] = 0.61078 * exp(17.269 * AT[i]/(237.3 + AT[i])) #saturated vapor pressure
    
    #interpolate lake area from lake volume, use bathymetric table *GSL.area and GSL.volume are inputs*
    LA[i] = interp.lin(LV[i], GSL.volume, GSL.area)
    
    #interpolate lake level from lake volume, use bathymetric table *GSL.level and GSL.volume are inputs*
    L.level[i] = interp.lin(LV[i], GSL.volume, GSL.level)
    
    #LST-dependent fractionation
    #re2H is the 2H ratio for water vapor in boundary layer air
    # re2H[i] = (rlw2H[i]/Alpha.2H[i]-rh.n[i]*f*ra2H) / (((1.-rh.n[i])/alphak) + (rh.n[i]*(1.0-f))) #Dee 2015
    # 
    # #re18O is the 18O ratio for water vapor in boundary layer air
    # re18O[i] = (rlw18O[i]/Alpha.18O[i]-rh.n[i]*f*ra18O) / (((1.-rh.n[i])/alphak) + (rh.n[i]*(1.0-f))) #Dee 2015
    # 
    # #equilibrium fractionations, temperature dependent
    # Alpha.2H[i] = exp(24844/(LST.k[i]^2)- 76.248/LST.k[i] + 0.05261) #Dee 2015
    # 
    # Alpha.18O[i] = exp(1137/(LST.k[i]^2)- 0.4156/LST.k[i] - 0.00207) #Dee 2015
    # 
    # #normalize relative humidity, using saline water vapor pressure coefficient
    # rh.n[i] = ifelse (rh > S.coeff[i], 1, rh/S.coeff[1])

    #saline water vapor pressure coefficient using a combination of data
    S.coeff[i] = 1-9.098e-07*sal[i]^2.267
    
    #use salinity table to interpolate salinity, *GSL.level and GSL.sali are inputs*
    sal[i] = interp.lin(L.level[i], GSL.level, GSL.sali)
    
    #Lake water isotope mass balance
    
    # rlw18O[i] = (Runoff[i - 1] * rRo.18O[i - 1] - Evap[i - 1] * re18O[i - 1] + LV[i - 1] * rlw18O[i - 1])/LV[i]
    # 
    # rlw2H[i] = (Runoff[i - 1] * rRo.2H[i - 1] - Evap[i - 1] * re2H[i - 1] + LV[i - 1] * rlw2H[i - 1])/LV[i]
    
    #Lake water mass balance
    LV[i] = LV[i-1] + Runoff[i-1] - Evap[i-1]
    
  }
  

  Evap[1] = E.rate[1]*LA[1]/1000
  
  #Finch and Calver 2008 #acceptible rates with a warm season bias
  E.rate[1] =  2.909* nsws * (LA[1]*1e6)^-0.05 * (S.coeff[1] - rh) * Vpfs[1]/(2.501-0.002361*LST[1])/1e6 * 183*1000
  #if S.coeff[i] - rh > 0, then evaporation happens, if not, condensation happens
  
  #fresh water saturation vapor pressure Teten's eq in kPa
  Vpfs[1] = 0.61078 * exp(17.269 * AT[1]/(237.3 + AT[1])) #saturated vapor pressure

  # rlw2H[1] = (Lw.d2H[1] * 1e-3) + 1 #lake water 2H ratio from delta value
  # rlw18O[1] = (Lw.d18O[1] * 1e-3) + 1 #lake water 18O ratio from delta value
  # 
  
  #####LST-dependent evaporation isotope calculations####

  #alpha.2H may need further correction for salinity, 18O does not (Gibson 2016)!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  # re2H[1] = (rlw2H[1]/Alpha.2H[1]-rh.n[1]*f*ra2H) / (((1.-rh.n[1])/alphak) + (rh.n[1] * (1 - f))) #Dee 2015
  # 
  # #re18O is the 18O ratio for water vapor in boundary layer air
  # re18O[1] = (rlw18O[1]/Alpha.18O[1]-rh.n[1]*f*ra18O) / (((1.-rh.n[1])/alphak) + (rh.n[1] * (1 - f))) #Dee 2015
  # 
  # #equilibrium fractionations, temperature dependent; Majoube, 1971
  # Alpha.2H[1] = exp(24844/(LST.k[1]^2)- 76.248/LST.k[1] + 0.05261)
  # 
  # Alpha.18O[1] = exp(1137/(LST.k[1]^2)- 0.4156/LST.k[1] - 0.00207)
  # 
  # #salinity normalized relative humidity
  # rh.n[1] = ifelse (rh > S.coeff[1], 1, rh/S.coeff[1])
  # 
  S.coeff[1] = 1-9.098e-07*sal[1]^2.267
  
  sal[1] <- interp.lin(L.level[1], GSL.level, GSL.sali)
  
  #####Water vapor isotopes initial values#####
  ##convert to water vapor ratios in air 
  # ra2H <- (air.d2H*1e-03) + 1 #convert to water vapor 2H ratio in air
  # 
  # ra18O <- (air.d18O*1e-03) + 1 #convert to water vapor 18O ratio in air
  # 
  # air.d2H = air.d18O*V.slope.m + V.intc.m #calculate vapor isotopes
  # 
  # #here V.intc and V.slope parameters are from *model input*
  # V.intc.m ~ dnorm(V.intc, 1 / 1 ^ 2) #with some uncertainty
  # 
  # V.slope.m ~ dnorm(V.slope, 1 / 0.1 ^ 2) #with some uncertainty
  # #lake water initial value
  # air.d18O ~ dnorm(air.d18O.int, 1/0.2^2) #allowed some variation
  # 
  # #an initial value that centers around estimates of modern seasonal water vapor 
  # air.d18O.int ~ dnorm((1-f.sea) * d18O.prec.mean + f.sea *d18O.V.mean, 1/0.5^2)
  # 
  # f.sea <- 0.5 #effective vapor for seasonal environments, Gibson 2016
  # 
  # #f: fraction of advected air over lake, stochastic, Tanganyika is 0.3
  # f ~ dnorm(0.5, 1/0.1^2)  #relatively dry place, also allow some variation
  # 
  # #fractionation coeff of kinetic fractionation (open ocean diffusivity, see Haese et al. 2013)
  # alphak = 0.994        
  # 
  # #####Lake water isotopes initial values#####
  # ##convert to lake water ratios
  # rlw2H[1] <- (Lw.d2H*1e-03) + 1 #convert to water vapor 2H ratio in air
  # 
  # rlw18O[1] <- (Lw.d18O*1e-03) + 1 #convert to water vapor 18O ratio in air
  # #lake water d18O and d2H using a line
  # #calculate Lw.d2H using d18O
  # Lw.d2H[1] = Lw.d18O[1]*Lw.slope.m + Lw.intc.m
  # 
  # #here Lw.intc[1:2] and Lw.slope[1:2] parameters are from *model input*
  # Lw.intc.m ~ dnorm(Lw.intc[1], 1 / Lw.intc[2] ^ 2) 
  # 
  # Lw.slope.m ~ dnorm(Lw.slope[1], 1 / Lw.slope[2] ^ 2) 
  # #lake water initial value
  # Lw.d18O[1] ~ dnorm(Lw.d18O.int, 1/0.2^2) #allowed some variation
  # 
  # #an initial value that centers around modern estimates, uninformative prior
  # Lw.d18O.int ~ dnorm(Lw.d18O.int.mean, 1/Lw.d18O.int.sd^2) T(-15,5)
  # 
  # Lw.d18O.int.mean = -5
  # 
  # Lw.d18O.int.sd = 4
  # 
  # #####Runoff isotopes time series#####
  # #runoff ratio from delta values
  # rRo.2H[1:t] = (Ro.d2H[1:t] * 1e-3) + 1 #runoff 2H ratio from delta value
  # 
  # rRo.18O[1:t] = (Ro.d18O[1:t] * 1e-3) + 1 #runoff 18O ratio from delta value
  # 
  # for (i in 2:t){
  #   #runoff d18O and d2H are correlated and evolve along MWL
  #   Ro.d2H[i] = Ro.d18O[i]*Ro.slope.m + Ro.intc.m
  #   
  #   Ro.d18O[i] = Ro.d18O[i - 1] + Ro.d18O.cps[i]
  #   
  #   Ro.d18O.cps[i] ~ dnorm(Ro.d18O.cps[i - 1] * Ro.cps.ac, Ro.d18O.pre)
  #   
  # }
  # Ro.d18O.cps[1] ~ dnorm(0, Ro.d18O.pre) #centered around 0, allowed some variation
  # 
  # Ro.cps.ac ~ dunif(0, 0.8)
  # 
  # Ro.d18O.pre ~ dgamma(Runoff.pre.shp, Runoff.pre.rate) # ~0.5 per mil error/10 years
  # Ro.d18O.pre.shp = 16
  # Ro.d18O.pre.rate = 4
  # 
  # #calculate Ro.d2H using d18O
  # Ro.d2H[1] = Ro.d18O[1]*Ro.slope.m + Ro.intc.m
  # 
  # #here MWL.intc[1:2] and MWL.slope[1:2] parameters are from *model input*
  # Ro.intc.m ~ dnorm(MWL.intc[1], 1 / MWL.intc[2] ^ 2) 
  # 
  # Ro.slope.m ~ dnorm(MWL.slope[1], 1 / MWL.slope[2] ^ 2) 
  # 
  # #Runoff initial value, uninformative prior
  # Ro.d18O[1] ~ dnorm(Ro.d18O.int, 1/0.2^2) #allowed some variation
  # 
  # #an initial value that centers around modern estimates
  # Ro.d18O.int ~ dnorm(Ro.d18O.int.mean, 1/Ro.d18O.int.sd^2) T(-30,0)
  # 
  # Ro.d18O.int.mean <- -15
  # 
  # Ro.d18O.int.sd <- 8
  # 
  #####Runoff amount#####
  #not an autocorrelated time series 
  for (i in 1:t){
    
    Runoff[i] ~ dlnorm(log(Runoff.int), Runoff.pre) #allow some variation
  
  }
 
  #an initial value that centers around modern estimates
  Runoff.int ~ dlnorm(Runoff.int.logmean, 1/Runoff.int.logsd^2) 
  
  Runoff.int.logmean = log(3.3) #mode input, precip + runoff + ground water (Mohammed 2011)
  
  Runoff.int.logsd = log(1.2)
  
  Runoff.pre ~ dgamma(Runoff.pre.shp, Runoff.pre.rate)
  Runoff.pre.shp = 30
  Runoff.pre.rate = 2
  
  #####LST time series####
  #normal distribution, in degrees C
  for (i in 2:t){
    LST.k[i] = 273.15 + LST[i]
    
    AT[i] = LST[i] + 10 #air temperature of the warm season is set at a constant offset
    
    LST[i] = LST[i - 1] + LST.cps[i]
    
    LST.cps[i] ~ dnorm(LST.cps[i - 1] * LST.cps.ac, LST.pre)
    
  }
  LST.cps[1] ~ dnorm(0, LST.pre) #allowed some variation
  
  LST.cps.ac ~ dunif(0.001, 0.8)
  
  AT[1] = LST[1] + 10
  # initiate the series with an reasonable prior
  LST.k[1] = 273.15 + LST[1]
  LST[1] ~ dnorm(LST.int, LST.pre) #allowed some variation
  
  #an uninformative initial value: 20+-5 degrees C with a warm season bias, Steenburgh et al 2000
  LST.int ~ dnorm(20, 1/5^2) T(10, 35)  
  
  LST.pre ~ dgamma(LST.pre.shp, LST.pre.rate) # ~0.1 degrees error/10 years
  LST.pre.shp = 200
  LST.pre.rate = 2

  #####starting values####
  
  nsws ~ dnorm(5.8, 1/0.2^2) #wind speed data from Steenburgh, 2000
  
  #relative humidity is set to center around 0.4 with a warm season bias (dry summer)
  rh ~ dbeta(rh.a, rh.b)
  
  rh.b ~ dlnorm(log(100),1/log(2)^2)
  
  rh.a ~ dlnorm(log(70),1/log(2)^2)
  
  #lake volume, use bathymetric table
  LV[1] = interp.lin(L.level[1], GSL.level, GSL.volume)
  
  #lake area, use bathymetric table
  LA[1] = interp.lin(L.level[1], GSL.level, GSL.area)
  
  #lake starting level: 1280 +-1 meters
  L.level[1] ~ dnorm(1280,1/1^2) T(1275,1286)#m
  
  # #NaCl, MgCl2, Na2SO4, KCl
  # salt.molar.mass = (58.44*0.7591 + 95.211*0.1092)
  # #NaCl, MgCl2, Na2SO4, KCl
  # salt.molarg = (0.7591/58.44 + 0.1092/95.211 + 0.0952/142.04 + 3.16/100/74.55)
  
}

#question is with dE, use eqn 3 in Gibbson 2016, need RH, vapor d, will link it to dL
#also consider salt water correction using aw: thermodynamic activity of salt water, which increases RH

#Kinetic isotopic separation (epsilon.k), double check
# epsilon.k.18O = 14.2*(1-RH/aw)
# epsilon.k.2H = 12.5*(1-RH/aw)
# 
# #Equilibrium isotopic separations (epsilon)
# epsilon.18O = 1000(Alpha.18O-1)
# epsilon.2H = 1000(Alpha.2H-1)
# 
# #equilibrium fractionations
# Alpha.18O = EXP((-7.685+6.7123*(10^3/LST.k[i])-1.6664*(10^6/LST.k[i]^2)+0.35041*(10^9/LST.k[i]^3))/1000)
# Alpha.2H =EXP((1158.8*(LST.k[i]^3/10^9)-1620.1*(LST.k[i]^2/10^6)+794.84*(LST.k[i]/10^3)-161.04+2.9992*(10^9/(LST.k[i]^3))/1000)
