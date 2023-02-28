####helper functions#####
MCMC.CI.bound <- function (MCMC.res, CI){
  require(KernSmooth)
  require(bayestestR)
  dim.MCMC <- dim(MCMC.res)
  #typically, the first element is # of iterations
  #the second element is the time series
  map.res <- rep(0, dim.MCMC[2])
  hdi.high <- rep(0, dim.MCMC[2])
  hdi.low <- rep(0, dim.MCMC[2])
  
  for(i in 1:dim.MCMC[2]){
    map.res[i] <- map_estimate(MCMC.res[,i], method = "KernSmooth")
    hdi <- hdi(MCMC.res[,i], ci = CI)
    hdi.low[i] <- hdi$CI_low
    hdi.high[i] <- hdi$CI_high
  }
  
  return (list(map.res, hdi.low, hdi.high, CI))
}
