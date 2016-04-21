immunity_parms = function(peakday = NA, par = NA){
  list(
    immune_fixparms = c(
      c = 2,
      K = 1e7 * 40,  #Dawes
      mu = 1/(4.5), #Dawes  is 1/4.5
      delta = 0.77, #from emery paper
      p = 1600, #from PNAS temperature paper
      alpha = 1, #one day latency
      initI = 0,
      initV = 0
    ),
    immune_fitparms = c(
      start_day = log10(2),
      beta = -10
      )
    ) 
}




