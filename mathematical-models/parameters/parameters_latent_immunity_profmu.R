immunity_parms = function(peakday = NA, par = NA){
  list(
    immune_fixparms = c(
      c = 2,
      K = 1e7 * 40,  #Dawes
      mu = 1/(4.5), #Dawes  is 1/4.5
      p = 1600, #from PNAS temperature paper
      alpha = 1,
      initI = 0,
      initV = 0
    ),
    immune_fitparms = c(
      beta = -10,
      start_day = log10(7),
      delta = 0
      )
    ) 
}




