immunity_parms = function(peakday = NA, par = NA){
  list(
    immune_fixparms = c(
      c = 2,
      K = 1e7 * 40,  #Dawes
      mu = 1/(4.5), #Dawes  is 1/4.5
      delta = 0.77, #from emery paper
      p = 1600, #from PNAS temperature paper
      initI = 0,
      initV = 0,
      alpha = 1,
      k = 0.01,
      beta = NA,
      start_day = NA
    ),
    immune_fitparms = c(
      KI = 5,
      death = -3,
      theta = log10(6.76)
    )
    ) 
}




