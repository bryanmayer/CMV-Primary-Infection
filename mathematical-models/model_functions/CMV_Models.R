## ------------------------------------------------------------------------
### This script was created by model_code_documentation/create_scripts.sh
### The following models are described in model_code_documentation/CMV_models.Rmd with corresponding pdf.

CMVModel_latent_linear = function(t, x, parms, AUCData = NULL){
  with(as.list(c(parms, x)), {
    
    if(as.logical(is.na(parms["lambda"]))) lambda = K * mu
    
    dS <- lambda - mu * S - beta * S * V
    dI0 <- beta * S * V - alpha * I0 - mu * I0
    dI <- alpha * I0  -  delta * I
    dV <- p * I - c * V 
    
    res <- c(dS, dI0, dI, dV)
    list(res)
  })
}

## ------------------------------------------------------------------------
CMVModel_latent_immunity_CTL = function(t, x, parms, AUCData = NULL){
  with(as.list(c(parms, x)), {
    if(as.logical(is.na(parms["lambda"]))) lambda = K * mu
    
    dS <- lambda - mu * S - beta * S * V
    dI0 <- beta * S * V - alpha * I0 - mu * I0
    dI <- alpha * I0  -  delta * I - k * I * Tcell
    dV <- p * I - c * V 
    dT <- theta * (I/(KI+I)) - death * Tcell
    
    res <- c(dS, dI0, dI, dV, dT)
    list(res)
  })
}

## ------------------------------------------------------------------------
CMVModel_latent_immunity_V = function(t, x, parms, AUCData = NULL){
  with(as.list(c(parms, x)), {
    if(as.logical(is.na(parms["lambda"]))) lambda = K * mu
    
    dS <- lambda - mu * S - beta * S * V
    dI0 <- beta * S * V - alpha * I0 - mu * I0
    dI <- alpha * I0  -  delta * I    
    dV <- p * I - c * V - k * V * Tcell
    dT <- theta * (V/(KT+V)) - death * Tcell
    
    res <- c(dS, dI0, dI, dV, dT)
    list(res)
  })
}

## ------------------------------------------------------------------------
find_peak = function(data, window = 1){
  #window is the log range around the peak, 0 would be peak
  
  peak = max(data$count)
  
  peak_range = which((peak - data$count) <= window) #find counts within window of measured peak
  peak_day = data$days2[min(peak_range)]
  return(peak_day)
}

