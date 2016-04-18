## ------------------------------------------------------------------------
### This script was created by model_code_documentation/create_scripts.sh
### The following functions are described in model_code_documentation/optim_functions.Rmd with corresponding pdf.

nloptr_call = function(fitted_parms, fixed_parms, fit_data, model_initial, CMV_model = CMVModel_latent_linear,
                       lowerBounds = NULL, upperBounds = NULL,  R0_test = T, AUCData = NULL, 
                       justParms = F, maxeval = 1000){
  #calls the base model by default
  
  #if either set of constraints are missing, runs unconstrained
  if(is.null(upperBounds) | is.null(lowerBounds)){
    fit = nloptr(
      x0 = fitted_parms,
      eval_f = optimize_fun,
      opts = list("algorithm"="NLOPT_LN_NELDERMEAD",
                  "xtol_rel"=1.0e-8, "maxeval" = maxeval),
      model_fun = CMV_model, inData = fit_data,
      inParms = fixed_parms,  init = model_initial, AUCData = AUCData,
      parmNames = names(fitted_parms), R0_test = R0_test)
  } 
  else{
    fit = nloptr(
      x0 = fitted_parms,
      eval_f = optimize_fun,
      lb = unname(lowerBounds),
      ub = unname(upperBounds),
      opts = list("algorithm"="NLOPT_LN_NELDERMEAD",
                  "xtol_rel"=1.0e-8, "maxeval" = maxeval),
      model_fun = CMV_model, inData = fit_data,
      inParms = fixed_parms,  init = model_initial, AUCData = AUCData,
      parmNames = names(fitted_parms), R0_test = R0_test)
  }
  parmsfit = 10^fit$sol
  names(parmsfit) = names(fitted_parms)
  
  outparms = c(parmsfit, fixed_parms)
  
  fit_output = as.data.frame(t(outparms))
  
  #default to include additional fit information
  if(!justParms){
    fit_output$mse = fit$obj
    fit_output$conv = fit$status
  }
  
  fit_output
}


## ------------------------------------------------------------------------
optimize_fun = function(fitParms, inParms, inData, init, model_fun = CMVModel_latent_linear, parmNames = NULL, 
                        AUCData = NULL, R0_test = F){
  inData = arrange(inData, days2)
  
  fitParms = 10 ^ fitParms
  if(length(parmNames) != 0) names(fitParms) = parmNames
  
  id_var = which(names(inParms) == "id")

  parms = c(fitParms, inParms)
  
  init["S"] = as.numeric(unname(parms["K"]))
  init["I"] = as.numeric(unname(parms["initI"]))
  init["V"] = as.numeric(unname(parms["initV"]))
  
  start_time = -round(unname(parms["start_day"]), 1)
  times = seq(start_time, max(inData$days2) + 50, 0.1)
  
  # this throws out parameters that have R0 > 100, this is the old R0 and is missing (1+mu/alpha) in denominator
  # so the actual constraint is > 82 since mu = 1/4.5 and alpha =1, that term =1.22
  # so when this R0 equals 100, the actual R0 is 100/1.22 - 81.97
  
  if(R0_test) if(with(as.list(c(fitParms, inParms)), (beta * K * p /(delta * c))) > 100) return(Inf)

  model_out = as.data.frame(lsoda(init, times, model_fun, parms, AUCData = AUCData))

  if (sum(model_out$V) == "NaN" | dim(model_out)[1] != length(times)) return(Inf)
  
  #if there is early oscillation
  if(model_out$time[which.max(model_out$V)] < 5) return(Inf)
  
  mse = cost_function(model_out, inData)
  if (mse == 0) return(Inf) # browser()
  return(mse)
}

## ------------------------------------------------------------------------
optimize_fun_test = function(fitOutput, inData, AUCData = NULL){
  
  modelIn = fitOutput$model
  
  if(modelIn == "CMVModel_latent_linear" | is.null(fitOutput$theta)) {
    init = with(fitOutput, c(S = K, I0 = 0, I = initI, V = initV))
    }else init = with(fitOutput, c(S = K, I0 = 0, I = initI, V = initV, Tcell = 0))

  times = seq(-round(fitOutput$start_day, 1), max(inData$days2) + 50, 0.1)
  
  out = as.data.frame(lsoda(init, times, get(modelIn), 
                            unlist(select(fitOutput, -model, -PatientID)), AUCData = AUCData))
  
  inData = arrange(inData, days2)
  
  mse = cost_function(out, inData)
  #print(mse)
  if (mse == 0) return(Inf) # browser()
  return(mse)
}


## ------------------------------------------------------------------------
cost_function = function(model, data, debug = F){
  if(debug) browser()
  data = arrange(data, days2)
  data = data %>% dplyr::mutate(days_model = days2) 
    
  sample_model = arrange(model[which(round(model$time, 1) %in% round(data$days_model, 1)), ], time)
  mse = sum((log10(sample_model$V) - data$count2)^2, na.rm = T) 
  if (mse == "NaN") browser()
  #print(mse)
  return(mse)
  }

