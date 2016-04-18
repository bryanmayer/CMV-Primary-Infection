## ------------------------------------------------------------------------
#### EDIT THIS IN THE DOCUMENTATION FILE (.Rmd) ONLY ####


boundary_set = function(parms, boundary = NULL, t_set = NA, par_set = NA, initV_upper = 2){
  if(length(parms) == 0 | is.null(boundary)) return(NULL)
  
  #remember these are all log10ed
  
  if(boundary %in% c("l", "lower", "L", "Lower")){
    set_values = c(
      beta = -14,
      delta = -4,
      theta = -3,
      KI = 0,
      KT = 2,
      death = -5,
      start_day = log10(0.1)
    )
    return_list = sapply(1:length(parms), function(i) which(names(set_values) == names(parms)[i]))
    
    return(unname(set_values[return_list]))
  }
  
  if(boundary %in% c("u", "upper", "U", "Upper")){
    set_values = c(
      beta = -6.5,
      delta = 3,
      theta = 2,
      KI = 8,
      KT = 10,
      death = 0,
      start_day = log10(25)
    )
    return_list = sapply(1:length(parms), function(i) which(names(set_values) == names(parms)[i]))
    
    return(unname(set_values[return_list]))
  }
  
  print("Incorrect boundary label given (use 'lower' or 'upper')")
  return(NULL)
}

## ------------------------------------------------------------------------

#this returns a list of a set of initial value draws for given set of parameters
#uses LHS package
#R0check is only for when beta and delta are fit together
draw_initial = function(fitparms, lower, upper,
                        total_draws = 50, save.out = T, R0check = F, file_name = NULL){
  totalFit = length(fitparms)
  
  if(totalFit != length(lower) | totalFit != length(upper)) {
    print(paste("mismatch, 1=lower, 2=upper:", which(!c(length(lower), length(upper)) %in% totalFit)))
    return(NULL)
  }
  
  
  lhs = randomLHS(total_draws, totalFit)
  
  starting_sets = llply(1:total_draws, function(i){
    tempdraws = (upper - lower) * lhs[i, ] + lower
    names(tempdraws) = names(fitparms)
    
    if(R0check){
      tempR0 = with(as.list(10^tempdraws), beta * 4e8 * 60 /(delta * 2))
      if(tempR0 < 0.5) return(NULL)  
    }

    tempdraws
  })
  
  starting_sets <- compact(starting_sets) #compact removes NULLs (remove nothing if R0check = F)
  total_kept = length(starting_sets) #should be same as total_draw when R0check = F
  
  if(save.out){
    starting_sets_save = ldply(1:total_kept, function(i) starting_sets[[i]])
    if(is.null(file_name)) file_name = Sys.Date()
    write.csv(starting_sets_save,
            paste(file_name,".csv", sep = ""),
            row.names = F)
  }
  
  return(starting_sets)

}

## ------------------------------------------------------------------------

#this tacks on parameters names that werent used in the model for rbinding in make_plot_data and main ldply
addnames_fit = function(output, parmnames){
  addnames = which(!parmnames %in% names(output))
  
  if(length(addnames) == 0) return(output)
  
  temp_output = cbind(output, t(rep(NA, length(addnames))))
  names(temp_output) = c(names(output), parmnames[addnames])
  
  return(temp_output)
  
}

