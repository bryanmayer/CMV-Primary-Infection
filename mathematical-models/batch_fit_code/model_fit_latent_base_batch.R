#example install.packages("dplyr", lib = "libraryFolder")

#The initial sets are independently generated (and saved) in each run

args<-(commandArgs(TRUE));
if(length(args)==0){
  print("No arguments supplied.")
  
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
  print(args)
}

library(plyr)
library(dplyr, lib = "libraryFolder")
library(nloptr, lib = "libraryFolder")
library(deSolve, lib = "libraryFolder")
library(pracma, lib = "libraryFolder")
library(lhs, lib = "libraryFolder")
library(doParallel)
registerDoParallel(cores = 8)

load("CMVprimaryEpisode.RData")
source("optim_functions.R")
source("CMV_Models.R")
source("parameters_latent_immunity_base.R")
source("parameter_boundary.R")


#baseparmsIn = read.csv("model_fit_growth_informed.csv") 
starttime <- Sys.time()


########## set up starting_sets #############
initial_immunity = immunity_parms(NA)
fixparms = initial_immunity$immune_fixparms

fitparms = initial_immunity$immune_fitparms
init_model = c(S = NA, I0 = 1, I = 0, V = NA)
all_parmnames = c(names(fixparms), names(fitparms))

model_fit_set = "CMVModel_latent_linear"

patient_set = unique(CMVPrimaryEpisodes$PatientID2)

model_fit_base = ldply(unique(CMVPrimaryEpisodes$PatientID2), function(id){
  
  subData = subset(CMVPrimaryEpisodes, PatientID2 == id)
  peakday = find_peak(subData, window = 1)  
  first_pos = subset(subData, days2 == 0)$count
  
  total_draws = 50
  lower = boundary_set(fitparms, "lower")
  upper = boundary_set(fitparms, "upper", initV_upper = first_pos)
  
  #this will break unless everything else is fixed
  beta_ind = which(names(fitparms) == "beta")
  # R0 >= 1
  lower[beta_ind] = with(as.list(c(fixparms)), 
                         log10((delta * c * (1 + mu/alpha)) / (K * p))
  )
  #R0 <= 50, make sure it doesnt go outside contraints, add a little noise if it's on the edge
  upper[beta_ind] = min(with(as.list(c(fixparms)), 
                             log10(50 * (delta * c * (1 + mu/alpha)) / (K * p))), -9 - runif(1)
  )
    
  starting_sets = draw_initial(fitparms, lower, upper, total_draws,
                               save.out = F)
  
  ldply(1:total_draws, function(draw){
    
    fitparms = starting_sets[[draw]]
    
    outdata = ldply(model_fit_set, function(model){
      
      AUCData = NULL
      
      fit = nloptr_call(fitparms, fixparms, subset(subData, days2 <= peakday), init_model, CMV_model = get(model), AUCData = AUCData,
                        lowerBounds = boundary_set(fitparms, "lower"), 
                        upperBounds = boundary_set(fitparms, "upper", initV_upper = first_pos))
      
      fit$model = model
      fit$PatientID2 = id
      fit$draw = draw
      fit$initR0 = with(as.list(c(fixparms, 10^fitparms)), beta * K * p /(delta * c * (1 + mu/alpha)))
      fit$R0 = with(fit, beta * K * p /(delta * c * (1 + mu/alpha)))
      
      return(fit[, c("PatientID2", "model", "draw", "initR0", "R0", all_parmnames, "mse", "conv")])
      
    }, .parallel = T)
    
    outdata
  })
  
}, .parallel = F)

file_out = paste("latent_initial_fit.csv", sep = "")

write.csv(model_fit_base, file_out, row.names = F)
Sys.time() - starttime
