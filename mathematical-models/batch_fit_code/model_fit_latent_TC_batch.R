#example install.packages("dplyr", lib = "libraryFolder")

#The initial sets are independently generated (and saved) in each run
#This fits beta and start_day across the entire data (base only fits the first phase)

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

if(!exists("p_in")) {p_set = 1:length(patient_set; p_in = 0)
} else if(p_in == 1) p_set = 1:7 else if(p_in == 2) p_set = 8:14

lower_start_day = log10(c(0.1, 1, seq(5, 45, 5)))
upper_start_day = log10(c(1, seq(5, 50, 5)))

model_fit_base = ldply(patient_set[p_set], function(id){
  
  subData = subset(CMVPrimaryEpisodes, PatientID2 == id)
  peakday = find_peak(subData, window = 1)  
  first_pos = subset(subData, days2 == 0)$count
  
  total_draws = 10
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
  
  ldply(1:length(upper_start_day), function(s_day){
    ldply(1:total_draws, function(draw){
      
      fitparms = starting_sets[[draw]]
      
      #pick between set bounds
      fitparms["start_day"] = runif(1, lower_start_day[s_day], upper_start_day[s_day])
      lb = boundary_set(fitparms, "lower")
      lb[1] = lower_start_day[s_day]
      ub = boundary_set(fitparms, "upper", initV_upper = first_pos)
      ub[1] = upper_start_day[s_day]
      
      outdata = ldply(model_fit_set, function(model){
        
        AUCData = NULL
        
        fit = nloptr_call(fitparms, fixparms, subData, init_model, CMV_model = get(model), AUCData = AUCData,
                          lowerBounds = lb, 
                          upperBounds = ub)
        
        fit$model = model
        fit$PatientID2 = id
        fit$draw = draw
        fit$initR0 = with(as.list(c(fixparms, 10^fitparms)), beta * K * p /(delta * c * (1 + mu/alpha)))
        fit$R0 = with(fit, beta * K * p /(delta * c * (1 + mu/alpha)))
        
        return(fit[, c("PatientID2", "model", "draw", "initR0", "R0", all_parmnames, "mse", "conv")])
        
      }, .parallel = F)
      
      outdata
      }, .parallel = T)
    })
  }, .parallel = F)

file_out = paste("model_latent_fit_profileSD_TC", p_in, ".csv", sep = "")

write.csv(model_fit_base, file_out, row.names = F)
Sys.time() - starttime
