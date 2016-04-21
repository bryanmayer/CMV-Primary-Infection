#example install.packages("dplyr", lib = "libraryFolder")

#this takes fixed initial values from base_model run and only fits the immunity parameters


library(plyr)
library(dplyr, lib = "libraryFolder")
library(readr, lib = "libraryFolder")
library(nloptr, lib = "libraryFolder")
library(deSolve, lib = "libraryFolder")
library(pracma, lib = "libraryFolder")
library(lhs, lib = "libraryFolder")
library(doParallel)
registerDoParallel(cores = 8)


load("CMVprimaryEpisode.RData")
load("initial_latent_values.RData")
source("optim_functions.R")
source("CMV_Models.R")
source("parameters_latent_immunityCTL.R")
source("parameter_boundary.R")

starttime <- Sys.time()

########## set up starting_sets #############
initial_immunity = immunity_parms(NA)
fixparms = initial_immunity$immune_fixparms

#from the prop input, assign delta#


fitparms = initial_immunity$immune_fitparms
init_model = c(S = NA, I0 = 1, I = NA, V = NA, Tcell = 0)
all_parmnames = c(names(fixparms), names(fitparms))


model_fit_set = c("CMVModel_latent_immunity_CTL")


#sample initial values
total_draws = 50
totalFit = length(fitparms)

lower = boundary_set(fitparms, "lower")
upper = boundary_set(fitparms, "upper")

lhs = randomLHS(total_draws, totalFit)
starting_sets = llply(1:total_draws, function(i){
  tempdraws = (upper - lower) * lhs[i, ] + lower
  names(tempdraws) = names(fitparms)
  tempdraws
})

patient_set = unique(CMVPrimaryEpisodes$PatientID2)

model_fit_immunity = ldply(patient_set, function(pid){

  subData = subset(CMVPrimaryEpisodes, PatientID2 == pid)
  initial_values = subset(base_model_fits, PatientID2 == pid)
  fixparms["beta"] = initial_values$beta
  fixparms["start_day"] = initial_values$start_day

  ldply(1:total_draws, function(draw){
    
    fitparms = starting_sets[[draw]]
    
    ldply(model_fit_set, function(model){
      
      fit = nloptr_call(fitparms, fixparms, subData, init_model, CMV_model = get(model),
                        lowerBounds = boundary_set(fitparms, "lower"), 
                        upperBounds = boundary_set(fitparms, "upper"))
      
      fit$model = model
      fit$PatientID2 = pid
      fit$R0 = initial_values$R0
      
      return(fit[, c("PatientID2", "model", "R0", all_parmnames, "mse", "conv")])
      
    }, .parallel = F)
    
  })
  
}, .parallel = T)

file_out = paste("latent_fit_ctl.csv", sep = "")

write.csv(model_fit_immunity, file_out, row.names = F)
Sys.time() - starttime
