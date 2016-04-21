#example install.packages("dplyr", lib = "libraryFolder")

#this takes in the mu_in parameter, which is an index for which mu value to fix
#if mu_in is not passed in, it should run all of them

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
library(readr, lib = "libraryFolder")
library(nloptr, lib = "libraryFolder")
library(deSolve, lib = "libraryFolder")
library(pracma, lib = "libraryFolder")
library(lhs, lib = "libraryFolder")
library(doParallel)
registerDoParallel(cores = 8)

load("CMVprimaryEpisode.RData")
source("optim_functions.R")
source("CMV_Models.R")
source("parameters_latent_immunity_profmu.R")
source("parameter_boundary.R")

starttime <- Sys.time()

#for LHS
total_draws = 20

mu_set = c(1e-4, 1e-3, 1/365, 1/30, 1/4.5, 1, 2)
if(!exists("mu_in")) mu_in = 1:length(mu_set)
print(mu_set[mu_in])

########## set up starting_sets #############
model_fit_set = "CMVModel_latent_linear"
init_model = c(S = NA, I0 = 1, I = 0, V = NA)

initial_immunity = immunity_parms(NA)
fixparms = initial_immunity$immune_fixparms
fitparms = initial_immunity$immune_fitparms
all_parmnames = c(names(fixparms), names(fitparms))

totalFit = length(fitparms)
lhs = randomLHS(total_draws, totalFit)

#manually setting these
beta_min = 10^-12
delta_min = 10^-3.5
beta_max = 10^-9
delta_max = 100

beta_i = which(names(fitparms) == "beta")
delta_i = which(names(fitparms) == "delta")
start_i = which("start_day" == names(fitparms))

lower = c(0, 0, 0)
lower[beta_i] = beta_min
lower[delta_i] = delta_min
lower[start_i] = 0.1

upper = c(0, 0, 0)
upper[beta_i] = beta_max
upper[delta_i] = delta_max
upper[start_i] = 50

starting_sets = llply(1:total_draws, function(i){
  tempdraws = (upper - lower) * lhs[i, ] + lower
  names(tempdraws) = names(fitparms)
  
  if(T){
    tempR0 = with(as.list(c(tempdraws, fixparms)), beta * K * p /(delta * c * (1 + mu/alpha)))
    if(tempR0 < 0.8) return(NULL)  
    if(tempR0 > 50) return(NULL)  
  }
  
  log10(tempdraws)
})

starting_sets <- compact(starting_sets) #compact removes NULLs (remove nothing if R0check = F)
print(length(starting_sets))

lower_start_day = log10(c(0.1, 1, seq(5, 45, 5)))
upper_start_day = log10(c(1, seq(5, 50, 5)))

model_fit_base = ldply(mu_set[mu_in], function(mu_fix){
  temp_mu =  mu_fix
  fixparms["mu"] = temp_mu
  
  ldply(1:length(starting_sets), function(draw){
    ldply(1:length(upper_start_day), function(s_day){

      fitparms = starting_sets[[draw]]
      
      #pick between set bounds
      fitparms["start_day"] = runif(1, lower_start_day[s_day], upper_start_day[s_day])
      lb = boundary_set(fitparms, "lower")
      lb[start_i] = lower_start_day[s_day]
      ub = boundary_set(fitparms, "upper")
      ub[start_i] = upper_start_day[s_day]
      
      ldply(unique(CMVPrimaryEpisodes$PatientID2), function(id){    
        subData = subset(CMVPrimaryEpisodes, PatientID2 == id)

        ldply(model_fit_set, function(model){

          fit = nloptr_call(fitparms, fixparms, subset(subData), init_model, CMV_model = get(model), 
                          lowerBounds = lb, 
                          upperBounds = ub)
        
          fit$model = model
          fit$PatientID2 = id
          fit$draw = draw
          fit$prop = fit$delta/temp_mu
          fit$initR0 = with(as.list(c(fixparms, 10^fitparms)), beta * K * p /(delta *  c * (1 + mu/alpha)))
          fit$R0 = with(fit, beta * K * p /(delta *  c * (1 + mu/alpha)))
          
          return(fit[, c("PatientID2", "model", "draw", "initR0", "prop", "R0", all_parmnames, "mse", "conv")])
        
          }, .parallel = F)
        }, .parallel = F) #Patient level
      }) #start_day
    }, .parallel = T) #initial values
  }, .parallel = F)

if(length(mu_in) > 1) mu_in = "all"

write_csv(model_fit_base, paste("latent_fit_profileSD_profmu_", mu_in, ".csv", sep = ""))
Sys.time() - starttime
