# Mathematical models of primary CMV

Besides containing the analysis files, this documentation directories contains all of the function, models, batch scripts, and model optimization results from the manuscript. The models presented in the article are all in the model_functions folder.  The main immune model presented is the CTL model. 

The results in the manuscript are in the analysis .Rmd and .pdf files. These R scripts are analysis of optimization results from fitting models to the infant episode data (done using the code in batch_fit_code/):
  * target_cell_model_analysis contains all of the analysis of the target cell model optimizations. R0 was calculated from the expansion phase analysis (the first section).
  * ctl_immune_model_analysis contains the analysis of the cytolytic immune model that was presented in the manuscript.
  * viral_immune_model_analysis contains an additional model of viral mediated immunity. This was not included in the manuscript.

# How to use this repo to simulate a model

##1) Use R 

Make sure appropriate packages are installed: **plyr**, **dplyr**, **deSolve**, **ggplot2**, **gridExtra**, **scales**, RColorBrewer, doParallel, pracma, reshape2, knitr.  Some of these libraries were utilized in very early iterations of the modeling and may not be necessary anymore, I bolded the ones I know are absolutely necessary.

##2) Source in necessary functions

The two most important scripts to simulating the models are processing_functions.R and CMV_Models.R from the "model_functions/" folder. optim_functions.R and scripts in the parameters/ folder are used for optimization.

```
#make sure to set the appropriate directory (here I'm in mathematical-models)
source("model_function/processing_functions.R")
source("model_function/CMV_Models.R")
```
processing_functions.R will give you error if you are missing packages.

##3) Assign your parameter values
Here is an example
```
parms = data.frame(
	model = "CMVModel_latent_immunity_CTL", #look in model_code_documentation/CMV_Models.pdf for a list
	PatientID2 = "test", #this is needed because of how I coded it
	initI = 0, #this is I in the equations
	initV = 0,
	K = 4e8, #this is S0
	mu = 1/4.5,
	delta = 0.77,
	p = 1600,
	c = 2,
	alpha = 1,
	beta = 4.20e-12,
	start_day = 5,
	k = 0.01,
	theta = 0.42,
	KI = 82.66,
	death = 7.83e-3 #this is gamma
	,stringsAsFactors = F)
```

I0 = 1 by default in my functions (this is the L population in the manuscript equations).

## 4A) Simple simulation plot using quickplot()

```
quickplot(parms)
```

quickplot will always plot in the time range (-start_day, 800).
If you so happen to have data of the form:

```
test_data = data.frame(days2 = 0:500, count = rnorm(501, 6))
```

quickplot will also layer that data onto the plot.

```
quickplot(parms, test_data)
```

But let's use some read data

```
#check directory if there an error
load("../data/CMVprimaryEpisode.RData")
infant_A_episode = subset(CMVPrimaryEpisodes, PatientID2 == "A")

quickplot(parms, infant_A_episode)
```

##4B) Get simulation data using make_plot_data() or make_highres_sim()

make_plot_data() will return output of the simulation.  Like quickplot, you can add a data set for it to plot over (it will be stored in the output with variable model = "data").  If data is provided, the simulation will only run over the length of that data (800 is default).  It can conduct multiple simulations if there are multiple models in the parm file.  The model simulations are stacked (tidy data for plotting) so they can be easily faceted.  make_plot_data does not match the parameter set's PatientID2 variable to the data variable.

count = log10(V)
```
model_sim = make_plot_data(parms)

ggplot(data = model_sim, aes(x = days_model, y = count)) + geom_point()

model_sim_withdata = make_plot_data(parms, infant_A_episode)

ggplot(data = model_sim_withdata, aes(x = days_model, y = count, colour = model)) + 
	geom_point()
```

There is also make_highres_sim() that will return a high resolution output of the simulation (0.01 time increments) but it only goes to 350 time currently.  It has a nice feature that it'll cycle through multiple patients and models (in the parm data.frame) and stack the results, if the parms data is setup accordingly.

```
model_sim = make_highres_sim(parms)

ggplot(data = model_sim, aes(x = days_model, y = count)) + geom_point()
```

if you set `first_time = T` in either function, it will return simulation output starting at -start_day instead of 0.

```
parms$start_day = 50

model_sim = make_plot_data(parms, first_time = T)

ggplot(data = model_sim, aes(x = days_model, y = count)) + geom_point()
```

##4C) Look at actual simulation results. Here we'll just look at episodes from infants A-D for simplicity using their best fits from the immune model optimization. The make_subject_plots takes this information, matches the PatientID2 variable and returns a data.frame with model simulation and the raw data (model = "data") by PatientID2. Plotting this is a little clunky right now, and the make_subject_plots will break unless return_data = T.
```
load("model_results/ctl_immune_model_fits.RData")
subset_episodes = subset(CMVPrimaryEpisodes, PatientID2 %in% LETTERS[1:4])
subset_fits = subset(immune_model_fits, PatientID2 %in% LETTERS[1:4])

simulation_data = make_subject_plots(fitOutput = subset_fits, allCMVData = subset_episodes, return_data = T)

ggplot(data = subset(simulation_data, model != 'data'),
  aes(x = days_model, y = count)) +
  geom_point(data = subset(simulation_data, model == 'data')) +
  geom_line(alpha = 0.7) +
  facet_wrap(~PatientID2)

```


##4D) Manually simulate the model using an ode solver

Note: you still need to source in CMV_Models.R:

```
init = c(S = parms$K, I0 = 1, I = 0, V = 0, Tcell = 0)
times = seq(-round(parms$start_day, 1), 500, 0.1)

model_sim = as.data.frame(lsoda(init, times, CMVModel_latent_immunity_CTL, unlist(select(parms, -model, -PatientID2))))

#you can use get(parms$model) instead of CMVModel_latent_immunity_CTL

with(model_sim, plot(time, log10(V)))

ggplot(data = model_sim, aes(x = time, y = log10(V))) + geom_point()
```

remember initV = 0 so log10(V)[1] does equal -infinity, but the plots are smart enough to clip.

