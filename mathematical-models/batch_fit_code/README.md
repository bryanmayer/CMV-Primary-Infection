The scripts in this folder are almost the exact same scripts I used to run batch jobs on the scicomp server at the Fred Hutch. I made the following changes:

* Changed PatientID to PatientID2 to reflect recoded infant IDs so the jobs could be rerun
* Fixed a typo in model_fit_latent_TC_batch.R where it wouldn't work without p_in being passed in.
* Removed my local absolute directory from combine_latent_mu_sets.R

Unless specified, all scripts were executed using a direct sbatch command in bash on the servers. The command was similar to that seen in both .sh files in this directory. Besides the packages listed in the scripts, the following files are also required for any of the scripts to run (given as directories in this repo). The three R files are described in detail in the mathematical-models/model_code_documentation/

* data/CMVprimaryEpisode.RData
* mathematical-models/model_functions/optim_functions.R
* mathematical-models/model_functions/CMV_Models.R
* mathematical-models/parameters/parameter_boundary.R

Specific file requirements to run an optimization are given below with each specific script. Each optimization has a corresponding parameter file list is used to inform the optimizer which parameters are fixed and which are to be fitted. The following optimizations were performed:

1) model_fit_latent_base_batch.R - Fits the expansion phase of the model to estimate beta, start_time, and R0.
  * uses mathematical-models/parameters/parameters_latent_immunity_base.R
  * generates the data that creates initial_latent_values.RData that is used later
2) model_fit_latent_TC_batch.R - Fits the entirety of the episode (for beta and start_time) without an immune response and with fixed delta and mu.
  * uses mathematical-models/parameters/parameters_latent_immunity_base.R
  * was parallelized on the server by passing in an iterator that splits the analysis up (p_in) with a bash script - sh_latent_fit.sh
3) model_fit_profileTC_batch.R - Fits the entirety of the episode without an immune response (for beta, delta, and start_time) for various levels of mu
  * uses mathematical-models/parameters/parameters_latent_immunity_profmu.R
  * was parallelized on the server by passing in an iterator for mu (mu_in) with a bash script - mu_latent_fit.sh
  * combine_latent_mu_sets.R is a script that can be used to combine the output files
4) latent_fit_immunityCTL_batch.R - Fits the entirety of the episode with a cytolytic immune response using the fits from 1) for beta and start_time
  * uses mathematical-models/parameters/parameters_latent_immunityCTL.R
  * uses results from expansion model fits in mathematical-models/model_results/initial_latent_values.RData
5) latent_fit_immunity_batch.R - Fits the entirety of the episode with a viral mediated immune response using the fits from 1) for beta and start_time
  * uses mathematical-models/parameters/parameters_latent_immunity.R
  * uses results from expansion model fits in mathematical-models/model_results/initial_latent_values.RData