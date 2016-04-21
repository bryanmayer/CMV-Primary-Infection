#!/bin/bash

for x in {1,2,3,4,5,6,7}; do
sbatch --cpus-per-task=8 --time=0-8 --wrap="R --no-save --no-restore '--args mu_in=$x' < latent_fit_profileTC_batch.R"
done