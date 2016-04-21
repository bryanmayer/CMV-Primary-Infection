#!/bin/bash

for x in {1,2}; do
sbatch --cpus-per-task=8 --time=0-8 --wrap="R --no-save --no-restore '--args p_in=$x' < model_fit_latent_TC_batch.R"
done