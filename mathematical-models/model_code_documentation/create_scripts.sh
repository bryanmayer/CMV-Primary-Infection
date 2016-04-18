#!/bin/bash
Rscript -e 'knitr::knit(input = "CMV_Models.Rmd", output = "../model_functions/CMV_Models.R", tangle = T)'
Rscript -e 'knitr::knit(input = "optim_functions.Rmd", output = "../model_functions/optim_functions.R", tangle = T)'
Rscript -e 'knitr::knit(input = "processing_functions.Rmd", output = "../model_functions/processing_functions.R", tangle = T)'
Rscript -e 'knitr::knit(input = "parameter_boundary.Rmd", output = "../parameters/parameter_boundary.R", tangle = T)'