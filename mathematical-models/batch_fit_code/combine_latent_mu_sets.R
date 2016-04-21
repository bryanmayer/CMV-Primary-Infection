library(plyr)
library(readr)

filedir = ""

write_csv(ldply(1:7, function(i) read_csv(paste(filedir, "latent_fit_profileSD_profmu_", i,".csv", sep = ""))),  
         paste(filedir, "latentTC_fit_profmu.csv", sep = ""))