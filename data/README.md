This folder contains the raw data used to create all of the analysis:

1) complete_episode_data.csv - The time series shedding data for all of the infants

2) final_episode_data.csv - The subsetted time series for the infants included in the main analysis with count2 variable, coded NA when count = 0.

3) CMVprimaryEpisode.Rdata - The R save space version of final_episode_data.csv. Loads a data.frame called CMVPrimaryEpisodes

4) primary_infection_notes.csv - Notes on visual inspection of primary infection data used to develop criteria for onset of oral shedding episode and for exclusion.

5) age_data.csv - the age data for each of the infants. This variable was calculated manually in the analysis presented in data-analysis/primary_episode_analysis.Rmd.

Data dictionary:

|**Variable name** | **Description** | **Type** | **Categories (if applicable)** |
|------------------|-----------------|------------------|-----------------|
|**PatientID2** | Recoded patient ID | chr | A - N |
|**Virus** | Virus ID | chr| "CMV"|
|**count** | log10 oral viral load | num | |
|**plasma_count** | log10 plasma viral load | num | |
|**gender** |  | num | 1 = ; 2 =  |
|**momhiv** | Maternal HIV status | chr | "neg", "pos" |
|**days** | Days after enrollment | num |  |
|**Use** |  indicator of whether infant has > 1 month data| chr | "N" "Y" |
|**days2** | Days after primary infection (oral shedding defintion) | num |  |
|**days_orig** | Days after primary infection (Gantt et al. defintion) | num |  |
|**days_dob** | Days after infant birth day| num |  |
|**last_swab_day** | Episode day of final observed swab (used for exclusion)| num |  |
|**count2** | log10 oral viral load (0 coded as NA) | num | |


Note: count2 only exists in final_episode_data.csv and CMVprimaryEpisode.Rdata
