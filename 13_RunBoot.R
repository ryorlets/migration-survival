# ------------------------------------------------------------------------------
# Run bootstrapping function
# ------------------------------------------------------------------------------

# # Load data
load('Step06_PreWeightsClean.Rdata')
d <- Step06_PreWeightsClean

# Load packages
pacman::p_load(foreach, 
               doParallel,
               tictoc,
               dplyr,
               ggplot2,
               parallelly
               )
####
# Run preceding files that create:
####
source('07_MICE.R')
# Function 1: Estimate sex-specific stabilized joint IPTW*IPCW
source('08_Weights.R')
# Function 2: Estimate sex-specific effects of migration on survival from pooled weighted logistic regression
source('09_Model.R')
# Function 3: Create sex-specific prediction data sets for causal survival curves
source('10_Prediction.R')
# Do not need to run functions in the empirical data set for the purposes of this file
# source('11_RunFunctions.R')
source('12_Boot.R')

#########
# Women
#########
# Use imputed data set because the bootstrap function needs to re-run the weighting,
# modeling, and prediction for each replicate

#tic()
# Apply bootstrap function that we wrote in 12_Boot.R to women

## feed it the pre-weights clean data instead
boot_women <- boot_curve(data = d %>% filter(sex == "Female"),
                         VISIT.MAX = 15,
                         #B = 1,
                         B = 1999,
                         # Allow it to use the number of cores available, 
                         # instead of manually re-setting every time
                         #ncores = 1)
                         ncores = parallelly::availableCores())
#toc()

# Save boot_women, which has all bootstrap replicates
saveRDS(boot_women, file = 'boot_women.Rds')

#########
# Men
#########

# Apply bootstrap function that we wrote in 12_Boot.R to men
boot_men <- boot_curve(data = d %>% filter(sex == "Male"),
                       VISIT.MAX = 15,
                       # B = 5, 
                       B = 1999,
                       ncores = parallelly::availableCores())
# Save boot_men, which has all bootstrap replicates
saveRDS(boot_men, file = 'boot_men.Rds')
