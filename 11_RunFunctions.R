# ------------------------------------------------------------------------------
# Run Functions
# ------------------------------------------------------------------------------

# NOTE: When running this file, check the restriction on the number of visits 
# in order to print the tables and curves that show the timeframe you want


# # Load data
load('Step06_PreWeightsClean.Rdata')
d <- Step06_PreWeightsClean

# Load packages
pacman::p_load(writexl,
               dplyr,
               readr,
               stringr,
               janitor,
               purrr)

# Source Functions
# Run preceding files that create the following:
source('07_MICE.R')
# Function 2: Estimate sex-specific stabilized joint IPTW*IPCW
source('08_Weights.R')
# Function 3: Estimate sex-specific effects of migration on survival from pooled 
# weighted logistic regression models
source('09_Model.R')
# Function 4: Create sex-specific prediction data sets for causal survival curves
source('10_Prediction.R')

# ------------------------------------------------------------------------------
# Function 1: Estimate sex-specific multiple imputation models via chained eqns
# ------------------------------------------------------------------------------
# From 07_MICE.R
# Women
#-------------------------------------------------------------------------------
women_imputed <- estimate_mice(data = d %>%
                                 filter(sex == "Female"))

# Men
#-------------------------------------------------------------------------------
men_imputed <- estimate_mice(data = d %>%
                         filter(sex == "Male"))

# Combined imputed data sets for men and women
data_imputed <- as.data.frame(rbind(men_imputed, women_imputed))
# saveRDS(data_imputed, file = 'data_imputed.Rds')

# ------------------------------------------------------------------------------
# Function 2: Estimate sex-specific stabilized joint IPTW*IPCW
# ------------------------------------------------------------------------------
# From 08_Weights.R
# Women
#-------------------------------------------------------------------------------
# Where i is the loop index
# Remember this returns a list
women_outcomes <- lapply(1:20, function(i){
  estimate_ipw(datasub = data_imputed %>%
                 filter(visit_num <= 15, 
                        sex == "Female", 
                        .imp == i))})

## Note: For women, only 1 woman remains at visit 16, so we subset to address 
## non-positivity Can also print this without any restrictions on the visits 
## so that we can see the full FUP for the supplement 
## datasub = data_imputed %>% filter(sex == 'Female')
#-------------------------------------------------------------------------------

# Men
#-------------------------------------------------------------------------------
men_outcomes <- lapply(1:20, function(i){
  estimate_ipw(datasub = data_imputed %>%
                  filter(visit_num <= 15,
                         sex == "Male",
                         .imp == i))})
                             
## Note: For men, we see data scarcity after 2.5 years, which would be visit 5)
# Can also print this without any restrictions on the visits so that we can see 
# full FUP for the supplement datasub = data_imputed %>% filter(sex == 'Male')

# Check function - check sex of all individuals in each data set
all(women_outcomes$sex == "Female")
all(men_outcomes$sex == "Male")

#-------------------------------------------------------------------------------
# Function 3: Fit sex-specific marginal structural models
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Women
#-------------------------------------------------------------------------------
# Apply function created in 09_Model.R
women_mod <- lapply(1:20, function(i){
  estimate_mod(data_mod = women_outcomes[[i]], 
               check = FALSE)})

# women_mod is what we want to use in the next function
# The next chunks of code are to create/save the averaged model output from 
# the models from the imputed data sets so we can report an MSM in the supplement
#-------

# if you want to present the msm as standalone supplementary results
women_averaged_model_output <- lapply(1:20, function(i) {
  broom::tidy(women_mod[[i]]$model)
}
) %>%
  bind_rows(.id = ".imp") %>%
  group_by(term) %>%
  summarise(estimate = mean(estimate))

# Save table of model output
#write_xlsx(women_averaged_model_output, 'Table_Women_AvgMSMOutput.xlsx')

# Men---------------------------------------------------------------------------
# Apply function created in 09_Model.R
men_mod <- lapply(1:20, function(i){
  estimate_mod(data_mod = men_outcomes[[i]], 
               check = FALSE)})

# men_mod is what we want to use in the next function
# The next chunks of code are to create/save the averaged model output from 
# the models from the imputed data sets so we can report an MSM in the supplement
#-------

men_averaged_model_output <- lapply(1:20, function(i) {
  broom::tidy(men_mod[[i]]$model)
}
) %>%
  bind_rows(.id = ".imp") %>%
  group_by(term) %>%
  summarise(estimate = mean(estimate))

# Save table of model output
write_xlsx(men_averaged_model_output, 'Table_Men_AvgMSMOutput.xlsx')

#-------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Function 4: Create sex-specific prediction data sets for causal survival curves
# ------------------------------------------------------------------------------
# From 10_Prediction.R -> This function cumulates the survival probabilities
# ------------------------------------------------------------------------------
# Women
women_predict <- bind_rows(lapply(1:20, function(i){
  estimate_predict(data = women_mod[[i]], 
                                  # If we want to restrict the maximum # of visits
                                  visit.max = 15, check = FALSE)}), .id = ".imp")

 #lapply(1:20, function(i) {
#  as.data.frame(women_mod[[i]]$model$fitted.values)
#}) %>%
#  bind_rows(.id = ".imp") %>%
#  group_by(.imp) %>%
#  summarise(n())

# Average across the 20 predicted data sets to get the average time-specific 
# sex-specific cumulative survival probability
women_averaged_plotpts <- women_predict %>%
  # Group by visit (we want an average for each time point)
  group_by(visit_num, exp_factor) %>%
  summarise(mean_survival = mean(mean_survival)) 
  
save(women_averaged_plotpts, file = "women_averaged_plotpts.Rdata")

# Save table of model output
write_xlsx(women_averaged_plotpts, 'Table_Women_AvgPlotPts.xlsx')

# Men
men_predict <- bind_rows(lapply(1:20, function(i){
  estimate_predict(data = men_mod[[i]], 
                   # If we want to restrict the maximum # of visits
                   visit.max = 15, check = FALSE)}), .id = ".imp")

# Average across the 20 predicted data sets to get the average time-specific 
# sex-specific cumulative survival probability
men_averaged_plotpts <- men_predict %>%
  # Group by visit (we want an average for each time point)
  group_by(visit_num, exp_factor) %>%
  summarise(mean_survival = mean(mean_survival)) 

save(men_averaged_plotpts, file = "men_averaged_plotpts.Rdata")

# When we average across these 20 predicted data sets for each men and women, we 
# get one curve for each men and women. It's an average at each time point, and 
# this is the point estimate used in the causal RD, RR calculations, and the point
# estimates for the cumulative probability at each time point (ie, the points on the plots)
# So after we average these, we have the final causal curves for the paper
# We need to next do the boot to do this 1999 times to get the CIs
