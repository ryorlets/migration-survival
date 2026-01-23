# Rachel Yorlets
# 07_MICE
# 18 Aug 2025
# The goal of this code is to write a function to impute any missing data 

# ------------------------------------------------------------------------------
# Function 1: Multiple imputation via chained equations for missing covariates
# ------------------------------------------------------------------------------
estimate_mice <- function(datasub, check = FALSE) {
  # Load packages
  pacman::p_load(mice,
                 dplyr,
                 lubridate,
                 reshape2,
                 RColorBrewer,
                 ggplot2,
                 tableone,
                 doBy,            #summaryBy()
                 writexl
  )
  # Make a numeric version of 'visit'
  datasub$visit_num <- as.numeric(datasub$visit)
  
  # Fit imputation model
  # Multiple imputation via chained eqns for missing variables
  #-------------------------------------------------------------------------------
  # The number of iterations is the quality of the imputations
  # The number of imputated data sets is about the statistical inference
  
  # Men
  i <- mice(data = datasub,
                # We want to set a seed because this is a stochastic process
                seed = 082022,
                # Number of imputed data sets
                m = 20,
                # maxit specifies the number of iterations. Default is 5
                maxit = 10,
                # predicted mean model
                method = 'pmm',
                formulas = list(
                  edu ~ exp + exp_lag + exp_hx + died_ever + visit_num + age + sei + cohabit + cd4_base,
                  sei ~ exp + exp_lag + exp_hx + died_ever + visit_num + age + edu + cohabit + cd4count + cd4_base,
                  cohabit ~ exp + exp_lag + exp_hx + died_ever + visit_num + age + edu + sei + cd4count + cd4_base,
                  cd4count ~ exp + exp_lag + exp_hx + died_ever + visit_num + age + edu + sei + cohabit + cd4_base,
                  cd4_base ~ exp + exp_lag + exp_hx + died_ever + visit_num + age + edu + sei + cohabit + cd4count
                )
  )
  complete(i, action = "long")
}

########################
# Notes
########################
#Level 2 variable: any baseline value, which does not change over time and is 
#clustered by a level, in this case, the IndividualId
#Level 1 variable: any variable value at time t

#Concept in imputation that you want to fill in data - always include an 
#'outcome_ever' variable in the model, and the exposure, and any other variables 
#needed to satisfy MCAR assumption, including things like time since ART start - 
#  everything in the weights
#######################
# Load packages
# pacman::p_load(mice,
#               dplyr,
#               lubridate,
#               reshape2,
#               RColorBrewer,
#               ggplot2,
#               tableone,
#               doBy,            #summaryBy()
#               writexl
#               )
# 
# # Load data
# load('Step06_PreWeightsClean.Rdata')
#   # View(Step06_PreWeightsClean)
# 
# # Rename data
# d <- Step06_PreWeightsClean
# 
# # Make a numeric version of 'visit'
# d$visit_num <- as.numeric(d$visit)
# 
# # Stratify by sex
# d_men <- subset(d, sex == 'Male')
# d_women <- subset(d, sex == 'Female')
# 
# # CD4
# # subset to remove rows in which cd4 is NA
# d_men_cd4countnotna <- subset(d_men, !is.na(cd4count))
# d_women_cd4countnotna <- subset(d_women, !is.na(cd4count))
# 
# ## men
# # Run the functions length, mean, and sd on the value of 'cd4' for each group, 
# # stratified by visit (note, we can add more groups if we did visit + othergroup)
# men_cd4 <- summaryBy(cd4count ~ visit_num, data = d_men_cd4countnotna, FUN = c(length, mean, sd))
# # Rename column with length to just n
# names(men_cd4)[names(men_cd4) == 'cd4count.length'] <- 'n'
# # Calculate standard error of the mean
# men_cd4$cd4count.se <- men_cd4$cd4count.sd / sqrt(men_cd4$n)
# # Calculate upper bound of 95% CI
# men_cd4$CI95 <- paste("(",round(men_cd4$cd4count.mean-(1.96*men_cd4$cd4count.se), digits = 1), round(men_cd4$cd4count.mean+(1.96*men_cd4$cd4count.se), digits = 1),")", sep = ", ")
# #men_cd4$cd4count.CLlower <- men_cd4$cd4count.mean-(1.96*men_cd4$cd4count.se)
# #men_cd4$cd4count.CLupper <- men_cd4$cd4count.mean+(1.96*men_cd4$cd4count.se)
# men_cd4_df <- men_cd4 %>%
#   filter(visit_num < 16) 
# 
# ## women
# # Run the functions length, mean, and sd on the value of 'cd4' for each group, 
# # stratified by visit (note, we can add more groups if we did visit + othergroup)
# women_cd4 <- summaryBy(cd4count ~ visit_num, data = d_women_cd4countnotna, FUN = c(length, mean, sd))
# # Rename column with length to just n
# names(women_cd4)[names(women_cd4) == 'cd4count.length'] <- 'n'
# # Calculate standard error of the mean
# women_cd4$cd4count.se <- women_cd4$cd4count.sd / sqrt(women_cd4$n)
# # Calculate upper bound of 95% CI
# women_cd4$CI95 <- paste("(",round(women_cd4$cd4count.mean-(1.96*women_cd4$cd4count.se), digits = 1), round(women_cd4$cd4count.mean+(1.96*women_cd4$cd4count.se), digits = 1),")")
# women_cd4$Interval <- paste(women_cd4$CI95, collapse = ", ")
# women_cd4_df <- women_cd4 %>%
#   filter(visit_num < 16) 
# 
# ci <- cbind(men_cd4_df, women_cd4_df)
# write_xlsx(ci, 'Observed CD4_Means_CIs.xlsx')
# 
# # Multiple imputation via chained eqns for missing variables
# #-------------------------------------------------------------------------------
# # The number of iterations is the quality of the imputations
# # The number of imputated data sets is about the statistical inference
# 
# # Men
# i_men <- mice(data = d_men,
#               # We want to set a seed because this is a stochastic process
#               seed = 082022,
#           # Number of imputed data sets
#           m = 20,
#           # maxit specifies the number of iterations. Default is 5
#           maxit = 10,
#           # predicted mean model
#           method = 'pmm',
#           formulas = list(
#             edu ~ exp + exp_lag + exp_hx + died_ever + visit_num + age + sei + cohabit + cd4_base,
#             sei ~ exp + exp_lag + exp_hx + died_ever + visit_num + age + edu + cohabit + cd4count + cd4_base,
#             cohabit ~ exp + exp_lag + exp_hx + died_ever + visit_num + age + edu + sei + cd4count + cd4_base,
#             cd4count ~ exp + exp_lag + exp_hx + died_ever + visit_num + age + edu + sei + cohabit + cd4_base,
#             cd4_base ~ exp + exp_lag + exp_hx + died_ever + visit_num + age + edu + sei + cohabit + cd4count
#           )
# )
# #i_men$loggedEvents
# #i_men$predictorMatrix
# #i_men$formulas
# #i_men$blocks
# 
# # Women
# i_women <- mice(data = d_women,
#               # We want to set a seed because this is a stochastic process
#               seed = 082022,
#               # Number of imputed data sets
#               m = 20,
#               # maxit specifies the number of iterations. Default is 5
#               maxit = 10,
#               # predicted mean model
#               method = 'pmm',
#               formulas = list(
#                 edu ~ exp + exp_lag + exp_hx + died_ever + visit_num + age + sei + cohabit + cd4_base,
#                 sei ~ exp + exp_lag + exp_hx + died_ever + visit_num + age + edu + cohabit + cd4count + cd4_base,
#                 cohabit ~ exp + exp_lag + exp_hx + died_ever + visit_num + age + edu + sei + cd4count + cd4_base,
#                 cd4count ~ exp + exp_lag + exp_hx + died_ever + visit_num + age + edu + sei + cohabit + cd4_base,
#                 cd4_base ~ exp + exp_lag + exp_hx + died_ever + visit_num + age + edu + sei + cohabit + cd4count
#               )
# )
# #i_women$loggedEvents
# #i_women$predictorMatrix
# #i_women$formulas
# #i_women$blocks
# 
# #-------------------------------------------------------------------------------
# # Diagnostics
# #-------------------------------------------------------------------------------
# 
# # Check average value within (0) observed non-missing data and (1-20) combined 
# # observed+imputed data sets, of which there are 20, per our imputation above
# #-------------------------------------------------------------------------------
# 
# # Men---------------------------------------------------------------------------
# # complete() allows us to look at the full data set with observed and imputed values
# i_men_complete <- complete(i_men, 
#                            action = 'long', 
#                            include = TRUE)
# # Check we have 21 rows for each imputed data set (20 imputed data sets, 1 observed)
# i_men_complete %>% count(.id)
# 
# # Women---------------------------------------------------------------------------
# # complete() allows us to look at the full data set with observed and imputed values
# i_women_complete <- complete(i_women, 
#                            action = 'long', 
#                            include = TRUE)
# # Check we have 21 rows for each imputed data set (20 imputed data sets, 1 observed)
# i_women_complete %>% count(.id)
# 
# #-------------------------------------------------------------------------------
# # We will print all of these checks below, after the code for each variable
# # Check edu
# #-------------------------------------------------------------------------------
# # Men---------------------------------------------------------------------------
# # isolate rows in original data where observed edu_max is missing
# i_men_comp_index_edu <- i_men_complete %>%
#   filter(.imp == 0 & is.na(edu)) %>%
#   pull(.id)
# 
# i_men_comp_edu <- i_men_complete %>%
#   mutate(was_edu_na = .id %in% i_men_comp_index_edu) %>%
#   group_by(.imp, was_edu_na) %>%
#   summarise(mean_edu = mean(edu, na.rm = TRUE)) %>%
#   filter(.imp == 0 | was_edu_na == TRUE)
# 
# # Women---------------------------------------------------------------------------
# # isolate rows in original data where observed edu_max is missing
# i_women_comp_index_edu <- i_women_complete %>%
#   filter(.imp == 0 & is.na(edu)) %>%
#   pull(.id)
# 
# i_women_comp_edu <- i_women_complete %>%
#   mutate(was_edu_na = .id %in% i_women_comp_index_edu) %>%
#   group_by(.imp, was_edu_na) %>%
#   summarise(mean_edu = mean(edu, na.rm = TRUE)) %>%
#   filter(.imp == 0 | was_edu_na == TRUE)
# 
# #----
# # For categorical variables, we can compare the proportion of values in each 
# # category, which is the same as the mean when we have a binary 0/1 variable
# 
# # Check sei
# #-------------------------------------------------------------------------------
# # Men---------------------------------------------------------------------------
# # isolate rows in original data where observed sei is missing
# i_men_comp_index_sei <- i_men_complete %>%
#   filter(.imp == 0 & is.na(sei)) %>%
#   pull(.id)
# 
# i_men_comp_sei <- i_men_complete %>%
#   mutate(was_sei_na = .id %in% i_men_comp_index_sei) %>%
#   group_by(.imp, was_sei_na) %>%
#   summarise(mean_sei = mean(sei, na.rm = TRUE)) %>%
#   filter(.imp == 0 | was_sei_na == TRUE)
# 
# # Women---------------------------------------------------------------------------
# # isolate rows in original data where observed SEIdx is missing
# i_women_comp_index_sei <- i_women_complete %>%
#   filter(.imp == 0 & is.na(sei)) %>%
#   pull(.id)
# 
# i_women_comp_sei <- i_women_complete %>%
#   mutate(was_sei_na = .id %in% i_women_comp_index_sei) %>%
#   group_by(.imp, was_sei_na) %>%
#   summarise(mean_sei = mean(sei, na.rm = TRUE)) %>%
#   filter(.imp == 0 | was_sei_na == TRUE)
# 
# # Check cohabit
# #-------------------------------------------------------------------------------
# # Men
# # isolate rows in original data where observed SEIdx is missing
# i_men_comp_index_cohab <- i_men_complete %>%
#   filter(.imp == 0 & is.na(cohabit)) %>%
#   pull(.id)
# 
# i_men_comp_cohab <- i_men_complete %>%
#   mutate(was_cohab_na = .id %in% i_men_comp_index_cohab) %>%
#   group_by(.imp, was_cohab_na) %>%
#   summarise(mean_cohabit = mean(cohabit, na.rm = TRUE)) %>%
#   filter(.imp == 0 | was_cohab_na == TRUE)
# 
# # Women
# # isolate rows in original data where observed SEIdx is missing
# i_women_comp_index_cohab <- i_women_complete %>%
#   filter(.imp == 0 & is.na(cohabit)) %>%
#   pull(.id)
# 
# i_women_comp_cohab <- i_women_complete %>%
#   mutate(was_cohab_na = .id %in% i_women_comp_index_cohab) %>%
#   group_by(.imp, was_cohab_na) %>%
#   summarise(mean_cohabit = mean(cohabit, na.rm = TRUE)) %>%
#   filter(.imp == 0 | was_cohab_na == TRUE)
# 
# # Check cd4_base
# #-------------------------------------------------------------------------------
# # Men
# # isolate rows in original data where observed baseline cd4 is missing
# i_men_comp_index_basecd4 <- i_men_complete %>%
#   filter(.imp == 0 & is.na(cd4_base)) %>%
#   pull(.id)
# 
# i_men_comp_basecd4 <- i_men_complete %>%
#   mutate(was_basecd4_na = .id %in% i_men_comp_index_basecd4) %>%
#   group_by(.imp, was_basecd4_na) %>%
#   summarise(mean_cd4_base = mean(cd4_base, na.rm = TRUE)) %>%
#   filter(.imp == 0 | was_basecd4_na == TRUE)
# 
# # Women
# # isolate rows in original data where observed baseline cd4 is missing
# i_women_comp_index_basecd4 <- i_women_complete %>%
#   filter(.imp == 0 & is.na(cd4_base)) %>%
#   pull(.id)
# 
# i_women_comp_basecd4 <- i_women_complete %>%
#   mutate(was_basecd4_na = .id %in% i_women_comp_index_basecd4) %>%
#   group_by(.imp, was_basecd4_na) %>%
#   summarise(mean_cd4_base = mean(cd4_base, na.rm = TRUE)) %>%
#   filter(.imp == 0 | was_basecd4_na == TRUE)
# 
# # Check cd4
# #-------------------------------------------------------------------------------
# # Men
# # isolate rows in original data where observed baseline cd4 is missing
# i_men_comp_index_cd4 <- i_men_complete %>%
#   filter(.imp == 0 & is.na(cd4count)) %>%
#   pull(.id)
# 
# i_men_comp_cd4 <- i_men_complete %>%
#   mutate(was_cd4_na = .id %in% i_men_comp_index_cd4) %>%
#   group_by(.imp, was_cd4_na, visit_num) %>%
#   summarise(mean_cd4 = mean(cd4count, na.rm = TRUE), 
#             n = n()) %>%
#   filter(.imp == 0 | was_cd4_na == TRUE)
# 
# # Women
# # isolate rows in original data where observed baseline cd4 is missing
# i_women_comp_index_cd4 <- i_women_complete %>%
#   filter(.imp == 0 & is.na(cd4count)) %>%
#   pull(.id)
# 
# i_women_comp_cd4 <- i_women_complete %>%
#   mutate(was_cd4_na = .id %in% i_women_comp_index_cd4) %>%
#   group_by(.imp, was_cd4_na, visit_num) %>%
#   summarise(mean_cd4 = mean(cd4count, na.rm = TRUE), 
#             n = n()) %>%
#   filter(.imp == 0 | was_cd4_na == TRUE)
# 
# # Print out tables of means (which, for binary variables, is the proportion of the higher value)
# # Men
# i_men_comp_all <- cbind(i_men_comp_edu[,c(1,3)], i_men_comp_sei[,3], i_men_comp_cohab[,3], i_men_comp_basecd4[,3])
# 
# # Men CD4
# # View(i_men_comp_cd4)
# plot_men <- ggplot(data = i_men_comp_cd4, aes(x = visit_num, y = mean_cd4)) +
#                      geom_point(aes(color = .imp > 0))
# 
# table(i_men_comp_cd4$mean_cd4, i_men_comp_cd4$visit_num)
# 
# # Women
# i_women_comp_all <- cbind(i_women_comp_edu[,c(1, 3)], i_women_comp_sei[,3], i_women_comp_cohab[,3], i_women_comp_basecd4[,3])
# 
# # Women CD4
# # View(i_women_comp_cd4[,3])
# plot_women <- ggplot(data = i_women_comp_cd4, aes(x = visit_num, y = mean_cd4)) +
#   geom_point(aes(color = .imp > 0))
# 
# 
# #-------------------------------------------------------------------------------
# # Following diagnostics described in van buren textbook
# #-------------------------------------------------------------------------------
# # Worm plot | plot(type = 3) specifies that we want a worm plot
# #-------------------------------------------------------------------------------
# plot(i_men, type = 3, variable = 'SEIdx')
# 
# 
# #-------------------------------------------------------------------------------
# # Check stripplot() to produce the numerical values per imputation
# # (ideal for continuous variables; only shows binomial dist for binary ones)
# #-------------------------------------------------------------------------------
# 
# # Men
# stripplot(i_men, sei + cd4_base ~ .imp)
# # can show individual plots like this: (i_men, sei ~ .imp)
# 
# # Women
# stripplot(i_women, sei + cd4_base ~ .imp)
# 
# #-------------------------------------------------------------------------------
# # Kernel plot (density plot) estimates the marginal distribution of the observed
# # data (blue) and m=20 estimates of the variable in each imputed data set
# #-------------------------------------------------------------------------------
# densityplot(i_men, .imp ~ edu + sei + cohabit + cd4_base)
# 
# densityplot(i_men, .imp ~ cd4count | visit_num, subset = visit_num < 15)
# 
# # Shows consistency across imputed data sets because lines are grouped
# 
# densityplot(i_women, .imp ~ edu + sei + cohabit + cd4_base)
# 
# densityplot(i_women, .imp ~ cd4count | visit_num, subset = visit_num < 15)
# 
# #-------------------------------------------------------------------------------
# # BW plot (Box and whisker) - for continuous variables being imputed
# # Allows us to check if the imputed data (see x-axis for imputation number)
# # have a similar distribution to the observed values, which are in blue at imputation number = 0
# #-------------------------------------------------------------------------------
# bwplot(i_men, cd4_base + sei ~ .imp)
# 
# bwplot(i_women, cd4_base + sei ~ .imp)
# 
# #-------------------------------------------------------------------------------
# # xy plot - for categorical variables
# #-------------------------------------------------------------------------------
# xyplot(i_men, edu + sei + cohabit ~ .imp)
# 
# xyplot(i_men, sei ~ .imp)
# 
# xyplot(i_men,
#        cohabit ~ .imp)
# 
# # Plot (Spaghetti plots, per jason)
# plot(x = i_men,
#      layout = c(2, 5))
# #-------------------------------------------------------------------------------
# # Worm plots
# #-------------------------------------------------------------------------------
# 
# # Custom function to run proportions for imputed v observed categorical variables
# # https://gist.github.com/NErler/0d00375da460dd33839b98faeee2fdab
# 
# #-------------------------------------------------------------------------------
# # Export imputed data sets
# # Could use 'include = TRUE' if we want to retain original variables (w/o imputation)
# men_imputed <- complete(i_men, action = "long")
# save(men_imputed, file = 'men_imputed.Rdata')
# 
# women_imputed <- complete(i_women, action = "long")
# save(women_imputed, file = 'women_imputed.Rdata')
# 
# # Combine imputed data sets
# data_imputed <- rbind(men_imputed, women_imputed)
# data_imputed <- as.data.frame(data_imputed)
# saveRDS(data_imputed, file = 'data_imputed.Rds')
# 
# # Note that this gives us a data frame with stacked imputed data sets
# # We have two new columns in this data set, and did not include unimputed variables
# # .imp indicates which of the imputed data sets the row belongs to
# # .id numbers the rows (I think)
