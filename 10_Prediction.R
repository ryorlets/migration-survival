# ------------------------------------------------------------------------------
# Function 4: Create sex-specific prediction data sets for causal survival curves
# This generates cumulative survival probabilities (not just time-specific probabilities)
# ------------------------------------------------------------------------------

# Never specify the value of the arguments in the function unless you want it to only ever take that value (be the default)
estimate_predict <- function(data, visit.max, check = FALSE) {
  # Load packages
  pacman::p_load(broom,     # tidy
                 dplyr
  )
  # We don't need to subset the data by sex because we have sex-specific data sets now

  ## -----------------------------------------------------------------------------------------------------------------------------
  # Initialize the new variables to create the prediction data set
  # Define the time variable, that is the maximum number of visits
  # visit.max already exists; it was created in Function 1 in this file
  
  # Treated
  # Create data frame with all time points under treatment
  data_treated <- data.frame(cbind(seq(1, visit.max), 1))
  colnames(data_treated) <- c('visit', 'exp_num')
  
  # Set variable classes to match those used in the model
  data_treated$visit <- as.factor(data_treated$visit)
  data_treated$exp_num <- as.numeric(data_treated$exp_num)
  
  ## -----------------------------------------------------------------------------------------------------------------------------
  # Untreated
  # Create data frame with all time points under no treatment
  data_untreated <- data.frame(cbind(seq(1, visit.max), 0))
  colnames(data_untreated) <- c('visit', 'exp_num')
  
  # Set variable classes to match those used in the model
  data_untreated$visit <- as.factor(data_untreated$visit)
  data_untreated$exp_num <- as.numeric(data_untreated$exp_num)
  
  ## -----------------------------------------------------------------------------------------------------------------------------
  # To the treated data set, calculate the predicted probability of survival (which is 1-probability of death)
  data_treated$p <- 1 - predict(data$model, newdata = data_treated, type = 'response')
  
  ## -----------------------------------------------------------------------------------------------------------------------------
  # To the untreated data set, calculate the predicted probability of survival (which is 1-probability of death)
  data_untreated$p <- 1 - predict(data$model, newdata = data_untreated, type = 'response')
  
  ## -----------------------------------------------------------------------------------------------------------------------------
  # Add a cumulative survival probability variable to the treated prediction data set
  data_treated <- data_treated %>%
    arrange(visit) %>%
    mutate(s = cumprod(p))
  
  ## -----------------------------------------------------------------------------------------------------------------------------
  # Add a cumulative survival probability variable to the untreated prediction data set
  data_untreated <- data_untreated %>%
    arrange(visit) %>%
    mutate(s = cumprod(p))
  
  ## -----------------------------------------------------------------------------------------------------------------------------
  # Create concatenated data set, only keep s, rand, and visit
  data_both <- dplyr::bind_rows(data_treated, data_untreated)
  
  # Order columns and retain only the columns with time (visit_rank_DoD), randomization/treatment (exp_outmig_pattern_bin_num), and the cumulative probability of survival (s)
  data_both <- data_both[ ,c('visit', 'exp_num', 's')]
  
  ## -----------------------------------------------------------------------------------------------------------------------------
  # Calculate the average cumulative survival probability at each time point
  data.results1 <- data_both %>%
    group_by(visit, exp_num) %>%
    dplyr::summarize(mean_survival = mean(s))
  
  ## -----------------------------------------------------------------------------------------------------------------------------
  # Add a variable that treats time as a numeric variable
  data.results1$visit_num <- as.numeric(as.character(data.results1$visit))
  
  # Add a row for each of Placebo and Treated where survival at time 0 is 1.
  data.results2 <- dplyr::bind_rows(c(visit_num = 0, exp_num = 0, mean_survival = 1), 
                                    c(visit_num = 0, exp_num = 1, mean_survival = 1), 
                                    data.results1)
  
  ## -----------------------------------------------------------------------------------------------------------------------------
  # Add a variable that treats randomization as a factor
  data.results2$exp_factor <- factor(data.results2$exp_num, labels = c("Resident", "Migrant"))
  
  #####################
  if (check == TRUE) {
    # Check maximum number of visits
    # print(visit.max)
    # View(data_untreated)
    # View(data_treated)
    # View(data_both)
    # View(data.results1)
    # View(data.results2)
    # View(data.results2)
    # Comparing the factor version of the exposure to the numeric version of the 
    # exposure in this version of the data set confirms that 0=Resident and 1=Migrant
    
  }
  # Return prediction data set
  data.results2
}