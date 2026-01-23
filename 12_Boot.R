# Function: draw a sample from the superpopulation and estimate a bootstrapped
#           confidence interval for the cumulative survival curve

# Subset data by sex and maximum number of visit for which we have adequate data

# We want to bootstrap, then impute, then run the stack of functions in 11_RunFunctions.R

# Set a default visit max here in the original function
boot_curve <- function(data, B, VISIT.MAX = 15, ncores = 1, trim = TRUE) {
  # B: number of bootstrap replicates
  # data: data set created before weights are calculated, and subset by sex and visit.max
  
  idlist <- unique(data$IndividualId)
  N <- length(idlist)
  
  registerDoParallel(ncores)
  # Loop from 1:number of bootstrap replicates I want to create, and put into a list
  boot_output <- foreach(i = seq_len(B)) %dopar% {
    # Step 1: Identify the persons we want to include in our random sample
    b <- data.frame(IndividualId = sample(idlist, size = N, replace = TRUE))
    
    # Step 2: Select all person-period observations associated with each person we 
    # identified in our random sample
    data_boot <- left_join(x = b, y = data, by = 'IndividualId', relationship = 'many-to-many')
    
    # Step 3: Estimate MICE
    data_imputed <- estimate_mice(data = data_boot)
    
    # Step 4: Estimate joint IPTW*IPCW in resampled data frame
    # Loop through the 1:20 data sets and call this function
    # When i==1, run the estimate_ipw() function for 1, and then it will proceed
    # to i==2, etc. and the lapply returns it as a list
    data_boot_ipw <- lapply(1:20, function(i) {
      estimate_ipw(datasub = data_imputed %>%
                     filter(visit_num <= VISIT.MAX & .imp == i))
      })
    
    # Step 5: Estimate pooled weighed logistic regression
    # Now it's running 1:20 of a list
    data_boot_ipw_mod <- lapply(1:20, function(i) {
      estimate_mod(data_mod = data_boot_ipw[[i]], 
                   check = FALSE)})
    
    # Calculate average coefficients only across imputed data sets
    # (Cannot use the SE raw from this output)
    averaged_msm_estimates <- lapply(1:20, function(i) {
      broom::tidy(data_boot_ipw_mod[[i]]$model)
    }
    ) %>%
      bind_rows(.id = ".imp") %>%
      group_by(term) %>%
      summarise(estimate = mean(estimate))
    
    # Step 6: Estimate prediction data sets
    data_pred <- bind_rows(
      lapply(1:20, function(i) {
        estimate_predict(data = data_boot_ipw_mod[[i]], 
                         # If we want to restrict the maximum # of visits
                         visit.max = VISIT.MAX, check = FALSE)
      }), 
      # We are re-introducing the .imp variable to make it explicit
      .id = ".imp"
    )
    
    # Average across the 20 predicted data sets to get the average time-specific 
    # sex-specific cumulative survival probability
    averaged_plotpts <- data_pred %>%
      # Group by visit (we want an average for each time point)
      ###### 10 Dec edit: add 'exp_num' per the 10_Prediction.R to group_by() here
      # note that we used exp_num because we needed a numeric version of the exposure that
      # we created in the treatment and no treatment prediction data sets
      group_by(exp_num, visit_num) %>%
      summarise(
        mean_survival = mean(mean_survival)
      )
    
    if(trim == TRUE){# Return the fitted model output
      data_boot_ipw_mod$data <- NULL
      data_boot_ipw_mod$y <- NULL
    }
    
    # Return a list
    # list(idlist = idlist, N = N, b = b, data_boot = data_boot)
    list(msm = averaged_msm_estimates, # averaged estimates in a data frame
         pred = averaged_plotpts #class dataframe
    )
  }
  # Ends the parallel session, closes all of the cores
  stopImplicitCluster()
  
  # Return output from boot_output
  boot_output
}


# test
# resample <- boot_curve(women_outcomes, 10)
