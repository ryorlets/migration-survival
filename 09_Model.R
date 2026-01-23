# ------------------------------------------------------------------------------
# Function 3: Estimate sex-specific effects of migration on survival from pooled 
# weighted logistic regression
# ------------------------------------------------------------------------------

####
# Function needs to take the output of the estimate_ipw() function
# Then set classes of variables
# Then fit the discrete-time hazard model and return logistic.fit.weighted

# By default, no checks will be run
estimate_mod <- function(data_mod, check = FALSE) {
  # Load packages
  pacman::p_load(broom,     # tidy
                 dplyr,
                 sandwich,  # sandwich estimator
                 ggplot2
                 )
  # No need to subset data because we will use subsetted sex-specific outcomes-only person-assessment data frame

## -----------------------------------------------------------------------------------------------------------------------------
  # Set classes of variables to be used in the model
  # Make the outcome of death a factor variable
  data_mod$death_assessed <- as.factor(data_mod$death_assessed)
  # Make the exposure of out-migration a numeric variable
  data_mod$exp_num <- as.numeric(as.character(data_mod$exp))
  # Model the time measurement of visit a factor variable
  data_mod$visit <- as.factor(data_mod$visit)

## -----------------------------------------------------------------------------------------------------------------------------
  # Fit the weighted pooled logistic regression
  logistic.fit.weighted <- glm(death_assessed ~ exp_num + visit + exp_num*visit, 
                                                family = quasibinomial(),
                                                #family = binomial(link = 'logit'), 
                                                data = data_mod,
                                                # Include stabilized joint IPTW*IPCW
                                                weights = sw)
  
  # Calculate the covariance matrix using 'HC0' (refers to the sandwich estimator)
  covmat <- vcovHC(logistic.fit.weighted, type = "HC0")
  
  #Calculate the standard error
  se <- sqrt(diag(covmat))
  
  # browser()
  if(length(se) != length(coef(logistic.fit.weighted)))
  message('SE vector is of length ', length(se), ' versus coef vector of length ', length(coef(logistic.fit.weighted)), '. cbind() warning is benign because the final SE is inestimable.')
  # Bind together model output
  model_output.logistic.fit.weighted <- cbind(
    # 1. exponentiated coefficients
    Estimate = exp(coef(logistic.fit.weighted)),
    # 2. robust standard errors
    `Robust SE` = se,
    # 3. 95% confidence intervals
    # Note that qnorm(0.975) approximately equals 1.96
    Lower = exp(coef(logistic.fit.weighted) - qnorm(0.975) * se),
    Upper = exp(coef(logistic.fit.weighted) + qnorm(0.975) * se)
  )
  
  # To prepare for creating the prediction data set, in which we need a unit of 
  # time, which will range from 1 to the max number of visits
  # Calculate the max visit number
  visit.max <- max(data_mod$visit_num)
  
  ###################################
  if (check == TRUE) {
    # Checks
    #class(datasub_mod$exp_outmig_pattern_bin)
    # factor
    #class(datasub_mod$exp_outmig_pattern_bin_num)
    # numeric
    # covmat
    # View(model_output.logistic.fit.weighted)
    # write_xlsx(model_output.logistic.fit.weighted, 'model_output.logistic.fit.weighted.17Feb2023.xlsx')
    
  }
  #################################
  # Prepare for creating the prediction data set
  # Create a data frame of the original model output
  model_output <- tidy(logistic.fit.weighted)
  # Now add the new visit_rank_DoD object to the data frame with the model output
  model_output <- cbind(model_output, visit.max)
  ####
  # Return the output from the weighted pooled logistic regression with the visit
  ####
  # model_output
  
  # Also return the clean model output with exponentiated coefficients and 
  # adjusted SEs
  # View(model_output.logistic.fit.weighted)
  
  # Coerce model_output into a data frame
  # Return second-third rows to focus on exposure and time
  model_output <- as.data.frame(model_output.logistic.fit.weighted)
  model_output <- knitr::kable(model_output[2:30, ], digits = 4)
  
  #model_output.logistic.fit.weighted <- as.data.frame(model_output.logistic.fit.weighted)
  
  ####
  # Also return the clean model output with exponentiated coefficients and 
  # adjusted SEs
  ####
  # We want to return both the raw model output to feed to the predict function 
  # and the clean model output with exponentiated coefficients and robust SE
  list(model = logistic.fit.weighted, 
       exp_coefs = model_output)
}

