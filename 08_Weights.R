#-------------------------------------------------------------------------------
# Function 2: Estimate sex-specific stabilized joint IPTW*IPCW
#-------------------------------------------------------------------------------
# By default, no checks will be run
 estimate_ipw <- function(datasub, check = FALSE) {
   # Load packages
   pacman::p_load(haven,
                  dplyr,
                  broom
                  )
   
   ## -----------------------------------------------------------------------------------------------------------------------------
   # 1.  IPTW Numerator
   # 1a. Fit the probability of being exposed
   num.fit.a <- glm(exp ~ exp_lag,
                  data = datasub,
                  family = binomial
                  )
   
   # Include a check to make sure we have subsetted by sex when we run this function
   if(length(unique(datasub$sex)) != 1){
     stop("Need to subset data to either 'Male' or 'Female' to use function")
   }
   
   # 1b. Calculate the probability using predict()
   datasub$num.pred.a <- predict(num.fit.a, newdata = datasub, type = 'response')
   
   ## -----------------------------------------------------------------------------------------------------------------------------
   # 2.  IPTW Denominator
   # 2a. Fit the probability of being exposed given prior exposure and current covariates and time period
   den.fit.a <- glm(exp ~ exp_lag + exp_hx + age + edu + sei + cohabit,
                  data = datasub,
                  family = binomial
                  )
   
   # 2b. Calculate the probability using predict()
   datasub$den.pred.a <- predict(den.fit.a, newdata = datasub, type = 'response')
   
   ## -----------------------------------------------------------------------------------------------------------------------------
   # Numerator
   datasub$num.a <- ifelse(datasub$exp == 1, datasub$num.pred.a, (1 - datasub$num.pred.a))
   # Denominator
   datasub$den.a <- ifelse(datasub$exp == 1, datasub$den.pred.a, (1 - datasub$den.pred.a))
   
   # Cumulate the weights for each person-assessment after the first visit
   datasub <- datasub %>%
     group_by(IndividualId) %>%
     arrange(IndividualId, visit) %>%
     mutate(
       num.cumpred.a = case_when(
         (visit == 1) ~ num.a,
         (visit != 1) ~ cumprod(num.a)
       ),
       den.cumpred.a = case_when(
         (visit == 1) ~ den.a,
         (visit != 1) ~ cumprod(den.a)
       )
     )
   
   ## -----------------------------------------------------------------------------------------------------------------------------
   # Calculate unstabilized IPTW
   datasub$iptw.w <- (1 / datasub$den.cumpred.a)
   
   ## -----------------------------------------------------------------------------------------------------------------------------
   # Calculate stabilized IPTW
   datasub$iptw.sw <- (datasub$num.cumpred.a / datasub$den.cumpred.a)
   
   ## -----------------------------------------------------------------------------------------------------------------------------
   # IPCW
   # Fit the probability of not being censored within the analytic sample (n=5052)
   num.fit.c <-
     glm((censored == 0) ~ exp + exp_lag,
         data = datasub,
         family = binomial
     )
   
   # Calculate the probability using predict()
   datasub$num.pred.c <-
     predict(num.fit.c, newdata = datasub, type = 'response')
   
   ## -----------------------------------------------------------------------------------------------------------------------------
   # IPCW
   # Fit the probability of not being censored given exposure, prior exposure, exposure in the year before ART initiation, and all current covariates plus the time period
   den.fit.c <-
     glm((censored == 0) ~ exp + exp_lag + exp_hx + age + edu + sei + cohabit,
         data = datasub,
         family = binomial
     )
   
   # Calculate the probability using predict()
   datasub$den.pred.c <-
     predict(den.fit.c, newdata = datasub, type = 'response')
   
   ## -----------------------------------------------------------------------------------------------------------------------------
   # Cumulate the weights for each person-assessment after the first visit
   datasub <- datasub %>%
     group_by(IndividualId) %>%
     arrange(IndividualId, visit) %>%
     mutate(
       num.cumpred.c = case_when(
         (visit == 1) ~ num.pred.c,
         (visit != 1) ~ cumprod(num.pred.c)
       ),
       den.cumpred.c = case_when(
         (visit == 1) ~ den.pred.c,
         (visit != 1) ~ cumprod(den.pred.c)
       )
     )
   
   ## -----------------------------------------------------------------------------------------------------------------------------
   # Calculate the unstabilized IPCW
   datasub$ipcw.w <- (1 / datasub$den.cumpred.c)
   
   ## -----------------------------------------------------------------------------------------------------------------------------
   # Calculate the stabilized IPCW
   datasub$ipcw.sw <- (datasub$num.cumpred.c / datasub$den.cumpred.c)
   
   ## -----------------------------------------------------------------------------------------------------------------------------
   # Calculate unstabilized joint IPTW*IPCW
   datasub$w <- datasub$iptw.w * datasub$ipcw.w

   ## ----------------------------------------------------------------------------------------------------------------------------
   # Calculate stabilized joint IPTW*IPCW
   datasub$sw <- datasub$iptw.sw * datasub$ipcw.sw
   
   ## ----------------------------------------------------------------------------------------------------------------------------
   # Subset data set to include only person-period assessments in which there is an outcome (death) measurement
   datasub_outcomes <- subset(datasub, !is.na(death_assessed))
   
###################################
   if (check == TRUE) {
     
     # 1. IPTW
     # Summarize fit numerator values
     #cat('summary of fit numerator values:', '\n\n')
     #print(summary(num.fit.a))
     #cat('\n\n')
     # num.fit.a.f_tidy <- tidy(num.fit.a.f, conf.int = TRUE)
     
     #cat('data with new IPTW numerator values:', '\n\n')
     #View(datasub[, c('IndividualId', 'visit_rank_DoD', 'exp_outmig_pattern_bin', 'lag_exp_outmig_pattern_bin', 'num.pred.a')])
     #cat('\n\n')
     
     # Summarize fit denominator values
     # Check coefficients from model (we removed suppression_ever because it is mostly missing)
     #print(summary(den.fit.a)) # provides AIC in output, but not BIC
     
     #View(datasub[, c('IndividualId',  'visit_rank_DoD', 'exp_outmig_pattern_bin', 'lag_exp_outmig_pattern_bin', 'num.a.f', 'den.a.f', 'num.cumpred.a.f', 'den.cumpred.a.f')])
     
     # View(datasub[, c('IndividualId',  'visit_rank_DoD', 'exp_outmig_pattern_bin', 'lag_exp_outmig_pattern_bin', 'num.a.f', 'den.a.f', 'num.cumpred.a.f', 'den.cumpred.a.f', 'iptw.w.f')])
     
     # 20a. Summary
     #summary(datasub$iptw.w)
     # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
     #    1       1       1     951       4 4326366
     
     # 20c. Ratio of each weight to the smallest weight: Joe suggests checking ratio of each weight to the smallest weight
     #datasub$ratio_iptw_w_to_smallest <- datasub$iptw.w.f / min(datasub$iptw.w.f)
     
     # 20d. Summary of ratios
     #summary(datasub$ratio_iptw_w_to_smallest)
     # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
     #    1       1       1     945       4 4299436
     # Note: because the minimum value of the ratio is 1, the ratio is the same value as the weight (and same as summary)
     
     # 20e. Check how many weights exceed 100
     #length(which(datasub$ratio_iptw_w_to_smallest > 100))
     # 4102
     #nrow(datasub)
     # 33442 person-assessment periods
     # 4102/33442 = 12%
     
     # print the values that are not the same
     # View(datasub[datasub$iptw.w.f != datasub$ratio_iptw_w_to_smallest, c('IndividualId', 'visit_rank_DoD', 'exp_outmig_pattern_bin', 'num.a.f', 'num.cumpred.a.f', 'den.a.f', 'den.cumpred.a.f', 'iptw.w.f', 'ratio_iptw_w_to_smallest')])
     # Here we see that the value of the weight and the value of the ratio is very very close, but even though the minimum is '1', the decimals make it not an exact replica, so the mean is slightly off when we compare the mean of the weights versus the mean of the ratio
     
     # Histograms
     # hist(datasub$ratio_iptw_w_to_smallest)
     # hist(log10(datasub$ratio_iptw_w_to_smallest))
     #
     # # 20g. Look at histogram of ratio of unstabilized treatment weights among person-assessments in which the binary exposure==1
     # datasub_A1 <- subset(datasub, exp_outmig_pattern_bin == 1)
     # nrow(datasub_A1) # 8324
     # datasub_A1$ratio_iptw_w_to_smallest <- datasub_A1$iptw_w/min(datasub_A1$iptw_w)
     # hist(log10(datasub_A1$ratio_iptw_w_to_smallest))
     #
     # # 20h. Look at histogram of ratio of unstabilized treatment weights among person-assessments in which the binary exposure==0
     # datasub_A0 <- subset(datasub, exp_outmig_pattern_bin == 0)
     # nrow(datasub_A0) # 25118
     # datasub_A0$ratio_iptw_w_to_smallest <- datasub_A0$iptw_w/min(datasub_A0$iptw_w)
     # hist(log10(datasub_A0$ratio_iptw_w_to_smallest))
     #
     # table(datasub$exp_outmig_pattern_bin, exclude = NULL)
     # #     0     1
     # # 25118  8324
     
     # View(datasub[datasub$iptw.w.f != datasub$ratio_iptw_w_to_smallest, c('IndividualId', 'visit_rank_DoD', 'exp_outmig_pattern_bin', 'num.a.f', 'num.cumpred.a.f', 'den.a.f', 'den.cumpred.a.f', 'iptw.w.f', 'iptw.sw.f', 'ratio_iptw_w_to_smallest')])
     
     ## -----------------------------------------------------------------------------------------------------------------------------
     # 20a. Summary
     # summary(datasub$iptw.sw.f)
     #   Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
     # 0.0493   0.6815   0.8910   1.0067   0.9810 579.0089
     
     # 20c. Ratio of each weight to the smallest weight: Joe suggests checking ratio of each weight to the smallest weight
     #datasub$ratio_iptw_sw_to_smallest <-
      # datasub$iptw.sw.f / min(datasub$iptw.sw.f)
     
     # 20d. Summary of ratios
     # summary(datasub$ratio_iptw_sw_to_smallest)
     # Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
     # 1.00    13.82    18.07    20.42    19.89 11742.22
     # Note: because the minimum value of the ratio is 1, the ratio is the same value as the weight (and same as the summary)
     
     # 20e. Check how many weights exceed 100
     #length(which(datasub$ratio_iptw_sw_to_smallest > 100))
     # 340
     #nrow(datasub)
     # 33442 person-assessment periods
     # 340/33442 = 1%
     
     # View(datasub[datasub$iptw.sw.f != datasub$ratio_iptw_sw_to_smallest, c('IndividualId', 'visit_rank_DoD', 'exp_outmig_pattern_bin', 'num.a.f', 'num.cumpred.a.f', 'den.a.f', 'den.cumpred.a.f', 'iptw.w.f', 'iptw.sw.f', 'ratio_iptw_w_to_smallest', 'ratio_iptw_sw_to_smallest')])
     
     # Histograms
     #hist(datasub$ratio_iptw_sw_to_smallest)
     #hist(log10(datasub$ratio_iptw_sw_to_smallest))
     #
     # # 20g. Look at histogram of ratio of unstabilized treatment weights among person-assessments in which the binary exposure==1
     # datasub_A1 <- subset(datasub, exp_outmig_pattern_bin == 1)
     # nrow(datasub_A1) # 8324
     # datasub_A1$ratio_iptw_sw_to_smallest <- datasub_A1$iptw_sw.f/min(datasub_A1$iptw_sw.f)
     # hist(log10(datasub_A1$ratio_iptw_sw_to_smallest))
     #
     # # 20h. Look at histogram of ratio of unstabilized treatment weights among person-assessments in which the binary exposure==0
     # datasub_A0 <- subset(datasub, exp_outmig_pattern_bin == 0)
     # nrow(datasub_A0) # 25118
     # datasub_A0$ratio_iptw_sw_to_smallest <- datasub_A0$iptw_sw.f/min(datasub_A0$iptw_sw.f)
     # hist(log10(datasub_A0$ratio_iptw_sw_to_smallest))
     #
     # table(datasub$exp_outmig_pattern_bin, exclude = NULL)
     # #     0     1
     # # 25118  8324
     
     # View(datasub[, c('IndividualId',  'visit_rank_DoD', 'censored', 'exp_outmig_pattern_bin', 'lag_exp_outmig_pattern_bin', 'num.pred.a.pattern.bin.f', 'num.cumpred.a.pattern.bin.f', 'num.pred.c.pattern.bin.f', 'num.cumpred.c.pattern.bin.f')])
     ######
     # IPCW
     # Check coefficients from denominator model
     # summary(den.fit.c)
     # Standard errors look okay
     
     # View(datasub[, c('IndividualId',  'visit_rank_DoD', 'exp_outmig_pattern_bin', 'lag_exp_outmig_pattern_bin', 'num.a.f', 'num.cumpred.a.f', 'den.a.f', 'den.cumpred.a.f', 'num.pred.c.f', 'num.cumpred.c.f',  'den.pred.c.f', 'den.cumpred.c.f')])
     
     # 28a. Summary
     #summary(datasub$ipcw.w.f)
     #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
     # 1.044   1.307   1.693   2.119   2.413  50.347
     
     # 28c. Ratio of each weight to the smallest weight
     # Joe suggests checking ratio of each weight to the smallest weight
     #datasub$ratio_ipcw_w_to_smallest <- datasub$ipcw.w.f / min(datasub$ipcw.w.f)
     
     # 28d. Summary of ratios
     #summary(datasub$ratio_ipcw_w_to_smallest)
     #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
     # 1.000   1.253   1.622   2.030   2.312  48.246
     # Since the minimum ratio is 1, this tells us that this is the same as the weights, and we see that the summary for the ratios is the same as the one above
     
     # 28e. Check how many weights exceed 100
     # length(which(datasub$ratio_ipcw_w_to_smallest > 100))
     # 0
     
     # Histograms
     # 28f. Look at histogram of ratio of unstabilized treatment weights among person-assessments in which the binary exposure==1
     # datasub_A1 <- subset(datasub, exp_outmig_pattern_bin == 1)
     # datasub_A1$ratio_w_ipcw_to_smallest <- datasub_A1$w_ipcw.pattern.bin.f/min(datasub_A1$w_ipcw.pattern.bin.f)
     # hist(log10(datasub_A1$ratio_w_ipcw_to_smallest))
     #
     # # 28g. Look at histogram of ratio of unstabilized treatment weights among person-assessments in which the binary exposure==0
     # datasub_A0 <- subset(datasub, exp_outmig_pattern_bin == 0)
     # datasub_A0$ratio_w_ipcw_to_smallest <- datasub_A0$w_ipcw.pattern.bin.f/min(datasub_A0$w_ipcw.pattern.bin.f)
     # hist(log10(datasub_A0$ratio_w_ipcw_to_smallest))
     
     ## -----------------------------------------------------------------------------------------------------------------------------
     # 28a. Summary
     # summary(datasub$ipcw.sw.f)
     #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
     # 0.3797  0.9082  0.9722  1.0073  1.0386 11.3701
     
     # 28c. Ratio of each weight to the smallest weight
     # Joe suggests checking ratio of each weight to the smallest weight
    # datasub$ratio_ipcw_sw_to_smallest <-
     #  datasub$ipcw.sw.f / min(datasub$ipcw.sw.f)
     
     # 28d. Summary of ratios
     # summary(datasub$ratio_ipcw_sw_to_smallest)
     #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
     # 1.000   2.392   2.560   2.653   2.735  29.942
     # Since the minimum ratio is 1, this tells us that this is the same as the weights, and we see that the summary for the ratios is the same as the one above
     
     # 28e. Check how many weights exceed 100
     # length(which(datasub$ratio_ipcw_sw_to_smallest > 100))
     # 0
     
     # Histograms
     # 28f. Look at histogram of ratio of unstabilized treatment weights among person-assessments in which the binary exposure==1
     # datasub_A1 <- subset(datasub, exp_outmig_pattern_bin == 1)
     # datasub_A1$ratio_w_ipcw_to_smallest <- datasub_A1$w_ipcw.pattern.bin.f/min(datasub_A1$w_ipcw.pattern.bin.f)
     # hist(log10(datasub_A1$ratio_w_ipcw_to_smallest))
     #
     # # 28g. Look at histogram of ratio of unstabilized treatment weights among person-assessments in which the binary exposure==0
     # datasub_A0 <- subset(datasub, exp_outmig_pattern_bin == 0)
     # datasub_A0$ratio_w_ipcw_to_smallest <- datasub_A0$w_ipcw.pattern.bin.f/min(datasub_A0$w_ipcw.pattern.bin.f)
     # hist(log10(datasub_A0$ratio_w_ipcw_to_smallest))
     
     # Evaluate quality of unstabilized joint IPTW*IPCW weights
     ## -----------------------------------------------------------------------------------------------------------------------------
     # 28a. Summary
     # summary(datasub$w.f)
     # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
     #    1       2       3    2843      10 9786176
     
     # 28c. Ratio of each weight to the smallest weight
     # Joe suggests checking ratio of each weight to the smallest weight
     # datasub$ratio_w_to_smallest <- datasub$w.f / min(datasub$w.f)
     
     # 28d. Summary of ratios
     # summary(datasub$ratio_w_to_smallest)
     # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
     #    1       2       3    2621       9 9020818
     # Since the minimum ratio is 1, this tells us that this is the same as the weights, and we see that the summary for the ratios is the same as the one above
     
     # 28e. Check how many weights exceed 100
     # length(which(datasub$ratio_w_to_smallest > 100))
     # 5104
     # nrow(datasub)
     #33442
     # 5104/33442 = 15% extreme values
     
     # Check quality of stabilized joint IPTW*IPCW
     ## -----------------------------------------------------------------------------------------------------------------------------
     # 28a. Summary
     # summary(datasub$sw.f)
     #   Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
     # 0.0230    0.6930    0.8743    1.1235    0.9683 2240.4815
     
     # 28c. Ratio of each weight to the smallest weight
     # Joe suggests checking ratio of each weight to the smallest weight
     # datasub$ratio_sw_to_smallest <- datasub$sw.f / min(datasub$sw.f)
     
     # 28d. Summary of ratios
     # summary(datasub$ratio_sw_to_smallest)
     # Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
     # 1.00    30.12    38.00    48.84    42.09 97391.88
     # Since the minimum ratio is 1, this tells us that this is the same as the weights, and we see that the summary for the ratios is the same as the one above
     
     # 28e. Check how many weights exceed 100
     # length(which(datasub$ratio_sw_to_smallest > 100))
     # 1112
     # nrow(datasub)
     #33442
     # 1112/33442 = 3% extreme values
     
     
     # cat('summary of sw ratios:', '\n\n')
     # print(summary(datasub$ratio_sw_to_smallest))
     # cat('\n\n')
     
     
     # cat('nrow in final data set:', nrow(datasub), '\n')
     
     ### Checks for final data set that has only person-period assessments in which an outcome was measured
     # length(unique(datasub$IndividualId))
     # 
     # nrow(datasub)
     #
     
     ## -----------------------------------------------------------------------------------------------------------------------------
     #datasub_outcomes <- subset(datasub, !is.na(death_assessed))
     #length(unique(datasub_outcomes$IndividualId))
     # 
     #nrow(datasub_outcomes)
     # 
     
     ## -----------------------------------------------------------------------------------------------------------------------------
     # Check number of individuals at each visit
     #table(datasub_outcomes$visit_rank_DoD, exclude = NULL)
     
     # Check number of individuals exposed and unexposed at each visit
     #table(datasub_outcomes$exp_outmig_pattern_bin, datasub_outcomes$visit_rank_DoD, exclude = NULL)
     
     # Check number of outcomes at each visit
     # table(datasub_outcomes$death_assessed, datasub_outcomes$visit_rank_DoD, exclude = NULL)
   }
   
   #return new data object that includes new weights variables
   datasub_outcomes
 }
 
 
 
 