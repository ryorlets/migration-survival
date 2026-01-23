# Load data saved from running 13_RunBoot.R
boot_women <- readRDS('boot_women.Rds')
boot_men <- readRDS('boot_men.Rds')

# Load packages
pacman::p_load(foreach, 
               doParallel,
               tictoc,
               dplyr,
               ggplot2,
               stringr     # str_wrap() for axis labels in ggplot
)

# Calculate the 95% CI for the MSM coefficients
####
# Template
# lapply() calls for a function for its second argument
# Here, we write an anonymous function because we do not need to execute another function, just extract the model output from a given bootstrap replicate
# boot_women[[1]] for example allows us to see the list of the msm (ie, the MSM model output, exponentiated) and the pred 
#boot_coef_1999 <- Reduce(rbind, lapply(boot_women, function(x) x[['msm']])) 

#boot_coef_1999_95CI <- 
#  boot_coef_1999 %>%
#  group_by(term) %>%
#  summarize(q025 = quantile(estimate, 0.025),
#            q975 = quantile(estimate, 0.975))

# Calculate the 95% CI for the survival curves
####
# Template
# boot_curve_1999 <- rbind(lapply(boot_women, function(x) x$pred)) %>%
#  group_by(exposure, t) %>%
#  mutate(q025 = quantile(cumsurv, 0.025),
#         q975 = quantile(cumsurv, 0.975))
####
# Jason suggests changing this to boot_women[1:999] to see if the estimates are different

#######
# Women
######
# Subset (reduce) and stack (rbind) the list boot_women to lists that include 
# the predictions from the causal curves only
women_boot_curve_1999 <- Reduce(rbind,
                                # First, get the number of each list within 
                                # boot_women (which is the number of bootstrap 
                                # replicates)
                                lapply(seq_len(length(boot_women)), 
                                       function(x) {
                                         # Extract the prediction list from the 
                                         # first boot_women list (which is the 
                                         # first bootstrap replicate)
                                         boot_women[[x]][['pred']] %>%
                                           # Label each returned data frame as 
                                           # the bootstrap replicate from which 
                                           # it was generated
                                           mutate(B = x,
                                                  Migrant = 1, 
                                                  Resident = 0)
                                       }))

# Check distribution of bootstrap replicate estimates of cumulative probability
# of survival for women
SupplBoxplot_BootstrapCumSurvProbWomen <- ggplot(data = women_boot_curve_1999,
                                 aes(x = factor(visit_num), y = mean_survival)) +
  geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.5) + 
  # Labels for axes
  labs(x = "Six-month study visits", 
       y = str_wrap("Cumulative survival probability", width = 10), 
       colour = '') +
  theme_light() +
  # Format legend
  theme(legend.position = "bottom",
        text = element_text(family = "Times New Roman"),
        axis.title.y = element_text(angle = 0, vjust = 0.5, 
                                    # Add space between y-axis label and y-axis                            
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = 
                                      margin(t = 10, r = 0, b = 0, l = 0))) +
  facet_wrap(vars(exp_num))

SupplBoxplot_BootstrapCumSurvProbWomen

# 95% CI around the bootstrap replicate causal curves 
women_boot_curve_1999_95CI <- women_boot_curve_1999 %>%
  group_by(exp_num, visit_num) %>%
  summarize(q025 = quantile(mean_survival, 0.025),
            q975 = quantile(mean_survival, 0.975))

# Calculation of RR and RD within bootstrap replicates
women_boot_curve_1999_RDRR <- women_boot_curve_1999 %>%
  # Select columns we need for calculating RDs and RRs
  select(B, visit_num, exp_num, mean_survival) %>%
  # seq_len gives you a sequence of integers from 1 to the n observations
  #mutate(rowid = seq_len(n())) %>%
  # transform from long to wide
  tidyr::pivot_wider(id_cols = c(B, visit_num), 
                     names_from = exp_num, 
                     values_from = mean_survival) %>%
  # Calculate risk of survival at each visit
  # Remember the probabilities are already survival probabilities (not death)
  mutate(Migrant = `1`, 
         Resident = `0`,
         rd = Migrant - Resident,
         rr = Migrant/Resident,
         # Create the visit number as a factor variable to number x-axis
         visit = as.factor(visit_num)) 

# make box plots of rd and rr distribution across bootstrap replicates
SupplBoxplot_BootstrapRRWomen <- ggplot(data = women_boot_curve_1999_RDRR,
          aes(group = visit_num, x = visit, y = rr)) +
  geom_boxplot(outlier.size = 0.5) + 
  # Labels for axes
  labs(x = "Six-month study visits", 
       y = str_wrap("Causal risk ratio", width = 10), 
       colour = '') +
  theme_light() +
  # Format legend
  theme(legend.position = "bottom",
        text = element_text(family = "Times New Roman"),
        axis.title.y = element_text(angle = 0, vjust = 0.5, 
                                    # Add space between y-axis label and y-axis                            
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = 
                                      margin(t = 10, r = 0, b = 0, l = 0))) 

# Check distribution of causal risk differences across bootstrap replicates
SupplBoxplot_BootstrapRDWomen <- ggplot(data = women_boot_curve_1999_RDRR,
                                        aes(group = visit_num, 
                                            x = visit, y = rd)) +
  geom_boxplot(outlier.size = 0.5) + 
  # Labels for axes
  labs(x = "Six-month study visits", 
       y = str_wrap("Causal risk difference", width = 10), 
       colour = '') +
  theme_light() +
  # Format legend
  theme(legend.position = "bottom",
        text = element_text(family = "Times New Roman"),
        axis.title.y = element_text(angle = 0, vjust = 0.5, 
                                    # Add space between y-axis label and y-axis                            
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = 
                                      margin(t = 10, r = 0, b = 0, l = 0))) 

# Calculate 95% CI around bootstrap replicate RRs and RDs
women_boot_curve_1999_RDRR_95CI <- women_boot_curve_1999_RDRR %>%
  group_by(visit_num) %>%
  # Calculate the quantiles around the RRs and RDs
  # Use summarize because we want to summarize the quantiles
  # Summarize means we can't retain the rd and rr columns, so we will merge back
  summarize(rd_q025 = quantile(rd, 0.025),
            rd_q975 = quantile(rd, 0.975),
            rr_q025 = quantile(rr, 0.025),
            rr_q975 = quantile(rr, 0.975)) 

#######
# Men
######
# Subset (reduce) and stack (rbind) the list boot_women to lists that include the predictions from the causal curves only
men_boot_curve_1999 <- Reduce(rbind,
                              # First, get the number of each list within boot_men (which is the number of bootstrap replicates)
                              lapply(seq_len(length(boot_men)), 
                                     function(x) {
                                       # Extract the prediction list from the first boot_men list (which is the first bootstrap replicate)
                                       boot_men[[x]][['pred']] %>%
                                         # Label each returned data frame as the bootstrap replicate from which it was generated
                                         mutate(B = x)
                                     }))

# Check distribution of bootstrap replicate estimates of cumulative probability
# of survival for women
SupplBoxplot_BootstrapCumSurvProbMen <- ggplot(data = men_boot_curve_1999,
                                                 aes(x = visit_num, y = mean_survival)) +
  geom_boxplot(outlier.size = 0.5) + 
  # Labels for axes
  labs(x = "Six-month study visits", 
       y = str_wrap("Cumulative survival probability", width = 10), 
       colour = '') +
  theme_light() +
  # Format legend
  theme(legend.position = "bottom",
        text = element_text(family = "Times New Roman"),
        axis.title.y = element_text(angle = 0, vjust = 0.5, 
                                    # Add space between y-axis label and y-axis                            
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = 
                                      margin(t = 10, r = 0, b = 0, l = 0))) +
  facet_wrap(vars(exp_num))

# 95% CI around the bootstrap replicate causal curves 
men_boot_curve_1999_95CI <- men_boot_curve_1999 %>%
  group_by(exp_num, visit_num) %>%
  summarize(q025 = quantile(mean_survival, 0.025),
            q975 = quantile(mean_survival, 0.975))

# Calculation of RR and RD within bootstrap replicates
men_boot_curve_1999_RDRR <- men_boot_curve_1999 %>%
  # Select columns we need for calculating RDs and RRs
  select(B, visit_num, exp_num, mean_survival) %>%
  # seq_len gives you a sequence of integers from 1 to the n observations
  #mutate(rowid = seq_len(n())) %>%
  # transform from long to wide
  tidyr::pivot_wider(id_cols = c(B, visit_num), 
                     names_from = exp_num, 
                     values_from = mean_survival) %>%
  # Calculate within each visit
  mutate(Migrant = `1`,
         Resident = `0`,
         rd = Migrant - Resident,
         rr = Migrant/Resident,
         visit = as.factor(visit_num)) 

# make box plots of rd and rr distribution across bootstrap replicates
SupplBoxplot_BootstrapRRMen <- ggplot(data = men_boot_curve_1999_RDRR,
                                        aes(group = visit_num, x = visit, y = rr)) +
  geom_boxplot(outlier.size = 0.5) + 
  # Labels for axes
  labs(x = "Six-month study visits", 
       y = str_wrap("Causal risk ratio", width = 10), 
       colour = '') +
  theme_light() +
  # Format legend
  theme(legend.position = "bottom",
        text = element_text(family = "Times New Roman"),
        axis.title.y = element_text(angle = 0, vjust = 0.5, 
                                    # Add space between y-axis label and y-axis                            
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = 
                                      margin(t = 10, r = 0, b = 0, l = 0))) 

# Check distribution of causal risk differences across bootstrap replicates
SupplBoxplot_BootstrapRDMen <- ggplot(data = men_boot_curve_1999_RDRR,
                                        aes(group = visit_num, 
                                            x = visit, y = rd)) +
  geom_boxplot(outlier.size = 0.5) + 
  # Labels for axes
  labs(x = "Six-month study visits", 
       y = str_wrap("Causal risk difference", width = 10), 
       colour = '') +
  theme_light() +
  # Format legend
  theme(legend.position = "bottom",
        text = element_text(family = "Times New Roman"),
        axis.title.y = element_text(angle = 0, vjust = 0.5, 
                                    # Add space between y-axis label and y-axis                            
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = 
                                      margin(t = 10, r = 0, b = 0, l = 0))) 

# Calculate the 95% CI around the bootstrap replicate RR and RD
men_boot_curve_1999_RDRR_95CI <- men_boot_curve_1999_RDRR %>%
  group_by(visit_num) %>%
  # Calculate the quantiles around the RRs and RDs
  # Use summarize because we want to summarize over all the bootstrap replicates
  summarize(rd_q025 = quantile(rd, 0.025),
            rd_q975 = quantile(rd, 0.975),
            rr_q025 = quantile(rr, 0.025),
            rr_q975 = quantile(rr, 0.975))
