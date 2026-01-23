# Will not save these data sets in new repo
# Load averaged_plotpts outputs for each sex from 11_RunFunctions.R
load("men_averaged_plotpts.Rdata")
load("women_averaged_plotpts.Rdata")
# Load data with bootstrap replicates saved from running 13_RunBoot.R
boot_women <- readRDS('boot_women.Rds')
boot_men <- readRDS('boot_men.Rds')
source('14_BootCI.R')

# Load packages
pacman::p_load(foreach, 
               doParallel,
               tictoc,
               dplyr,
               ggplot2,
               spatstat.utils, # paren()
               stringr,
               scales          # package to tailor elements of ggplot
)

# Tables

####
# Women
####

# Select columns from data frame with predicted mean survival probabilities
women_predict_clean <- women_averaged_plotpts %>%
  select(visit_num, exp_factor, mean_survival)

# Reshape women_predict_clean from long to wide to prepare for merge
women_predict_wide <- women_predict_clean %>%
  tidyr::pivot_wider(id_cols = visit_num, 
                     names_from = exp_factor, 
                     values_from = mean_survival) %>%
  # Move the Resident column to come after the Migrant column so that the 
  # rr and rd columns make sense afterwards
  relocate(Resident, .after = Migrant) %>%
  # Calculate rr and rd
  mutate(RD = Migrant - Resident,
         RR = Migrant/Resident)

# Merge causal RR and RD with bootstrapped 95% CI around them
women_table <- merge(women_predict_wide, women_boot_curve_1999_RDRR_95CI, 
                     by = 'visit_num')

# Format the table
women_table_clean <- women_table %>%
  # Round estimates to two decimal places
  round(., 2) %>%
  # Put CI into one column and put parentheses around the values
  mutate(RDCI = paren(paste(rd_q025, rd_q975, sep = ','), type = '('),
         RRCI = paren(paste(rr_q025, rr_q975, sep = ','), type = '(')) %>%
  # Drop columns with only upper and lower bounds because now we have columns
  # with the full CIs
  select(-c(rd_q025, rd_q975, rr_q025, rr_q975)) %>%
  # Re-order columns
  relocate(RD, .before = RDCI) %>%
  relocate(RR, .before = RRCI) %>%
  # Rename new=old
  rename(Visit = visit_num)

women_table_clean <- print(women_table_clean, showAllLevels = TRUE, includeNA = TRUE)
write_xlsx(as.data.frame(women_table_clean), 'Table_WomenPredicted_02Apr2025.xlsx')


####
# Men
####

# Select columns from data frame with predicted mean survival probabilities
men_predict_clean <- men_averaged_plotpts %>%
  select(visit_num, exp_factor, mean_survival)

# Reshape women_predict_clean from long to wide to prepare for merge
men_predict_wide <- men_predict_clean %>%
  tidyr::pivot_wider(id_cols = visit_num, 
                     names_from = exp_factor, 
                     values_from = mean_survival) %>%
  # Move the Resident column to come after the Migrant column so that the 
  # rr and rd columns make sense afterwards
  relocate(Resident, .after = Migrant) %>%
  # Calculate rr and rd
  mutate(RD = Migrant - Resident,
         RR = Migrant/Resident)

# Merge the causal RR and RD with the bootstrap replicate 95% CI
men_table <- merge(men_predict_wide, men_boot_curve_1999_RDRR_95CI, 
                   by = 'visit_num')

# Format the table
men_table_clean <- men_table %>%
  # Round estimates to two decimal places
  round(., 2) %>%
  # Put CI into one column and put parentheses around the values
  mutate(RDCI = paren(paste(rd_q025, rd_q975, sep = ', '), type = '('),
         RRCI = paren(paste(rr_q025, rr_q975, sep = ', '), type = '(')) %>%
  # Drop columns with only upper and lower bounds because now we have columns
  # with the full CIs
  select(-c(rd_q025, rd_q975, rr_q025, rr_q975)) %>%
  # Re-order columns
  relocate(RD, .before = RDCI) %>%
  relocate(RR, .before = RRCI) %>%
  # Rename new=old
  rename(Visit = visit_num)

men_table_clean <- print(men_table_clean, showAllLevels = TRUE, includeNA = TRUE)
write_xlsx(as.data.frame(men_table_clean), 'Table_MenPredicted_02Apr2025.xlsx')

