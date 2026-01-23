# File: 06_PreWeightsClean.R
# Author: Rachel Yorlets
# Update date: 10 Apr 2025
# Purpose: Clean fully merged surveillance + viral + cd4 data restricted to 
# analytic sample in order to prepare for calculating stabilized 
# person-period-specific inverse probability of treatment and censoring weights

# Use 'window' baseline CD4 as the cd4 for visit==1 

## -----------------------------------------------------------------------------
pacman::p_load(haven, 
               dplyr,
               lubridate,
               Hmisc,       # Lag()
               broom,       # tidy
               tableone,    # CreateTableOne
               DescTools,   # use LOCF function
               finalfit,    # missing_plot()
               ggplot2,
               mgcv         # gam() for generalized additive model
               )

## -----------------------------------------------------------------------------
load('Step05_OutcomeSTROBE.Rdata')
# Origin: Created by running 05_OutcomeSTROBE coding file
# Description: person-assessments for participants included in the analytic data 
# set (n=5052) with censored== assigned
# Rows: Rows that contain assessments (includes EarliestARTInitDate and 
# 6-month visits) only
# Columns: IndividualId, EndDate_PT, VisitDate_DoD, visit_rank_DoD, Died_Assess, 
# DoD, time-varying covariates summarized across 6-month time period, censored

## -----------------------------------------------------------------------------
df2 <- Step05_OutcomeSTROBE
#-------------------------------------------------------------------------------
# Check number of rows and person-period assessments in analytic data set
#-------------------------------------------------------------------------------
length(unique(df2$IndividualId))
# 5052 people in the analytic data set
nrow(df2)
# 44708 person-period assessments (ie, for all 5052 persons in the analytic 
# sample after post-censoring and post-event rows have been removed)
# View(df2[ ,c('IndividualId', 'visit_rank_DoD', 'exp_outmig_pattern_bin')])
#-------------------------------------------------------------------------------
# Rename visit_rank_DoD
#-------------------------------------------------------------------------------
df2 <- df2 %>%
  rename(visit = visit_rank_DoD)
df2$visit <- as.factor(df2$visit)
#-------------------------------------------------------------------------------
# Review distribution of sex (effect measure modifier) in analytic data set
#-------------------------------------------------------------------------------
df2_m <- subset(df2, sex == 'Male')
length(unique(df2_m$IndividualId))
# 1322 persons in analytic sample were men

df2_f <- subset(df2, sex == 'Female')
length(unique(df2_f$IndividualId))
# nrow(df1_f)
# 3730 persons in analytic sample were women (74% women)
#-------------------------------------------------------------------------------
# Exposure 
#-------------------------------------------------------------------------------
# Exposure - format of time-varying exposure
#-------------------------------------------------------------------------------
# Create migration exposure as a binary variable for weights
df2 <- df2 %>%
  group_by(IndividualId) %>%
  # Remember outmig_pattern is polytomous numeric variable with levels 0, 1, 2
  mutate(exp = ifelse(exp > 0, 1, exp))

df2$exp <- as.factor(df2$exp)

# Re-run sex-subsets so that we can look at sex-specific exposure distributions
df2_m <- subset(df2, sex == 'Male')
df2_f <- subset(df2, sex == 'Female')

# Evaluate exposure by visit by sex
# Among women
table(df2_f$exp, df2_f$visit)
#     1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
# 0 2903 2799 2668 2578 2440 2271 2014 1701 1425 1169  977  790  630  441  242   70
# 1  827  761  779  778  803  779  721  662  564  477  391  295  224  157   84   22

length(unique(df2_f$IndividualId))
# 3730
nrow(df2_f)
# 33442

# Among men
table(df2_m$exp, df2_m$visit)
#    1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
# 0 954 901 859 819 784 723 640 544 437 379 342 283 206 138  78  22
# 1 368 321 312 304 289 275 246 218 197 172 138 120  93  63  30  11

length(unique(df2_m$IndividualId))
# 1322
nrow(df2_m)
# 11266
#--------------------------------------------------------------------------------
# Exposure - lagged time-varying exposure (a covariate)
#--------------------------------------------------------------------------------
# Because we will regress exposure on the exposure that happened in the prior 
# assessment, create a lagged version of the exposure variable 
df2 <- df2 %>%
  arrange(IndividualId, visit) %>%
  # We use lag() from dplyr, and set the default to 0, meaning that all baseline 
  # person-assessments for which there is no lagged exposure by definition get 
  # assigned a '0'
  # We do this because participants enter the study 'treatment-naive' in theory
  # Problem/complexity with a prevalent user design is that you know they aren't 
  # all unexposed at baseline
  mutate(exp_lag = lag(exp, default = '0')) 

# Assign it as a factor variable
df2$exp_lag <- as.factor(df2$exp_lag)
# View(df2[ ,c('IndividualId', 'visit', 'exp', 'exp_lag')])
#-------------------------------------------------------------------------------
# Exposure - history in year before ART initiation, polytomous 3-level variable
#-------------------------------------------------------------------------------
df2$exp_hx <- as.factor(df2$exp_hx)
#-------------------------------------------------------------------------------
# Exposure - Create 'ever exposed' binary variable
#-------------------------------------------------------------------------------
df2 <- df2 %>%
  arrange(IndividualId, visit) %>%
  group_by(IndividualId) %>%
  # if an individual had one or more person-period assessments in which they
  # had out-migrated, they have a history of migration, or 'mig_ever == 1'
  mutate(exp_ever = ifelse(length(which(exp == 1)) >= 1, 1, 0))

# View(df2[, c('IndividualId', 'sex', 'visit', 'exp', 'exp_ever')]) 
#-------------------------------------------------------------------------------
# Check baseline missingness before fitting weights in next file
#-------------------------------------------------------------------------------
df2_baseline <- subset(df2, visit == 1)
# nrow(df2_baseline)
# 5052

# Check how many people were missing any non-CD4 covariate at baseline
length(which(is.na(df2_baseline$exp_hx) 
             | is.na(df2_baseline$children) 
             | is.na(df2_baseline$age)
             | is.na(df2_baseline$sex) 
             | is.na(df2_baseline$sei)
             | is.na(df2_baseline$edu) 
             | is.na(df2_baseline$cohabit)))
# 880
880/length(unique(df2_baseline$IndividualId))
# 17.4% of participants are missing one or more covariates at baseline

# We can also do it this way
with(df2_baseline, length(which(is.na(exp_hx) 
                                | is.na(children) 
                                | is.na(age)
                                | is.na(sex) 
                                | is.na(sei) 
                                | is.na(edu) 
                                | is.na(cohabit))))
# 880
#-------------------------------------------------------------------------------
# Among men
with(df2_baseline %>%
       filter(sex == 'Male'),
     length(which(is.na(exp_hx) | is.na(children) 
                  | is.na(age)| is.na(sex) 
                  | is.na(sei) | is.na(edu) 
                  | is.na(cohabit)))
)
# 223
#-------------------------------------------------------------------------------
# Among women
with(df2_baseline %>%
       filter(sex == 'Female'),
     length(which(is.na(exp_hx) | is.na(children) 
                  | is.na(age)| is.na(sex) 
                  | is.na(sei) | is.na(edu) 
                  | is.na(cohabit)))
)
# 657
#-------------------------------------------------------------------------------
# Check missingness at any time point, including baseline
#-------------------------------------------------------------------------------
length(which(is.na(df2$exp_hx) 
             | is.na(df2$children) 
             | is.na(df2$age)
             | is.na(df2$sex) 
             | is.na(df2$sei) 
             | is.na(df2$edu) 
             | is.na(df2$cohabit)))
# 5200 person-period assessments are missing 
#-------------------------------------------------------------------------------
# Check missingness of individual variables
#-------------------------------------------------------------------------------
# Exposure (out-migration pattern as a binary variable)
sum(is.na(df2$exp))
# 0, which is expected because there is no missingness in Resident, which was 
# used to create the exposure

# Exposure: nonres_duration_days (number of days spent away from home)
sum(is.na(df2$nonres_duration_mos))
# 0

# Exposure history in year before ART initiation, 3-level polytomous variable
sum(is.na(df2$exp_hx))
# 0

# Children
sum(is.na(df2$children))
# 0

# Age
sum(is.na(df2$age))
# 0

# Sex (numeric version where 1==male, 2==female)
sum(is.na(df2$sex))
# 0

# SEIdx
length(which(is.na(df2$sei)))
# 2453 missing

length(which(is.na(df2$sei) & (df2$visit == 1)))
# 368 are missing baseline SEIdx
# 368/5052 = 7.2% missingness at baseline

# among women
length(which(is.na(df2$sei) & (df2$visit == 1) & df2$sex == 'Female'))
# 272
272/(length(which(df2$sex == 'Female' & (df2$visit == 1))))
# 272 women are missing SEI at baseline, which is 7.3% of women at baseline

# among men
length(which(is.na(df2$sei) & df2$sex == 'Male' & (df2$visit == 1)))
# 96
96/(length(which(df2$sex == 'Male' & df2$visit == 1)))
# 96 men are missing SEI at baseline, which is 7.3% of men at baseline

# Joe raised that taking the median over all values, including future values of 
# SEI means that you could technically be conditioning on a future value of a 
# covariate, which induces bias, so we won't use the median

# We will impute missing SEI

# Education
length(which(is.na(df2$edu)))
# 4051
# df2$edu_max[is.na(df2$edu_max)] <- 'missing'
# We will impute missing education

# Cohabitation
length(which(is.na(df2$cohabit)))
# 1249
# We will impute missing cohabitation

# Suppression ever
length(which(is.na(df2$suppression)))
# 44483
table(df2$suppression, exclude = NULL)
#   0       1 NA 
# 108     117   44483 
# Since we have 44708 assessment-periods, we see that almost all (44483) are NA
# We likely wouldn't want to use suppression_ever in the weights anyway 
# because it could be collinear and a mediator for people who die of HIV/AIDS

# cd4
# Remember that in 02_DerivingSample.Rmd, we identified the cd4 count that is closest
# to the earliest known ART initiation date, allowing a window of 12 months before ART 
# start and 6 months after.
# Now, we need to assign that baseline CD4 value as the cd4 value for visit==1 
# because we don't need to impute that value later.
df3 <- df2 %>%
  group_by(IndividualId) %>%
  mutate(cd4count = ifelse(visit == 1, cd4_base, cd4))

#------------------------------------------------------------------------------
# Save data set with cleaned covariates, ready for weights calculations
Step06_PreWeightsClean <- df3
save(Step06_PreWeightsClean, file = 'Step06_PreWeightsClean.Rdata')
