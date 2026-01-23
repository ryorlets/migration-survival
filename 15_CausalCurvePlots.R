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
               paletteer,
               spatstat.utils, # paren()
               stringr,
               scales          # package to tailor elements of ggplot
)

# Plots
# Women

# Merge together women data sets for the causal survival curves
  # women_boot_1999_95CI, which has the bootstrapped 95% confidence intervals 
  # women_predict, which has the time-specific predicted point estimates

# First, clean women_predict to have the same columns as women_boot_1999_95CI
# Remember now that we have imputed data sets, we have 'women_averaged_plotpts'
# instead of 'women_predict' for our sex-specific, time-specific, 
# exposure-specific cumulative survival probabilities
women_predict_curve <- women_averaged_plotpts %>%
  # Select columns we need for calculating RDs and RRs
  select(visit_num, exp_factor, mean_survival) 

# Merge women_predict_curve and women_boot_curve_1999_95CI
## First, make sure that variable names match up for merge
women_predict_curve_clean <- women_predict_curve %>%
  mutate(exp_num = ifelse(exp_factor == "Resident", 0, 1))

## Now merge
women_curve <- merge(women_predict_curve_clean, women_boot_curve_1999_95CI, 
                     by = c("visit_num", "exp_num")) %>%
                     arrange(visit_num)

# Plot the causal cumulative survival curve data
women_plot <- ggplot(women_curve, 
                     aes(x = visit_num, y = mean_survival)) + 
              geom_line(aes(colour = exp_factor)) +
              geom_point(aes(colour = exp_factor))+
              geom_errorbar(aes(ymin = q025, ymax = q975, color = exp_factor),
                            width = 0.2,
                            position = position_dodge(width = 0.15)) +
              # Labels for axes
              labs(x = "Six-month study visits", 
                   y = str_wrap("Cumulative survival probability", width = 10), 
                   colour = '') +
              # Set limits and intervals for x- and y-axes
              scale_x_continuous(limits = c(0, 10.5), breaks = seq(from=0, to=10.5, by=1)) +
              # the next line eliminates labels beyond visit 10
              #scale_x_continuous(expand = c(0,0)) +
              scale_y_continuous(limits = c(0.7, 1), breaks = seq(0.7, 1, 0.05)) +
              # plot title
              # ggtitle("Figure 3. Women's survival weighted for time-varying confounding and selection bias") +
              # Specify theme for plot
              theme_light() +
              #theme(text = element_text(family = "Times New Roman")) +
              # Format legend
              theme(legend.position = "bottom",
                    text = element_text(family = "serif"),
                    axis.title.y = element_text(angle = 0, vjust = 0.5, 
                    # Add space between y-axis label and y-axis                            
                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
                    axis.title.x = element_text(margin = 
                                                  margin(t = 10, r = 0, b = 0, l = 0))) + 
              # Specify colors
              scale_colour_manual(values = c('#ED8B00FF', '#00496FFF'))
              # scale_colour_manual(values = c('#D72000FF', '#FFAD0AFF'))
                #+
              # If we want to facet wrap by sex (using full data set, not sex-specific data)
              # facet_wrap(~ Sex) +
              # format labels that go with facet_wrap
              # theme(
              #  strip.text.x = element_text(
              #    color = 'black',
              #    # Increase font size of facet wrap labels
              #    size = 11.5),
              #  strip.background = element_rect(fill = 'white')
              #)

women_plot

#png(filename = "PerProtocolSurv.png", width = 2*1060, height = 2*1024, units = 'px', res = 72*5)
#p4
#dev.off()

women_plot15 <- ggplot(women_curve, 
                       aes(x = visit_num, y = mean_survival)) + 
  geom_line(aes(colour = exp_factor)) +
  geom_point(aes(colour = exp_factor))+
  geom_errorbar(aes(ymin = q025, ymax = q975, color = exp_factor),
                width = 0.2,
                position = position_dodge(width = 0.15)) +
  # Labels for axes
  labs(x = "Six-month study visits", 
       y = str_wrap("Cumulative survival probability", width = 10), 
       colour = '') +
  # Set limits and intervals for x- and y-axes
  scale_x_continuous(limits = c(0, 15.5), breaks = seq(from=0, to=15.5, by=1)) +
  # the next line eliminates labels beyond visit 10
  #scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0.7, 1), breaks = seq(0.7, 1, 0.05)) +
  # Specify theme for plot
  theme_light() +
  # Format legend
  theme(legend.position = "bottom",
        text = element_text(family = "serif"),
        axis.title.y = element_text(angle = 0, vjust = 0.5, 
                                    # Add space between y-axis label and y-axis                            
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = 
                                      margin(t = 10, r = 0, b = 0, l = 0))) + 
  # Specify colors
  scale_colour_manual(values = c('#ED8B00FF', '#00496FFF'))
# If we want to facet wrap by sex (using full data set, not sex-specific data)
# facet_wrap(~ Sex) +
# format labels that go with facet_wrap
# theme(
#  strip.text.x = element_text(
#    color = 'black',
#    # Increase font size of facet wrap labels
#    size = 11.5),
#  strip.background = element_rect(fill = 'white'))
women_plot15

# Men

# Merge together men data sets for the causal survival curves
# men_boot_1999_95CI, which has the bootstrapped 95% confidence intervals 
# men_predict, which has the time-specific predicted point estimates

# First, clean men_predict to have the same columns as men_boot_1999_95CI
men_predict_curve <- men_averaged_plotpts %>%
  # Select columns we need for calculating RDs and RRs
  select(visit_num, exp_factor, mean_survival) 

# Merge men_predict_curve and men_boot_curve_1999_95CI
## First, make sure that variable names match up for merge
men_predict_curve_clean <- men_predict_curve %>%
  mutate(exp_num = ifelse(exp_factor == "Resident", 0, 1))

## Now merge
men_curve <- merge(men_predict_curve_clean, men_boot_curve_1999_95CI, 
                     by = c("visit_num", "exp_num")) %>%
                     arrange(visit_num)

men_plot <- ggplot(men_curve, 
                     aes(x = visit_num, y = mean_survival)) + 
            geom_line(aes(colour = exp_factor)) +
            geom_point(aes(colour = exp_factor))+
            geom_errorbar(aes(ymin = q025, ymax = q975, color = exp_factor),
                width = 0.2,
                position = position_dodge(width = 0.15)) +
            # Labels for axes
            labs(x = "Six-month study visits", 
            y = str_wrap("Cumulative survival probability", width = 10), 
            colour = '') +
            # Set limits and intervals for x- and y-axes
            scale_x_continuous(limits = c(0, 2.5), breaks = seq(from=0, to=2.5, by=1)) +
            scale_y_continuous(limits = c(0.7, 1), breaks = seq(0.7, 1, 0.05)) +
            # plot title
            #ggtitle("Figure 3. Men's survival weighted for time-varying confounding and selection bias") +
            theme_light() +
            theme(legend.position = "bottom",
                  text = element_text(family = "Times New Roman"),
                  axis.title.y = element_text(angle = 0, vjust = 0.5, 
                                              # Add space between y-axis label and y-axis                            
                                              margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.title.x = element_text(margin = 
                                                margin(t = 10, r = 0, b = 0, l = 0))) +
            scale_colour_manual(values = c('#ED8B00FF', '#00496FFF'))

men_plot

# Causal survival curve with full follow-up time, ignoring data scarcity
men_plot15 <- ggplot(men_curve, 
                       aes(x = visit_num, y = mean_survival)) + 
  geom_line(aes(colour = exp_factor)) +
  geom_point(aes(colour = exp_factor))+
  geom_errorbar(aes(ymin = q025, ymax = q975, color = exp_factor),
                width = 0.2,
                position = position_dodge(width = 0.15)) +
  # Labels for axes
  labs(x = "Six-month study visits", 
       y = str_wrap("Cumulative survival probability", width = 10), 
       colour = '') +
  # Set limits and intervals for x- and y-axes
  scale_x_continuous(limits = c(0, 15.5), breaks = seq(from=0, to=15.5, by=1)) +
  # the next line eliminates labels beyond visit 10
  #scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0.0, 1), breaks = seq(0.0, 1, 0.10)) +
  # Specify theme for plot
  theme_light() +
  # Format legend
  theme(legend.position = "bottom",
        text = element_text(family = "serif"),
        axis.title.y = element_text(angle = 0, vjust = 0.5, 
                                    # Add space between y-axis label and y-axis                            
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = 
                                      margin(t = 10, r = 0, b = 0, l = 0))) + 
  # Specify colors
  scale_colour_manual(values = c('#ED8B00FF', '#00496FFF'))

men_plot15


