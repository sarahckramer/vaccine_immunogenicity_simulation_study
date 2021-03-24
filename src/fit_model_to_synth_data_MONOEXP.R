# ---------------------------------------------------------------------------------------------------------------------
# Code to fit nonlinear mixed effects models to synthetic antibody kinetic data, both accounting and not accounting
# for seasonality in the degree of antibody response to vaccination
# ---------------------------------------------------------------------------------------------------------------------

# Setup ---------------------------------------------------------------------------------------------------------------

dev.off()

# Load libraries:
library(reshape2)
library(nlme)
library(ggplot2)
library(gridExtra)

# Set participant #'s/timepoints to test:
n_participants <- c(25, 50, 100, 250, 500, 1000, 2000, 5000)
select_timepoints <- c(379, 730, 1095, 2190)

# # Alternate timepoints to check:
# select_timepoints <- list(c(379, 730, 1095, 3650),
#                           c(379, 730, 1095),
#                           c(379, 730, 2190),
#                           c(379, 548, 730),
#                           c(379, 548, 1095),
#                           c(379, 730),
#                           c(379, 1095))

# Output plots/results as the code runs?:
noisy <- TRUE

# Read in functions ---------------------------------------------------------------------------------------------------

source('src/functions_ab_titers_over_time.R')
source('src/functions_assess_output.R')

# Read in/format data -------------------------------------------------------------------------------------------------

# Read in noise-laden "data":
ab_titers <- read.csv('data/synth_ab_kinetics/obs_data_MONO.csv')
ab_titers_ORIG <- ab_titers

# Read in month of vaccination:
vacc_month <- read.csv('data/synth_ab_kinetics/vacc_months_MONO.csv', header = FALSE)

# Loop through participant numbers and fit ----------------------------------------------------------------------------

for (n in n_participants) {
  print(n)
  
  # Set seed so that same subjects chosen each time:
  set.seed(3970396)
  
  # Reformat data frame:
  ab_titers <- format_synth_data(ab_titers_ORIG, vacc_month, n, select_timepoints)
  
  # # Plot data:
  # p1 <- ggplot(data = ab_titers, aes(x = time, y = value, color = subject)) + geom_line() + geom_point() +
  #   # geom_vline(xintercept = 0, lty = 2) +
  #   theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Ab Titers') +
  #   scale_x_continuous(breaks = seq(0, max(ab_titers$time), by = 30)) +
  #   scale_y_continuous(breaks = seq(0, max(ab_titers$value), by = 1.0))#, trans = 'log')
  # print(p1)
  
  # Try to first fit w/o seasonality:
  m1 <- try(
    nlme(log(value) ~ calculate_ab_titers_LOG_postOnly(time, log_beta, log_r_1, mono = TRUE),
         data = ab_titers,
         fixed = log_beta + log_r_1 ~ 1,
         random = pdDiag(log_beta + log_r_1 ~ 1),
         start = c(log_beta = log(18.0), log_r_1 = log(log(2)/365)))
  )
  
  # Check whether random effects associated with season of vaccination:
  if (class(m1)[1] != 'try-error' & noisy) {
    rand_effects <- ranef(m1)
    rand_effects$subject <- rownames(rand_effects)
    rand_effects <- merge(rand_effects, unique(ab_titers[, c('subject', 'vacc_month')]), by = 'subject')
    plot(rand_effects) # generally, clear seasonal patterns in beta estimates
  }
  
  # Now fit seasonal model:
  m2 <- try(
    nlme(log(value) ~ calculate_ab_titers_LOG_postOnly_seasonal(time, vacc_month, log_beta0,
                                                                logit_beta1, phi_hat, log_r_1,
                                                                mono = TRUE),
         data = ab_titers,
         fixed = log_beta0 + logit_beta1 + phi_hat + log_r_1 ~ 1,
         random = pdDiag(log_beta0 + log_r_1 ~ 1),
         start = c(log_beta0 = log(18.0), logit_beta1 = qlogis(0.1), phi_hat = qlogis(1/12),
                   log_r_1 = log(log(2)/365)))
  )
  
  # If fit fails, try selecting another subset of subjects:
  counter = 0
  while (class(m2)[1] == 'try-error' & counter < 10) {
    counter <- counter + 1
    print(counter)
    
    ab_titers <- format_synth_data(ab_titers_ORIG, vacc_month, n, select_timepoints)
    
    m2 <- try(
      nlme(log(value) ~ calculate_ab_titers_LOG_postOnly_seasonal(time, vacc_month, log_beta0,
                                                                  logit_beta1, phi_hat, log_r_1,
                                                                  mono = TRUE),
           data = ab_titers,
           fixed = log_beta0 + logit_beta1 + phi_hat + log_r_1 ~ 1,
           random = pdDiag(log_beta0 + log_r_1 ~ 1),
           start = c(log_beta0 = log(18.0), logit_beta1 = qlogis(0.1), phi_hat = qlogis(1/12),
                     log_r_1 = log(log(2)/365)))
    )
  }
  
  # If more than 10 extra trys, give up and move to next iteration:
  if (class(m2)[1] == 'try-error') next
  
  # Review model diagnostics:
  if (noisy) {
    print(plot(m2))
    print(qqnorm(m2, abline = c(0, 1)))
    print(pairs(m2))
    if (class(m1)[1] != 'try-error' & counter == 0) print(anova(m1, m2))
  }
  
  # Do we need both random effects?:
  m3 <- update(m2, random = pdDiag(log_beta0 ~ 1))
  m4 <- update(m2, random = pdDiag(log_r_1 ~ 1))
  if (noisy) {
    print(anova(m2, m3))
    print(anova(m2, m4))
  }
  
  # Extract parameter fits and random effect sds:
  results.df <- get_param_est(m2, seasonal = TRUE, mono = TRUE)
  print(results.df)
  
  if (!dir.exists('results/')) {
    dir.create('results/')
  }
  if (!dir.exists('results/ab_kinetics_res/')) {
    dir.create('results/ab_kinetics_res/')
  }
  
  write.csv(results.df, file = paste0('results/ab_kinetics_res/res_n', n, '_MONO.csv'))
  
  # Clean up:
  rm(m1, m2, m3, m4)
}

# Clean up ------------------------------------------------------------------------------------------------------------

rm(list = ls())
# dev.off()
