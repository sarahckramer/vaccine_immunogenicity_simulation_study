# ---------------------------------------------------------------------------------------------------------------------
# Code to fit nonlinear mixed effects models to synthetic antibody kinetic data, both accounting and not accounting
# for seasonality in the degree of antibody response to vaccination
# ---------------------------------------------------------------------------------------------------------------------

# Setup ---------------------------------------------------------------------------------------------------------------

# Load libraries:
library(reshape2)
library(nlme)
library(ggplot2)
library(gridExtra)

# Specify date (for results outputs):
ymd <- '20210317'

# Set participant #'s/timepoints to test:
n_participants <- c(25, 50, 100, 250, 500, 1000, 2000, 5000)
select_timepoints <- c(379, 730, 1095, 2190)

# Read in functions ---------------------------------------------------------------------------------------------------

source('src/functions_ab_titers_over_time.R')
source('src/functions_assess_output.R')

# Read in/format data -------------------------------------------------------------------------------------------------

# Read in noise-laden "data":
ab_titers <- read.csv('data/prelim_check_20210223/obs_data.csv')
ab_titers_ORIG <- ab_titers

# Read in month of vaccination:
vacc_month <- read.csv('data/prelim_check_20210223/vacc_months.csv', header = FALSE)

# Loop through participant numbers and fit ----------------------------------------------------------------------------

for (n in n_participants) {
  print(n)
  
  # Set seed so that same subjects chosen each time:
  set.seed(3970395)
  
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
    nlme(log(value) ~ calculate_ab_titers_LOG_postOnly(time, log_beta, log_r_1, log_r_2, logit_rho),
         data = ab_titers,
         fixed = log_beta + log_r_1 + log_r_2 + logit_rho ~ 1,
         random = pdDiag(log_beta + log_r_1 + log_r_2 ~ 1),
         groups = ~subject,
         start = c(log_beta = log(18.0), log_r_1 = log(log(2)/30),
                   log_r_2 = log(log(2)/3650), logit_rho = qlogis(0.75)))
  )
  
  # Check whether random effects associated with season of vaccination:
  if (class(m1)[1] != 'try-error') {
    rand_effects <- ranef(m1)
    rand_effects$subject <- rownames(rand_effects)
    rand_effects <- merge(rand_effects, unique(ab_titers[, c('subject', 'vacc_month')]), by = 'subject')
    plot(rand_effects) # generally, clear seasonal patterns in beta estimates
  }
  
  # Now fit seasonal model:
  m2 <- try(
    nlme(log(value) ~ calculate_ab_titers_LOG_postOnly_seasonal(time, vacc_month, log_beta0, logit_beta1,
                                                                phi_hat, log_r_1, log_r_2, logit_rho),
         data = ab_titers,
         fixed = log_beta0 + logit_beta1 + phi_hat + log_r_1 + log_r_2 + logit_rho ~ 1,
         random = pdDiag(log_beta0 + log_r_1 + log_r_2 ~ 1),
         groups = ~subject,
         start = c(log_beta0 = log(18.0), logit_beta1 = qlogis(0.1), phi_hat = qlogis(1/12),
                   log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650), logit_rho = qlogis(0.75)))
  )
  
  # If fit fails, try selecting another subset of subjects:
  counter = 0
  while (class(m2)[1] == 'try-error' & counter < 10) {
    counter <- counter + 1
    print(counter)
    
    ab_titers <- format_synth_data(ab_titers_ORIG, vacc_month, n, select_timepoints)
    
    m2 <- try(
      nlme(log(value) ~ calculate_ab_titers_LOG_postOnly_seasonal(time, vacc_month, log_beta0, logit_beta1,
                                                                  phi_hat, log_r_1, log_r_2, logit_rho),
           data = ab_titers,
           fixed = log_beta0 + logit_beta1 + phi_hat + log_r_1 + log_r_2 + logit_rho ~ 1,
           random = pdDiag(log_beta0 + log_r_1 + log_r_2 ~ 1),
           groups = ~subject,
           start = c(log_beta0 = log(18.0), logit_beta1 = qlogis(0.1), phi_hat = qlogis(1/12),
                     log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650), logit_rho = qlogis(0.75)))
    )
  }
  
  # If more than 10 extra trys, give up and move to next iteration:
  if (class(m2)[1] == 'try-error') next
  
  # Review model diagnostics:
  plot(m2)
  qqnorm(m2, abline = c(0, 1))
  pairs(m2)
  if (class(m1)[1] != 'try-error' & counter == 0) anova(m1, m2)
  
  # Do we need random effect of r2?:
  m3 <- update(m2, random = pdDiag(log_beta0 + log_r_1 ~ 1))
  anova(m2, m3)
  
  # Get and output results ----------------------------------------------------------------------------------------------
  
  # Extract parameter fits and random effect sds:
  results.df <- get_param_est(m2, seasonal = TRUE)
  print(results.df)
  
  if (!dir.exists(paste0('results/PRELIM_nlme_res_', ymd, '/'))) {
    dir.create(paste0('results/PRELIM_nlme_res_', ymd, '/'))
  }
  write.csv(results.df, file = paste0('results/PRELIM_nlme_res_', ymd, '/res_n', n_participants, '.csv'),
            row.names = FALSE)
  
  # Clean up:
  rm(m1, m2, m3)
}

# Potentially useful references ---------------------------------------------------------------------------------------

# https://stackoverflow.com/questions/26449969/backward-selection-in-lme-singularity-in-backsolve-occured
# https://stackoverflow.com/questions/50505290/singularity-in-backsolve-at-level-0-block-1-in-lme-model

# Clean up ------------------------------------------------------------------------------------------------------------

rm(list = ls())
dev.off()
