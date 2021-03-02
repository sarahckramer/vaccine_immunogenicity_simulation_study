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
ymd <- '20210223'

# Set timepoints to select/test:
select_timepoints <- list(c(379, 730, 1095, 2190),
                          c(379, 730, 1095, 3650),
                          c(379, 730, 1095),
                          c(379, 730, 2190),
                          c(379, 548, 730),
                          c(379, 548, 1095),
                          c(379, 730),
                          c(379, 1095))

# Read in functions ---------------------------------------------------------------------------------------------------

source('src/functions_ab_titers_over_time.R')
source('src/functions_assess_output.R')

# Read in/format data -------------------------------------------------------------------------------------------------

# Read in noise-laden "data":
ab_titers <- read.csv('data/prelim_check_20210225/obs_data_MONO.csv')
ab_titers_ORIG <- ab_titers

# Read in month of vaccination:
vacc_month <- read.csv('data/prelim_check_20210225/vacc_months_MONO.csv', header = FALSE)

# Set seed so that same subjects chosen each time:
set.seed(3970395)

# Reformat data frame:
ab_titers <- format_synth_data(ab_titers_ORIG, vacc_month, 5000, select_timepoints[[1]])

# # Plot data:
# p1 <- ggplot(data = ab_titers, aes(x = time, y = value, color = subject)) + geom_line() + geom_point() +
#   # geom_vline(xintercept = 0, lty = 2) +
#   theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Ab Titers') +
#   scale_x_continuous(breaks = seq(0, max(ab_titers$time), by = 30)) +
#   scale_y_continuous(breaks = seq(0, max(ab_titers$value), by = 1.0))#, trans = 'log')
# print(p1)

# Fit w/o seasonality -------------------------------------------------------------------------------------------------

# First fit without random effects:
m1 <- nls(log(value) ~ calculate_ab_titers_LOG_postOnly_MONO(time, log_beta, log_r),
          data = ab_titers,
          start = c(log_beta = log(18.0), log_r = log(log(2)/365)))

# Next try nlsList (separate fits for each person):
m2 <- nlsList(log(value) ~ calculate_ab_titers_LOG_postOnly_MONO(time, log_beta, log_r) | subject,
              data = ab_titers,
              start = c(log_beta = log(18.0), log_r = log(log(2)/365)))
# plot(intervals(m2))
# plot(m2, subject ~ resid(.), abline = 0 )
pairs(m2)

# Fit nlme:
m3 <- nlme(log(value) ~ calculate_ab_titers_LOG_postOnly_MONO(time, log_beta, log_r),
           data = ab_titers,
           fixed = log_beta + log_r ~ 1,
           random = pdDiag(log_beta + log_r ~ 1),
           groups = ~subject,
           start = c(log_beta = log(18.0), log_r = log(log(2)/365)))
plot(m3)
pairs(m3)
qqnorm(m3, abline = c(0, 1))
# qqnorm(m3, ~ ranef(.))

# Check whether random effects associated with season of vaccination:
rand_effects <- ranef(m3)
rand_effects$subject <- rownames(rand_effects)
rand_effects <- merge(rand_effects, unique(ab_titers[, c('subject', 'vacc_month')]), by = 'subject')
plot(rand_effects) # clear seasonal patterns in beta estimates

# Fit model including seasonality -------------------------------------------------------------------------------------

# First without random effects:
m4 <- nls(log(value) ~ calculate_ab_titers_LOG_postOnly_seasonal_MONO(time, vacc_month, log_beta0, logit_beta1,
                                                                      phi_hat, log_r),
          data = ab_titers,
          start = c(log_beta0 = log(18.0), logit_beta1 = qlogis(0.15), phi_hat = qlogis(2/12),
                    log_r = log(log(2)/60)))

# m5 <- nlsList(log(value) ~ calculate_ab_titers_LOG_postOnly_seasonal_MONO(time, vacc_month, log_beta0, logit_beta1,
#                                                                           phi_hat, log_r) | subject,
#               data = ab_titers,
#               start = c(log_beta0 = log(18.0), logit_beta1 = qlogis(0.15), phi_hat = qlogis(2/12),
#                         log_r = log(log(2)/365)))
# plot(intervals(m5))
# pairs(m5)

m6 <- nlme(log(value) ~ calculate_ab_titers_LOG_postOnly_seasonal_MONO(time, vacc_month, log_beta0,
                                                                       logit_beta1, phi_hat, log_r),
           data = ab_titers,
           fixed = log_beta0 + logit_beta1 + phi_hat + log_r ~ 1,
           random = pdDiag(log_beta0 + log_r ~ 1),
           groups = ~subject,
           start = c(log_beta0 = log(18.0), logit_beta1 = qlogis(0.1), phi_hat = qlogis(1/12),
                     log_r = log(log(2)/365)))

plot(m6)
pairs(m6)
qqnorm(m6, abline = c(0, 1))

# rand_effects <- ranef(m6)
# rand_effects$subject <- rownames(rand_effects)
# rand_effects <- merge(rand_effects, unique(ab_titers[, c('subject', 'vacc_month')]), by = 'subject')
# plot(rand_effects)

# Optimize/fine-tune model fit ----------------------------------------------------------------------------------------

# Do we need both random effects?:
m7 <- update(m6, random = pdDiag(log_beta0 ~ 1))
m8 <- update(m6, random = pdDiag(log_r ~ 1))
anova(m6, m7)
anova(m6, m8)
rm(m7, m8)

# Get and output results ----------------------------------------------------------------------------------------------

# Now extract parameter fits and random effect sds:
results.df <- get_param_est(m6, seasonal = TRUE, mono = TRUE)
print(results.df)

# if (!dir.exists(paste0('results/PRELIM_nlme_res_', ymd, '/'))) {
#   dir.create(paste0('results/PRELIM_nlme_res_', ymd, '/'))
# }
# 
# write.csv(results.df, file = paste0('results/PRELIM_nlme_res_', ymd, '/res_n', n_participants, '_t',
#                                     sample_interval, '_MONO.csv'),
#           row.names = FALSE)

# Potentially useful references ---------------------------------------------------------------------------------------

# https://stackoverflow.com/questions/26449969/backward-selection-in-lme-singularity-in-backsolve-occured
# https://stackoverflow.com/questions/50505290/singularity-in-backsolve-at-level-0-block-1-in-lme-model

# Clean up ------------------------------------------------------------------------------------------------------------

rm(list = ls())
dev.off()
