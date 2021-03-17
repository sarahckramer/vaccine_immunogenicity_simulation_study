# ---------------------------------------------------------------------------------------------------------------------
# Compare model estimates to true parameter values across multiple sample sizes (# of participants included)
# (For now, this plots results for the biexponential model only)
# ---------------------------------------------------------------------------------------------------------------------

# Setup ---------------------------------------------------------------------------------------------------------------

dev.off()

# Load libraries:
library(ggplot2)

# Read in/format results ----------------------------------------------------------------------------------------------

# Read in:
n_participants <- c(25, 50, 100, 250, 500, 1000, 2000, 5000)

param_est <- NULL
for (i in n_participants) {
  
  filename <- paste0('results/ab_kinetics_res/res_n', i, '.csv')
  if (file.exists(filename)) {
    dat <- read.csv(filename)
    names(dat)[1] <- 'param'
    dat$subjects <- i
    param_est <- rbind(param_est, dat)
  } else {
    print(paste0('res_n', i, ' is missing.'))
  }
  
}
rm(i, dat, filename)

# Format:
param_est$param <- factor(param_est$param)

# Compare to truth ----------------------------------------------------------------------------------------------------

# Specify true values:
sd_truth <- 0.20
param_est$truth <- NA

param_est$truth[param_est$param == 'beta0'] <- 18.0
param_est$truth[param_est$param == 'beta1'] <- 0.2
param_est$truth[param_est$param == 'phi'] <- 1.0
param_est$truth[param_est$param == 'r_long'] <- log(2)/3650.0
param_est$truth[param_est$param == 'r_short'] <- log(2)/30.0 - param_est$truth[param_est$param == 'r_long']
param_est$truth[param_est$param == 'rho'] <- 0.70

# Calculate absolute and relative errors (parameter values):
param_est$abs_err_val <- param_est$Value - param_est$truth
param_est$rel_err_val <- param_est$abs_err_val / param_est$truth

# Calculate absolute and relative errors (random effect sds):
param_est$abs_err_val_sd <- param_est$sd_random - sd_truth
param_est$rel_err_val_sd <- param_est$abs_err_val_sd / sd_truth

# Plot results --------------------------------------------------------------------------------------------------------

# Plot estimates (w/ ranges) vs. truth for each subject/interval combination:
p1 <- ggplot(data = param_est) +
  geom_pointrange(aes(x = factor(subjects), y = Value, ymin = lower, ymax = upper)) +
  geom_hline(aes(yintercept = truth)) +
  facet_wrap(~ param, scales = 'free_y') +
  theme_classic() +
  labs(x = '# of Subjects', y = 'Parameter Value')
print(p1)
# Estimates become more accurate/have smaller confidence intervals as # of participants increases

# Plot estimated random effects sds for each sample size and parameter:
p2 <- ggplot(data = param_est[!(param_est$param %in% c('beta1', 'phi', 'rho')), ]) +
  geom_point(aes(x = factor(subjects), y = sd_random)) +
  geom_hline(yintercept = sd_truth) +
  facet_wrap(~ param) +
  theme_bw() +
  labs(x = '# of Subjects', y = 'SD of Random Effects')
print(p2)
# Fits well for beta0; consistently greatly underestimated for r_short and r_long

# Look at distribution of error for each parameter (value and sd) over all combinations:
p3 <- ggplot(data = param_est) +
  geom_histogram(aes(x = rel_err_val, y = 0.01 * ..density..),
                 binwidth = 0.01, col = 'white', fill = 'steelblue2') +
  geom_vline(xintercept = 0, lty = 2) +
  facet_wrap(~ param, scales = 'free_x') +
  theme_classic() +
  labs(x = 'Relative Error', y = '')
p4 <- ggplot(data = param_est[!(param_est$param %in% c('beta1', 'phi', 'rho')), ]) +
  geom_histogram(aes(x = rel_err_val_sd, y = 0.01 * ..density..),
                 binwidth = 0.01, col = 'white', fill = 'coral') +
  geom_vline(xintercept = 0, lty = 2) +
  facet_wrap(~ param, scales = 'free_x') +
  theme_classic() +
  labs(x = 'Relative Error', y = '')
print(p3)
print(p4)
# For fixed effects: r1 and r2 have most error, and beta0 the least; rho consistently slightly overestimated
# For random effects: most accurate for beta0, greatly underestimated for r_short and r_long

# Clean up:
rm(list = ls())
