### Compare model estimates to true parameter values ###
# For now, this is by number of subjects/timepoints included

# Specify date:
ymd <- '20210203'

# Read in results:
n_participants <- c(50, 100, 250, 500, 1000)
timepoint_intervals <- c(1, 2, 3, 4)

param_est <- NULL
for (i in n_participants) {
  for (j in timepoint_intervals) {
    
    filename <- paste0('results/PRELIM_nlme_res_', ymd, '/res_n', i,'_t', j, '.csv')
    if (file.exists(filename)) {
      dat <- read.csv(filename)
      dat$param <- c('beta', 'rho', 'r1', 'r2')
      dat$subjects <- i
      dat$interval <- j
      param_est <- rbind(param_est, dat)
    }
    
  }
}
rm(i, j, dat, filename)

# Format results:
param_est$param <- factor(param_est$param)
param_est$interval <- factor(param_est$interval)
levels(param_est$interval) <- c('30 Days', '60 Days', '120 Days', 'Realistic')

# Specify true values:
sd_truth <- 0.10
param_est$truth <- NA

# param_est$truth[param_est$param == 'alpha'] <- 8.0
# param_est$truth[param_est$param == 'm'] <- log(2)/42.0
param_est$truth[param_est$param == 'beta'] <- 18.0
param_est$truth[param_est$param == 'rho'] <- 0.70
param_est$truth[param_est$param == 'r1'] <- log(2)/30.0 - log(2)/3650.0
param_est$truth[param_est$param == 'r2'] <- log(2)/3650.0

# Calculate absolute and relative errors (parameter values):
param_est$abs_err_val <- param_est$Value - param_est$truth
param_est$rel_err_val <- param_est$abs_err_val / param_est$truth

# Calculate absolute and relative errors (random effect sds):
param_est$abs_err_val_sd <- param_est$sd_random - sd_truth
param_est$rel_err_val_sd <- param_est$abs_err_val_sd / sd_truth

# Plot estimates (w/ ranges) vs. truth for each subject/interval combination:
p1 <- ggplot(data = param_est) +
  geom_pointrange(aes(x = factor(subjects), y = Value, ymin = lower, ymax = upper)) +
  geom_hline(aes(yintercept = truth)) +
  facet_grid(param ~ interval, scales = 'free_y') +
  theme_classic() +
  labs(x = '# of Subjects', y = 'Parameter Value')
print(p1)
# Estimates become more accurate with more participants, but not necessarily with more timepoints
    # Exception is r1, which is more accurately and more exactly estimated with more data points
# Estimates also have smaller confidence intervals as # of participants increases
# Estimates for most parameters appear quite accurate (rho is consistently underestimated, but not by much)

# Plot estimated sds for each combination and parameter:
p2 <- ggplot(data = param_est[param_est$param != 'rho', ]) +
  geom_point(aes(x = factor(subjects), y = sd_random)) +
  geom_hline(yintercept = sd_truth) +
  facet_grid(param ~ interval) +
  theme_bw() +
  labs(x = '# of Subjects', y = 'SD of Random Effects')
print(p2)
# Seems to have most trouble with r2, but also some with r1, esp. when fewer timepoints
# Tends to underestimate variability in r1; can go over or under with r2

# Look at distribution of error for each parameter (value and sd) over all combinations:
p3 <- ggplot(data = param_est) +
  geom_histogram(aes(x = rel_err_val, y = 0.01 * ..density..),
                 binwidth = 0.01, col = 'white', fill = 'steelblue2') +
  geom_vline(xintercept = 0, lty = 2) +
  facet_wrap(~ param, scales = 'free_x') +
  theme_classic() +
  labs(x = 'Relative Error', y = '')
p4 <- ggplot(data = param_est[param_est$param != 'rho', ]) +
  geom_histogram(aes(x = rel_err_val_sd, y = 0.01 * ..density..),
                 binwidth = 0.01, col = 'white', fill = 'coral') +
  geom_vline(xintercept = 0, lty = 2) +
  facet_wrap(~ param, scales = 'free_x') +
  theme_classic() +
  labs(x = 'Relative Error', y = '')
print(p3)
print(p4)
# For fixed effects: r1 and r2 have most error; rho consistently underestimated
# For random effects: sd most accurate for beta/m, least for r1/r2

# Look at error in each parameter (value and sd) by # of subjects, and by interval, separately:
p5 <- ggplot(data = param_est) +
  geom_histogram(aes(x = rel_err_val, y = 0.01 * ..density..),
                 binwidth = 0.01, col = 'white', fill = 'steelblue2') +
  geom_vline(xintercept = 0, lty = 2) +
  facet_wrap(~ subjects, scales = 'free_x', nrow = 1) +
  theme_classic() +
  labs(x = 'Relative Error', y = '')
p6 <- ggplot(data = param_est) +
  geom_histogram(aes(x = rel_err_val, y = 0.01 * ..density..),
                 binwidth = 0.01, col = 'white', fill = 'steelblue2') +
  geom_vline(xintercept = 0, lty = 2) +
  facet_wrap(~ interval, scales = 'free_x', nrow = 1) +
  theme_classic() +
  labs(x = 'Relative Error', y = '')
grid.arrange(p5, p6)
# Difficult to tell w/o looking at parameters separately, but seems interval has more influence than subjects
# This is not consistent with p1; but here we are looking at all parameters together

# p7 <- ggplot(data = param_est) +
#   geom_histogram(aes(x = rel_err_val_sd, y = 0.01 * ..density..),
#                  binwidth = 0.01, col = 'white', fill = 'coral') +
#   geom_vline(xintercept = 0, lty = 2) +
#   facet_wrap(~ subjects, scales = 'free_x', nrow = 1) +
#   theme_classic() +
#   labs(x = 'Relative Error', y = '')
# p8 <- ggplot(data = param_est) +
#   geom_histogram(aes(x = rel_err_val_sd, y = 0.01 * ..density..),
#                  binwidth = 0.01, col = 'white', fill = 'coral') +
#   geom_vline(xintercept = 0, lty = 2) +
#   facet_wrap(~ interval, scales = 'free_x', nrow = 1) +
#   theme_classic() +
#   labs(x = 'Relative Error', y = '')
# grid.arrange(p7, p8)
# # Difficult to tell much from this
