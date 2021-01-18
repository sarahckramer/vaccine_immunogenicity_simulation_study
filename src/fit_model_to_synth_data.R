
# # Install needed packages and update renv lock file:
# install.packages('reshape2')
# install.packages('nlme')
# install.packages('lme4')
# install.packages('ggplot2')
# renv::snapshot()

# Load libraries:
library(reshape2)
library(nlme)
library(ggplot2)
library(gridExtra)

# Set seed so that same 100 subjects chosen each time:
set.seed(3970395)

# Read in noise-laden "data":
ab_titers <- read.csv('data/prelim_check_20210115/obs_data.csv')
# ab_titers <- read.csv('data/prelim_check_20210115/truth.csv')

# Reformat data frame:
ab_titers$time <- 0:(dim(ab_titers)[1] - 1)
ab_titers$time <- as.numeric(as.character(ab_titers$time))
ab_titers <- melt(ab_titers, id.vars = 'time')
names(ab_titers)[2] <- 'subject'

# Select subset of subjects?:
ab_titers <- ab_titers[ab_titers$subject %in% levels(ab_titers$subject)[sample(1:1000, 100)], ]
ab_titers$subject <- factor(ab_titers$subject)

# Select key timepoints:
ab_titers_ORIG <- ab_titers # store "full" data
ab_titers <- ab_titers_ORIG[ab_titers_ORIG$time %in% c(seq(0, 360, by = 30), 379, seq(390, 720, by = 30)), ]
# ab_titers <- ab_titers_ORIG[ab_titers_ORIG$time %in% c(seq(0, 360, by = 30), 379, seq(390, 3650, by = 30)), ] # allows nlsList to converge

# Plot data:
p1 <- ggplot(data = ab_titers, aes(x = time, y = value, color = subject)) + geom_line() + geom_point() +
  geom_vline(xintercept = 365, lty = 2) +
  theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Ab Titers') +
  scale_x_continuous(breaks = seq(0, max(ab_titers$time), by = 30)) +
  scale_y_continuous(breaks = seq(0, max(ab_titers$value), by = 1.0))#, trans = 'log')
print(p1)

# Read in functions:
source('functions_ab_titers_over_time.R')
source('functions_assess_output.R')

# First fit without random effects:
m1 <- nls(log(value) ~ calculate_ab_titers_LOG(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2),
          data = ab_titers,
          start = c(log_alpha = log(5.0), log_m = log(log(2)/30), log_beta = log(18.0), logit_rho = logit(0.75),
                    log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650)))
plot(m1)

# Next try nlsList (separate fits for each person):
m2 <- nlsList(log(value) ~ calculate_ab_titers_LOG(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2) | subject,
              data = ab_titers,
              start = c(log_alpha = log(5.0), log_m = log(log(2)/30), log_beta = log(18.0), logit_rho = logit(0.75),
                        log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650)))
plot(intervals(m2)) # almost half not fit, but m seems most variable based on this plot; also alpha and beta
# plot(m2, subject ~ resid(.), abline = 0 )
pairs(m2)
# warnings for almost half of them! this seems mostly fixed by setting r2 constant, or else using a longer timeframe for data

# Fit nlme:
# m3 <- nlme(log(value) ~ calculate_ab_titers_LOG(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2),
#            data = ab_titers,
#            fixed = log_alpha + log_m + log_beta + logit_rho + log_r_1 + log_r_2 ~ 1,
#            random = log_alpha + log_m + log_beta + log_r_1 + log_r_2 ~ 1,
#            groups = ~subject,
#            start = c(log_alpha = log(5.0), log_m = log(log(2)/30), log_beta = log(18.0), logit_rho = logit(0.75),
#                      log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650)))
# # might be able to get estimates if allowed to run for several hours, but iterations seem not to converge quite often anyway

m4 <- nlme(log(value) ~ calculate_ab_titers_LOG(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2),
           data = ab_titers,
           fixed = log_alpha + log_m + log_beta + logit_rho + log_r_1 + log_r_2 ~ 1,
           random = pdDiag(log_alpha + log_m + log_beta + log_r_1 + log_r_2 ~ 1),
           groups = ~subject,
           start = c(log_alpha = log(5.0), log_m = log(log(2)/30), log_beta = log(18.0), logit_rho = logit(0.75),
                     log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650)))
# however, making small change to how random effects are specified allows fitting!
# won't converge anymore if every 60 days instead of every 30 days; have to remove r1 and r2, or m and r2 (keeping m seems more important; r_1 and beta seem correlated)
plot(m4)
pairs(m4) # no obvious strong correlations
qqnorm(m4, abline = c(0, 1))
qqnorm(m4, ~ ranef(.))
viz_model_fit(m4, ab_titers, logscale = TRUE)

# Try to account for any correlation structure?:
plot(ACF(m4, maxLag = 10), alpha = 0.05) # seems to be some autocorrelation?
m5 <- update(m4, correlation = corAR1(abs(ACF(m4)[2, 2])))
anova(m4, m5) # doesn't significantly improve fit
# viz_model_fit(m5, ab_titers, logscale = TRUE)
rm(m5)

# Now extract parameter fits and random effect sds:
results.df <- get_param_est(m4)
print(results.df)





