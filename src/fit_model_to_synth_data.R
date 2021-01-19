
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

# Set number of data points/subjects:
n_participants <- 50 # 50, 100, 250, 500, 1000
sample_interval <- 1 # 1-4; selecting from list below
select_timepoints <- list(c(seq(0, 360, by = 30), 379, seq(379+21, 1090, by = 30)),
                          c(seq(0, 360, by = 60), 379, seq(379+51, 1090, by = 60)),
                          c(0, 180, 360, 379, 550, 730, 1095),
                          c(0, 360, 379, 730, 1095, 2190))

# Set seed so that same subjects chosen each time:
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
ab_titers <- ab_titers[ab_titers$subject %in% levels(ab_titers$subject)[sample(1:1000, n_participants)], ]
ab_titers$subject <- factor(ab_titers$subject)

# Select key timepoints:
ab_titers <- ab_titers[ab_titers$time %in% select_timepoints[[sample_interval]], ]
# need at least 2 points pre-vaccine, at least 4 post?

# Plot data:
p1 <- ggplot(data = ab_titers, aes(x = time, y = value, color = subject)) + geom_line() + geom_point() +
  geom_vline(xintercept = 365, lty = 2) +
  theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Ab Titers') +
  scale_x_continuous(breaks = seq(0, max(ab_titers$time), by = 30)) +
  scale_y_continuous(breaks = seq(0, max(ab_titers$value), by = 1.0))#, trans = 'log')
print(p1)

# Read in functions:
source('src/functions_ab_titers_over_time.R')
source('src/functions_assess_output.R')

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

# Fit nlme:
m3 <- nlme(log(value) ~ calculate_ab_titers_LOG(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2),
           data = ab_titers,
           fixed = log_alpha + log_m + log_beta + logit_rho + log_r_1 + log_r_2 ~ 1,
           random = pdDiag(log_alpha + log_m + log_beta + log_r_1 + log_r_2 ~ 1),
           groups = ~subject,
           start = c(log_alpha = log(5.0), log_m = log(log(2)/30), log_beta = log(18.0), logit_rho = logit(0.75),
                     log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650)))
plot(m3) # evidence of some increase in variance with values higher than about 1.5
plot(m3, resid(.) ~ log(value), abline = 0)
plot(m3, resid(.) ~ time, abline = 0)
pairs(m3) # r_1 and r_2 might be correlated (moreso when more data points) - explore whether both random effects are needed
qqnorm(m3, abline = c(0, 1)) # looks better when more data points are included
qqnorm(m3, ~ ranef(.))
# viz_model_fit(m3, ab_titers, logscale = TRUE)

# Do we need random effect of r2?
m4 <- update(m3, random = pdDiag(log_alpha + log_m + log_beta + log_r_1 ~ 1))
anova(m3, m4)
rm(m4)

# Try to account for any correlation structure?:
plot(ACF(m3, maxLag = 10), alpha = 0.05) # seems to be some autocorrelation?
m4 <- update(m3, correlation = corAR1(abs(ACF(m3)[2, 2])))
anova(m3, m4) # doesn't significantly improve fit (close for 250/2)
# viz_model_fit(m4, ab_titers, logscale = TRUE)
rm(m4)

# What about variance functions?:
m4 <- update(m3, weights = varPower(form = ~log(value)))
anova(m3, m4)
plot(m4, resid(.) ~ log(value), abline = 0)
# seems to improve for some (more likely to make a difference when more timepoints/participants), but don't see much visual evidence

# Now extract parameter fits and random effect sds:
results.df <- get_param_est(m3)
print(results.df)

if (!dir.exists('results/PRELIM_nlme_res_20210118/')) {
  dir.create('results/PRELIM_nlme_res_20210118/')
}

write.csv(results.df, file = paste0('results/PRELIM_nlme_res_20210118/res_n', n_participants, '_t', sample_interval, '.csv'), row.names = FALSE)

# Clean up:
rm(list = ls())
dev.off()
