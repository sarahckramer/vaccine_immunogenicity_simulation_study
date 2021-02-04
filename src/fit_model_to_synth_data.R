
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

# Specify date (for results outputs):
ymd <- '20210203'

# Set number of data points/subjects:
n_participants <- 50 # 50, 100, 250, 500, 1000
sample_interval <- 1 # 1-4; selecting from list below
select_timepoints <- list(c(379, seq(379+21, 1090, by = 30)),
                          c(seq(379+51, 1090, by = 60)),
                          c(379, 550, 730, 1095),
                          c(379, 730, 1095, 2190))

# select_timepoints <- list(c(365, 379, seq(379+21, 1090, by = 30)),
#                           c(365, seq(379+51, 1090, by = 60)),
#                           c(365, 379, 550, 730, 1095),
#                           c(365, 379, 730, 1095, 2190))
# # alternative to also include timepoint right before vacc; but here we're assuming negligible titers at this point, so might not need it

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

# Set time of vaccine = time 0:
ab_titers$time <- ab_titers$time - 365

# Plot data:
p1 <- ggplot(data = ab_titers, aes(x = time, y = value, color = subject)) + geom_line() + geom_point() +
  # geom_vline(xintercept = 0, lty = 2) +
  theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Ab Titers') +
  scale_x_continuous(breaks = seq(0, max(ab_titers$time), by = 30)) +
  scale_y_continuous(breaks = seq(0, max(ab_titers$value), by = 1.0))#, trans = 'log')
print(p1)

# Read in functions:
source('src/functions_ab_titers_over_time.R')
source('src/functions_assess_output.R')

# First fit without random effects:
m1 <- nls(log(value) ~ calculate_ab_titers_LOG_postOnly(time, log_beta, logit_rho, log_r_1, log_r_2),
          data = ab_titers,
          start = c(log_beta = log(18.0), logit_rho = qlogis(0.75), log_r_1 = log(log(2)/30),
                    log_r_2 = log(log(2)/3650)))
summary(m1)

# Next try nlsList (separate fits for each person):
m2 <- nlsList(log(value) ~ calculate_ab_titers_LOG_postOnly(time, log_beta, logit_rho, log_r_1, log_r_2) | subject,
              data = ab_titers,
              start = c(log_beta = log(18.0), logit_rho = qlogis(0.75), log_r_1 = log(log(2)/30),
                        log_r_2 = log(log(2)/3650)))
plot(intervals(m2)) # almost half not fit, but m seems most variable based on this plot; also alpha and beta
# plot(m2, subject ~ resid(.), abline = 0 )
pairs(m2)

# Fit nlme:
m3 <- nlme(log(value) ~ calculate_ab_titers_LOG_postOnly(time, log_beta, logit_rho, log_r_1, log_r_2),
           data = ab_titers,
           fixed = log_beta + logit_rho + log_r_1 + log_r_2 ~ 1,
           random = pdDiag(log_beta + log_r_1 + log_r_2 ~ 1),
           groups = ~subject,
           start = c(log_beta = log(18.0), logit_rho = qlogis(0.75), log_r_1 = log(log(2)/30),
                     log_r_2 = log(log(2)/3650)))
m3.alt <- nlme(m2, random = pdDiag(log_beta + log_r_1 + log_r_2 ~ 1)) # fit is very similar
plot(m3.alt) # evidence of some increase in variance with values higher than about 1.5
# plot(m3, resid(.) ~ log(value), abline = 0)
# plot(m3, resid(.) ~ time, abline = 0)
pairs(m3.alt) # r_1 and r_2 might be correlated (moreso when more data points) - explore whether both random effects are needed
qqnorm(m3.alt, abline = c(0, 1)) # looks better when more data points are included
# qqnorm(m3, ~ ranef(.))
# viz_model_fit(m3, ab_titers, logscale = TRUE)

# Do we need random effect of r2?
m4 <- update(m3.alt, random = pdDiag(log_beta + log_r_1 ~ 1))
anova(m3.alt, m4)
rm(m4)

# # Try to account for any correlation structure?:
# plot(ACF(m3, maxLag = 10), alpha = 0.05) # seems to be some autocorrelation?
# m4 <- update(m3, correlation = corAR1(abs(ACF(m3)[2, 2])))
# anova(m3, m4) # doesn't significantly improve fit (close for 250/2)
# # viz_model_fit(m4, ab_titers, logscale = TRUE)
# rm(m4)
# 
# # What about variance functions?:
# m4 <- update(m3, weights = varPower(form = ~log(value)))
# anova(m3, m4)
# plot(m4, resid(.) ~ log(value), abline = 0)
# # seems to improve for some (more likely to make a difference when more timepoints/participants), but don't see much visual evidence

# Now extract parameter fits and random effect sds:
results.df <- get_param_est(m3.alt)
print(results.df)

# if (!dir.exists(paste0('results/PRELIM_nlme_res_', ymd, '/'))) {
#   dir.create(paste0('results/PRELIM_nlme_res_', ymd, '/'))
# }
# 
# write.csv(results.df, file = paste0('results/PRELIM_nlme_res_', ymd, '/res_n', n_participants, '_t', sample_interval, '.csv'), row.names = FALSE)

# # saemix?:
# ab_titers$log_value <- log(ab_titers$value)
# library(saemix)
# ab_mix <- saemixData(name.data = ab_titers, name.group = 'subject', name.predictors = 'time', name.response = 'log_value', units = list(x = 'time', y = 'log_value'))
# ab_mod <- function(psi, id, xidep) {
#   time <- xidep[, 1]
# 
#   log_alpha <- psi[id, 1]
#   log_m <- psi[id, 2]
#   log_beta <- psi[id, 3]
#   logit_rho <- qlogis(psi[id, 4])
#   log_r_1 <- psi[id, 5]
#   log_r_2 <- psi[id, 6]
# 
#   log_value <- calculate_ab_titers_LOG(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2)
#   return(log_value)
# }
# saemix_model <- saemixModel(model = ab_mod,
#                             psi0 = matrix(c(log(5.0), log(log(2)/30), log(18.0), 0.75, log(log(2)/30), log(log(2)/3650)),
#                                           ncol = 6,
#                                           byrow = TRUE,
#                                           dimnames = list(NULL, c('log_alpha', 'log_m', 'log_beta', 'logit_rho', 'log_r_1', 'log_r_2'))),
#                             transform.par = c(0, 0, 0, 3, 0, 0),
#                             fixed.estim = c(1, 1, 1, 1, 1, 1),
#                             covariance.model = matrix(c(1, 1, 1, 0, 1, 1,
#                                                         1, 1, 1, 0, 1, 1,
#                                                         1, 1, 1, 0, 1, 1,
#                                                         0, 0, 0, 0, 0, 0,
#                                                         1, 1, 1, 0, 1, 1,
#                                                         1, 1, 1, 0, 1, 1),
#                                                       ncol = 6, byrow = TRUE))
# opt <- list(seed = 94352514, save = FALSE, save.graphs = FALSE)
# m4 <- saemix(saemix_model, ab_mix, opt)
# summary(m4)
# # similar estimates of parameter values and sds

# # brms?:
# ab_titers$logvalue <- log(ab_titers$value)
# library(brms)
# f1 <- logvalue ~ log((exp(logalpha) * exp(-exp(logm) * time)) * (time < (365 + 14)) +
#                        ((exp(logalpha) * exp(-exp(logm) * time)) + exp(logbeta) *
#                           ((exp(logitrho) / (exp(logitrho) + 1)) * (exp(-exp(logr1) * (time - (365 + 14)))) +
#                              (1 - exp(logitrho) / (exp(logitrho) + 1)) * (exp(-exp(logr2) * (time - (365 + 14)))))) * (time >= (365 + 14)))
# prior_1 <- c(set_prior("normal(log(8.0), 0.1)", nlpar = 'logalpha'),
#              set_prior("normal(log(log(2)/42), 0.1)", nlpar = 'logm'),
#              set_prior("normal(log(18.0), 0.1)", nlpar = 'logbeta'),
#              set_prior("normal(qlogis(0.75), 0.05)", nlpar = 'logitrho'),
#              set_prior("normal(log(log(2)/30), 0.1)", nlpar = 'logr1'),
#              set_prior("normal(log(log(2)/3650), 0.1)", nlpar = 'logr2'))
# form = bf(f1, nl = TRUE) + list(logitrho~1, logalpha~(1|2|subject), logm~(1|2|subject),
#                                 logbeta~(1|2|subject), logr1~(1|2|subject), logr2~(1|2|subject))
# m4 <- brm(form, data = ab_titers, prior = prior_1)
# summary(m4)
# # right now at least, really slow; seems to maybe have same issues as saemix with getting NAs from function?

# Clean up:
rm(list = ls())
dev.off()
