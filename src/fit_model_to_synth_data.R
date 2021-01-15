
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

# Function to visualize model fit:
viz_model_fit <- function(m, dat, logscale = F) {
  print(summary(m))
  
  if (logscale == T) {
    dat$fitted <- exp(fitted(m))
  } else {
    dat$fitted <- fitted(m)
  }
  
  dat$residuals_fixed <- m$residuals[, 1]
  dat$residuals_subject <- m$residuals[, 2]
  
  p1 <- ggplot(data = dat) + geom_line(aes(x = time, y = value, color = subject)) + geom_point(aes(x = time, y = fitted, color = subject)) +
    theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Ab Titers') + 
    scale_x_continuous(breaks = seq(0, 120, by = 15)) + scale_y_continuous(breaks = seq(0, 0.4, by = 0.05))#, trans = 'log')
  
  p21 <- ggplot(data = dat) + geom_point(aes(x = time, y = residuals_fixed, color = subject)) +
    theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Residuals (Fixed)') +
    scale_x_continuous(breaks = seq(0, 120, by = 15))
  p22 <- ggplot(data = dat) + geom_point(aes(x = time, y = residuals_subject, color = subject)) +
    theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Residuals (Subject)') +
    scale_x_continuous(breaks = seq(0, 120, by = 15))
  
  grid.arrange(p1, p21, p22, layout_matrix = rbind(c(1, 1), c(1, 1), c(2, 3)))
}

# # Set seed so that same 100 subjects chosen each time:
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
# ab_titers <- ab_titers_ORIG[ab_titers_ORIG$time %in% c(0, 15, 30, 45, 60, 75, 90, 105, 120), ]
ab_titers <- ab_titers_ORIG[ab_titers_ORIG$time %in% c(seq(0, 360, by = 30), 379, seq(390, 720, by = 30)), ]
# ab_titers <- ab_titers_ORIG[ab_titers_ORIG$time %in% c(seq(0, 360, by = 30), 379, seq(390, 3650, by = 30)), ] # allows nlsList to converge

# # And specify step function where vaccine occurs:
# ab_titers$h_1 <- ifelse(ab_titers$time >= (365 + 14), 1, 0)

# Plot data:
p1 <- ggplot(data = ab_titers, aes(x = time, y = value, color = subject)) + geom_line() +
  geom_vline(xintercept = 365, lty = 2) +
  theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Ab Titers') +
  scale_x_continuous(breaks = seq(0, max(ab_titers$time), by = 30)) +
  scale_y_continuous(breaks = seq(0, max(ab_titers$value), by = 1.0))#, trans = 'log')
print(p1)

# Now begin trying to fit:
calculate_ab_titers <- function(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2) {
  # if (log_r_2 < -10) {
  #   log_r_2 <- -10
  # } # tried to restrict r_2 to prevent issue with nlsList, but no luck
  
  alpha = exp(log_alpha)
  m = exp(log_m)
  beta = exp(log_beta)
  rho = exp(logit_rho) / (exp(logit_rho) + 1)
  r_1 = exp(log_r_1)
  r_2 = exp(log_r_2)
  
  v1 <- alpha * exp(-m * time)
  v2 <- alpha * exp(-m * time) + beta * (rho * exp(-r_1 * (time - (365 + 14))) + (1 - rho) * exp(-r_2 * (time - (365 + 14))))
  value <- v1 * (time < (365 + 14)) + v2 * (time >= (365 + 14))
  return(value)
}

# # alternative way to write same function:
# calculate_ab_titers <- function(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2) {
#   v1 <- exp(log_alpha) * exp(-exp(log_m) * time)
#   v2 <- exp(log_alpha) * exp(-exp(log_m) * time) +
#     exp(log_beta) * ((exp(logit_rho) / (exp(logit_rho) + 1)) * exp(-exp(log_r_1) * (time - (365 + 14))) +
#                        (1 - (exp(logit_rho) / (exp(logit_rho) + 1))) * exp(-exp(log_r_2) * (time - (365 + 14))))
#   value <- v1 * (time < (365 + 14)) + v2 * (time >= (365 + 14))
#   return(value)
# }

logit <- function(p) {
  return(log(p / (1 - p)))
}

m1 <- nls(value ~ calculate_ab_titers(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2),
          data = ab_titers,
          start = c(log_alpha = log(5.0), log_m = log(log(2)/30), log_beta = log(18.0), logit_rho = logit(0.75),
                    log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650)))
exp(summary(m1)$coefficients['log_alpha', ])
exp(summary(m1)$coefficients['log_beta', ])
log(2)/exp(summary(m1)$coefficients['log_m', ])
log(2)/exp(summary(m1)$coefficients['log_r_1', ])
log(2)/exp(summary(m1)$coefficients['log_r_2', ])
exp(summary(m1)$coefficients['logit_rho', ]) / (exp(summary(m1)$coefficients['logit_rho', ]) + 1)
# good! seems to converge/fit well for now

# m2 <- nls(value ~ calculate_ab_titers(time, log(alpha), log(m), log(beta), logit(rho), log(r_1), log(r_2)),
#           data = ab_titers,
#           start = c(alpha = 5.0, m = log(2)/30, beta = 18.0, rho = 0.75, r_1 = log(2)/30, r_2 = log(2)/3650))
# summary(m2)$coefficients
# same estimates - use m1 and estimate the logs/logits of the values directly

plot(m1)
plot(m1, subject ~ resid(.), abline = 0)

# Here a selfStart function is possible, if desired

# Move on to fit each "group" (individual) separately:
m2 <- nlsList(value ~ calculate_ab_titers(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2) | subject,
              data = ab_titers,
              start = c(log_alpha = log(5.0), log_m = log(log(2)/30), log_beta = log(18.0), logit_rho = logit(0.75),
                        log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650)))
plot(intervals(m2))
plot(m2, subject ~ resid(.), abline = 0 )
# warnings for almost half of them!
# Writing the model in a different way does not help, nor does setting more realistic start values
# So maybe we need some constraint saying that r2 << r1? Also possible that a selfStart function would help

# Holding r_2 constant fixes the problem in most or even all cases; holding rho doesn't, even though no random effect there:
calculate_ab_titers_log_params_r2const <- function(time, log_alpha, log_m, log_beta, logit_rho, log_r_1) {
  v1 <- exp(log_alpha) * exp(-exp(log_m) * time)
  v2 <- exp(log_alpha) * exp(-exp(log_m) * time) +
    exp(log_beta) * ((exp(logit_rho) / (exp(logit_rho) + 1)) * exp(-exp(log_r_1) * (time - (365 + 14))) +
                       (1 - (exp(logit_rho) / (exp(logit_rho) + 1))) * exp(-exp(log(2)/3650) * (time - (365 + 14))))
  value <- v1 * (time < (365 + 14)) + v2 * (time >= (365 + 14))
  return(value)
}
m3 <- nlsList(value ~ calculate_ab_titers_log_params_r2const(time, log_alpha, log_m, log_beta, logit_rho, log_r_1) | subject,
              data = ab_titers,
              start = c(log_alpha = log(8.0), log_m = log(log(2)/40), log_beta = log(18.0), logit_rho = logit(0.75),
                        log_r_1 = log(log(2)/20)))
plot(intervals(m3))
# This makes me wonder whether a longer time period of time would work better, or if less noisy data would fit better; or maybe r_2 just too low, or not variable enough?
# Using "truth" instead of "obs" (so, data with no error added) actually has it fail for all 100 participants
# Same issue occurs if we include all "participants" instead of just 100
# Including all 10 years of "data" allows m2 to converge! So just not enough info about r_2 to fit?










# m2 <- nlme(value ~ alpha * exp(-m * time) +
#              h_1 * beta * (rho * exp(-r_1 * (time - 379)) +
#                              (1 - rho) * exp(-r_2 * (time - 379))),
#            data = ab_titers,
#            fixed = alpha + m + beta + rho + r_1 + r_2 ~ 1,
#            random = alpha + m + beta + r_1 + r_2 ~ 1,
#            groups = ~subject,
#            start = c(5.0, log(2)/30, 18.0, 0.75, log(2)/30, log(2)/3650))

# for log scale: can use glmmPQL in MASS package instead (https://www.juliapilowsky.com/2018/10/19/a-practical-guide-to-mixed-models-in-r/)

# FIND A WAY TO CHECK AGAINST OBSERVED SDs, TOO
# TRY GEE; LME4; BRMS; RSTAN; ETC.
# ALTERNATIVE MODELS?
# REORGANIZE ROW ORDER
# IS IT IMPORTANT THAT TIME BE NUMERIC RATHER THAN INTEGER?
# consider weights/correlation


