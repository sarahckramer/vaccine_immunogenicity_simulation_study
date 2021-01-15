
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
ab_titers <- ab_titers[ab_titers$time %in% c(0, 30, 60, 90, 120, 180, 240, 300, 365), ]

# Plot "data":
p1 <- ggplot(data = ab_titers, aes(x = time, y = value, color = subject)) + geom_line() +
  theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Ab Titers') + 
  scale_x_continuous(breaks = seq(0, 365, by = 15)) + scale_y_continuous(breaks = seq(0, max(ab_titers$value), by = 0.5))#, trans = 'log')
print(p1)

### Explore possibility of fitting ONLY before vaccination ###
m1 <- nlme(log(value) ~ log(A_m) - m * time,
           data = ab_titers,
           fixed = A_m + m ~ 1,
           random = A_m + m ~ 1,
           groups = ~subject,
           method = 'ML',
           start = c(5.0, log(2)/30))
viz_model_fit(m1, ab_titers, logscale = TRUE)

# Clean up:
rm(ab_titers, p1, m1)

### Move on to trying to fit including data post-first primary dose ###

# First, get increased time range from data:
# ab_titers <- ab_titers_ORIG[ab_titers_ORIG$time %in% c(0, 15, 30, 45, 60, 75, 90, 105, 120), ]
ab_titers <- ab_titers_ORIG[ab_titers_ORIG$time %in% c(seq(0, 360, by = 30), 379, seq(390, 720, by = 30)), ]

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
calculate_ab_titers <- function(time, alpha, m, beta, rho, r_1, r_2) {
  v1 <- alpha * exp(-m * time)
  v2 <- alpha * exp(-m * time) + beta * (rho * exp(-r_1 * (time - (365 + 14))) + (1 - rho) * exp(-r_2 * (time - (365 + 14))))
  value <- v1 * (time < (365 + 14)) + v2 * (time >= (365 + 14))
  return(value)
}

m1 <- nls(value ~ calculate_ab_titers(time, alpha, m, beta, rho, r_1, r_2),
          data = ab_titers,
          start = c(alpha = 5.0, m = log(2)/30, beta = 18.0, rho = 0.75, r_1 = log(2)/30, r_2 = log(2)/3650))
summary(m1)
log(2)/summary(m1)$coefficients['m', 1]
log(2)/summary(m1)$coefficients['r_1', 1]
log(2)/summary(m1)$coefficients['r_2', 1]
# good! seems to converge for now

plot(m1)
plot(m1, subject ~ resid(.), abline = 0)

# Force parameters to remain positive:
calculate_ab_titers_log_params <- function(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2) {
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

logit <- function(p) {
  return(log(p / (1 - p)))
}

m2 <- nls(value ~ calculate_ab_titers_log_params(time, log(alpha), log(m), log(beta), logit(rho), log(r_1), log(r_2)),
          data = ab_titers,
          start = c(alpha = 5.0, m = log(2)/30, beta = 18.0, rho = 0.75, r_1 = log(2)/30, r_2 = log(2)/3650))
# same results

# # Fit on log scale instead?:
# m3 <- nls(log(value) ~ log(calculate_ab_titers_log_params(time, log(alpha), log(m), log(beta), logit(rho), log(r_1), log(r_2))),
#           data = ab_titers,
#           start = c(alpha = 5.0, m = log(2)/30, beta = 18.0, rho = 0.75, r_1 = log(2)/30, r_2 = log(2)/3650))
# # very similar results; I think this is how we'll want to fit it eventually

# Here a selfStart function is possible, if desired

# Move on to fit each "group" (individual) separately:
m3 <- nlsList(value ~ calculate_ab_titers_log_params(time, log(alpha), log(m), log(beta), logit(rho), log(r_1), log(r_2)) | subject,
              data = ab_titers,
              start = c(alpha = 5.0, m = log(2)/30, beta = 18.0, rho = 0.75, r_1 = log(2)/30, r_2 = log(2)/3650))
plot(intervals(m3))
# warning: NAs produced by log - so not forcing positive as well as it should?

m4 <- nlsList(value ~ calculate_ab_titers(time, alpha, m, beta, rho, r_1, r_2) | subject,
              data = ab_titers,
              start = c(alpha = 5.0, m = log(2)/30, beta = 18.0, rho = 0.75, r_1 = log(2)/30, r_2 = log(2)/3650))
plot(intervals(m4))
# but this does the same thing, it just allows the values to be negative, rather than returning an error...

m5 <- nlsList(value ~ calculate_ab_titers_log_params(time, alpha, m, beta, rho, r_1, r_2) | subject,
              data = ab_titers,
              start = c(alpha = log(5.0), m = log(log(2)/30), beta = log(18.0), rho = logit(0.75), r_1 = log(log(2)/30), r_2 = log(log(2)/3650)))
plot(exp(intervals(m5))) # ignore rho here
# several errors

plot(m3, subject ~ resid(.), abline = 0 )
plot(m4, subject ~ resid(.), abline = 0 )
plot(m5, subject ~ resid(.), abline = 0 )

# Would coding the model a bit differently help?:
calculate_ab_titers_log_params <- function(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2) {
  v1 <- exp(log_alpha) * exp(-exp(log_m) * time)
  v2 <- exp(log_alpha) * exp(-exp(log_m) * time) +
    exp(log_beta) * ((exp(logit_rho) / (exp(logit_rho) + 1)) * exp(-exp(log_r_1) * (time - (365 + 14))) +
                       (1 - (exp(logit_rho) / (exp(logit_rho) + 1))) * exp(-exp(log_r_2) * (time - (365 + 14))))
  value <- v1 * (time < (365 + 14)) + v2 * (time >= (365 + 14))
  return(value)
}

m3 <- nlsList(value ~ calculate_ab_titers_log_params(time, log(alpha), log(m), log(beta), logit(rho), log(r_1), log(r_2)) | subject,
              data = ab_titers,
              start = c(alpha = 5.0, m = log(2)/30, beta = 18.0, rho = 0.75, r_1 = log(2)/30, r_2 = log(2)/3650))
plot(intervals(m3))
# nope, same issue

m4 <- nlsList(value ~ exp(log_alpha) * exp(-exp(log_m) * time) +
                h_1 * (exp(log_beta) * ((exp(logit_rho) / (exp(logit_rho) + 1)) * exp(-exp(log_r_1) * (time - (365 + 14))) +
                                          (1 - (exp(logit_rho) / (exp(logit_rho) + 1))) * exp(-exp(log_r_2) * (time - (365 + 14))))) | subject,
              data = ab_titers,
              start = c(log_alpha = log(5.0), log_m = log(log(2)/30), log_beta = log(18.0), logit_rho = logit(0.75),
                        log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650)))
plot(exp(intervals(m4)))
# even including the model explicitly gives these issues

# What about just changing the start values to be better?:
m3 <- nlsList(value ~ calculate_ab_titers_log_params(time, log(alpha), log(m), log(beta), logit(rho), log(r_1), log(r_2)) | subject,
              data = ab_titers,
              start = c(alpha = 8.0, m = log(2)/40, beta = 18.0, rho = 0.75, r_1 = log(2)/60, r_2 = log(2)/3650))
plot(intervals(m3))
# doesn't necessarily improve; sometimes makes it worse!

# I think we do need to be estimating the log/logit of the values and not the values themselves, though, to keep things unconstrained:
m3 <- nlsList(value ~ calculate_ab_titers_log_params(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2) | subject,
              data = ab_titers,
              start = c(log_alpha = log(8.0), log_m = log(log(2)/40), log_beta = log(18.0), logit_rho = logit(0.75),
                        log_r_1 = log(log(2)/20), log_r_2 = log(log(2)/165)))
plot(intervals(m3))
# even if I make the starting values more realistic, I still can't get fits for everyone...

# What if I hold r_2 constant?:
calculate_ab_titers_log_params_r2const <- function(time, log_alpha, log_m, log_beta, logit_rho, log_r_1) {
  v1 <- exp(log_alpha) * exp(-exp(log_m) * time)
  v2 <- exp(log_alpha) * exp(-exp(log_m) * time) +
    exp(log_beta) * ((exp(logit_rho) / (exp(logit_rho) + 1)) * exp(-exp(log_r_1) * (time - (365 + 14))) +
                       (1 - (exp(logit_rho) / (exp(logit_rho) + 1))) * exp(-exp(log(2)/3650) * (time - (365 + 14))))
  value <- v1 * (time < (365 + 14)) + v2 * (time >= (365 + 14))
  return(value)
}
m4 <- nlsList(value ~ calculate_ab_titers_log_params_r2const(time, log_alpha, log_m, log_beta, logit_rho, log_r_1) | subject,
              data = ab_titers,
              start = c(log_alpha = log(8.0), log_m = log(log(2)/40), log_beta = log(18.0), logit_rho = logit(0.75),
                        log_r_1 = log(log(2)/20)))
plot(intervals(m4))
# now it succeeds for all but 3?

# At least look at the parameter values for individuals that can vs. cannot be fit:
na3 <- names(which(is.na(summary(m3)$coefficients[, 1, 'log_r_1'])))
na4 <- names(which(is.na(summary(m4)$coefficients[, 1, 'log_r_1'])))
all(na4 %in% na3) # TRUE

# don't actually have the individual parameter values saved, but can see how the "data" look:
p2 <- ggplot(data = ab_titers, aes(x = time, y = value, color = subject)) + geom_line() +
  geom_vline(xintercept = 365, lty = 2) +
  geom_line(data = ab_titers[ab_titers$subject %in% na4, ], aes(x = time, y = value, group = subject), col = 'black', lwd = 1.0) +
  theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Ab Titers') +
  scale_x_continuous(breaks = seq(0, max(ab_titers$time), by = 30)) +
  scale_y_continuous(breaks = seq(0, max(ab_titers$value), by = 1.0))#, trans = 'log')
print(p2) # anomalously high - maybe a lower r_1?
p3 <- ggplot(data = ab_titers, aes(x = time, y = value, color = subject)) + geom_line() +
  geom_vline(xintercept = 365, lty = 2) +
  geom_line(data = ab_titers[ab_titers$subject %in% na3, ], aes(x = time, y = value, group = subject), col = 'black') +
  theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Ab Titers') +
  scale_x_continuous(breaks = seq(0, max(ab_titers$time), by = 30)) +
  scale_y_continuous(breaks = seq(0, max(ab_titers$value), by = 1.0))#, trans = 'log')
print(p3) # unclear what the pattern here might be

# I guess we don't have rho as a random effect, so we can actually hold that constant:
calculate_ab_titers_log_params_rhoconst <- function(time, log_alpha, log_m, log_beta, log_r_1, log_r_2) {
  v1 <- exp(log_alpha) * exp(-exp(log_m) * time)
  v2 <- exp(log_alpha) * exp(-exp(log_m) * time) +
    exp(log_beta) * (0.75 * exp(-exp(log_r_1) * (time - (365 + 14))) +
                       0.25 * exp(-exp(log_r_2) * (time - (365 + 14))))
  value <- v1 * (time < (365 + 14)) + v2 * (time >= (365 + 14))
  return(value)
}
m5 <- nlsList(value ~ calculate_ab_titers_log_params_rhoconst(time, log_alpha, log_m, log_beta, log_r_1, log_r_2) | subject,
              data = ab_titers,
              start = c(log_alpha = log(8.0), log_m = log(log(2)/20), log_beta = log(18.0),
                        log_r_1 = log(log(2)/20), log_r_2 = log(log(2)/3650)))
plot(intervals(m5)) # this actually seems to get worse

# So maybe we need some constraint saying that r2 << r1? Also possible that a selfStart function would help








# m2 <- nlme(value ~ alpha * exp(-m * time) +
#              h_1 * beta * (rho * exp(-r_1 * (time - 379)) +
#                              (1 - rho) * exp(-r_2 * (time - 379))),
#            data = ab_titers,
#            fixed = alpha + m + beta + rho + r_1 + r_2 ~ 1,
#            random = alpha + m + beta + r_1 + r_2 ~ 1,
#            groups = ~subject,
#            start = c(5.0, log(2)/30, 18.0, 0.75, log(2)/30, log(2)/3650))

# # interestingly, seems to do better if we set m = r, but still only if very high-dimension (every 1-5 days) data:
# # for log scale: can use glmmPQL in MASS package instead (https://www.juliapilowsky.com/2018/10/19/a-practical-guide-to-mixed-models-in-r/)

# FIND A WAY TO CHECK AGAINST OBSERVED SDs, TOO
# TRY GEE; LME4; BRMS; RSTAN; ETC.
# ALTERNATIVE MODELS?
# REORGANIZE ROW ORDER
# IS IT IMPORTANT THAT TIME BE NUMERIC RATHER THAN INTEGER?
# consider weights/correlation


