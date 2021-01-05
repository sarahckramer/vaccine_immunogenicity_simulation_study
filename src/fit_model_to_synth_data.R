
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
# set.seed(3970395)

# Read in noise-laden "data":
ab_titers <- read.csv('data/prelim_check_20210104/obs_data.csv')

# Reformat data frame:
ab_titers$time <- 0:(dim(ab_titers)[1] - 1)
ab_titers$time <- as.numeric(as.character(ab_titers$time))
ab_titers <- melt(ab_titers, id.vars = 'time')
names(ab_titers)[2] <- 'subject'

# Select subset of subjects?
ab_titers <- ab_titers[ab_titers$subject %in% levels(ab_titers$subject)[sample(1:1000, 100)], ]
ab_titers$subject <- factor(ab_titers$subject)

# Select key timepoints:
ab_titers_ORIG <- ab_titers # store "full" data
ab_titers <- ab_titers[ab_titers$time %in% c(0, 15, 30, 45, 60), ]

# Plot "data":
p1 <- ggplot(data = ab_titers, aes(x = time, y = value, color = subject)) + geom_line() +
  theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Ab Titers') + 
  scale_x_continuous(breaks = seq(0, 60, by = 15)) + scale_y_continuous(breaks = seq(0, 0.4, by = 0.05))#, trans = 'log')
print(p1)

### Explore possibility of fitting ONLY before vaccination ###
m1 <- nlme(log(value) ~ log(A_m) - m * time,
           data = ab_titers,
           fixed = A_m + m ~ 1,
           random = A_m + m ~ 1,
           groups = ~subject,
           method = 'ML',
           start = c(0.14, log(2)/30))
m2 <- nlme(log(value) ~ log(A_m) - m * time,
           data = ab_titers,
           fixed = A_m + m ~ 1,
           random = A_m + m ~ 1,
           groups = ~subject,
           method = 'ML',
           weights = varPower(),
           start = c(0.14, log(2)/30))
m3 <- nlme(log(value) ~ log(A_m) - m * time,
           data = ab_titers,
           fixed = A_m + m ~ 1,
           random = A_m + m ~ 1,
           groups = ~subject,
           method = 'ML',
           correlation = corCAR1(),
           start = c(0.14, log(2)/30))
m4 <- nlme(log(value) ~ log(A_m) - m * time,
           data = ab_titers,
           fixed = A_m + m ~ 1,
           random = A_m + m ~ 1,
           groups = ~subject,
           method = 'ML',
           weights = varPower(),
           correlation = corCAR1(),
           start = c(0.14, log(2)/30))
anova.lme(m1, m2) # sig difference - second model has higher logLik; note that sometimes not sig - depends on random sample
anova.lme(m1, m3) # no sig difference - adding correlation structure does not seem to help
anova.lme(m2, m4) # no sig difference
# if doing all on log-scale, no significant differences anywhere

viz_model_fit(m1, ab_titers, logscale = TRUE)
viz_model_fit(m2, ab_titers, logscale = TRUE)
# that said, it seems like the residuals show a pattern when "weights" are used, and don't otherwise; m1 also seems to better capture extreme values?

# use m1/m2 moving forward
rm(m3, m4)

# Comparing log to not log scale: residuals look better
# but note that this is much simpler before incorporating vaccination b/c easy way to apply log to right side of the equation

# # With how few data points will this work?:
# ab_titers <- ab_titers[ab_titers$time %in% c(0, 30, 60), ]
# m5 <- nlme(log(value) ~ log(A_m) - m * time,
#            data = ab_titers,
#            fixed = A_m + m ~ 1,
#            random = A_m + m ~ 1,
#            groups = ~subject,
#            method = 'ML',
#            weights = varPower(),
#            start = c(0.14, log(2)/30))
# viz_model_fit(m5, ab_titers, logscale = T)
# # this sometimes converges, sometimes not - must depend on randomly-selected 100 participants
# 
# ab_titers <- ab_titers[ab_titers$time %in% c(0, 60), ]
# m6 <- nlme(log(value) ~ log(A_m) - m * time,
#            data = ab_titers,
#            fixed = A_m + m ~ 1,
#            random = A_m + m ~ 1,
#            groups = ~subject,
#            method = 'ML',
#            weights = varPower(),
#            start = c(0.14, log(2)/30))
# viz_model_fit(m6, ab_titers, logscale = T)
# # same here - technically linear when log-transformed, so I guess that makes sense, but often doesn't converge - and we are unlikely to get more than 1 data point from this period...

# Clean up:
rm(ab_titers, p1, m1, m2)

### Move on to trying to fit including data post-first primary dose ###

# First, get increased time range from data:
# ab_titers <- ab_titers_ORIG[ab_titers_ORIG$time %in% c(0, 15, 30, 45, 60, 75, 90, 105, 120), ]
ab_titers <- ab_titers_ORIG[ab_titers_ORIG$time %in% seq(0, 120, by = 1), ]

# And specify step function where vaccine occurs:
ab_titers$h_1 <- ifelse(ab_titers$time >= 60, 1, 0)

# Plot data:
p1 <- ggplot(data = ab_titers, aes(x = time, y = value, color = subject)) + geom_line() +
  theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Ab Titers') + 
  scale_x_continuous(breaks = seq(0, 120, by = 15)) + scale_y_continuous(breaks = seq(0, 0.5, by = 0.05))#, trans = 'log')
print(p1)

# Now begin trying to fit:
m1 <- nlme(value ~ A_m * exp(-m * time) +
             h_1 * (beta_1 * ((rho / (r - c_s)) * (exp(-c_s * (time - 60)) - exp(-r * (time - 60))) +
                                ((1 - rho)/(r - c_l)) * (exp(-c_l * (time - 60)) - exp(-r * (time - 60))))),
           data = ab_titers,
           fixed = A_m + m + beta_1 + rho + r + c_s + c_l ~ 1,
           # random = A_m + m + beta_1 + rho + r + c_s + c_l ~ 1,
           random = beta_1 ~ 1,
           groups = ~subject,
           method = 'ML',
           weights = varPower(),
           start = c(0.20, log(2)/46, 0.035, 0.98, log(2)/20, log(2)/5.5, log(2)/1650))
# no convergence - "singularity in backsolve"
# or: Error in nlme.formula(value ~ A_m * exp(-m * time) + h_1 * (beta_1 * ((rho/(r -  : step halving factor reduced below minimum in PNLS step
# doesn't converge even with all data points included - not sure if the issue is too few data points, rapid increase at time t=tau_1, both, or something else
# and note that this is with only using some of the random effects

# interestingly, seems to do better if we set m = r, but still only if very high-dimension (every 1-5 days) data:
m2 <- nlme(value ~ A_m * exp(-r * time) +
             h_1 * (beta_1 * ((rho / (r - c_s)) * (exp(-c_s * (time - 60)) - exp(-r * (time - 60))) +
                                ((1 - rho)/(r - c_l)) * (exp(-c_l * (time - 60)) - exp(-r * (time - 60))))),
           data = ab_titers,
           fixed = A_m + r + beta_1 + rho + c_s + c_l ~ 1,
           # random = A_m + m + beta_1 + rho + r + c_s + c_l ~ 1,
           random = A_m + beta_1 + r ~ 1,
           groups = ~subject,
           method = 'ML',
           weights = varPower(),
           start = c(0.20, log(2)/20, 0.035, 0.98, log(2)/5.5, log(2)/1650))
viz_model_fit(m2, ab_titers)
# however, even this is fitting rho as greater than 1, which shouldn't be allowed

# Re-parameterize rho so that forced between 0 and 1:
m3 <- nlme(value ~ A_m * exp(-r * time) +
             h_1 * (beta_1 * (((exp(rho_hat) / (exp(rho_hat) + 1)) / (r - c_s)) * (exp(-c_s * (time - 60)) - exp(-r * (time - 60))) +
                                ((1 - (exp(rho_hat) / (exp(rho_hat) + 1)))/(r - c_l)) * (exp(-c_l * (time - 60)) - exp(-r * (time - 60))))),
           data = ab_titers,
           fixed = A_m + r + beta_1 + rho_hat + c_s + c_l ~ 1,
           # random = A_m + m + beta_1 + rho + r + c_s + c_l ~ 1,
           random = A_m ~ 1,
           groups = ~subject,
           method = 'ML',
           weights = varPower(),
           start = c(0.20, log(2)/20, 0.035, log(0.98 / (1 - 0.98)), log(2)/5.5, log(2)/1650))
# won't converge - "singularity"

# Try to fit on log scale instead:
m4 <- nlme(log(value) ~ log(A_m * exp(-r * time) +
                              h_1 * (beta_1 * ((rho / (r - c_s)) * (exp(-c_s * (time - 60)) - exp(-r * (time - 60))) +
                                                 ((1 - rho)/(r - c_l)) * (exp(-c_l * (time - 60)) - exp(-r * (time - 60)))))),
           data = ab_titers,
           fixed = A_m + r + beta_1 + rho + c_s + c_l ~ 1,
           # random = A_m + m + beta_1 + rho + r + c_s + c_l ~ 1,
           random = A_m + beta_1 + r ~ 1,
           groups = ~subject,
           method = 'ML',
           weights = varPower(),
           start = c(0.20, log(2)/20, 0.035, 0.98, log(2)/5.5, log(2)/1650))
viz_model_fit(m4, ab_titers, logscale = TRUE)
# still fits rho > 1; but not at least that residuals look much better - I imagine this is the way to go

m5 <- nlme(log(value) ~ log(A_m * exp(-r * time) +
                              h_1 * (beta_1 * (((exp(rho_hat) / (exp(rho_hat) + 1)) / (r - c_s)) * (exp(-c_s * (time - 60)) - exp(-r * (time - 60))) +
                                                 ((1 - (exp(rho_hat) / (exp(rho_hat) + 1)))/(r - c_l)) * (exp(-c_l * (time - 60)) - exp(-r * (time - 60)))))),
           data = ab_titers,
           fixed = A_m + r + beta_1 + rho_hat + c_s + c_l ~ 1,
           # random = A_m + m + beta_1 + rho + r + c_s + c_l ~ 1,
           random = A_m ~ 1,
           groups = ~subject,
           method = 'ML',
           weights = varPower(),
           start = c(0.20, log(2)/20, 0.035, log(0.98 / (1 - 0.98)), log(2)/5.5, log(2)/1650))
# does not converge
# alternative: can use glmmPQL in MASS package instead (https://www.juliapilowsky.com/2018/10/19/a-practical-guide-to-mixed-models-in-r/)

# Finally, try writing out an explicit piecewise function first:
# calculate_ab_titers <- function(time, A_m, m, beta_1, rho, r, c_s, c_l) {
#   v1 <- A_m * exp(-m * time)
#   v2 <- A_m * exp(-m * time) + beta_1 * ((rho / (r - c_s)) * (exp(-c_s * (time - 60)) - exp(-r * (time - 60))) +
#                                            ((1 - rho)/(r - c_l)) * (exp(-c_l * (time - 60)) - exp(-r * (time - 60))))
#   value <- v1 * (time < 60) + v2 * (time >= 60)
#   return(value)
# }

calculate_ab_titers <- function(time, A_m, r, beta_1, rho, c_s, c_l) {
  v1 <- A_m * exp(-r * time)
  v2 <- A_m * exp(-r * time) + beta_1 * ((rho / (r - c_s)) * (exp(-c_s * (time - 60)) - exp(-r * (time - 60))) +
                                           ((1 - rho)/(r - c_l)) * (exp(-c_l * (time - 60)) - exp(-r * (time - 60))))
  value <- v1 * (time < 60) + v2 * (time >= 60)
  return(value)
}

m6 <- nlme(value ~ calculate_ab_titers(time, A_m, r, beta_1, rho, c_s, c_l),
           data = ab_titers,
           fixed = A_m + r + beta_1 + rho + c_s + c_l ~ 1,
           # random = A_m + m + beta_1 + rho + r + c_s + c_l ~ 1,
           random = A_m + beta_1 + r ~ 1,
           groups = ~subject,
           method = 'ML',
           weights = varPower(),
           start = c(0.20, log(2)/20, 0.035, 0.98, log(2)/5.5, log(2)/1650))
viz_model_fit(m6, ab_titers)
# again, rho > 1; also seems to yield same results as model coded with h_1




# FIND A WAY TO CHECK AGAINST OBSERVED SDs, TOO
# TRY GEE; LME4; BRMS; RSTAN; ETC.
# ALTERNATIVE MODELS?
# REORGANIZE ROW ORDER
# IS IT IMPORTANT THAT TIME BE NUMERIC RATHER THAN INTEGER?





# # Try a simpler model?:
# m3 <- nlme(value ~ A_m * exp(-m * time) + h_1 * (alpha_1 * exp(-r * (time - 75))),
#            data = ab_titers,
#            fixed = A_m + m + alpha_1 + r ~ 1,
#            random = A_m + alpha_1 ~ 1,
#            groups = ~subject,
#            # method = 'REML',
#            start = c(0.20, log(2)/46, 1.0, log(2)/45))
# # but also seems way too high at t=60... why? Oh! b/c here it's the Ab titers that are generated in boosts, not the ASCs!
# # still not working
# # works if only use alpha as random effect (or only alpha and one other param)
# # with more data points, works with all random effects included! (although still throws warning m)
# ab_titers$fitted <- fitted(m3)
# ggplot(data = ab_titers) + geom_line(aes(x = time, y = value, color = subject)) + geom_point(aes(x = time, y = fitted, color = subject)) + theme(legend.position = 'None')
# # this one will fit fine if we remove the maternal part and only model times 75-120 (decreasing )

# See https://dpmartin42.github.io/posts/Piecewise-growth for simply plotting different slopes over time (non-mechanistic)
