
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
ab_titers <- ab_titers_ORIG[ab_titers_ORIG$time %in% c(seq(0, 360, by = 30), 379, seq(390, 720, by = 30)), ]
# ab_titers <- ab_titers_ORIG[ab_titers_ORIG$time %in% c(seq(0, 360, by = 30), 379, seq(390, 3650, by = 30)), ] # allows nlsList to converge

# And specify step function where vaccine occurs:
ab_titers$h_1 <- ifelse(ab_titers$time >= (365 + 14), 1, 0)

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
plot(intervals(m3)) # beta seems to be fit with most variability
# Including all 10 years of "data" allows m2 to converge!

pairs(m3) # alpha/m, r_1/rho pos correlated; rho/beta neg correlated?

# Fit nlme:
m4 <- nlme(value ~ calculate_ab_titers(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2),
           data = ab_titers,
           fixed = log_alpha + log_m + log_beta + logit_rho + log_r_1 + log_r_2 ~ 1,
           random = log_alpha + log_m + log_beta + log_r_1 + log_r_2 ~ 1,
           groups = ~subject,
           start = c(log_alpha = log(5.0), log_m = log(log(2)/30), log_beta = log(18.0), logit_rho = logit(0.75),
                     log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650)))

# Try with r2 held constant:
m5 <- nlme(value ~ calculate_ab_titers_log_params_r2const(time, log_alpha, log_m, log_beta, logit_rho, log_r_1),
           data = ab_titers,
           fixed = log_alpha + log_m + log_beta + logit_rho + log_r_1 ~ 1,
           random = log_alpha + log_beta ~ 1,
           groups = ~subject,
           start = c(log_alpha = log(5.0), log_m = log(log(2)/30), log_beta = log(18.0), logit_rho = logit(0.75),
                     log_r_1 = log(log(2)/30)))
# Error: Cholesky-Zerlegung kann nicht bestimmt werden: Führender Minor der Ordnung 1 ist nicht positiv definit
# converges if only alpha and beta allowed to vary randomly
summary(m1)
summary(m5) # standard errors actually seem to be mostly larger here - expect them to be smaller (b/c of r2?)
pairs(m5) # suggests perfect correlation between alpha and beta? - indicates overparameterized?

m6 <- update(m5, random = log_beta ~ 1)
anova(m5, m6) # statistically equivalent
summary(m6) # reasonable estimates, and good estimate of stddev of beta

# What if we start it out with the nlsList results?:
m6 <- nlme(m3)
# Error in solve.default(pdMatrix(a, factor = TRUE)) : system is computationally singular: reciprocal condition number = 1.21573e-39

# But also try "full" model with fewer random effects:
m7 <- nlme(value ~ calculate_ab_titers(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2),
           data = ab_titers,
           fixed = log_alpha + log_m + log_beta + logit_rho + log_r_1 + log_r_2 ~ 1,
           random = log_alpha + log_beta ~ 1,
           groups = ~subject,
           start = c(log_alpha = log(5.0), log_m = log(log(2)/30), log_beta = log(18.0), logit_rho = logit(0.75),
                     log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650)))
# alpha and beta only; otherwise:
# Error in nlme.formula(value ~ calculate_ab_titers(time, log_alpha, log_m, : Singularität in backsolve auf Stufe 0, Block 1
summary(m7)
pairs(m7)
m8 <- update(m7, random = log_beta ~ 1)
anova(m7, m8) # here including alpha is sig. better

# What about log scale?:
m9 <- nlme(log(value) ~ log(calculate_ab_titers(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2)),
           data = ab_titers,
           fixed = log_alpha + log_m + log_beta + logit_rho + log_r_1 + log_r_2 ~ 1,
           random = log_alpha + log_beta ~ 1,
           groups = ~subject,
           start = c(log_alpha = log(5.0), log_m = log(log(2)/30), log_beta = log(18.0), logit_rho = logit(0.75),
                     log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650)))
m7$coefficients$fixed
m9$coefficients$fixed
# similar estimates; a bit slower for r2 (~20 years rather than "true" 10)

m10 <- update(m9, random = log_alpha + log_beta + log_r_1 ~ 1)
anova(m9, m10) # no difference though, and r_1 seems correlated with other two
pairs(m10)

m11 <- update(m10, random = log_alpha + log_beta + log_r_1 + log_r_2 ~ 1)
m12 <- update(m10, random = log_alpha + log_m + log_beta + log_r_1 ~ 1)
anova(m9, m12) # m12 has negative AIC?

# m13 <- update(m11, random = log_alpha + log_m + log_beta + log_r_1 + log_r_2 ~ 1)
# anova(m9, m13); anova(m11, m13)
# # m13 slow to update - but no clear errors, so maybe this would also fit with time?

# Can this also be fit from scratch?:
m14 <- nlme(log(value) ~ log(calculate_ab_titers(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2)),
           data = ab_titers,
           fixed = log_alpha + log_m + log_beta + logit_rho + log_r_1 + log_r_2 ~ 1,
           random = log_alpha + log_m + log_beta + log_r_1 ~ 1,
           groups = ~subject,
           start = c(log_alpha = log(5.0), log_m = log(log(2)/30), log_beta = log(18.0), logit_rho = logit(0.75),
                     log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650)))
# note: done here w/o random effect for r2; this is the same as m12

# Create fxn that itself returns the log:
calculate_ab_titers_LOG <- function(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2) {
  alpha = exp(log_alpha)
  m = exp(log_m)
  beta = exp(log_beta)
  rho = exp(logit_rho) / (exp(logit_rho) + 1)
  r_1 = exp(log_r_1)
  r_2 = exp(log_r_2)
  
  v1 <- alpha * exp(-m * time)
  v2 <- alpha * exp(-m * time) + beta * (rho * exp(-r_1 * (time - (365 + 14))) + (1 - rho) * exp(-r_2 * (time - (365 + 14))))
  value <- v1 * (time < (365 + 14)) + v2 * (time >= (365 + 14))
  return(log(value))
}

m14 <- nlme(log(value) ~ calculate_ab_titers_LOG(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2),
            data = ab_titers,
            fixed = log_alpha + log_m + log_beta + logit_rho + log_r_1 + log_r_2 ~ 1,
            random = log_alpha + log_m + log_beta + log_r_1 ~ 1,
            groups = ~subject,
            start = c(log_alpha = log(5.0), log_m = log(log(2)/30), log_beta = log(18.0), logit_rho = logit(0.75),
                      log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650)))
# same fits as m12

# Check model assumptions:
plot(m14) # no pattern around zero, but fits seem to have much more negative range than positive
# fitting on log scale and including all these random effects definitely looks better! also looked funny when we held r2 constant, but could be b/c not log scale yet
# plot(augPred(m14, level = 0:1))
qqnorm(m14, abline = c(0, 1)) # looks mostly fine; again, other models don't look good

# Continue to explore random effects structure:
m15 <- update(m14, random = pdDiag(log_alpha + log_m + log_beta + log_r_1 ~ 1))
anova(m14, m15) # not sig different

# Can get log_r_2 in there if I use this!:
m16 <- update(m14, random = pdDiag(log_alpha + log_m + log_beta + log_r_1 + log_r_2 ~ 1))
anova(m14, m16) # not actually sig better, but still
qqnorm(m16, ~ ranef(.))

# for varClasses: options are varExp, varPower, varConstPower, varConstProp, varIdent, varFixed, varComb
m17 <- update(m16, weights = varPower())
anova(m16, m17)
# all either make no difference or make things worse, but can certainly explore once we have the real data

# for corClasses: corAR1, corARMA, corCAR1, corCompSymm (rest seem to be spatial?: corExp, corGaus, corLin, corRatio, corSpher, corSymm)
plot(ACF(m16, maxLag = 10), alpha = 0.05) # most actually seem to exceed the 0.05 confidence level

m17 <- update(m16, correlation = corCAR1(0.134))
anova(m16, m17)
# again, none seem to improve the model - try again with real data

# What about gnls?:
m17 <- gnls(model = log(value) ~ calculate_ab_titers_LOG(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2),
            data = ab_titers,
            params = log_alpha + log_m + log_beta + logit_rho + log_r_1 + log_r_2 ~ 1,
            start = c(log_alpha = log(5.0), log_m = log(log(2)/30), log_beta = log(18.0), logit_rho = logit(0.75),
                      log_r_1 = log(log(2)/30), log_r_2 = log(log(2)/3650)))

plot(m17) # variance seems to decrease at higher values
plot(ACF(m17, form = ~ 1 | subject), alpha = 0.05)

m18 <- update(m17, corr = corAR1(0.814, form = ~1 | subject))
plot(ACF(m18, form = ~ 1 | subject), alpha = 0.05)
anova(m17, m18) # m18 is much better, but ACF looks very similar?

m19 <- update(m18, weights = varPower())
anova(m18, m19) # also sig better

anova(m16, m19, test = FALSE) # but nlme still seems to fit much better - will have to explore for real data, but stick to nlme for now

# m20 <- nlme(value ~ alpha * exp(-m * time) +
#              h_1 * beta * (rho * exp(-r_1 * (time - 379)) +
#                              (1 - rho) * exp(-r_2 * (time - 379))),
#            data = ab_titers,
#            fixed = alpha + m + beta + rho + r_1 + r_2 ~ 1,
#            random = alpha + m + beta + r_1 + r_2 ~ 1,
#            groups = ~subject,
#            start = c(5.0, log(2)/30, 18.0, 0.75, log(2)/30, log(2)/3650))
# # interesting that this fits, but m4 above has trouble converging - only difference seems to be that m4 reparameterizes everything to be positive
# m20$coefficients # some of these are negative, which shouldn't be permitted
