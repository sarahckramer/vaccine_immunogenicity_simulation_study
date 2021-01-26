### Functions to observe and analyze model fit ###

viz_model_fit <- function(m, dat, logscale = F) {
  # Function to plot the fitted values and residuals of an nlme object
  # param m: The fitted model object
  # param dat: The data frame used for model fitting
  # param logscale: Was the model fit on a log scale? (Boolean)
  
  if (logscale == T) {
    dat$fitted <- exp(fitted(m))
  } else {
    dat$fitted <- fitted(m)
  }
  
  dat$residuals_fixed <- m$residuals[, 1]
  dat$residuals_subject <- m$residuals[, 2]
  
  p1 <- ggplot(data = dat) + geom_line(aes(x = time, y = value, color = subject)) + geom_point(aes(x = time, y = fitted, color = subject)) +
    theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Ab Titers') + 
    scale_x_continuous(breaks = seq(0, max(dat$time), by = 60)) + scale_y_continuous(breaks = seq(0, max(dat$value), by = 1.0))#, trans = 'log')
  
  p21 <- ggplot(data = dat) + geom_point(aes(x = time, y = residuals_fixed, color = subject)) +
    theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Residuals (Fixed)') +
    scale_x_continuous(breaks = seq(0, max(dat$time), by = 60))
  p22 <- ggplot(data = dat) + geom_point(aes(x = time, y = residuals_subject, color = subject)) +
    theme_classic() + theme(legend.position = 'None') + labs(x = 'Time (Days)', y = 'Residuals (Subject)') +
    scale_x_continuous(breaks = seq(0, max(dat$time), by = 60))
  
  grid.arrange(p1, p21, p22, layout_matrix = rbind(c(1, 1), c(1, 1), c(2, 3)))
}


get_param_est <- function(m) {
  # Function to format and return estimates from an nlme model fit
  # param m: The fitted model object
  # returns: A data frame containing parameter estimates, 95% confidence intervals, and st. devs. of random effects
  
  # Extract results from model:
  res.df <- as.data.frame(summary(m)$tTable)
  res.sd <- as.numeric(VarCorr(m)[c('log_alpha', 'log_m', 'log_beta', 'log_r_1', 'log_r_2'), 'StdDev'])
  
  # Remove columns not of interest:
  res.df <- res.df[, 1:2]
  
  # Convert point estimates to scale of interest:
  res.df[c(1:3, 5:6), 1] <- exp(res.df[c(1:3, 5:6), 1]) # exponentiate all except rho
  res.df[5, 1] <- res.df[5, 1] + res.df[6, 1] # although we fit r1+r2, the original model was run using short- and long-term decay rates
  
  # Calculate confidence intervals for fixed effects (use delta method):
  res.df$lower <- res.df$Value - 1.96 * res.df$Value * res.df$Std.Error
  res.df$upper <- res.df$Value + 1.96 * res.df$Value * res.df$Std.Error
  res.df$lower[4] <- exp(res.df[4, 1]) / (exp(res.df[4, 1]) + 1) - 1.96 * (exp(res.df[4, 1]) / (exp(res.df[4, 1]) + 1)**2) * res.df[4, 2]
  res.df$upper[4] <- exp(res.df[4, 1]) / (exp(res.df[4, 1]) + 1) + 1.96 * (exp(res.df[4, 1]) / (exp(res.df[4, 1]) + 1)**2) * res.df[4, 2]
  
  # And convert rho back to natural scale:
  res.df[4, 1] <- exp(res.df[4, 1]) / (exp(res.df[4, 1]) + 1) # inverse logit of rho
  
  # Remove standard errors:
  res.df <- res.df[, c(1, 3:4)]
  
  # Convert rates to half lives:
  res.df[c(2, 5:6), ] <- log(2) / res.df[c(2, 5:6), ] # convert rates to half-lives
  res.df[c(2, 5:6), 2:3] <- res.df[c(2, 5:6), 3:2] # correct upper and lower bounds for half-lives
  
  # Add in st. dev. of random effects:
  res.df$sd_random <- NA
  res.df$sd_random[c(1:3, 5:6)] <- res.sd
  
  # Correct rownames to no longer say "log":
  rownames(res.df) <- c('alpha', 'm', 'beta', 'rho', 'r_short', 'r_long')
  
  # Return results:
  return(res.df)
}
