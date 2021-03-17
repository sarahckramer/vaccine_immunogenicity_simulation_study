# ---------------------------------------------------------------------------------------------------------------------
# Functions to observe and analyze model fit
# ---------------------------------------------------------------------------------------------------------------------

get_param_est <- function(m, seasonal = FALSE, mono = FALSE, r2const = FALSE) {
  # Function to format and return estimates from an nlme model fit
  # param m: The fitted model object
  # param seasonal: Is seasonality included in the model?
  # param mono: Is model being fit mono-exponential (as opposed to biexponential)?
  # param r2const: Is r2 held constant in fitting procedure?
  # returns: A data frame containing parameter estimates, 95% confidence intervals, and st. devs. of random effects
  
  # Specify random effects:
  if (seasonal) {
    random.parms <- c('log_beta0', 'log_r_1', 'log_r_2')
    log.parms <- c('log_beta0', 'log_r_1', 'log_r_2')
    logit.parms <- c('logit_beta1', 'phi_hat', 'logit_rho')
    
    if (r2const) {
      random.parms <- c('log_beta0', 'log_r_1')
      log.parms <- c('log_beta0', 'log_r_1')
    }
    
    if (mono) {
      random.parms <- c('log_beta0', 'log_r')
      log.parms <- c('log_beta0', 'log_r')
      logit.parms <- c('logit_beta1', 'phi_hat')
    }
    
  } else {
    random.parms <- c('log_beta', 'log_r_1', 'log_r_2')
    log.parms <- c('log_beta', 'log_r_1', 'log_r_2')
    logit.parms <- 'logit_rho'
    
    if (mono) {
      random.parms <- c('log_beta', 'log_r')
      log.parms <- c('log_beta', 'log_r')
      logit.parms <- NULL
    }
  }
  
  # Extract results from model:
  res.df <- as.data.frame(summary(m)$tTable)
  res.sd <- as.numeric(VarCorr(m)[random.parms, 'StdDev'])
  
  # Remove columns not of interest:
  res.df <- res.df[, 1:2]
  
  # Convert point estimates to scale of interest:
  res.df[log.parms, 1] <- exp(res.df[log.parms, 1]) # exponentiate all params on log scale
  
  # Calculate confidence intervals for fixed effects (use delta method):
  res.df[log.parms, 'lower'] <- res.df[log.parms, 1] - 1.96 * res.df[log.parms, 1] * res.df[log.parms, 2]
  res.df[log.parms, 'upper'] <- res.df[log.parms, 1] + 1.96 * res.df[log.parms, 1] * res.df[log.parms, 2]
  
  res.df[logit.parms, 'lower'] <- plogis(res.df[logit.parms, 1]) - 1.96 *
    (exp(res.df[logit.parms, 1]) / (exp(res.df[logit.parms, 1]) + 1) ** 2) * res.df[logit.parms, 2]
  res.df[logit.parms, 'upper'] <- plogis(res.df[logit.parms, 1]) + 1.96 *
    (exp(res.df[logit.parms, 1]) / (exp(res.df[logit.parms, 1]) + 1) ** 2) * res.df[logit.parms, 2]
  
  # And convert logit-transformed params back to natural scale:
  res.df[logit.parms, 1] <- plogis(res.df[logit.parms, 1]) # inverse logit of all params on logit scale
  
  # Finish conversion for phi:
  if (seasonal) {
    res.df['phi_hat', c(1, 3:4)] <- res.df['phi_hat', c(1, 3:4)] * 12
  }
  
  # Remove standard errors:
  res.df <- res.df[, c(1, 3:4)]
  
  # Add in st. dev. of random effects:
  res.df[random.parms, 'sd_random'] <- res.sd
  
  # Correct rownames to no longer say "log":
  if (seasonal) {
    if (r2const) {
      rownames(res.df) <- c('beta0', 'beta1', 'phi', 'r_short', 'rho')
    } else if (mono) {
      rownames(res.df) <- c('beta0', 'beta1', 'phi', 'r')
    } else {
      rownames(res.df) <- c('beta0', 'beta1', 'phi', 'r_short', 'r_long', 'rho')
    }
  } else {
    if (mono) {
      rownames(res.df) <- c('beta', 'r')
    } else {
      rownames(res.df) <- c('beta', 'r_short', 'r_long', 'rho')
    }
  }
  
  # Return results:
  return(res.df)
}
