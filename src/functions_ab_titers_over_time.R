# ---------------------------------------------------------------------------------------------------------------------
# Functions describing antibody titers over time
# ---------------------------------------------------------------------------------------------------------------------

format_synth_data <- function(dat, vacc_dat, n_participants, t_samples) {
  # Formats synthetic antibody titer data to be used in model fitting
  # param dat: Data frame containing all timepoints (rows) for all participants (cols)
  # param vacc_dat: Data frame containing a single column with the month of vaccination for each participant
  # param n_participants: The number of participants to sample from the data
  # param t_samples: The timepoints at which sampling is performed
  # returns: A formatted data frame for use in model fitting
  
  dat$time <- 0:(dim(dat)[1] - 1)
  dat$time <- as.numeric(as.character(dat$time))
  dat <- melt(dat, id.vars = 'time')
  names(dat)[2] <- 'subject'
  dat <- dcast(subject ~ time, data = dat, value.var = 'value')
  dat$vacc_month <- vacc_dat$V1
  dat <- melt(dat, id.vars = c('subject', 'vacc_month'))
  names(dat)[3] <- 'time'
  dat$time <- as.numeric(as.character(dat$time))
  dat$vacc_month <- dat$vacc_month - 12
  dat <- dat[order(dat$subject), ]
  
  # Select subset of subjects:
  dat <- dat[dat$subject %in% levels(dat$subject)[sample(1:length(levels(dat$subject)), n_participants)], ]
  dat$subject <- factor(dat$subject)
  
  # Select key timepoints:
  dat <- dat[dat$time %in% t_samples, ]
  
  # Set time of vaccine = time 0:
  dat$time <- dat$time - 365
  
  # Convert to grouped data object:
  dat <- groupedData(value ~ time | subject, data = dat, outer = ~vacc_month)
  
  # Return formatted "data":
  return(dat)
}


calculate_ab_titers_LOG_postOnly_seasonal <- function(time, v_time, log_beta0, logit_beta1, phi_hat, log_r_1, log_r_2 = log(log(2) / 3650.0), logit_rho = Inf, mono = FALSE, r2const = FALSE) {
  # Calculates log of antibody titers over time using non-mechanistic, bi-exponential model, assuming maternal Ab
  # negligible at vaccine timepoint
  # param time: Time in days, with day of vaccine = 0
  # param v_time: Month of initial vaccination (1-12 = Jan-Dec)
  # param log_beta0: Natural log of the average boost in antibodies 2 weeks after vaccination
  # param log_beta1: Logit of the magnitude of seasonal variation in antibody boost
  # param phi_hat: A transformation of the month of maximal vaccine impact, such that phi (below) varies between 0 and 12
  # param log_r_1: Natural log of the rate of antibody decay (short-lived)
  # param log_r_2: Natural log of the rate of antibody decay (long-lived)
  # param logit_rho: Logit of the proportion of antibodies decaying at faster rate r_1
  # param mono: Should the data be fit to a mono- rather than bi-exponential model?
  # param r2const: Should the parameter r_2 be held constant, rather than fit?
  # returns: The log of simulated antibody titers over time

  log_r_1[log_r_1 > 500] <- 500 # otherwise r_1 is Inf
  log_r_2[log_r_2 > 500] <- 500 # otherwise r_2 is Inf
  
  beta0 = exp(log_beta0)
  beta1 = plogis(logit_beta1)
  phi = 12 / (1 + exp(-phi_hat)) # or: plogis(phi_hat) * 12
  r_1 = exp(log_r_1)
  r_2 = exp(log_r_2)
  rho = plogis(logit_rho)
  
  if (mono) {
    rho = 1.0
    r_2 = 0.0
  }

  beta = beta0 * (1 + beta1 * cos((2 * pi / 12) * (v_time - phi)))
  if (any(is.na(beta))) {
    print('NAs in beta!')
  }

  value <- beta * (rho * exp(-(r_1 + r_2) * (time - 14)) + (1 - rho) * exp(-r_2 * (time - 14)))
  if (any(is.na(value))) {
    print('NAs in output!')
  }
  value[value == 0] <- 1e-323 # any lower will give NAs

  return(log(value))
}


calculate_ab_titers_LOG_postOnly <- function(time, log_beta, log_r_1, log_r_2 = log(log(2) / 3650.0), logit_rho = Inf, mono = FALSE, r2const = FALSE) {
  # Calculates log of antibody titers over time using non-mechanistic, bi-exponential model, assuming maternal Ab
  # negligible at vaccine timepoint
  # param time: Time in days, with day of vaccine = 0
  # param log_beta: Natural log of the boost in antibodies 2 weeks after vaccination
  # param log_r_1: Natural log of the rate of antibody decay (short-lived)
  # param log_r_2: Natural log of the rate of antibody decay (long-lived)
  # param logit_rho: Logit of the proportion of antibodies decaying at faster rate r_1
  # param mono: Should the data be fit to a mono- rather than bi-exponential model?
  # param r2const: Should the parameter r_2 be held constant, rather than fit?
  # returns: The log of simulated antibody titers over time
  
  log_r_1[log_r_1 > 500] <- 500 # otherwise r_1 is Inf
  log_r_2[log_r_2 > 500] <- 500 # otherwise r_2 is Inf
  
  beta = exp(log_beta)
  r_1 = exp(log_r_1)
  r_2 = exp(log_r_2)
  rho = plogis(logit_rho)
  
  if (mono) {
    rho = 1.0
    r_2 = 0.0
  }
  
  value <- beta * (rho * exp(-(r_1 + r_2) * (time - 14)) + (1 - rho) * exp(-r_2 * (time - 14)))
  
  if (any(is.na(value))) {
    print('NAs in output!')
  }
  
  value[value == 0] <- 1e-323 # any lower will give NAs
  
  return(log(value))
}
