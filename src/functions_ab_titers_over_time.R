### Functions describing antibody titers over time ###

calculate_ab_titers_LOG_postOnly <- function(time, log_beta, logit_rho, log_r_1, log_r_2) {
  # Calculates log of antibody titers over time using non-mechanistic, bi-exponential model, assuming maternal Ab
  # negligible at vaccine timepoint
  # param time: Time in days, with day of vaccine = 0
  # param log_beta: Natural log of the boost in antibodies 2 weeks after vaccination
  # param logit_rho: Logit of the proportion of antibodies decaying at faster rate r_1
  # param log_r_1: Natural log of the rate of antibody decay (short-lived)
  # param log_r_2: Natural log of the rate of antibody decay (long-lived)
  # returns: The log of simulated antibody titers over time
  
  beta = exp(log_beta)
  rho = plogis(logit_rho)
  r_1 = exp(log_r_1)
  r_2 = exp(log_r_2)
  
  value <- beta * (rho * exp(-(r_1 + r_2) * (time - 14)) + (1 - rho) * exp(-r_2 * (time - 14)))
  
  if (any(is.na(value))) {
    print('NAs in output!')
  }
  
  value[value == 0] <- 1e-323 # any lower will give NAs
  
  return(log(value))
}


calculate_ab_titers_LOG <- function(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2) {
  # Calculates log of antibody titers over time using non-mechanistic, bi-exponential model
  # param time: Time in days
  # param log_alpha: Natural log of the initial amount of maternal antibodies
  # param log_m: Natural log of the rate of decay of maternal antibodies
  # param log_beta: Natural log of the boost in antibodies 2 weeks after vaccination
  # param logit_rho: Logit of the proportion of antibodies decaying at faster rate r_1
  # param log_r_1: Natural log of the rate of antibody decay (short-lived)
  # param log_r_2: Natural log of the rate of antibody decay (long-lived)
  # returns: The log of simulated antibody titers over time
  
  alpha = exp(log_alpha)
  m = exp(log_m)
  beta = exp(log_beta)
  rho = plogis(logit_rho)
  r_1 = exp(log_r_1)
  r_2 = exp(log_r_2)
  
  v1 <- alpha * exp(-m * time)
  v2 <- alpha * exp(-m * time) + beta * (rho * exp(-(r_1 + r_2) * (time - (365 + 14))) + (1 - rho) * exp(-r_2 * (time - (365 + 14))))
  
  v1 <- v1 * (time < (365 + 14))
  v2 <- v2 * (time >= (365 + 14))
  
  v1[is.na(v1)] <- 0
  v2[is.na(v2)] <- 0
  
  value <- v1 + v2
  value[value == 0] <- 1e-323 # any lower will give NAs
  
  return(log(value))
}
