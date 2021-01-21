### Functions describing antibody titers over time ###

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
  rho = exp(logit_rho) / (exp(logit_rho) + 1)
  r_1 = exp(log_r_1)
  r_2 = exp(log_r_2)
  
  v1 <- alpha * exp(-m * time)
  v2 <- alpha * exp(-m * time) + beta * (rho * exp(-r_1 * (time - (365 + 14))) + (1 - rho) * exp(-r_2 * (time - (365 + 14))))
  value <- v1 * (time < (365 + 14)) + v2 * (time >= (365 + 14))
  
  # value[is.na(value)] <- 1e-16 # ? Seems necessary for saemix, but impairs nlme
  
  return(log(value))
}


calculate_ab_titers <- function(time, log_alpha, log_m, log_beta, logit_rho, log_r_1, log_r_2) {
  # Calculates antibody titers over time using non-mechanistic, bi-exponential model
  # param time: Time in days
  # param log_alpha: Natural log of the initial amount of maternal antibodies
  # param log_m: Natural log of the rate of decay of maternal antibodies
  # param log_beta: Natural log of the boost in antibodies 2 weeks after vaccination
  # param logit_rho: Logit of the proportion of antibodies decaying at faster rate r_1
  # param log_r_1: Natural log of the rate of antibody decay (short-lived)
  # param log_r_2: Natural log of the rate of antibody decay (long-lived)
  # returns: Simulated antibody titers over time
  
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


calculate_ab_titers_log_params_r2const <- function(time, log_alpha, log_m, log_beta, logit_rho, log_r_1) {
  # Calculates antibody titers over time using non-mechanistic, bi-exponential model, holding r_2 constant (r2 = 3650.0)
  # param time: Time in days
  # param log_alpha: Natural log of the initial amount of maternal antibodies
  # param log_m: Natural log of the rate of decay of maternal antibodies
  # param log_beta: Natural log of the boost in antibodies 2 weeks after vaccination
  # param logit_rho: Logit of the proportion of antibodies decaying at faster rate r_1
  # param log_r_1: Natural log of the rate of antibody decay (short-lived)
  # returns: Simulated antibody titers over time
  
  v1 <- exp(log_alpha) * exp(-exp(log_m) * time)
  v2 <- exp(log_alpha) * exp(-exp(log_m) * time) +
    exp(log_beta) * ((exp(logit_rho) / (exp(logit_rho) + 1)) * exp(-exp(log_r_1) * (time - (365 + 14))) +
                       (1 - (exp(logit_rho) / (exp(logit_rho) + 1))) * exp(-exp(log(2)/3650) * (time - (365 + 14))))
  value <- v1 * (time < (365 + 14)) + v2 * (time >= (365 + 14))
  return(value)
}


logit <- function(p) {
  # Returns the logit of input parameter, p
  return(log(p / (1 - p)))
}
