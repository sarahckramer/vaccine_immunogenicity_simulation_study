### Draws 10 sets of fixed/median parameter values using Latin hypercube sampling to use in generating synthetic data ###

# Ensure reproducibility:
set.seed(40764301)

# Read in libraries:
library(tgp)

# Set parameter bounds:
lower <- c(10.0, 0.05, 0.00, 0.50, 10.0, 365.0)
upper <- c(25.0, 0.50, 12.0, 0.95, 60.0, 7300.0)
# order: beta0, beta1, phi, rho, log(2)/r1, log(2)/r2
param_bound <- cbind(lower, upper)

# Draw using LHS:
input_params <- lhs(10, param_bound)

# Round phi to integers:
input_params[,3] <- floor(input_params[, 3])

# Output results:
if (!dir.exists('data/')) {
  dir.create('data/')
}
if (!dir.exists('data/inputs/')) {
  dir.create('data/inputs/')
}
write.table(input_params, file = 'data/inputs/lhs_params.txt', sep = '\t', row.names = FALSE, col.names = FALSE)
