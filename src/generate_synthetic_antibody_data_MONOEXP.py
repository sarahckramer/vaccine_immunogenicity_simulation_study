# Run model with random effects by participant #

# Import necessary functions:
import os
import pandas as pd
import matplotlib.pyplot as plt
from functions_python import *

#######################################################################################################################

# Set date:
ymd = '20210225'

# Set global parameters:
N_pop = 5000  # number of "participants"
response_delay = 14  # 2 week delay in Ab response
beta_1 = np.float64(0.20)  # extent of seasonal variability in beta
phi = np.float64(1.0)  # maximum in January

# Set parameter medians:
maternal_antibodies_median = 8.0
beta_0s_median = np.array([18.0])  # median avg. titer; single dose for now
half_life_maternal_median = 42.0  # days
half_life_ab_median = 365.0  # days

# Get standard deviation from desired coefficient of variations (var/mean of distribution):
sd = calculate_sigma_from_cv(0.2)

# Set start, end, and vaccination timepoints:
tm_start = 0
tm_end = 3652  # 730  # 10 years
vacc_timepoints = np.array([365])  # vaccinate between 12 and 15 months; can be second dose at 4-6 yrs

#######################################################################################################################

# Select birth month (1-12) for each participant:
rng = np.random.default_rng()
birth_months = rng.integers(low=1, high=13, size=N_pop)

# Convert to month of first vaccination based on vacc_timepoints:
vacc_months = birth_months + np.floor(vacc_timepoints[0] / 30)

# Set values of random effects for each participant:
maternal_antibodies = generate_random_effects(maternal_antibodies_median, sd, N_pop)
half_life_maternal = generate_random_effects(half_life_maternal_median, sd, N_pop)
half_life_ab = generate_random_effects(half_life_ab_median, sd, N_pop)

# Check sd is correct:
print(np.std(maternal_antibodies) / np.mean(maternal_antibodies))
print(np.std(half_life_maternal) / np.mean(half_life_maternal))
print(np.std(half_life_ab) / np.mean(half_life_ab))
# print(np.std(get_rate_from_half_life(half_life_long)))  # want: 0.010?

# And get distribution for each value of beta_0, as well:
beta_0s = np.zeros([len(beta_0s_median), N_pop])
for i in range(len(beta_0s_median)):
    beta_0s[i] = generate_random_effects(beta_0s_median[i], sd, N_pop)
    print(np.std(beta_0s[i]) / np.mean(beta_0s[i]))  # check sd

# Calculate beta for each participant (random "median" value + seasonal effect):
betas = beta_0s * (1 + beta_1 * np.cos((2 * np.pi / 12) * (vacc_months - phi)))

#######################################################################################################################

# Run!
sim_titers = calculate_Ab_titers_monoexp(tm_start, tm_end, vacc_timepoints, maternal_antibodies, betas,
                                         half_life_maternal, half_life_ab, N_pop, response_delay)

# # Plot simulated "data":
# if not os.path.isdir('results/PRELIM_plotTiters_' + ymd + '/'):
#     os.mkdir('results/PRELIM_plotTiters_' + ymd + '/')
#
# plt.figure()
# plt.plot(sim_titers)
# plt.yscale('log')
# plt.xlabel('Time (Days)')
# plt.ylabel('Ab Titers')
# plt.tight_layout()
# plt.savefig('results/PRELIM_plotTiters_' + ymd + '/ab_titers_over_time.png', dpi=300)

# Write "true" values to file:
if not os.path.isdir('data/'):
    os.mkdir('data/')
if not os.path.isdir('data/prelim_check_' + ymd + '/'):
    os.mkdir('data/prelim_check_' + ymd + '/')

true_vals = pd.DataFrame(sim_titers)
true_vals.to_csv('data/prelim_check_' + ymd + '/truth_MONO.csv', na_rep='NA', index=False)
del true_vals

# Add random noise:
sim_titers = add_random_noise(sim_titers, 0.3)

# Write synthetic data to file:
obs_vals = pd.DataFrame(sim_titers)
obs_vals.to_csv('data/prelim_check_' + ymd + '/obs_data_MONO.csv', na_rep='NA', index=False)

# Write month of vaccine to file:
vacc_months = pd.DataFrame(vacc_months)
vacc_months.to_csv('data/prelim_check_' + ymd + '/vacc_months_MONO.csv', na_rep='NA', index=False, header=False)