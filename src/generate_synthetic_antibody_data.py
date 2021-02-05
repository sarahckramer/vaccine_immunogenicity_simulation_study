# Run model with random effects by participant #

# Import necessary functions:
import os
import pandas as pd
import matplotlib.pyplot as plt
from functions_python import *

#######################################################################################################################

# Set date:
ymd = '20210202'

# Set global parameters:
N_pop = 1000  # number of "participants"
response_delay = 14  # 2 week delay in Ab response
prop_short = np.float64(0.70)  # fit this as fixed for now
beta_1 = np.float64(0.20)  # extent of seasonal variability in beta
phi = np.float64(1.0)  # maximum in January

# Set parameter medians:
maternal_antibodies_median = 8.0
beta_0s_median = np.array([18.0])  # median avg. titer; single dose for now
half_life_maternal_median = 42.0  # days
half_life_short_median = 30.0  # days
half_life_long_median = 3650.0  # days

# Get standard deviation from desired coefficient of variations (var/mean of distribution):
sd = calculate_sigma_from_cv(0.1)

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
half_life_short = generate_random_effects(half_life_short_median, sd, N_pop)
half_life_long = generate_random_effects(half_life_long_median, sd, N_pop)

# Check sd is correct:
print(np.std(maternal_antibodies) / np.mean(maternal_antibodies))
print(np.std(half_life_maternal) / np.mean(half_life_maternal))
print(np.std(half_life_short) / np.mean(half_life_short))
print(np.std(half_life_long) / np.mean(half_life_long))
# print(np.std(get_rate_from_half_life(half_life_long)))  # want: 0.010?

# # To keep prop_short between 0-1, use logit-normal instead:
# prop_short = generate_random_effects(prop_short_median, sd, N_pop, True)
# # to get distribution with same sd/cv, need to experiment - no analytic solution
# N_pop = 1000000
# a = generate_random_effects(prop_short, 1.1805, N_pop, True)
# print(np.std(a) / np.mean(a))  # this depends on both the median and the sd
# plt.hist(a, bins=100)  # but values may need to be tighter here, to keep values realistic

# And get distribution for each value of beta_0, as well:
beta_0s = np.zeros([len(beta_0s_median), N_pop])
for i in range(len(beta_0s_median)):
    beta_0s[i] = generate_random_effects(beta_0s_median[i], sd, N_pop)
    print(np.std(beta_0s[i]) / np.mean(beta_0s[i]))  # check sd

# Calculate beta for each participant (random "median" value + seasonal effect):
betas = beta_0s * (1 + beta_1 * np.cos((2 * np.pi / 12) * (vacc_months - phi)))

#######################################################################################################################

# Run!
sim_titers = calculate_Ab_titers_biexp(tm_start, tm_end, vacc_timepoints, maternal_antibodies, betas,
                                       half_life_maternal, half_life_short, half_life_long, prop_short,
                                       N_pop, response_delay)

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

# Note: For visualization purposes, plotted data were an earlier run with only 100 "participants;" 1000 were used to
# generate the first round of synthetic data for fitting

# Write "true" values to file:
if not os.path.isdir('data/'):
    os.mkdir('data/')
if not os.path.isdir('data/prelim_check_' + ymd + '/'):
    os.mkdir('data/prelim_check_' + ymd + '/')

true_vals = pd.DataFrame(sim_titers)
true_vals.to_csv('data/prelim_check_' + ymd + '/truth.csv', na_rep='NA', index=False)
del true_vals

# Add random noise:
sim_titers = add_random_noise(sim_titers, 0.1)

# # Plot noise-laden "data":
# plt.figure()
# plt.plot(sim_titers)
# plt.yscale('log')
# plt.xlabel('Time (Days)')
# plt.ylabel('Ab Titers')
# plt.tight_layout()
# plt.savefig('results/PRELIM_plotTiters_' + ymd + '/ab_titers_over_time_NOISE.png', dpi=300)

# Write synthetic data to file:
obs_vals = pd.DataFrame(sim_titers)
obs_vals.to_csv('data/prelim_check_' + ymd + '/obs_data.csv', na_rep='NA', index=False)

# Write month of vaccine to file:
vacc_months = pd.DataFrame(vacc_months)
vacc_months.to_csv('data/prelim_check_' + ymd + '/vacc_months.csv', na_rep='NA', index=False, header=False)
