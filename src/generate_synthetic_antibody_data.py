# Run model with random effects by participant #

# Import necessary functions:
import os
import pandas as pd
from functions_python import *

#######################################################################################################################

# Set date:
ymd = '20210210'

# Set global parameters:
N_pop = 1000  # number of "participants"
response_delay = 14  # 2 week delay in Ab response
maternal_antibodies_median = 8.0
half_life_maternal_median = 42.0  # days

# Set start, end, and vaccination timepoints:
tm_start = 0
tm_end = 3652  # 730  # 10 years
vacc_timepoints = np.array([365])  # vaccinate between 12 and 15 months; can be second dose at 4-6 yrs

#######################################################################################################################

# Read in input parameter values:
input_params = pd.read_csv('data/inputs/lhs_params.txt', header=None, sep='\t')
input_params = input_params.to_numpy(dtype=np.float64)

#######################################################################################################################

# Select birth month (1-12) for each participant:
rng = np.random.default_rng()
birth_months = rng.integers(low=1, high=13, size=N_pop)

# Convert to month of first vaccination based on vacc_timepoints:
vacc_months = birth_months + np.floor(vacc_timepoints[0] / 30)

# Write month of vaccine to file:
if not os.path.isdir('data/'):
    os.mkdir('data/')
if not os.path.isdir('data/prelim_check_' + ymd + '/'):
    os.mkdir('data/prelim_check_' + ymd + '/')

vacc_months_pd = pd.DataFrame(vacc_months)
vacc_months_pd.to_csv('data/prelim_check_' + ymd + '/vacc_months.csv', na_rep='NA', index=False, header=False)
del vacc_months_pd

#######################################################################################################################

# LOOP THROUGH PARAMETER SETS
for param_set in range(input_params.shape[0]):
    print(param_set)

    # Get fixed effects for each parameter:
    beta_0s_median = np.array([input_params[param_set, 0]])
    beta_1 = input_params[param_set, 1]
    phi = input_params[param_set, 2]
    prop_short = input_params[param_set, 3]
    half_life_short_median = input_params[param_set, 4]
    half_life_long_median = input_params[param_set, 5]

    ###################################################################################################################

    # LOOP THROUGH COEFFICIENTS OF VARIATION
    for desired_cv in [0.1, 0.2]:

        # Get standard deviation from desired coefficient of variations (var/mean of distribution):
        sd = calculate_sigma_from_cv(desired_cv)

        # Set values of random effects for each participant:
        maternal_antibodies = generate_random_effects(maternal_antibodies_median, sd, N_pop)
        half_life_maternal = generate_random_effects(half_life_maternal_median, sd, N_pop)
        half_life_short = generate_random_effects(half_life_short_median, sd, N_pop)
        half_life_long = generate_random_effects(half_life_long_median, sd, N_pop)

        # Check sd is correct:
        print('True: ' + str(desired_cv))
        print(np.std(maternal_antibodies) / np.mean(maternal_antibodies))
        print(np.std(half_life_maternal) / np.mean(half_life_maternal))
        print(np.std(half_life_short) / np.mean(half_life_short))
        print(np.std(half_life_long) / np.mean(half_life_long))

        # And get distribution for each value of beta_0, as well:
        beta_0s = np.zeros([len(beta_0s_median), N_pop])
        for i in range(len(beta_0s_median)):
            beta_0s[i] = generate_random_effects(beta_0s_median[i], sd, N_pop)
            print(np.std(beta_0s[i]) / np.mean(beta_0s[i]))  # check sd
        print()

        # Calculate beta for each participant (random "median" value + seasonal effect):
        betas = beta_0s * (1 + beta_1 * np.cos((2 * np.pi / 12) * (vacc_months - phi)))

        ###############################################################################################################

        # Run!
        sim_titers = calculate_Ab_titers_biexp(tm_start, tm_end, vacc_timepoints, maternal_antibodies, betas,
                                               half_life_maternal, half_life_short, half_life_long, prop_short,
                                               N_pop, response_delay)

        # Write "true" values to file:
        true_vals = pd.DataFrame(sim_titers)
        true_vals.to_csv('data/prelim_check_' + ymd + '/truth_' + str(param_set) + '_cv' + str(desired_cv) + '.csv',
                         na_rep='NA', index=False)
        del true_vals

        ###############################################################################################################

        # LOOP THROUGH NOISE AMOUNTS
        for noise in [0.1, 0.2]:
            print(noise)

            # Add random noise:
            noisy_titers = add_random_noise(sim_titers, noise)

            # Write synthetic data to file:
            obs_vals = pd.DataFrame(noisy_titers)
            obs_vals.to_csv('data/prelim_check_' + ymd + '/obs_data_' + str(param_set) + '_cv' + str(desired_cv) +
                            '_noise' + str(noise) + '.csv',
                            na_rep='NA', index=False)
