# Run model with random effects by participant #

# Import necessary functions:
import pandas as pd
from functions_python import *
from functions_python_plot import *

#######################################################################################################################

# Set global parameters:
N_pop = 1000  # number of "participants"
prop_short = np.float64(0.95)  # fit this as fixed for now

# Set parameter medians:
maternal_antibodies_median = 0.2
betas_median = np.array([10.0])  # single dose for now
# betas_median = np.array([0.035, 0.2, 1.0, 10.0])
half_life_maternal_median = 46.0  # days
half_life_short_median = 30.0  # days
half_life_long_median = 365.0  # days

# Get standard deviation from desired coefficient of variations (var/mean of distribution):
sd = calculate_sigma_from_cv(0.1)

# Set start, end, and vaccination timepoints:
tm_start = 0
tm_end = 730  # 2 years; tm_end = 1825  # 5 years
# vacc_timepoints = np.array([60, 120, 180, 365])  # vaccinate at 2, 4, 6 months, + booster at 12
vacc_timepoints = np.array([365])  # vaccinate between 12 and 15 months; can be second dose at 4-6 yrs

#######################################################################################################################

# Set values for each participant:
maternal_antibodies = generate_random_effects(maternal_antibodies_median, sd, N_pop)
half_life_maternal = generate_random_effects(half_life_maternal_median, sd, N_pop)
half_life_short = generate_random_effects(half_life_short_median, sd, N_pop)
half_life_long = generate_random_effects(half_life_long_median, sd, N_pop)

# Check sd is correct:
print(np.std(maternal_antibodies) / np.mean(maternal_antibodies))
print(np.std(half_life_maternal) / np.mean(half_life_maternal))
print(np.std(half_life_short) / np.mean(half_life_short))
print(np.std(half_life_long) / np.mean(half_life_long))

# # To keep prop_short between 0-1, use logit-normal instead:
# prop_short = generate_random_effects(prop_short_median, sd, N_pop, True)
# # to get distribution with same sd/cv, need to experiment - no analytic solution
# N_pop = 1000000
# a = generate_random_effects(prop_short, 1.1805, N_pop, True)
# print(np.std(a) / np.mean(a))  # this depends on both the median and the sd
# plt.hist(a, bins=100)  # but values may need to be tighter here, to keep values realistic

# And get distribution for each value of beta, as well:
betas = np.zeros([len(betas_median), N_pop])
for i in range(len(betas_median)):
    betas[i] = generate_random_effects(betas_median[i], sd, N_pop)
    print(np.std(betas[i]) / np.mean(betas[i]))  # check sd
#######################################################################################################################

# Run!
sim_titers = calculate_Ab_titers_biexp(tm_start, tm_end, vacc_timepoints, maternal_antibodies, betas,
                                       half_life_maternal, half_life_short, half_life_long, prop_short,
                                       N_pop)

# # Plot simulated "data":
# plot_synth_ab_titers(sim_titers, save_path='results/PRELIM_plotTiters_20210112/',
#                      save_filename='ab_titers_over_time.png')
# # Note: For visualization purposes, plotted data were an earlier run with only 100 "participants;" 1000 were used to
# # generate the first round of synthetic data for fitting

# Write "true" values to file:
if not os.path.isdir('data/'):
    os.mkdir('data/')
if not os.path.isdir('data/prelim_check_20210112/'):
    os.mkdir('data/prelim_check_20210112/')

true_vals = pd.DataFrame(sim_titers)
true_vals.to_csv('data/prelim_check_20210112/truth.csv', na_rep='NA', index=False)
del true_vals

# Add random noise:
sim_titers = add_random_noise(sim_titers, 0.1)

# # Plot noise-laden "data":
# plot_synth_ab_titers(sim_titers, save_path='results/PRELIM_plotTiters_20210112/',
#                      save_filename='ab_titers_over_time_NOISE.png')

# Write synthetic data to file:
obs_vals = pd.DataFrame(sim_titers)
obs_vals.to_csv('data/prelim_check_20210112/obs_data.csv', na_rep='NA', index=False)
