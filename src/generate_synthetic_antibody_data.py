# Run model with random effects by participant #

# Import necessary functions:
import pandas as pd
from functions_python import *
from functions_python_plot import *

#######################################################################################################################

# Set global parameters:
N_pop = 1000  # number of "participants"
response_delay = 14  # 2 week delay in Ab response
prop_short = np.float64(0.70)  # fit this as fixed for now

# Set parameter medians:
maternal_antibodies_median = 8.0
betas_median = np.array([18.0])  # single dose for now
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
# print(np.std(get_rate_from_half_life(half_life_long)))  # want: 0.010?

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
                                       N_pop, response_delay)

# # Plot:
# plt.figure(figsize=(20, 9))
# plt.plot(sim_titers)
# # plt.yscale('log')
# plt.hlines(y=(10, 5), xmin=tm_start, xmax=tm_end)
# plt.vlines(x=(60, 120, 180, 548, 730, 1095, 2191), ymin=0, ymax=30)  # 2/4/6mo; 6m, 1y, 2y, 5y post vaccination
# plt.tight_layout()

# # Check for realism:
# print('\nCheck Realistic:')
# print(np.mean(sim_titers[vacc_timepoints[0]]))  # want: 11-18
# print(len(np.where(sim_titers[vacc_timepoints[0]] >= 5.0)[0]) / N_pop)  # want: 85-95% (or 100?)
# print(len(np.where(sim_titers[vacc_timepoints[0] + 365] >= 5.0)[0]) / N_pop)  # want: ~95%
# print(len(np.where(sim_titers[vacc_timepoints[0] + 365*5] >= 5.0)[0]) / N_pop)  # want: ~88-100%
# print(len(np.where(sim_titers[vacc_timepoints[0] + 365*9] >= 5.0)[0]) / N_pop)  # less than 100%?
# # Difficult to tell what's "realistic," as there seems to be a lot of boosting from natural exposure in the literature

# # Plot simulated "data":
# # plot_synth_ab_titers(sim_titers, save_path='results/PRELIM_plotTiters_20210113/',
# #                      save_filename='ab_titers_alpha' + str(maternal_antibodies_median) +
# #                                    '_m' + str(half_life_maternal_median) +
# #                                    '_beta' + str(betas_median) +
# #                                    '_r1' + str(half_life_short_median) +
# #                                    '_r2' + str(half_life_long_median) +
# #                                    '_rho' + str(prop_short) + '.png')
# plot_synth_ab_titers(sim_titers, save_path='results/PRELIM_plotTiters_20210115/',
#                      save_filename='ab_titers_over_time.png')

# Note: For visualization purposes, plotted data were an earlier run with only 100 "participants;" 1000 were used to
# generate the first round of synthetic data for fitting

# Write "true" values to file:
if not os.path.isdir('data/'):
    os.mkdir('data/')
if not os.path.isdir('data/prelim_check_20210115/'):
    os.mkdir('data/prelim_check_20210115/')

true_vals = pd.DataFrame(sim_titers)
true_vals.to_csv('data/prelim_check_20210115/truth.csv', na_rep='NA', index=False)
del true_vals

# Add random noise:
sim_titers = add_random_noise(sim_titers, 0.1)

# # Plot noise-laden "data":
# plot_synth_ab_titers(sim_titers, save_path='results/PRELIM_plotTiters_20210115/',
#                      save_filename='ab_titers_over_time_NOISE.png')

# Write synthetic data to file:
obs_vals = pd.DataFrame(sim_titers)
obs_vals.to_csv('data/prelim_check_20210115/obs_data.csv', na_rep='NA', index=False)
