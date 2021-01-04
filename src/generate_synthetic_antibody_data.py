
# Run model with random effects by participant #

# Import necessary functions:
import pandas as pd
from functions_python import *
from functions_python_plot import *

# def propagate_Ab_model(tmStrt, tmEnd, tmStep, A0, rho, beta, beta_s, beta_l, r, d_s, d_l, tau):
#     dBs = rho * beta - c_s * B_s
# but if the system can be solved... we can actually just plug all desired t into the one equation...
# Although using the differential equations might make it easier to incorporate effects like age, season of birth, etc.?
# Also could run using full ODEs and plot out not only Ab titers but also levels of ASCs

#######################################################################################################################

# Set global parameters:
N_pop = 1000  # number of "participants"

# Set parameter means/sd:
maternal_antibodies_median = 0.2
betas_median = np.array([0.035, 0.2, 1.0, 10.0])
half_life_maternal_median = 46.0  # days
half_life_igg_median = 20.3  # days
half_life_short_median = 5.5  # days
half_life_long_median = 1648  # days
prop_short_median = 0.98
sd = 0.20

# Set start, end, and vaccination timepoints:
tm_start = 0
tm_end = 730  # 2 years; tm_end = 1825  # 5 years
vacc_timepoints = np.array([60, 120, 180, 365])  # vaccinate at 2, 4, 6 months, + booster at 12

#######################################################################################################################

# Set values for each participant:
maternal_antibodies = generate_random_effects(maternal_antibodies_median, sd, N_pop)
half_life_maternal = generate_random_effects(half_life_maternal_median, sd, N_pop)
half_life_igg = generate_random_effects(half_life_igg_median, sd, N_pop)
half_life_short = generate_random_effects(half_life_short_median, sd, N_pop)
half_life_long = generate_random_effects(half_life_long_median, sd, N_pop)

# To keep prop_short between 0-1, use logit-normal instead:
prop_short = generate_random_effects(prop_short_median, sd, N_pop, True)

# And get distribution for each value of beta, as well:
betas = np.zeros([len(betas_median), N_pop])
for i in range(len(betas_median)):
    betas[i] = generate_random_effects(betas_median[i], sd, N_pop)

#######################################################################################################################

# Run!
sim_titers = calculate_Ab_titers(tm_start, tm_end, vacc_timepoints, maternal_antibodies, betas,
                                 half_life_maternal, half_life_igg, half_life_short, half_life_long,
                                 prop_short, N_pop)

# Plot simulated "data":
# plot_synth_ab_titers(sim_titers, save_path='results/PRELIM_plotTiters_20210104/',
#                      save_filename='ab_titers_over_time.png')
# Note: For visualization purposes, plotted data were an earlier run with only 100 "participants;" 1000 were used to
# generate the first round of synthetic data for fitting

# Write "true" values to file:
if not os.path.isdir('data/'):
    os.mkdir('data/')
if not os.path.isdir('data/prelim_check_20210104/'):
    os.mkdir('data/prelim_check_20210104/')

true_vals = pd.DataFrame(sim_titers)
true_vals.to_csv('data/prelim_check_20210104/truth.csv', na_rep='NA', index=False)
del true_vals

# Add random noise:
sim_titers = add_random_noise(sim_titers)

# Plot noise-laden "data":
# plot_synth_ab_titers(sim_titers, save_path='results/PRELIM_plotTiters_20210104/',
#                      save_filename='ab_titers_over_time_NOISE.png')

# Write synthetic data to file:
obs_vals = pd.DataFrame(sim_titers)
obs_vals.to_csv('data/prelim_check_20210104/obs_data.csv', na_rep='NA', index=False)
