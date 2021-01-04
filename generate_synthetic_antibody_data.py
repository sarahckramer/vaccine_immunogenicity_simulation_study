import os
import numpy as np
import matplotlib.pyplot as plt

# Run model with random effects by participant #

# def propagate_Ab_model(tmStrt, tmEnd, tmStep, A0, rho, beta, beta_s, beta_l, r, d_s, d_l, tau):
#     dBs = rho * beta - c_s * B_s
# but if the system can be solved... we can actually just plug all desired t into the one equation...
# Although using the differential equations might make it easier to incorporate effects like age, season of birth, etc.?
# Also could run using full ODEs and plot out not only Ab titers but also levels of ASCs

#######################################################################################################################

# Set up folders to store any outputs:
save_path = 'results/PRELIM_randomEffectsGen/'
if not os.path.isdir(save_path):
    os.mkdir(save_path)

#######################################################################################################################

# Set global parameters:
N = 100  # number of "participants"

# Set parameter means/sd:
maternal_antibodies_median = 0.2
betas_median = np.array([0.035, 0.2, 1.0, 10.0])
half_life_maternal_median = 46.0  # days
half_life_igg_median = 20.3  # days
half_life_short_median = 5.5  # days
half_life_long_median = 1648  # days
prop_short_median = 0.98
sd = 0.20

#######################################################################################################################


# Set values for each participant:
def generate_random_effects(median_val, sd_val, pop_size, logit=False):
    if logit:
        out = np.random.normal(loc=np.log(median_val / (1 - median_val)), scale=sd_val, size=pop_size)
        out = np.exp(out) / (1 + np.exp(out))
        return out
    else:
        log_val = np.log(median_val)
        # out = np.random.lognormal(mean=log_val, sigma=sd_perc_val*np.abs(log_val), size=pop_size)
        out = np.random.lognormal(mean=log_val, sigma=sd_val, size=pop_size)
        # out = np.exp(np.random.normal(loc=log_val, scale=sd_val, size=pop_size))
        return out


maternal_antibodies = generate_random_effects(maternal_antibodies_median, sd, N)
half_life_maternal = generate_random_effects(half_life_maternal_median, sd, N)
half_life_igg = generate_random_effects(half_life_igg_median, sd, N)
half_life_short = generate_random_effects(half_life_short_median, sd, N)
half_life_long = generate_random_effects(half_life_long_median, sd, N)

# To keep prop_short between 0-1, use logit-normal instead:
prop_short = generate_random_effects(prop_short_median, sd, N, True)

# And get distribution for each value of beta, as well:
betas = np.zeros([len(betas_median), N])
for i in range(len(betas_median)):
    betas[i] = generate_random_effects(betas_median[i], sd, N)
# print(np.where(betas == 0))

# Plot histograms to check distributions:
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=[12, 9.5])
axs[0, 0].hist(maternal_antibodies, bins=int(N / 5))
axs[0, 0].title.set_text('Maternal Ab Levels')
axs[0, 1].hist(prop_short, bins=int(N / 5))
axs[0, 1].title.set_text('Prop. Short-Lived ASCs')
axs[1, 0].hist(half_life_maternal, bins=int(N / 5))
axs[1, 0].title.set_text('Half-Life (Maternal Ab)')
axs[1, 1].hist(half_life_long, bins=int(N / 5))
axs[1, 1].title.set_text('Half-Life (Long-Lived ASCs)')
axs[2, 0].hist(half_life_short, bins=int(N / 5))
axs[2, 0].title.set_text('Half-Life (Short-Lived ASCs)')
axs[2, 1].hist(half_life_igg, bins=int(N / 5))
axs[2, 1].title.set_text('Half-Life (IGG)')
for i in range(3):
    for j in range(2):
        axs[i, j].set_xlabel('Value')
        axs[i, j].set_ylabel('Count')
fig.tight_layout()
plt.savefig(save_path + 'hist_random_effects.png', dpi=300)

# Check for expected medians:
f = open(save_path + 'check_medians.txt', 'w')

print('Diff. between observed and expected medians:', file=f)
print(np.median(maternal_antibodies) - maternal_antibodies_median, file=f)
print(np.median(half_life_maternal) - half_life_maternal_median, file=f)
print(np.median(half_life_igg) - half_life_igg_median, file=f)
print(np.median(half_life_short) - half_life_short_median, file=f)
print(np.median(half_life_long) - half_life_long_median, file=f)
print(np.median(prop_short) - prop_short_median, file=f)
for i in range(len(betas_median)):
    print(np.median(betas[i]) - betas_median[i], file=f)
print(file=f)

print('Diff. between observed and expected medians, relative to expected:', file=f)
print((np.median(maternal_antibodies) - maternal_antibodies_median) / maternal_antibodies_median, file=f)
print((np.median(half_life_maternal) - half_life_maternal_median) / half_life_maternal_median, file=f)
print((np.median(half_life_igg) - half_life_igg_median) / half_life_igg_median, file=f)
print((np.median(half_life_short) - half_life_short_median) / half_life_short_median, file=f)
print((np.median(half_life_long) - half_life_long_median) / half_life_long_median, file=f)
print((np.median(prop_short) - prop_short_median) / prop_short_median, file=f)
for i in range(len(betas_median)):
    print((np.median(betas[i]) - betas_median[i]) / betas_median[i], file=f)
print(file=f)

f.close()

# Check standard deviations? Value of "sd" is for "underlying normal distribution," not actual distribution:
f = open(save_path + 'check_sds.txt', 'w')

print('Observed standard deviations:', file=f)
print(np.std(maternal_antibodies), file=f)
print(np.std(half_life_maternal), file=f)
print(np.std(half_life_igg), file=f)
print(np.std(half_life_short), file=f)
print(np.std(half_life_long), file=f)
print(np.std(prop_short), file=f)
for i in range(len(betas_median)):
    print(np.std(betas[i]), file=f)
print(file=f)

print('Observed standard deviations, as proportion of expected median:', file=f)
print(np.std(maternal_antibodies) / maternal_antibodies_median, file=f)
print(np.std(half_life_maternal) / half_life_maternal_median, file=f)
print(np.std(half_life_igg) / half_life_igg_median, file=f)
print(np.std(half_life_short) / half_life_short_median, file=f)
print(np.std(half_life_long) / half_life_long_median, file=f)
print(np.std(prop_short) / prop_short_median, file=f)  # MAKE 20% somehow??
for i in range(len(betas_median)):
    print(np.std(betas[i]) / betas_median[i], file=f)
print(file=f)

f.close()

#######################################################################################################################

# Set start, end, and vaccination timepoints:
tm_start = 0
tm_end = 730  # 2 years; tm_end = 1825  # 5 years
vacc_timepoints = np.array([60, 120, 180, 365])  # vaccinate at 2, 4, 6 months, + booster at 12


# Function to calculate rates from half-lives:
def get_rate_from_half_life(d):
    return np.log(2) / d


# Function to calculate antibody titers at each time point for each participant:
def calculate_Ab_titers(t_start, t_end, tau, A_m, beta, d_m, d_a, d_s, d_l, rho):
    # Start by calculating rates of decay from half-lives:
    m = get_rate_from_half_life(d_m)
    r = get_rate_from_half_life(d_a)
    c_s = get_rate_from_half_life(d_s)
    c_l = get_rate_from_half_life(d_l)

    # Store Ab titers at each time point:
    Ab_titers = np.zeros([len(range(t_start, t_end + 1)), N])

    # Loop through each time and calculate:
    t = t_start
    while t <= t_end:
        maternal_ab = A_m * np.exp(-m * t)
        A_t = maternal_ab

        # add antibody response to each dose that has been given so far
        if t >= min(tau):
            for i in np.where(tau <= t)[0]:
                tau_i = tau[i]

                vaccine_ab = beta[i] * ((rho / (r - c_s)) *
                                        (np.exp(-c_s * (t - tau_i)) - np.exp(-r * (t - tau_i))) +
                                        ((1 - rho) / (r - c_l)) *
                                        (np.exp(-c_l * (t - tau_i)) - np.exp(-r * (t - tau_i))))
                A_t = A_t + vaccine_ab

        # append titers at this timepoint to list
        # Ab_titers.append(A_t)
        Ab_titers[t] = A_t

        # move forward one day
        t += 1

    # Return titers over time:
    return Ab_titers


# Run!
sim_titers = calculate_Ab_titers(tm_start, tm_end, vacc_timepoints, maternal_antibodies, betas,
                                 half_life_maternal, half_life_igg, half_life_short, half_life_long,
                                 prop_short)

# Plot simulated "data":
plt.figure()
plt.plot(sim_titers)
plt.yscale('log')
plt.xlabel('Time (Days)')
plt.ylabel('Ab Titers')
plt.tight_layout()
plt.savefig(save_path + 'ab_titers_over_time.png', dpi=300)
