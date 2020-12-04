import numpy as np
import matplotlib.pyplot as plt

# First goal is simply to simulate with no random effects or effect of season/schedule - same params for everyone


# def propagate_Ab_model(tmStrt, tmEnd, tmStep, A0, rho, beta, beta_s, beta_l, r, d_s, d_l, tau):
#     dBs = rho * beta - c_s * B_s
# but if the system can be solved... we can actually just plug all desired t into the one equation...
# although using the differential equations might make it easier to incorporate effects like age, season of birth, etc.?


# Set global parameters:
N = 100  # number of "participants"

# Set parameter means/sd:
# For now, single value:
maternal_antibodies = 0.2
betas = np.array([0.035, 0.2, 1.0, 10.0])
half_life_maternal = 46.0  # days
half_life_igg = 20.3
half_life_short = 5.5
half_life_long = 1648
prop_short = 0.98

# # betas = np.array([0.25, 1.0, 4.0, 12.9])
# half_life_maternal = 46.0  # days
# half_life_igg = 30.9
# half_life_short = 9.9
# half_life_long = 2287
# prop_short = 0.957

tm_start = 1
# tm_end = 1825  # 5 years
tm_end = 730  # 2 years
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
    Ab_titers = []

    # Loop through each time and calculate:
    t = t_start
    while t <= t_end:
        maternal_ab = A_m * np.exp(-m * t)
        A_t = maternal_ab

        # add antibody response to each dose that has been given so far
        if t >= min(tau):
            for i in np.where(tau <= t)[0]:
                # i = i
                # print(i)

                tau_i = tau[i]
                # print(tau_i)

                vaccine_ab = beta[i] * ((rho / (r - c_s)) *
                                        (np.exp(-c_s * (t - tau_i)) - np.exp(-r * (t - tau_i))) +
                                        ((1 - rho) / (r - c_l)) *
                                        (np.exp(-c_l * (t - tau_i)) - np.exp(-r * (t - tau_i))))
                A_t = A_t + vaccine_ab

        # append titers at this timepoint to list
        Ab_titers.append(A_t)

        # move forward one day
        t += 1

    # Return titers over time:
    return Ab_titers


# Run!
# values = [0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99]
# for prop_short in values:
sim_titers = calculate_Ab_titers(tm_start, tm_end, vacc_timepoints, maternal_antibodies, betas,
                                 half_life_maternal, half_life_igg, half_life_short, half_life_long,
                                 prop_short)
plt.plot(sim_titers)
plt.yscale('log')
# plt.legend(labels=values)
# plt.savefig(fname='results/by_half_life_propShort.png')
