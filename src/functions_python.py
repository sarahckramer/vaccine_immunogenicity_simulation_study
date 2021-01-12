# All (python) functions needed to generate synthetic data #
import numpy as np


def generate_random_effects(median_val, sd_val, pop_size, logit=False):
    """
    Generates individual-level parameter values as lognormal (or logit-normal) distributions around some median value

    :param median_val: Median value of parameter of interest
    :param sd_val: Standard deviation to be used by underlying normal distribution
    :param pop_size: Size of the synthetic population
    :param logit: Whether or not to generate logit-normal- rather than lognormal-distributed values (default False)
    :return: Individual-level parameter values
    """

    if logit:
        out = np.random.normal(loc=np.log(median_val / (1 - median_val)), scale=sd_val, size=pop_size)
        out = np.exp(out) / (1 + np.exp(out))
        return out
    else:
        log_val = np.log(median_val)
        out = np.random.lognormal(mean=log_val, sigma=sd_val, size=pop_size)
        return out


def get_rate_from_half_life(d):
    """
    Calculates rates of decay from half-lives

    :param d: Half-life
    :return: Corresponding decay rate
    """
    return np.log(2) / d


def calculate_sigma_from_cv(desired_cv):
    """
    Calculates standard deviation of a lognormal distribution corresponding to the input coefficient of variation

    :param desired_cv: The intended coefficient of variation
    :return: Standard deviation for the underlying normal distribution of a lognormal distribution
    """
    sigma = np.sqrt(np.log(desired_cv ** 2 + 1))
    return sigma


def calculate_Ab_titers_biexp(t_start, t_end, tau, alpha, beta, d_m, d_s, d_l, rho, N):
    """
    Calculates synthetic antibody titers at each time point for each participant according to a non-mechanistic
    model accounting for biexponential antibody decay rates

    :param t_start: First timepoint at which titers are calculated
    :param t_end: Final timepoint at which titers are calculated (inclusive)
    :param tau: Timepoints at which vaccine is given (1-D array)
    :param alpha: Maternal antibody titers for each participant (1-D array)
    :param beta: Individual-level values for beta at each vaccination timepoint (1 or 2-D array)
    :param d_m: Half-life of maternal antibodies for each participant (1-D array)
    :param d_s: Half-life of IgG antibodies for each participant (short-lived) (1-D array)
    :param d_l: Half-life of IgG antibodies for each participant ("long-lived") (1-D array)
    :param rho: Proportion of antibodies decaying at rate d_s (1-D array OR float64)
    :param N: Size of the synthetic population
    :return: Antibody titers for each time point (rows) and participant (columns)
    """
    # Start by calculating rates of decay from half-lives:
    m = get_rate_from_half_life(d_m)
    r1 = get_rate_from_half_life(d_s)
    r2 = get_rate_from_half_life(d_l)

    # Store Ab titers at each time point:
    Ab_titers = np.zeros([len(range(t_start, t_end + 1)), N])

    # Loop through each timepoint and calculate:
    t = t_start
    while t <= t_end:
        maternal_ab = alpha * np.exp(-m * t)
        A_t = maternal_ab

        # add antibody response to each dose that has been given so far
        if t >= min(tau):
            for i in np.where(tau <= t)[0]:
                tau_i = tau[i]

                vaccine_ab = beta[i] * (rho * (np.exp(-r1 * (t - tau_i))) +
                                        (1 - rho) * (np.exp(-r2 * (t - tau_i))))
                A_t = A_t + vaccine_ab

        # append titers at this timepoint to list
        Ab_titers[t] = A_t

        # move forward one day
        t += 1

    # Return titers over time:
    return Ab_titers


def calculate_Ab_titers_WhiteM3(t_start, t_end, tau, A_m, beta, d_m, d_a, d_s, d_l, rho, N):
    """
    Calculates synthetic antibody titers at each time point for each participant according to model #3
    in White et al. (2014)

    :param t_start: First timepoint at which titers are calculated
    :param t_end: Final timepoint at which titers are calculated (inclusive)
    :param tau: Timepoints at which vaccine is given (1-D array)
    :param A_m: Maternal antibody titers for each participant (1-D array)
    :param beta: Individual-level values for beta at each vaccination timepoint (1- or 2-D array)
    :param d_m: Half-life of maternal antibodies for each participant (1-D array)
    :param d_a: Half-life of IgG antibodies for each participant (1-D array)
    :param d_s: Half-life of short-lived antibody secreting cells for each participant (1-D array)
    :param d_l: Half-life of long-lived antibody secreting cells for each participant (1-D array)
    :param rho: Proportion of antibody secreting cells that are short-lived for each participant (1-D array)
    :param N: Size of the synthetic population
    :return: Antibody titers for each time point (rows) and participant (columns)
    """

    # Start by calculating rates of decay from half-lives:
    m = get_rate_from_half_life(d_m)
    r = get_rate_from_half_life(d_a)
    c_s = get_rate_from_half_life(d_s)
    c_l = get_rate_from_half_life(d_l)

    # Store Ab titers at each time point:
    Ab_titers = np.zeros([len(range(t_start, t_end + 1)), N])

    # Loop through each timepoint and calculate:
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


def add_random_noise(synth_titers, cv_natural):
    """
    Adds normally-distributed noise to "true" antibody titer values

    :param synth_titers: "True" (synthetic) antibody titers for all timepoints (rows) and individuals (columns)
    :param cv_natural: Desired coefficient of variation on the natural scale
    :return: Noise-laden synthetic antibody titers for each timepoint (rows) and individual (columns)
    """

    # Calculate desired sd:
    sd = calculate_sigma_from_cv(cv_natural)

    # Draw from normal distribution around log of "true" values:
    # noisy_titers_orig = np.random.normal(loc=synth_titers, scale=0.1*synth_titers)
    noisy_titers = np.random.lognormal(mean=np.log(synth_titers), sigma=sd)
    # print(np.std(np.reshape((noisy_titers - synth_titers) / synth_titers, [731000, ])))
    # print(np.std(np.reshape((noisy_titers_orig - synth_titers) / synth_titers, [731000, ])))

    # # Ensure that no values are <0:
    # neg_indices = np.where(noisy_titers < 0)
    # while len(neg_indices[0] > 0):
    #     # noisy_titers[neg_indices] = np.random.normal(loc=synth_titers[neg_indices], scale=0.1*synth_titers[neg_indices])
    #     noisy_titers[neg_indices] = np.random.lognormal(mean=np.log(synth_titers[neg_indices]), sigma=sd)
    #     neg_indices = np.where(noisy_titers < 0)

    # Return noise-laden values:
    return noisy_titers
