# Python functions for plotting input and output values #
import os
import matplotlib.pyplot as plt


def plot_hist_parameter_vals(mat_ab, prop_s, hl_mat, hl_long, hl_short, hl_igg, N, save_path, save_filename):
    """
    Plots histograms of input parameters to check/explore distributions

    :param mat_ab: Maternal antibody titers for each participant (1-D array)
    :param prop_s: Proportion of antibody secreting cells that are short-lived for each participant (1-D array)
    :param hl_mat: Half-life of maternal antibodies for each participant (1-D array)
    :param hl_long: Half-life of long-lived antibody secreting cells for each participant (1-D array)
    :param hl_short: Half-life of short-lived antibody secreting cells for each participant (1-D array)
    :param hl_igg: Half-life of IgG antibodies for each participant (1-D array)
    :param N: Size of the synthetic population
    :param save_path: Directory in which to save plot
    :param save_filename: Name of file to save plot as
    """

    # First ensure that save_path exists:
    if not os.path.isdir(save_path):
        os.mkdir(save_path)

    # Plots:
    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=[12, 9.5])
    axs[0, 0].hist(mat_ab, bins=int(N / 5))
    axs[0, 0].title.set_text('Maternal Ab Levels')
    axs[0, 1].hist(prop_s, bins=int(N / 5))
    axs[0, 1].title.set_text('Prop. Short-Lived ASCs')
    axs[1, 0].hist(hl_mat, bins=int(N / 5))
    axs[1, 0].title.set_text('Half-Life (Maternal Ab)')
    axs[1, 1].hist(hl_long, bins=int(N / 5))
    axs[1, 1].title.set_text('Half-Life (Long-Lived ASCs)')
    axs[2, 0].hist(hl_short, bins=int(N / 5))
    axs[2, 0].title.set_text('Half-Life (Short-Lived ASCs)')
    axs[2, 1].hist(hl_igg, bins=int(N / 5))
    axs[2, 1].title.set_text('Half-Life (IGG)')
    for i in range(3):
        for j in range(2):
            axs[i, j].set_xlabel('Value')
            axs[i, j].set_ylabel('Count')
    fig.tight_layout()
    plt.savefig(save_path + save_filename, dpi=300)


# Plot simulated "data":
def plot_synth_ab_titers(synth_titers, save_path, save_filename):
    """
    Plots simulated data by participant

    :param synth_titers: Generated antibody titers over time (rows) for each participant (columns)
    :param save_path: Directory in which to save plot
    :param save_filename: Name of file to save plot as
    """

    # First ensure that save_path exists:
    if not os.path.isdir(save_path):
        os.mkdir(save_path)

    # Plot:
    plt.figure()
    plt.plot(synth_titers)
    plt.yscale('log')
    plt.xlabel('Time (Days)')
    plt.ylabel('Ab Titers')
    plt.tight_layout()
    plt.savefig(save_path + save_filename, dpi=300)
