import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import argparse                 # Will be used mainly to select if we're creating nominal, errors or errors+corrections systems
from scipy.signal import find_peaks


"""
This script will calculate the beta-beating of the lattice before and
after applying the corrections.
"""

# We'll parse the IP 
parser = argparse.ArgumentParser()
# -- Choose which system to simulate
parser.add_argument(
        "-ip", "--interaction_point",
        help="Which interaction point to plot",
        dest="ip",
    )
parsed_args = parser.parse_args()
ip = parsed_args.ip

nominal = "outputs/nominal/nominal_measurements_twiss.parquet"
apj = "_measurements_twiss.parquet"

before_path = "outputs/errors/E"
after_path = "outputs/corrections/C"

errors_path_nominal = before_path + nominal
errors_path_apj = before_path + apj 

corrections_path_nominal = after_path + nominal
corrections_path_apj = after_path + apj 


# For FCC
left_arc = [15000, 20000]
right_arc = [23500, 30000]


def calc_beta_beating(nom_path, apj_path):
    """ This function takes the nominal and experimental twiss files and calculates the beta beating """

    nom = pd.read_parquet(nom_path)
    
    nom_s = nom["S"]
    nom_beta_x = nom["BETX"] 
    nom_beta_y = nom["BETY"]

    apj = pd.read_parquet(apj_path)

    apj_s = apj["S"]
    apj_beta_x = apj["BETX"] 
    apj_beta_y = apj["BETY"]


    beta_beating_x = (apj_beta_x - nom_beta_x) / nom_beta_x
    beta_beating_y = (apj_beta_y - nom_beta_y) / nom_beta_y


    left_arc_mask = (nom_s < left_arc[1]) & (nom_s > left_arc[0])
    right_arc_mask = (nom_s < right_arc[1]) & (nom_s > right_arc[0])
    arcs_mask = left_arc_mask | right_arc_mask


    # Find local maxima
    x_arcs = beta_beating_x[arcs_mask]
    y_arcs = beta_beating_y[arcs_mask]
    peaks_x, _ = find_peaks(x_arcs.to_numpy())
    peaks_y, _ = find_peaks(y_arcs.to_numpy())

    # Recover original Pandas indices
    peaks_x = x_arcs.index[peaks_x]
    peaks_y = y_arcs.index[peaks_y]

    # Calculate mean of maxima
    meanbb_x_arcs = x_arcs.loc[peaks_x].mean()
    meanbb_y_arcs = y_arcs.loc[peaks_y].mean()

    s_ip = (nom["S"].loc[nom["NAME"] == ip])

    return beta_beating_x, beta_beating_y, nom_s, s_ip, (meanbb_x_arcs, meanbb_y_arcs)


bx_bef, by_bef, sb, s_ip, mean_err = calc_beta_beating(nominal, errors_path_apj)
bx_aft, by_aft, sa, s_ip, mean_corr = calc_beta_beating(nominal, corrections_path_apj)

print()
print(mean_err[0], mean_err[1], mean_corr[0], mean_corr[1])


def plot(s_ip):
    # Publication style
    sns.set_theme(
        style="ticks",
        context="paper",
    )

    plt.rcParams.update({
        "font.size": 14,
        "axes.labelsize": 14,
        "axes.titlesize": 14,
        "legend.fontsize": 14,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "lines.linewidth": 1.8,
        "axes.linewidth": 0.8,
    })


    # Create figure with 2 rows, 1 column, sharing the x-axis
    fig, (ax_top, ax_bottom) = plt.subplots(
        2, 1,
        figsize=(7.0, 8.0),          # same width, reduced height for 2 subplots
        sharex=True,
        constrained_layout=True,
    )
    s_ip = s_ip.values[0]
    ax_top.set_title(r"$\beta$-Beating en eje X")
    ax_top.plot(sb, bx_bef, label = "Before corrections")
    ax_top.plot(sa, bx_aft, label = "After corrections")
    ax_top.axvline(s_ip, color = "black", lw = 0.5)
    ax_top.axhline(0.0, color = "black", lw = 0.5)
    ax_top.grid()
    ax_top.legend()
    ax_top.set_xlabel("s [m]")
    ax_top.set_ylabel(r"$\beta$-beating")


    ax_bottom.set_title(r"$\beta$-Beating en eje Y")
    ax_bottom.plot(sb, by_bef, label = "Before corrections")
    ax_bottom.plot(sa, by_aft, label = "After corrections")
    ax_bottom.axvline(s_ip, color = "black", lw = 0.5)
    ax_bottom.axhline(0.0, color = "black", lw = 0.5)
    ax_bottom.grid()
    ax_bottom.legend()
    ax_bottom.set_xlabel("s [m]")
    ax_bottom.set_ylabel(r"$\beta$-beating")

    # combined_bef = np.column_stack((sb, bx_bef, by_bef))
    # np.savetxt('beta_beating_before.csv', combined_bef, delimiter=',', header='s,bbx,bby', comments='')
    #
    # combined_aft = np.column_stack((sa, bx_aft, by_aft))
    # np.savetxt('beta_beating_after.csv', combined_aft, delimiter=',', header='s,bbx,bby', comments='')

    # ------------------------
    # Polish and export
    # ------------------------
    sns.despine()   # clean axis spines

    # Create a single legend from the top axis (both subplots share the same labels)
    handles, labels = ax_top.get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", ncol=1, frameon=False)

    # Save the figure
    plt.show()
    plt.close()


# plot(s_ip)
