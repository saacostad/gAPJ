import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

"""
This script will calculate the beta-beating of the lattice before and
after applying the corrections.
"""

nominal = "outputs/nominal/nominal_measurements_twiss.parquet"
apj = "_measurements_twiss.parquet"

before_path = "outputs/errors/E"
after_path = "outputs/corrections/C"

errors_path_nominal = before_path + nominal
errors_path_apj = before_path + apj 

corrections_path_nominal = after_path + nominal
corrections_path_apj = after_path + apj 


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

    return beta_beating_x, beta_beating_y, nom_s


bx_bef, by_bef, sb = calc_beta_beating(nominal, errors_path_apj)
bx_aft, by_aft, sa = calc_beta_beating(nominal, corrections_path_apj)



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

ax_top.set_title(r"$\beta$-Beating en eje X")
ax_top.plot(sb, bx_bef, label = "Before corrections")
ax_top.plot(sa, bx_aft, label = "After corrections")
ax_top.grid()
ax_top.legend()
ax_top.set_xlabel("s [m]")
ax_top.set_ylabel(r"$\beta$-beating")


ax_bottom.set_title(r"$\beta$-Beating en eje Y")
ax_bottom.plot(sb, by_bef, label = "Before corrections")
ax_bottom.plot(sa, by_aft, label = "After corrections")
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



