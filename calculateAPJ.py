"""
This script calculates the APJ variables from the trackone file
TODO: I'm planning to add a Xtrack output file parser too
"""

import tfs as tfs    # To read the twiss files
from turn_by_turn.trackone import read_tbt  # To read the trackone
import pandas as pd     # To get rid of the tfs that are ugly

# These two always go
import numpy as np  
import matplotlib.pyplot as plt

# self-made modules
from modules.APJ.ActionPhaseJump import calculate_APJ   
from modules.APJ.orbit_tools import avermax_2017, avermax_2022
from modules.data_tools.simulateSystem_parser import create_parser_args, parse  # To parse the code's parameters

# Parceros
import tomllib                  # This is to parse the config file
import argparse

print("""
\tA C T I O N   A N D   P H A S E   J U M P
\tVersion alpha
""" )



# ---------------------------------------------
#           FUNCTION TO SAVE DATA
# ---------------------------------------------
def save_APJ_var(plane, s, NAMES, var, name):
    """ Creates the APJ_var.sdds file in the same format as the OG APJ """

    global save_path 

    with open(f"{save_path}{name}", 'w') as f:
        
        for row in range(s.shape[0]):

            f.write(f'{plane} "{NAMES[row]}" {s[row]} {var[row]}\n')

    print(f"----> saved {save_path}{name}")


# --------------------------------
#       PARSING CONFIG 
# --------------------------------
# -- Parse command line

print("Parsing command line arguments...")

parser = argparse.ArgumentParser()
create_parser_args(parser)
parsed_args = parser.parse_args()

system = parse(parsed_args)

# Where to find the data to work with
dic_key = 'Nominal' if system == 'N' else 'Errors' if system == 'E' else 'Corrections' if system == 'C' else 'Invalid'

system_config = None        # Configuration dict for special case use
nominal_config = None        # Configuration dict for special case use

# Read the config file and save it 
with open("configuration.toml", "rb") as f:

    # We only need the nominal and system's data
    general_config = tomllib.load(f)
    
    nominal_config = general_config["NominalSystem"]
    system_config = general_config[f"{dic_key}System"]

print(f"Read from NominalSystem and {dic_key}System entry")

reference_bpm = system_config["ref_bpm"]      # The s_e from where we'll select the avermax orbits

# Paths
trackone_path = system_config["main_output_path"]+ "/" + system_config["track_path"] + ".parquet"     # The original trackone path
twiss_path = nominal_config["main_output_path"] + "/" +  nominal_config["measure_path"] + ".parquet"  # The twiss files to use    
save_path = system_config["main_output_path"] + "/" + system_config["APJ_path"]                       # Where to save the APJ files

# The arcs of IP2 so we calculate the avermax traj
arc =system_config["left_arc"] 
rightArc = system_config["right_arc"]

plot = system_config["plot_APJ"]

_treshold = system_config["avermax_TH"]      # Treshold to use for avermax calc

AVM_alg = system_config["avermax_algorythm"]

# --------------------------------
#       LECTURA DEL TRACKONE
# --------------------------------

# Messages to check the state of the script
print(f"Reading turn-by-turn data from {trackone_path}...")

"""
Something like this if we were reading original trackone files
# Se lee el TrackOne
# to_data = read_tbt(trackone_path)         
# Extraemos los datos
particle = to_data.matrices[0]

Xs = pd.DataFrame(particle.X)
Ys = pd.DataFrame(particle.Y)
"""

# Read the .parquet data
data = pd.read_parquet(trackone_path, engine="pyarrow")

# Create 2 different dataframes for each axis
Xs = data[data["PLANE"] == "x"].drop(columns="PLANE")
Ys = data[data["PLANE"] == "y"].drop(columns="PLANE")

# We'll re-index
Xs = Xs.set_index("NAME")
Ys = Ys.set_index("NAME")

# TODO: here I would have to add the noise

# -----------------------------------
#       HACEMOS LA TRAYECTORIA DIFF 
# -----------------------------------

print("Fitering betatron motion...")

# Calculamos el promedio en cada BPM
Xs_mean = Xs.mean(axis = 1).to_numpy()
Ys_mean = Ys.mean(axis = 1).to_numpy()
 

# TODO: without errors it works amazingly good if we add the mean
# HACK: well there are still jumps but not that big. At least this will work to generalize the APJ
Xs_new = Xs - Xs_mean[:, np.newaxis]
Ys_new = Ys - Ys_mean[:, np.newaxis]


Xs = Xs_new
Ys = Ys_new

# Obtenemos los elementos en las optics de momento
optics_elements_names = Xs.index

# --------------------------------
#       LECTURA DEL TWISS
# --------------------------------

print(f"Reading twiss data from {twiss_path}")

twiss_data = pd.read_parquet(twiss_path)

# Filtro el twiss para que tenga los mismos elementos que el trackone
twiss_data = (
    twiss_data
    .set_index("NAME")
    .loc[Xs.index]
)


# ----------------------------------------------
#       REGRESIÓN DE LAS ÓRBITAS 
#  Por cada turn, vamos a ajustar acción y fase 
# ----------------------------------------------


if AVM_alg == 2017: 
    avermax_x, avermax_y = avermax_2017(twiss_data, Xs, Ys, 
                                   arc, reference_bpm.lower(), _treshold,
                                   log = True)
elif AVM_alg == 2022: 
    avermax_x, avermax_y = avermax_2022(twiss_data, Xs, Ys, 
                                   arc, reference_bpm.lower(), _treshold,
                                   log = True)


# Guardamos otros datos importantes
S = np.array(list(twiss_data["S"]))
NAMES = list(twiss_data.index)

# S = np.array(list(twiss_data["S"]))
# NAMES = list(twiss_data.index)

zeros_list = np.zeros_like(S)
ones_list = np.ones_like(S)

Jx, Jy, deltax, deltay  = calculate_APJ(twiss_data, avermax_x, avermax_y, last_APJ = True, arcs=[arc, rightArc])




print(f"-> Saving APJ files in {save_path}...")

save_APJ_var(0, S, NAMES, Jx, "HAction.sdds")
save_APJ_var(0, S, NAMES, deltax, "HPhase.sdds")
save_APJ_var(1, S, NAMES, deltay, "VPhase.sdds")
save_APJ_var(1, S, NAMES, Jy, "VAction.sdds")


print(""" 
\t\tF i n i s h e d   A P J   c a l c u l a t i o n
\t\t                                  Version alpha""" )


# Plot if desired
if plot:
    print("\n--> Plotting APJ variables\n")
    plt.plot(S, Jx)
    plt.title("Jx")
    plt.show()

    plt.plot(S, Jy)
    plt.title("Jy")
    plt.show()

    plt.plot(S, deltax)
    plt.title("deltax")
    plt.show()

    plt.plot(S, deltay)
    plt.title("deltay")
    plt.show()


