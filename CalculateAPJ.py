import tfs as tfs    # To read the twiss files
from turn_by_turn.trackone import read_tbt  # Para leer el trackone
import pandas as pd     # Para poder usar todo eso 

import numpy as np      # Porque siempre va
import matplotlib.pyplot as plt

from modules.ActionPhaseJump import calculate_APJ   
from modules.orbit_tools import avermax   


print(""" \n\n \t\tACTION AND PHASE JUMP
\t\t  Santiago's version \n """ )



# ---------------------------------------------
#           FUNCTION TO SAVE DATA
# ---------------------------------------------
def save_APJ_var(plane, s, NAMES, var, name):
    """ Creates the APJ_var.sdds file in the same format as the OG APJ """

    global save_path 

    with open(f"{save_path}/{name}", 'w') as f:
        
        for row in range(s.shape[0]):

            f.write(f'{plane} "{NAMES[row]}" {s[row]} {var[row]}\n')

    print(f"saved {save_path}/{name}")


# --------------------------------
#       ADJUSTABLES 
# --------------------------------

# TODO: if this script actually works, then I'll have to parse all of this 
# TODO: tho i'd prefer to create a config file, would be funnier

# BPM: BPMSW.1L1.B1
reference_bpm = "QC3L.1"      # The s_e from where we'll select the avermax orbits
# reference_bpm = "BPMSW.1L1.B1"      # The s_e from where we'll select the avermax orbits

folder = "sin_errores"

trackone_path = f"data/{folder}/trackone"          # The original trackone path
twiss_path = f"data/{folder}/twiss.optics"       # The twiss files to use    
avermax_path = f"data/{folder}/avermax.sdds.new"        # Avermax only to make some tests 


save_path = f"data/output"


# The arcs of IP2 so we calculate the avermax traj
arc = (15000, 20000)
rightArc = (23500, 30000)


# --------------------------------
#       LECTURA DEL TRACKONE
# --------------------------------

# Messages to check the state of the script
print(f"Reading turn-by-turn data from {trackone_path}")

# Se lee el TrackOne
to_data = read_tbt(trackone_path)

# Extraemos los datos
particle = to_data.matrices[0]

Xs = pd.DataFrame(particle.X)
Ys = pd.DataFrame(particle.Y)

# -----------------------------------
#       HACEMOS LA TRAYECTORIA DIFF 
# -----------------------------------

print("Fitering betatron motion...")

# Calculamos el promedio en cada BPM
Xs_mean = np.asarray(np.mean(Xs, axis = 1))
Ys_mean = np.asarray(np.mean(Ys, axis = 1))


# TODO: without errors it works amazingly good if we add the mean
# HACK: well there are still jumps but not that big. At least this will work to generalize the APJ
Xs_new = Xs - Xs_mean[:, np.newaxis]
Ys_new = Ys - Ys_mean[:, np.newaxis]


Xs = Xs_new
Ys = Ys_new

# Obtenemos los elementos en las optics de momento
optics_elements_names = Xs.loc[::, 0].index

# --------------------------------
#       LECTURA DEL TWISS
# --------------------------------

print(f"\nReading twiss data from {twiss_path}")

twiss_data = pd.DataFrame(tfs.read(twiss_path))

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


avermax_x, avermax_y = avermax(twiss_data, Xs, Ys, 
                               arc = arc, ref_bpm = reference_bpm,
                               log = True)



# Guardamos otros datos importantes
S = np.array(list(twiss_data["S"]))[1:]
NAMES = list(twiss_data.index)[1:]

# S = np.array(list(twiss_data["S"]))
# NAMES = list(twiss_data.index)

zeros_list = np.zeros_like(S)
ones_list = np.ones_like(S)

Jx, Jy, deltax, deltay  = calculate_APJ(twiss_data, avermax_x, avermax_y, last_APJ = True, arcs=[arc, rightArc])




print("\n\t\t SAVING THE APJ FILES")

save_APJ_var(0, S, NAMES, Jx, "HmaxAPJaction_nofilt.sdds")
save_APJ_var(0, S, NAMES, deltax, "HmaxAPJphase_nofilt.sdds")
save_APJ_var(1, S, NAMES, deltay, "VmaxAPJphase_nofilt.sdds")
save_APJ_var(1, S, NAMES, Jy, "VmaxAPJaction_nofilt.sdds")

print("\n\nFinished APJ calculation")

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


