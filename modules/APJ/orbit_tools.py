"""
This script contains functions for orbit calculations:

    adjust_turn: adjusts the betatron equations using the nominal lattice functions on J and delta 
    avermax: creates the avermax orbit
"""


import numpy as np      # Porque siempre va
from scipy.optimize import least_squares as ls      # Para hacer el ajuste del avermax
import pandas as pd     # Para poder usar todo eso 

from modules.APJ.ActionPhaseJump import calculate_APJ


# --------------------------------------------------
#       MAPPING PHASE
# --------------------------------------------------
def phase_map(phase):
    """ This function maps a phase from -pi to pi """
    
    return (phase + np.pi) % (2 * np.pi) - np.pi

#--------------------------------------------------------------------------------------
#               AVERMAX FUNCTION
#--------------------------------------------------------------------------------------

def avermax(twiss_data, Xs, Ys, arc, ref_bpm, treshold, log = False):
    """ This function takes the twiss data, orbits (x and y) and other parameters to calculate the avermax trajectory """

    if log: 
        print(f"""
STARTING AVERMAX TRAJECTORY CALCULATION
-> Total turns read from trackone: {len(Xs.columns)}
-> Using a treshold: {treshold%np.pi:.2f}pi
-> Reference BPM: {ref_bpm}
        """)

    # We get the reference phases
    ref_bpm_phase_x = phase_map(2.*np.pi*((twiss_data.loc[ref_bpm, 'MUX']) % (1.)) - (np.pi / 2))
    ref_bpm_phase_y = phase_map(2.*np.pi*((twiss_data.loc[ref_bpm, 'MUY']) % (1.)) - (np.pi / 2)) 
   
    # The reference for the max min
    neg_ref_bpm_phase_x = phase_map(2.*np.pi*((twiss_data.loc[ref_bpm, 'MUX']) % (1.)) + (np.pi / 2))
    neg_ref_bpm_phase_y = phase_map(2.*np.pi*((twiss_data.loc[ref_bpm, 'MUY']) % (1.)) + (np.pi / 2))

    # We'll create a mask on the arc region we're interested in 
    arc_mask = twiss_data['S'].between(arc[0], arc[1])
    
    # Select the respective columns
    twiss_filtered = twiss_data[arc_mask]
    Xs_filtered = Xs[arc_mask]
    Ys_filtered = Ys[arc_mask]

    # This is the list where we'll save the turns
    x_turns = list()
    y_turns = list()


    count_x = 0   # To count how many orbits we'll be using in the avermax
    count_y = 0   
    for col in Xs_filtered.columns:
        

        # Get the z positions on the respective ARC
        x = Xs_filtered[col]
        y = Ys_filtered[col] 
        
        # Now, we get the action and phase there
        _, _, deltax, deltay = calculate_APJ(twiss_filtered, x, y, mean=True)
        deltax = deltax
        deltay = deltay

        # Save the phase differences and map them from -pi to pi for the max turns
        Delta_delta_x = phase_map(deltax - ref_bpm_phase_x)
        Delta_delta_y = phase_map(deltay - ref_bpm_phase_y)


        # Save the phase differences and map them from -pi to pi for the min turns
        neg_Delta_delta_x = phase_map(deltax - neg_ref_bpm_phase_x)
        neg_Delta_delta_y = phase_map(deltay - neg_ref_bpm_phase_y)
        

        # We want the delta to be close psi_s_e +- \pi/2 
        if abs(Delta_delta_x) < treshold:

            # If both conditions are satisfied, then we can add the turn to the respective list 
            x_turns.append(Xs[col].to_numpy())
            count_x += 1
        
        # Now we also select the negative maxima
        elif abs(neg_Delta_delta_x) < treshold:

            # But here, we'll flip the sign
            x_turns.append(-Xs[col].to_numpy())
            count_x += 1

        # -- Same but for the y axis

        if  abs(neg_Delta_delta_y) < treshold:
            y_turns.append(-Ys[col].to_numpy())
            count_y += 1

        # Now we also select the negative maxima
        elif abs(Delta_delta_y) < treshold:
            
            # But here, we'll flip the sign
            y_turns.append(Ys[col].to_numpy())
            count_y += 1
    
    # Lastly, we get the mean of the turns and multiply by 1000
    avermax_x = np.mean(np.array(x_turns), axis = 0)
    avermax_y = np.mean(np.array(y_turns), axis = 0)
        
    
    if log:
        print(f"-> Used {count_x} X orbits from a total of {(count_x / len(Xs.columns))*100:.2f}% for the avermax trajectory")
        print(f"-> Used {count_y} Y orbits from a total of {(count_y / len(Xs.columns))*100:.2f}% for the avermax trajectory")
        print(f"\nAVERMAX trajectory calculation finished.")


    return avermax_x, avermax_y









