""" 
This script will calculate the integrals of the beta and psi functions inside 
an element of the accelerator
"""
import pandas as pd
import tfs 
import numpy as np
from scipy.integrate import quad     # To integarte
from cpymad.madx import Madx





def calculate_integral(madx, sequence_name, optics_class):
    """ Performs the calculation of the integrals """
    
    print(f"--> performing integrals calculation...")

    optics_lower = [opt.lower() for opt in optics_class]  # A list with the optics keywords in lower just in case
    seq = madx.sequence[sequence_name]      # Create a new sequence that we can actually pass


    print("     \\__ Getting elements information....")
    # -- Will get information about the elements
    seq_elements = seq.elements     # So the for loop doesn't call it each time

    element_info = []
    for elem in seq_elements:

        # Get the keyword
        keyword = elem.base_type.name

        # If the element is not of our interests, skip it 
        if keyword.lower() not in optics_lower: continue

        # If we're dealing with an optics element, then we get it's values
        length = getattr(elem, "l", 0.0)
        centre = getattr(elem, "at", 0.0)
        K = getattr(elem, 'k1l', 0.0)

        # We'll append this data to the list 
        element_info.append({
                'name': elem.name,
                's_start': centre - length/2.,
                's_end': centre + length/2.,
                'length': length,
                's': centre,
                'k': K
            })

    elem_df = pd.DataFrame(element_info)    # Convert this to a dataframe for better use
    
    print("     \\__ Performing entrance elements twiss")

    madx.input("SELECT, FLAG=interpolate, CLEAR;")
    for optic in optics_class:
        madx.input(f"SELECT, FLAG=interpolate, CLASS={optic}, AT=0.0;")
    madx.input("SELECT, FLAG=twiss, CLEAR;")
    madx.input("SELECT, FLAG=twiss, COLUMN=name, s, betx, bety, mux, muy, alfx, alfy;")

    twiss = madx.twiss().dframe() 
    twiss = twiss[twiss["keyword"] == "quadrupole"][["name", "s", "betx", "bety", "mux", "muy", "alfx", "alfy"]]   # This is kinda the best filtering I can do, but it's more than enough
   



    print("     \\__ Performing integrals")
    # # -- Now, we'll perform the integrals
    results = []
    
    # We'll iterate over each quadrupole
    for _, row, in elem_df.iterrows():
    
        # Get the info of the main elements
        name = row['name']
        KK = row['k']
        K = np.sqrt(KK)
        
        # TODO: calculate the integrals with the analytic formula: 1. retrieve the betz and alfz 
        # from the twiss dataframe (the names end with :1). 2. Program the beta function.

        # Save in the dictionary
        results.append({
                'NAME': row['name'].upper(),
                'S': row['s'],
                'BETX': int_betx,   # Save the integrals here
                'BETY': int_bety,
                'MUX': np.mean(subset["mux"]),
                'MUY': np.mean(subset["muy"]),
                'L': row['length']
            })
    
    # Convert this to a dataframe
    res_df = pd.DataFrame(results)
    tfs.write('integrals_twiss.tfs', res_df)




