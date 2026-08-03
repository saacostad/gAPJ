""" 
This script will calculate the integrals of the beta and psi functions inside 
an element of the accelerator
"""
from os.path import join

import pandas as pd
import tfs 
import numpy as np
from scipy.integrate import quad     # To integarte
from cpymad.madx import Madx



def get_beta_function(sign, KK, beta0, alfa0):
    """ This function returns a ready-to-use function representing beta(s) inside of a quadrupole.
    It takes the sign (direction of the beam * plane * charge of the particle, either -1 or 1 depending on the convention used),
    the element's magnetic gradient K and entrance beta_0 and alfa_0 """
    
    K = np.sqrt(KK)             # Helper variable
    case = np.sign(sign * KK)   # Get if we're focusing or defocusing
    
    # Focusing case 
    if case > 0:

        def beta(s):
            ks = K*s 
            sin = np.sin(ks)
            cos = np.cos(ks)
            return beta0*cos**2 + 2.*alfa0*(sin*cos/K) + ((1.+alfa0**2)/beta0)*(sin**2)/KK
        return beta

    else: 
        # Defocusin case: change sin/cos for sinh/cosh
        def beta(s):
            ks = K*s 
            sin = np.sinh(ks)
            cos = np.cosh(ks)
            return beta0*cos**2 + 2.*alfa0*(sin*cos/K) + ((1.+alfa0**2)/beta0)*(sin**2)/KK
        return beta





def calculate_integral(madx, sign, optics_class):
    """ Performs the calculation of the integrals, takes: 
        madx: the cpymad.Madx session in order to calculate the twiss at the entrance,
        sign: the sign of the particle's charge * direction of the beam
        optics_class: the list of the keywords of the magnetic elements to get"""
    
    print(f"--> performing integrals calculation...")
    
    print("     \\__ Performing entrance elements twiss")

    madx.input("SELECT, FLAG=interpolate, CLEAR;")
    madx.input(f"SELECT, FLAG=interpolate, CLASS={", ".join(optics_class)}, AT=0.0;")
    madx.input("SELECT, FLAG=twiss, CLEAR;")
    madx.input("SELECT, FLAG=twiss, COLUMN=name, s, betx, bety, mux, muy, alfx, alfy;")

    twiss = madx.twiss().dframe() # Perform twiss
    
    # Get only the elements we're interested in
    twiss = twiss[twiss['keyword'].isin(optics_class)][["name", "s", "betx", "bety", "mux", "muy", "alfx", "alfy"]]   # This is kinda the best filtering I can do, but it's more than enough
   

    print("     \\__ Performing integrals")
    # # -- Now, we'll perform the integrals
    results = []
    
    # Convert this to a dataframe
    res_df = pd.DataFrame(results)
    tfs.write('integrals_twiss.tfs', res_df)




