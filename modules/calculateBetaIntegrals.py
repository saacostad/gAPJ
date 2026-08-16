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



def calculate_integrals(twiss, beam_params):
    return




