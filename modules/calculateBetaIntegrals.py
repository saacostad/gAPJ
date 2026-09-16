""" 
This script will calculate the integrals of the beta and psi functions inside 
an element of the accelerator
"""
from os.path import join

import pandas as pd
import numpy as np
from scipy.integrate import quad     # To integarte



def get_beta_function(_sign, sign_K, KK, beta0, alfa0):
    """ This function returns a ready-to-use function representing beta(s) inside of a quadrupole.
    It takes the sign (direction of the beam * plane * charge of the particle, either -1 or 1 depending on the convention used),
    the element's magnetic gradient K and entrance beta_0 and alfa_0 """
    
    # In the case we get weird quadrupoles, just return 0
    if KK == 0.0:
        return lambda x: 0.0

    K = np.sqrt(abs(KK))         # Helper variable
    case = np.sign(_sign * sign_K)   # Get if we're focusing or defocusing
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
            sinh = np.sinh(ks)
            cosh = np.cosh(ks)
            return beta0*cosh**2 + 2.*alfa0*(sinh*cosh/K) + ((1.+alfa0**2)/beta0)*(sinh**2)/KK
        return beta



def calculate_integrals(twiss, beam_params):
    
    # First, we calculate the general sign
    sign = np.sign(beam_params["dir"] * beam_params["charge"])      # TODO: check this bitch because I've got no clue if I gotta add a - or not
    # Here I will create the new pandas dataframe 
    results = list()

    # And here I'll make the calculation to populate it
    for row in twiss.itertuples():




        # TODO: I was passing K^2 but maybe it is actually K
        betax_func = get_beta_function(sign, np.sign(row.K1), np.abs(row.K1), row.BETX, row.ALFX)
        betay_func = get_beta_function(-sign, np.sign(row.K1), np.abs(row.K1), row.BETY, row.ALFY)

        integralx = quad(betax_func, 0, row.L)[0]
        integraly = quad(betay_func, 0, row.L)[0]

        results.append({
            "NAME": row.NAME,
            "S": row.S,
            "BETX": row.BETX,
            "BETY": row.BETY,
            "MUX": row.MUX,
            "MUY": row.MUY,
            "IBX": integralx,
            "IBY": integraly,
            "L": row.L
            })

        if row.NAME in ["mqxa.3l1", "mqxb.b2l1", "mqxb.a2l1", "mqxa.1l1", "mqxa.1r1", "mqxb.a2r1"]:
            print("")
            print("QP: ", row.NAME)
            print("k - ", row.K1)
            print("l - ", row.L)
            print("betax - ", row.BETX)
            print("betay - ", row.BETY)
            print("INT BETX = ", integralx)
            print("INT BETY = ", integraly)

    # Lastly, we create the dataframe
    integrals_dataframe = pd.DataFrame(results)

    return integrals_dataframe




