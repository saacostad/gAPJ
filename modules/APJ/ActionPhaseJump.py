""" 
This module contains all the functions needed to calculate the Action and Phase Jumps
"""

import numpy as np 
import warnings         # So I get prettier warnings

# Custom warning to make everything cleaner
def custom_warning(message, category, filename, lineno, file=None, line=None):
    print(f" -W- {filename.split('/')[-1]}{lineno} → {message}")

warnings.showwarning = custom_warning



# -----     CALCULATION USING THE 2009 FORMULAS
def calculate_J(BETA, PSI, Z):
    """
    Compute local action J for each adjacent BPM pair.
    """
    # Convert phase from turns to radians
    psi_rad = 2.0 * np.pi * PSI
    
    # Shift arrays for pair (i, i+1)
    beta_i = np.sqrt(BETA[:-1])          # sqrt(beta_i)
    beta_j = np.sqrt(BETA[1:])           # sqrt(beta_{i+1})
    z_i = Z[:-1]
    z_j = Z[1:]
    psi_i = psi_rad[:-1]
    psi_j = psi_rad[1:]
    dpsi = psi_j - psi_i
    
    # Avoid singularities 
    sin_dpsi = np.sin(dpsi)
    mask = np.abs(sin_dpsi) < 1e-12
    # sin_dpsi[mask] = np.sign(sin_dpsi[mask]) * 1e-12 
    

    # Action formula
    J = ((z_i / beta_i)**2 + (z_j / beta_j)**2 - (2.0 * np.cos(dpsi) * z_i * z_j / (beta_i * beta_j))) / (2.0 * sin_dpsi**2)
    
    return J


def calculate_delta(BETA, PSI, Z):
    """
    Compute local phase delta for each adjacent BPM pair.
    Returns: array delta of length N-1 (radians)
    """
    psi_rad = 2.*np.pi*PSI
    
    beta_i = np.sqrt(BETA[:-1])
    beta_j = np.sqrt(BETA[1:])
    z_i = Z[:-1]
    z_j = Z[1:]
    psi_i = psi_rad[:-1]
    psi_j = psi_rad[1:]
    
    numerator = (z_i / beta_i) * np.sin(psi_j) - (z_j / beta_j) * np.sin(psi_i)
    denominator = (z_i / beta_i) * np.cos(psi_j) - (z_j / beta_j) * np.cos(psi_i)    

    delta = np.arctan2(numerator, denominator)
    return delta




# -----     CALCULATION USING A MORE ROBUST FORMULA
def calculate_J_and_P(BETA, PSI, Z):
    """ Action J and Phase P calculation from the inversion formula """

    # Convert phase from turns to radians
    psi_rad = 2.0 * np.pi * PSI
    
    # -- Shift arrays for pair (i, i+1)

    # Get the beta functions
    beta_1 = np.sqrt(np.roll(BETA, -1))         
    beta_2 = np.sqrt(np.roll(BETA, 1))           

    # Get the reduces coordinates
    z_1 = np.roll(Z, -1) / beta_1               # This is the z_red on the original code
    z_2 = np.roll(Z, 1) / beta_2
    
    # Get the phases
    psi_1 = np.roll(psi_rad, -1)
    psi_2 = np.roll(psi_rad, 1)
    dpsi = psi_1 - psi_2

    # Avoid singularities 
    sin_dpsi = np.sin(dpsi)
    mask = np.abs(sin_dpsi) < 1e-12
    sin_dpsi[mask] = np.sign(sin_dpsi[mask]) * 1e-12
    
    # TODO: Find a way so the denominator cannot be 0
    A = (z_1 * np.sin(psi_2) - z_2 * np.sin(psi_1)) / (sin_dpsi)
    B = (z_2 * np.cos(psi_1) - z_1 * np.cos(psi_2)) / (sin_dpsi)

    # J will simply be (A^2 + B^2)/2
    J = (np.power(A, 2) + np.power(B, 2)) / 2.

    # Phase will be atan2(-A, B)
    P = np.arctan2(-A, B)

    # We'll make the first and last entries be equal to the second ones
    J[0] = J[1]
    J[-1] = J[-2]
    P[0] = P[1]
    P[-1] = P[-2]

    return J, P
    

# TODO: this function could give an stimation of the jumps
def calculate_jump(s, apj_var, arc_left, arc_right):
    """ This function takes the left and right side arcs, one of the APJ variables
    and the S position of the measurements and calculates the jump """

    return 




def calculate_APJ(TWISS, X, Y, mean = False, last_APJ = False, arcs = None):
    """
    This function is in charge of calculating the Action and Phase Jumps given

    TWISS: the pandas dataframe that contains the beta and mu functions
    Z: the orbit's turn data
    """

    if last_APJ: 
        print(f"""STARTING APJ CALCULATION
        """)

    BETX = np.asarray(TWISS["BETX"].to_numpy())
    BETY = np.asarray(TWISS["BETY"].to_numpy())
    MUX = np.asarray(TWISS["MUX"].to_numpy())
    MUY = np.asarray(TWISS["MUY"].to_numpy())

    
    # This block is for the 2009 calculation formula
    # Jx = calculate_J(BETX, MUX, np.asarray(X))
    # Jy = calculate_J(BETY, MUY, np.asarray(Y))
    # px = calculate_delta(BETX, MUX, np.asarray(X))
    # py = calculate_delta(BETY, MUY, np.asarray(Y))


    # This block is for the new calculation formula
    Jx, px = calculate_J_and_P(BETX, MUX, np.asarray(X))
    Jy, py = calculate_J_and_P(BETY, MUY, np.asarray(Y))


    if mean: 
        return np.mean(Jx), np.mean(Jy), np.angle(np.mean(np.exp(1j * px))), np.angle(np.mean(np.exp(1j * py))) 

    if last_APJ: 
        print(f"-> Finished APJ calculation")
        
        if arcs != None: 
            print("-> Stimated APJ jumps: TODO")

    return Jx, Jy, px, py


