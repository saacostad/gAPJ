"""
This script constains everything needed for the calculation of the QGCs 
according to the general model
"""

import numpy as np 
from scipy.integrate import quad


def get_lattice_functions(_sign, sign_K, KK, beta0, alfa0, psi0):
    """ This function returns a ready-to-use function representing beta(s) and psi(s) inside of a quadrupole.
    It takes the sign (direction of the beam * plane * charge of the particle, either -1 or 1 depending on the convention used),
    the element's magnetic gradient K and entrance beta_0 and alfa_0 """
    
    # In the case we get weird quadrupoles, just return 0
    if KK == 0.0:
        return lambda x: 0.0

    K = np.sqrt(abs(KK))            # Helper variable
    case = np.sign(_sign * sign_K)  # Get if we're focusing or defocusing
    
    beta = lambda x: None           # Create a function holder

    # Focusing case 
    if case > 0:

        def beta(s):
            ks = K*s 
            sin = np.sin(ks)
            cos = np.cos(ks)
            return beta0*cos**2 - 2.*alfa0*(sin*cos/K) + ((1.+alfa0**2)/beta0)*(sin**2)/KK

    else: 
        # Defocusin case: change sin/cos for sinh/cosh
        def beta(s):
            ks = K*s 
            sinh = np.sinh(ks)
            cosh = np.cosh(ks)
            return beta0*cosh**2 - 2.*alfa0*(sinh*cosh/K) + ((1.+alfa0**2)/beta0)*(sinh**2)/KK
    
    # Psi only depends on the beta function, let it be focusing or defocusing case
    def psi(s):
        return psi0 + quad(lambda x: 1. / beta(x), 0, s)[0]

    
    return beta, psi





# TODO: here is where we'll change the integrals
def CreatePropagationElements(B, P, A, K, L, sign):
    """ We'll create the M propagation matrix and the u and v propagation vectors elements.
    -- Inputs:
    B: NOMINAL ENTRANCE beta functions of the correction quadrupoles
    P: Same for Psi functions
    A: Same but alpha
    K: Same but the quadrupole gradient strength
    sign: Sign for wheter we are focusing or defocusing (beams direction * axis)
    """

    # CREATE THE LATTICE FUNCTIONS
    # _______________________________________

    # Create a list for the functions
    betas = list()
    psis = list()                       # HACK: psis is gonna be use in the case we slice the quadrupoles
    int_betas = list()                  # Create a list for the beta integrals which we'll be using
    for i in range(len(B)):
        # Retreive the respective functions
        beta, psi = get_lattice_functions(sign, np.sign(K[i]), np.abs(K[i]), B[i], A[i], P[i])
        int_beta = quad(beta, 0, L[i])[0]

        betas.append(beta)
        psis.append(psi)
        int_betas.append(int_beta)
    

    # CREATE THE U AND V VECTORS
    # _______________________________________

    temp_r = int_betas / np.sqrt(B)         # Common term for the propagation vectors 
    temp_r = np.sqrt(int_betas)
    u = np.sin(P) * temp_r
    v = np.cos(P) * temp_r


    # CREATE THE PROPAGATION MATRIX
    # _______________________________________   

    # -- Prepare the matrix
    n = len(P)  # Number of Quadrupoles
    hatM = np.zeros((n, n))    # Create a square matrix first
    ii, jj = np.tril_indices(n, k=-1) # Take only the indexes of the lower triangular section


    # -- Create the M element evaluator function
    def evaluate_M_element(i, j):
        # return np.sqrt(B[j] / B[i]) * int_betas[i] * np.sin(P[j] - P[i])
        return np.sqrt(int_betas[i] * int_betas[j]) * np.sin(P[j] - P[i])


    for i in ii:
        for j in jj: 
            hatM[i, j] = evaluate_M_element(i, j)    # Fill the gaps with what we need 

    return hatM, u, v






def CreateQ(hatM, K, p):
    """
    This function creates the Q matrix given that 

    \\hat Q = \\hat{\\Delta K} ( 1 - p \\hat{M} \\hat{\\Delta K} )^{-1}

    with 

    \\hat{\\Delta K} = \\Diag(\\Delta K_i)
    M_{ij} = \\sin(\\psi_i - \\psi_j) \\sqrt(\\beta_j \\beta_i) for j < i else 0

    where grad is the "order of the inverse" of the term in parenthesis, given by truncating 
    the Neumman series of p \\hat{M} \\hat{\\Delta K} when 

    ( 1 - p \\hat{M} \\hat{\\Delta K} )^{-1} = \\sum_{k = 1}^{grad} (p \\hat{M} \\hat{\\Delta K} )^k
    """

    # Matrix form of \Delta K
    hatK = np.diag(K)
   
    # Now we create the whole matrix in parenthesis and we invert it 
    #---------------------------------------------------------------------------
    
    n = len(K) 
    Q = hatK @ np.linalg.inv(np.eye(n) - p * hatM @ hatK)

    return Q

    

def CreateConstants(Q, u, v):
    """
    This function takes the system vectors/matrices and calculates the coefficients 
    of the expansion we care about. Returns a tuple with 

    \\vec v^T \\hat Q \\vec u
    \\vec v^T \\hat Q \\vec v
    \\vec v^u \\hat Q \\vec u
    \\vec v^u \\hat Q \\vec v
    """

    return v @ Q @ u, v @ Q @ v, u @ Q @ u, u @ Q @ v


"""
MAIN FUNCTION TO CREATE THE SYSTEM
"""
def createSystem_base2(M, u, v, Kerrors, delta0, axis):
    """ This function takes the Betas, Phis and errors (and the \\delta_0 value) and creates and evaluates
    the "Simulated" side of the system of equations.

    Inputs: Kerros, Betas, Phis -> np.array of len QP   
            delta0 ->   constant \\delta_0   
            axis: {'X', 'Y'} -> decides the p value
    """

    p = 1.0 if axis == 'X' else -1.0

    # Create the matricial system
    Q = CreateQ(M, Kerrors, p)

    # Create the big constants
    vQu, vQv, uQu, uQv = CreateConstants(Q, u, v)

    # The cos(\delta_0) for the calculations
    cosd = np.cos(delta0)
    sind = np.sin(delta0)
    
    # Return the two cosntants
    return cosd*vQu - sind*vQv, -cosd*uQu + sind*uQv

    
    
