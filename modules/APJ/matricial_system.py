"""
This script constains everything needed for the calculation of the QGCs 
according to the general model
"""

import numpy as np 
from scipy.integrate import quad

# TODO: would be better to import the function or to do something with the beta function
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
            return beta0*cos**2 - 2.*alfa0*(sin*cos/K) + ((1.+alfa0**2)/beta0)*(sin**2)/KK
        return beta

    else: 
        # Defocusin case: change sin/cos for sinh/cosh
        def beta(s):
            ks = K*s 
            sinh = np.sinh(ks)
            cosh = np.cosh(ks)
            return beta0*cosh**2 - 2.*alfa0*(sinh*cosh/K) + ((1.+alfa0**2)/beta0)*(sinh**2)/KK
        return beta



# TODO: here is where we'll change the integrals
def CreateM(B, P, A, K, L, sign):
    """ We'll create the M propagation matrix.
    -- Inputs:
    B: NOMINAL ENTRANCE beta functions of the correction quadrupoles
    P: Same for Psi functions
    A: Same but alpha
    K: Same but the quadrupole gradient strength
    sign: Sign for wheter we are focusing or defocusing (beams direction * axis)
    """
    
    # -- Prepare the matrix
    n = len(P)  # Number of Quadrupoles
    hatM = np.zeros((n, n))    # Create a square matrix first
    ii, jj = np.tril_indices(n, k=-1) # Take only the indexes of the lower triangular section
    

    #-----------------------------------------  
    #       Create auxiliary functions
    #-----------------------------------------  

    # -- Create an array of the beta and psi functions for every corrections quadrupole
    betas = list()
    psis= list()

    for i in range(n):
        beta_func = get_beta_function(  sign, np.sign(K[i]), np.abs(K[i]), B[i], A[i]  )
        betas.append(beta_func)      # Create beta function

        # Create psi function
        def make_psi(i, beta):
            def psi(s):
                return P[i] + quad(lambda x: 1 / beta(x), 0, s)[0] 

            return psi

        psis.append(make_psi(i, beta_func))

    
    # -- We create the integrands 
    def create_int_with_sin(beta, psi):

        def integrand(s):
            return np.sqrt(beta(s)) * np.sin(psi(s))

        return integrand

    def create_int_with_cos(beta, psi):

        def integrand(s):
            return np.sqrt(beta(s)) * np.cos(psi(s))

        return integrand


    # # -- Create a function that evaluates the M matrix components
    # def evaluate_M_element(i, j):
    #     return -quad(lambda s: create_int_with_cos(betas[i], psis[i])(s), 0, L[i])[0] * quad(lambda s: create_int_with_sin(betas[j], psis[j])(s), 0, L[j])[0]
    #     + quad(lambda s: create_int_with_cos(betas[j], psis[j])(s), 0, L[j])[0] * quad(lambda s: create_int_with_sin(betas[i], psis[i])(s), 0, L[i])[0]
   

    # -- Create a function that evaluates the M matrix components
    def evaluate_M_element(i, j):
        return -np.sqrt(B[i])*np.cos(P[i])* quad(lambda s: create_int_with_sin(betas[j], psis[j])(s), 0, L[j])[0]
        +quad(lambda s: create_int_with_cos(betas[j], psis[j])(s), 0, L[j])[0] * np.sqrt(B[i])*np.cos(P[i]) 

    #-----------------------------------------  
    #       CREATE THE MATRIX
    #-----------------------------------------  

    for i in ii:
        for j in jj: 
            hatM[i, j] = evaluate_M_element(i, j)    # Fill the gaps with what we need 

    return hatM





def NeummanInverse(T, grad = 8):
    """
    This function creates the Neuman Truncated Inverse of

    (1 - p\\hat M \\hat {\\Delta K})

    according to the formula

    (1 - p\\hat M \\hat {\\Delta K})^{-1} = sum_{k = 1}^{grad} (p \\hat M \\hat \\Delta K)^k 
    """
    n = T.shape[0]
    I = np.eye(n, dtype=T.dtype)
    
    S = I.copy()
    Tk = I.copy()
    
    for _ in range(1, grad):
        Tk = Tk @ T
        S += Tk
        
    return S


def CreateQ(hatM, K, p, grad = 8):
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
    if grad == hatM.shape[0]:
        Q = hatK @ np.linalg.inv(np.eye(n) - p * hatM @ hatK)
    else:
        # print("calculating inverse with Neumman")
        Q = hatK @ NeummanInverse(p * hatM @ hatK, grad)

    return Q

    

def CreateSystem(hatM, B, P, K, p = 1.0, grad = 8):
    """
    This function takes the system's \\eta's, \\phi's and quadrupole magnetic errors \\Delta K
    and returns a tuple with the matrix \\hat Q, and the vectors \\vec u and \\vec v following that

    \\vec u = {\\sqrt(\\beta_i) \\sin(\\psi_i)}
    \\vec v = {\\sqrt(\\beta_i) \\cos(\\psi_i)}

    \\hat Q = \\hat{\\Delta K} ( 1 - p \\hat{M} \\hat{\\Delta K} )^{-1}
    """

    u = np.sqrt(B) * np.sin(P)
    v = np.sqrt(B) * np.cos(P)

    Q = CreateQ(hatM, K, p, grad)

    return Q, u, v


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


def createFirstOrderMatrix(ux, vx, uy, vy):
    """ This function creates the "first order matrix", that is, approximation the system to first order 
    and writting \\hat G \\vec Delta k = \\vec a, \\hat G would be the first order matrix """
    
    rows = [ux * vx, vx * vx,ux * ux, ux * vx,
            uy * vy, vy * vy,uy * uy, uy * vy,]

    return np.array(rows)




"""
HERE I'LL TRY TO CREATE THE SYSTEM WITH ONLY 2 BASES SIN(PSI) AND COS(PSI)
"""
def createSystem_base2(M, Kerrors, Betas, Phis, delta0, axis, grad = 8):
    """ This function takes the Betas, Phis and errors (and the \\delta_0 value) and creates and evaluates
    the "Simulated" side of the system of equations.

    Inputs: Kerros, Betas, Phis -> np.array of len QP   
            delta0 ->   constant \\delta_0   
            axis: {'X', 'Y'} -> decides the p value
    """

    p = 1.0 if axis == 'X' else -1.0

    # Create the matricial system
    Q, u, v = CreateSystem(M, Betas, Phis, Kerrors, p, grad)

    # Create the big constants
    vQu, vQv, uQu, uQv = CreateConstants(Q, u, v)

    # The cos(\delta_0) for the calculations
    cosd = np.cos(delta0)
    sind = np.sin(delta0)
    
    # Return the two cosntants
    return cosd*vQu - sind*vQv, -cosd*uQu + sind*uQv

    
    
