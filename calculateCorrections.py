"""
This script is in charge of calculating the APJ corrections according to the general formalism
"""

# The must-be packages
import numpy as np 
import pandas as pd
from scipy.optimize import least_squares

# The math package for calculation
from modules.APJ.matricial_system import createSystem_base2 as createSystem   # This is the left hand side of the equation

# Packages for command input
import argparse                 # Will be used mainly to select if we're creating nominal, errors or errors+corrections systems
import tomllib                  # This is to parse the config file
from modules.data_tools.simulateSystem_parser import create_parser, parse_system, parse_commands  # To parse the code's parameters

print("""
      \t\tC O R R E C T I O N S   C A L C U L A T I O N
      \t\t              version alpha                   

Using the APJ general formalism""")

# ----------------------------
# Parse config file arguments
# ----------------------------

print("Parsing system...")

parser = argparse.ArgumentParser()
create_parser(parser)
parsed_args = parser.parse_args()

lattice_config = None       # Configuration dict for general use
nomi_config = None
errs_config = None
cors_config = None

# Read the config file and save it 
print(f"Parsing config file arguments in {parsed_args.config_file}...")
with open(parsed_args.config_file, "rb") as f:

    general_config = tomllib.load(f)

    lattice_config = general_config[f"LatticeFiles"]
    nomi_config = general_config[f"NominalSystem"]
    errs_config = general_config[f"ErrorsSystem"]               # We need both systems: to know where the files are
    cors_config = general_config[f"CorrectionsSystem"]          # And to then know where to output

print(f"Read from LatticeFiles, ErrorsSystem and CorrectionsSystem entries")

# For the command line parsing
arg = parsed_args

# Path of the output files we'll be dealing with
out_path = errs_config["main_output_path"] + "/APJ" if not arg.sequence_path else arg.sequence_path
MUXpath = out_path + "/HAction.sdds"
MUYpath = out_path + "/VAction.sdds"
PHASEXpath = out_path + "/HPhase.sdds"
PHASEYpath = out_path + "/VPhase.sdds"

# Path of the integrals file from where to get the lattice functions
nomi_out_path = nomi_config["main_output_path"] if not arg.twiss_path else arg.twiss_path
integrals_path = nomi_out_path + "/" + nomi_config["integrals_path"] + ".parquet"

# Correction settings
cors_config["left_arc"][0] = float(arg.left_arc_start) if arg.left_arc_start else cors_config["left_arc"][0] 
cors_config["left_arc"][1] = float(arg.left_arc_end) if arg.left_arc_end else cors_config["left_arc"][1] 
cors_config["right_arc"][0] = float(arg.right_arc_start) if arg.right_arc_start else cors_config["right_arc"][0] 
cors_config["right_arc"][1] = float(arg.right_arc_end) if arg.right_arc_end else cors_config["right_arc"][1]
cors_config["modifications_path"] = arg.modifications_path if arg.modifications_path  else cors_config["modifications_path"]
ip = arg.ip if arg.ip  else cors_config["ip"]


# Number of correction quadrupoles / 2
N = int(len(cors_config["correction_quadrupoles"]) / 2)

# We'll create a function to calculate the APJ parameters easily
def get_APJ_parameter(path, axis, left_arc, right_arc, _ip = ip):
    """ Given a sdds file with the APJ calculations, it formats it on the region of interest, filters the axis to deal with, 
    and calculates the average of this parameter on the arc. 

    Returns the average value of the important parameter """
    
    # Which columns we'll be looking for
    filter = 0 if axis == 'X' else 1
    
    # Values list
    left_values = []
    right_values = []
    ip_value = 0.0

    with open(path, 'r') as file:

        # Iterate the file
        lines = file.readlines()

        for line in lines:
            
            # Separate each column
            data = line.split()
            
            # Get the data from the line
            axis_val = int(data[0])
            value = float(data[-1])
            s = float(data[2])

            # Check if it belongs to any of the interest regions/filters
            if axis_val == filter:
                if left_arc[0] <= s <= left_arc[1]:
                        left_values.append(value)
                elif right_arc[0] <= s <= right_arc[1]: 
                    right_values.append(value)
                elif data[1] == f'"{_ip}"': 
                    # TODO: this whole function sucks because I'm using column index instead of a good data structure
                    ip_value = value

       

    # Return the mean of the values encountered
    return np.mean(left_values), np.mean(right_values), ip_value




def get_observed_system(mxp, myp, pxp, pyp, where = "arcs"):
    """ Given the 4 APJ .sdds paths, this function gets the values of each one of the APJ variables and,
    according to the theory, creates the right hand side vector to be solved by the system of equations 
    
    where: "arcs" for left-right arc analysis, "left" for left arc - ip analysis, "right" for ip - right arc analysis

    OUTPUT: np.array with the constants of RHS  |   value of delta_0x   |   value of delta_0y"""

    # Obtenemos acciones y fases para eje x
    J0x, J1x, JxIP = get_APJ_parameter(mxp, 'X', leftArc, rightArc)
    P0x, P1x, PxIP = get_APJ_parameter(pxp, 'X', leftArc, rightArc)

    # Igualmente para el eje y
    J0y, J1y, JyIP = get_APJ_parameter(myp, 'Y', leftArc, rightArc)
    P0y, P1y, PyIP = get_APJ_parameter(pyp, 'Y', leftArc, rightArc)



    def calculate_S_contribution(J0, J1, P0, P1):
        """ Calculates the \\sin(\\psi_s) contribution according to the system of equations """
        return np.sqrt(J1/J0)*np.cos(P1) - np.cos(P0)

    def calculate_C_contribution(J0, J1, P0, P1):
        """ Calculates the \\cos(\\psi_s) contribution according to the system of equations """
        return -np.sqrt(J1/J0)*np.sin(P1) + np.sin(P0)
    
    if where == "arcs":
        # We calculate the RHS constants for the arcs
        _J0x, _J1x, _P0x, _P1x = J0x, J1x, P0x, P1x
        _J0y, _J1y, _P0y, _P1y = J0y, J1y, P0y, P1y
    elif where == "left":
        # Or for the left/right regions
        _J0x, _J1x, _P0x, _P1x = J0x, JxIP, P0x, PxIP
        _J0y, _J1y, _P0y, _P1y = J0y, JyIP, P0y, PyIP
    elif where == "right":
        _J0x, _J1x, _P0x, _P1x = JxIP, J1x, PxIP, P1x
        _J0y, _J1y, _P0y, _P1y = JyIP, J1y, PyIP, P1y

    SinCont_X = calculate_S_contribution(_J0x, _J1x, _P0x, _P1x)
    cosCont_X = calculate_C_contribution(_J0x, _J1x, _P0x, _P1x)
    SinCont_Y = calculate_S_contribution(_J0y, _J1y, _P0y, _P1y)
    cosCont_Y = calculate_C_contribution(_J0y, _J1y, _P0y, _P1y)

    # Return the RHS vector
    return np.array([SinCont_X, cosCont_X, SinCont_Y, cosCont_Y]), _P0x, _P0y



# TODO: make this function read the original twiss files
def get_quadrupoles_lattice_functions(path, QPlist):
    """ Function that reads the integrals file generated and filters out the beta and phi functions
    for the quadrupoles that are being used at the moment """
    
    # Read the original data
    int_data = pd.read_parquet(path) 
   
    # Select only the desired rows
    df = int_data[np.isin(int_data["NAME"], QPlist)]

    df["MUX"] *= 2.*np.pi
    df["MUY"] *= 2.*np.pi

    # Return given dataframe
    return df



"""
=================================================================
        MAIN EXECUTION OF THE SCRIPT
=================================================================
"""
if __name__ == '__main__':
    
    print("\nPreparing system...")

    leftArc, rightArc, QUADRUPOLES_SELECTION = cors_config["left_arc"], cors_config["right_arc"], cors_config["correction_quadrupoles"]

    print(f"  -> Using arcs on s = ({leftArc[0]}, {leftArc[1]}) U ({rightArc[0]}, {rightArc[1]})")
    print(f"  -> Correction quadrupoles = {QUADRUPOLES_SELECTION}")
    
    print("\nCreating system...")

    print("  -> Creating RHS of system of equations")
    # Initial guess for the errors
    ERR_init = np.zeros(len(QUADRUPOLES_SELECTION))

    # First, we get the right hand side vector of the system 
    RHS, delta0_x, delta0_y = get_observed_system(MUXpath, MUYpath, PHASEXpath, PHASEYpath)

    # TODO: make this optional in case we only want arcs analysis
    RHS_L, delta0_x_L, delta0_y_L = get_observed_system(MUXpath, MUYpath, PHASEXpath, PHASEYpath, where = "left")
    RHS_R, delta0_x_R, delta0_y_R = get_observed_system(MUXpath, MUYpath, PHASEXpath, PHASEYpath, where = "right")

    print("  -> Creating LHS of system of equations")
    print(f"  \\__ Getting correction quadrupoles optical parameters from {integrals_path}")

    # In order to create the left hand side, we need to retreive the lattice functions of the quadrupoles of interest
    latticeDF = get_quadrupoles_lattice_functions(integrals_path, QUADRUPOLES_SELECTION)

    # We will now create simple lists of the lattice functions for easier access
    # TODO: make this so it takes a left/right quadrupoles vector 
    BETX = latticeDF['IBX'].to_numpy()      # We're not using BETX/Y because we actually want the integrals of them
    BETY = latticeDF['IBY'].to_numpy()
    MUX = latticeDF['MUX'].to_numpy()
    MUY = latticeDF['MUY'].to_numpy()
    L = latticeDF['L'].to_numpy()
    
    

    BETX_L, BETX_R = BETX[:N], BETX[-N:]
    BETY_L, BETY_R = BETY[:N], BETY[-N:]
    MUX_L, MUX_R = MUX[:N], MUX[-N:]
    MUY_L, MUY_R = MUY[:N], MUY[-N:]    # We'll create the residual function to use with Least_Squares()
    
    save_Q = True

    def residual(K):
        global Qx, Qy

        # We create the constants for both axis
        Sx, Cx, Qx = createSystem(K, BETX, MUX, delta0_x, axis = 'X')
        Sy, Cy, Qy = createSystem(K, BETY, MUY, delta0_y, axis = 'Y')

        Sx_L, Cx_L, _ = createSystem(K[:N], BETX_L, MUX_L, delta0_x_L, axis = 'X')
        Sy_L, Cy_L, _ = createSystem(K[:N], BETY_L, MUY_L, delta0_y_L, axis = 'Y')
        
        Sx_R, Cx_R, _ = createSystem(K[-N:], BETX_R, MUX_R, delta0_x_R, axis = 'X')
        Sy_R, Cy_R, _ = createSystem(K[-N:], BETY_R, MUY_R, delta0_y_R, axis = 'Y') 
        
        # Return the residual
        # return np.array([Sx, Cx, -Sy, -Cy]) - RHS
        return np.array([Sx, Cx, -Sy, -Cy, Sx_L, Cx_L, -Sy_L, -Cy_L, Sx_R, Cx_R, -Sy_R, -Cy_R]) - np.concatenate((RHS, RHS_L, RHS_R))


    print("\nSolving the system...")

    """ CALCULATE THE ERRORS STIMATIONS """
    ERR_estimations = least_squares(residual, ERR_init, ftol = 1e-15).x
   
    # --> PRINT Q MATRICES TO CHECK VALIDITY OF GENERAL FORMULATION
    # print("-"*30, "\n   Q matrices of the system")
    # with np.printoptions(precision=1, suppress=True):
    #     print("---> Para X:")
    #     print(Qx)
    #     print("---> Para Y:")
    #     print(Qy)


    print("="*25)
    print("\nErrors estimation: \n")
    
    for i in range(len(ERR_estimations)):
        # print(f"\t{QUADRUPOLES_SELECTION[i]}:  {ERR_estimations[i]:.2g} \t\t residue of {residual(ERR_estimations)[i]:.2g}")
        print(f"\t{QUADRUPOLES_SELECTION[i]}:  {ERR_estimations[i]:.2g}")

    # Write to the file
    # Build lookup dictionary
    err_dict = dict(zip(QUADRUPOLES_SELECTION, ERR_estimations))
    
    print(f"\n  -> Writting errors + corrections to {cors_config["modifications_path"]}")

    # Read and modify file
    with open(cors_config["modifications_path"], "r") as f:
        lines = f.readlines()

    new_lines = []

    for line in lines:
        # Get the quadrupole name
        quad_name = line.split("\t")[0]

        if quad_name in err_dict:

            # Get the correction, original value and new value
            err = err_dict[quad_name]
            or_val = float(line.split("\t")[1])
            new_val = err + or_val
                
            print(f"  \\__ {quad_name}: {or_val:.3g} + {err:.2g} = {new_val:.3g}")

            # Remove trailing semicolon/newline, append new term, add semicolon back
            parts = line.strip().split("\t")    # Get the cols of the file
            parts[1] = str(new_val)
            line = "\t".join(parts) 
            line = line + "\n"

        new_lines.append(line)

    # Write back
    with open(cors_config["modifications_path"], "w") as f:
        f.writelines(new_lines)


    print("Finished writing to file")
