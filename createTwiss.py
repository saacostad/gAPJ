"""
This script is in charge of reading the .madx sequence file 
for a given accelerator and create the necessary twiss files 
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Packages for command input
import argparse                 # Will be used mainly to select if we're creating nominal, errors or errors+corrections systems
import tomllib                  # This is to parse the config file

# Libraries for madx
from cpymad.madx import Madx
from modules.data_tools.read_beam_parameters import read_beam_parameters  # To parse the parameters .txt
from modules.data_tools.simulateSystem_parser import create_parser_args, parse  # To parse the code's parameters

# ---------------------------
# MAIN FUNCTION
# ---------------------------
def create_twiss(parameters_path, sequence_path, 
                 create_measurement_twiss, create_quads_data,
                 twiss_path, quadrupoles_parameters_path,
                 errors_path = None):

    # -- Create a mad-x connection
    # Read beam parameters
    beam_params = read_beam_parameters(parameters_path)

    # Create madx connection
    madx = Madx() 
    # Let mad-X log on terminal 
    madx.option(echo=False, info=False)

    # We define the beam to use 
    beam_input = "BEAM, " + ", ".join(f"{key}={value}" for key, value in beam_params.items())
    madx.input(beam_input)


    if create_quads_data or create_measurement_twiss:
        # We call the sequence
        madx.input(f'CALL, FILE="{sequence_path}";')
    else: 
        print("No option was given.")


    if create_measurement_twiss:

        # We perform the twiss
        # HACK: We only want certain twiss parameters, not all the table
        madx.input(f'USE, SEQUENCE="{sequence_name}";')     # Select the lattice to use

        # If we have errors or corrections, apply them
        if errors_path != None:
            madx.input(f'CALL, FILE="{errors_path}";')

        madx.input(f'SELECT, FLAG=twiss, CLEAR;')           # Tell the twiss to clear everything
    
        # We will select all the measuring elements for our twiss
        for measure in measure_class:
            madx.input(f'SELECT, FLAG=twiss, CLASS={measure}, COLUMN={", ".join(twiss_parameters)};')

        # Call the twiss function
        madx.twiss(file=twiss_path)


    if create_quads_data:
        # -- CREATE QUADRUPOLE STRENGTHS AND LENGTHS PATHS
        madx.input(f'USE, SEQUENCE="{sequence_name}";')     # Select the lattice to use

        # If we have errors or corrections, apply them
        if errors_path != None:
            madx.input(f'CALL, FILE="{errors_path}";')

        madx.input(f'SELECT, FLAG=twiss, CLEAR;')           # Tell the twiss to clear everything
        
        # We'll select all the optics elements we want
        for optic in optics_class:
            madx.input(f'SELECT, FLAG=twiss, CLASS={optic}, COLUMN={", ".join(quad_parameters)};')

        madx.input('SELECT, FLAG=twiss, PATTERN="^IP.*";')                # We'll also select the IPs here

        # Call the twiss function
        madx.twiss(file=quadrupoles_parameters_path)




# ---------------------------
# MAIN SCRIPT
# ---------------------------

# Print welcome message
print("""
\n
C r e a t i n g   t w i s s   f i l e
                        version alpha
""")

# -- Parse command line

print("Parsing command line arguments...")

parser = argparse.ArgumentParser()
create_parser_args(parser)
parsed_args = parser.parse_args()

system = parse(parsed_args)

# Where to find the data to work with
dic_key = 'Nominal' if system == 'N' else 'Errors' if system == 'E' else 'Corrections' if system == 'C' else 'Invalid'

# ----------------------------
# Parse config file arguments
# ----------------------------

print("Parsing command line arguments...")
print(f"Reading from LatticeFiles and {dic_key}System entry")

lattice_config = None       # Configuration dict for general use
system_config = None        # Configuration dict for special case use

# Read the config file and save it 
with open("configuration.toml", "rb") as f:

    general_config = tomllib.load(f)

    lattice_config = general_config[f"LatticeFiles"]
    system_config = general_config[f"{dic_key}System"]


# -- General options: where are the lattices, what elements to get, etc...

# -- Paths to lattice files 
input_main_path = lattice_config["main_input_path"]
sequence_path = input_main_path + lattice_config["sequence_path"]   # Where the sequence is located
parameters_path = input_main_path + lattice_config["params_path"]   # Parameters such as the beam energy and so


# -- Definitions
sequence_name = lattice_config["sequence_name"]             # The name of the sequence
measure_class = lattice_config["measure_class"]             # Which elements will be used as measurement points
optics_class = lattice_config["optics_class"]               # Which elements will be used as measurement points
twiss_parameters = lattice_config["measure_parameters"]     # What parameters we want to save from the twiss file
quad_parameters = lattice_config["optics_parameters"]       # What parameters we need for the quads (optics + integrals)


# -- We will check the different systems now and create the corresponding twiss
match system:
    case "N":
        # -- Nominal twiss options
        create_measurement_twiss = system_config["measure_twiss"] 
        create_quads_data = system_config["optics_twiss"]

        # -- Paths to nominal output files
        main_out = system_config["main_output_path"]
        twiss_path = main_out + system_config["measure_path"]                # Output file for the nominal twiss
        quadrupoles_path = main_out + system_config["optics_path"]           # Output file for the quads strengths and lengths

        create_twiss(parameters_path, sequence_path, create_measurement_twiss, create_quads_data, twiss_path, quadrupoles_path)



