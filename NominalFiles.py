"""
This script is in charge of reading the .madx sequence file 
for a given accelerator and create the necessary twiss files 
that will be used as the nominal parameters
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Libraries for madx
from cpymad.madx import Madx
from modules.data_tools.read_beam_parameters import read_parameters  # To parse the parameters .txt

# ---------------------------
# INPUT FILES
# ---------------------------

# -- Options
create_measurement_twiss = True 
create_quads_data = True

# -- Paths to input files 
input_main_path = "inputs/lattices/higgs/"
sequence_path = input_main_path + "fccee_h.madx"               # Where the sequence is located
parameters_path = input_main_path + "h_parameter_table.txt"    # Parameters such as the beam energy and so

# -- Paths to output files
twiss_path = "outputs/nominal_measurements_twiss.tfs"                                    # Output file for the nominal twiss
quadrupoles_parameters_path = "outputs/nominal_quads_twiss.tfs"           # Output file for the quads strengths and lengths

# -- Definitions
sequence_name = "fccee_p_ring"                                              # The name of the sequence
measure_class = "marker"                                    # Which elements will be used as measurement points
twiss_parameters = ["name", "s", "betx", "bety", "mux", "muy"]   # What parameters we want to save from the twiss file
quad_parameters = ["name", "s", "betx", "bety", "mux", "muy", "K1L"]                   # What parameters we need for the quads (optics + integrals)


# ---------------------------
# MAIN SCRIPT
# ---------------------------

# -- Create a mad-x connection

# Read beam parameters
beam_params = read_parameters(parameters_path)

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
    madx.input(f'SELECT, FLAG=twiss, CLEAR;')           # Tell the twiss to clear everything
    madx.input(f'SELECT, FLAG=twiss, CLASS={measure_class}, COLUMN={", ".join(twiss_parameters)};')

    # Call the twiss function
    twiss = madx.twiss(file=twiss_path)


if create_quads_data:
    # -- CREATE QUADRUPOLE STRENGTHS AND LENGTHS PATHS
    madx.input(f'USE, SEQUENCE="{sequence_name}";')     # Select the lattice to use
    madx.input(f'SELECT, FLAG=twiss, CLEAR;')           # Tell the twiss to clear everything
    madx.input(f'SELECT, FLAG=twiss, CLASS=QUADRUPOLE, COLUMN={", ".join(quad_parameters)};')
    madx.input('SELECT, FLAG=twiss, PATTERN="^IP.*";')                # We'll also select the IPs here

    # Call the twiss function
    quads = madx.twiss(file=quadrupoles_parameters_path)


