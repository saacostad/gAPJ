""" 
This script performs the TBT simulation given errors or errors+correction files.
It also performs the twiss if it is given as a parameter
"""

# We always use these packages
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt

# Packages for command input
import argparse                 # Will be used mainly to select if we're creating nominal, errors or errors+corrections systems
import tomllib                  # This is to parse the config file

# Packages for madx handling
from cpymad.madx import Madx
from modules.data_tools.read_beam_parameters import read_beam_parameters  # To parse the parameters .txt
from modules.data_tools.simulateSystem_parser import create_parser_args, parse  # To parse the code's parameters

# Print welcome message
print("""
\t\t System  Simulator
\t\t   Version alpha
\n
""")

"""
VARIABLES AND PARAMETERS FOR THE SYSTEM
Here we'll define all the possible variables and parameters for our systems before simulating.
Some variables will be accesed via a command line parser, other will be accessed from a configuration file 
"""



"""
PARSING THE COMMAND LINE ENTRIES
This will parse the following options:
    1. System: [nominal, errors, corrections] ->    select which system to handle. Nominal will only make the respective twiss, 
                                                    errors and corrections will do the particle tracking
"""

print("Parsing command line arguments...")

parser = argparse.ArgumentParser()
create_parser_args(parser)
parsed_args = parser.parse_args()

system = parse(parsed_args)


print(
    f"System: {'nominal' if system == 'N' else 'errors' if system == 'E' else 'corrections' if system == 'C' else 'invalid'}"
)

