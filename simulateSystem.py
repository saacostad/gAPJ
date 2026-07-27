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

# Packages for madx handling
from cpymad.madx import Madx
from modules.data_tools.read_beam_parameters import read_parameters  # To parse the parameters .txt
from modules.data_tools.simulateSystem_parser import create_parser_args, parse  # To parse the code's parameters

"""
VARIABLES AND PARAMETERS FOR THE SYSTEM
Here we'll define all the possible variables and parameters for our systems before simulating.
Some variables will be accesed via a command line parser, other will be accessed from a configuration file 
"""

system = "N"        # Default system [N: nominal, E: errors, C: corrections]


"""
PARSING THE COMMAND LINE ENTRIES
This will parse the following options:
    1. System: [nominal, errors, corrections] ->    select which system to handle. Nominal will only make the respective twiss, 
                                                    errors and corrections will do the particle tracking
"""

parser = argparse.ArgumentParser()
create_parser_args(parser)
parsed_args = parser.parse_args()

system = parse(parsed_args)


