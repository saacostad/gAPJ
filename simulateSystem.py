"""
This script is in charge of reading the .madx sequence file 
for a given accelerator and create the necessary twiss files 
"""

import numpy as np                  # Used for random number generation
import tfs                          # I'll read the measurement twiss to generate the PTC observers 

# Packages for command input
import argparse                 # Will be used mainly to select if we're creating nominal, errors or errors+corrections systems
import tomllib                  # This is to parse the config file

# Libraries for madx
from cpymad.madx import Madx
from modules.data_tools.calc_integrals import calculate_integral
from modules.data_tools.read_beam_parameters import read_beam_parameters  # To parse the parameters .txt
from modules.data_tools.simulateSystem_parser import create_parser_args, parse  # To parse the code's parameters

# I will import this one to get rid of a troublesome file 
import os



# ---------------------------
# MAIN FUNCTIONS
# ---------------------------
def simulate_system(parameters_path, sequence_path, 
                    create_measurement_twiss, create_quads_data,
                    twiss_path, quadrupoles_parameters_path, track_path,
                    debug, errors_path = None, 
                    track_flag = False, tracking_config = None,
                    make_integrals = False):

    # -- Create a mad-x connection
    # Read beam parameters
    print("Parsing beam parameters...")
    beam_params = read_beam_parameters(parameters_path)
    
    print("\n\t -- Executing MAD-X -- ")
    # Create madx connection
    madx = Madx(stdout=False, stderr=True) 
    # Let mad-X log on terminal 
    debug_flag = "" if debug else "-"
    madx.input(f"OPTION, {debug_flag}ECHO, {debug_flag}INFO, {debug_flag}WARN, {debug_flag}DEBUG, {debug_flag}TWISS_PRINT;")
    
    # We define the beam to use 
    beam_input = "BEAM, " + ", ".join(f"{key}={value}" for key, value in beam_params.items())
    print(f"\nSelected beam:\n{beam_input.lower()}")
    madx.input(beam_input)


    if create_quads_data or create_measurement_twiss:
        # We call the sequence
        print(f"\nCalling sequence {sequence_name} in {sequence_path}")
        madx.input(f'CALL, FILE="{sequence_path}";')
    else: 
        print("No option was given.")
        return
    

    madx.input(f'USE, SEQUENCE="{sequence_name}";')     # Select the lattice to use

    # If we have errors or corrections, apply them
    if errors_path != None:
        print(f"Adding modifications in {errors_path}")
        madx.input(f'CALL, FILE="{errors_path}";')


    #   ---------------------------
    #       CREATE TWISS FILES 
    #   ---------------------------

    if create_measurement_twiss:
        
        print(f"\nCreating twiss for measuring elements...")
        # We perform the twiss
        # HACK: We only want certain twiss parameters, not all the table

        madx.input(f'SELECT, FLAG=twiss, CLEAR;')           # Tell the twiss to clear everything
    
        # We will select all the measuring elements for our twiss
        for measure in measure_class:
            madx.input(f'SELECT, FLAG=twiss, CLASS={measure}, COLUMN={", ".join(twiss_parameters)};')

        # Call the twiss function
        madx.twiss(file=twiss_path)


    if create_quads_data:
        print(f"\nCreating twiss for optics elements...")
        # -- CREATE QUADRUPOLE STRENGTHS AND LENGTHS PATHS

        madx.input(f'SELECT, FLAG=twiss, CLEAR;')           # Tell the twiss to clear everything
        
        # We'll select all the optics elements we want
        for optic in optics_class:
            madx.input(f'SELECT, FLAG=twiss, CLASS={optic}, COLUMN={", ".join(quad_parameters)};')

        madx.input('SELECT, FLAG=twiss, PATTERN="^IP.*";')                # We'll also select the IPs here

        # Call the twiss function
        madx.twiss(file=quadrupoles_parameters_path)



        # -- Calculate the integrals 
        if make_integrals: 
            
            # -- Sign calculation: here we'll check which sign the integrals will have (beam direction * particle's charge)
            particle = beam_params['PARTICLE']
            direction = beam_params['BV']

            # This function will calculate the integrals to not make this script too messy
            calculate_integral(madx, 1, optics_class)
            # TODO: add the sign to actually make this work



            





    #   -------------------------------
    #       CREATE TRACKING TABLE
    #   -------------------------------

    if track_flag:
        
        print(f"\nCreating track on measuring elements...")
        print(f"-> Preparing tracking.") 

        x = tracking_config["x"]
        y = tracking_config["y"]
        px = tracking_config["px"]
        py = tracking_config["py"]

        if tracking_config["random_positions"]:
            print("-> Created random injection positions")
            x = np.random.normal(0.0, tracking_config["x_sigma"])
            y = np.random.normal(0.0, tracking_config["y_sigma"])
        
        if tracking_config["random_velocities"]:
            print("-> Created random injection velocities")
            px = np.random.normal(0.0, tracking_config["px_sigma"])
            py = np.random.normal(0.0, tracking_config["py_sigma"])   

        madx.input(f'SELECT, FLAG=track, CLEAR;')        # Tell the tracker to clear everything
    
        # We will select all the measuring elements for our twiss
        for measure in measure_class:
            madx.input(f'SELECT, FLAG=track, CLASS={measure}')
        
        print(f"-> Performing tracking for {tracking_config["N"]} turns...")

        # -- Tracking

        # Set debug level
        madx.input(f'PTC_SETSWITCH, DEBUGLEVEL={2 if debug else 0}, MAPDUMP={1 if debug else 0};')
        madx.input(f'PTC_CREATE_UNIVERSE, SYMPRINT={not debug};')

        madx.input(f'PTC_CREATE_LAYOUT, MODEL={tracking_config["model"]}, METHOD={tracking_config["method"]}, NST={tracking_config["nst"]};')
        madx.input(f'PTC_START, X={x}, PX={px}, Y={y}, PY={py};')

        # - Place the PTC observers 
        measure_elements = tfs.read(twiss_path, index="NAME").index.tolist()

        for element in measure_elements:
            madx.input(f'PTC_OBSERVE, PLACE={element};')

        madx.input(f'PTC_TRACK, TURNS={tracking_config["N"]}, DUMP=true, ONETABLE=true, FILE="{track_path}", ELEMENT_BY_ELEMENT=true;')
        madx.input(f'PTC_TRACK_END;')

        print("")

        # Get rid of that ugly file 
        os.system("rm internal_mag_pot.txt")

        # Call the twiss function
        madx.twiss(file=twiss_path)


# ---------------------------
# MAIN SCRIPT
# ---------------------------

# Print welcome message
print("""
\n
C r e a t i n g   s y s t e m   s i m u l a t i o n 
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

print("Parsing config file arguments...")

lattice_config = None       # Configuration dict for general use
system_config = None        # Configuration dict for special case use
tracking_config = None      # Configuration for tracking
beam_config = None          # Configuration for beam

# Read the config file and save it 
with open("configuration.toml", "rb") as f:

    general_config = tomllib.load(f)

    lattice_config = general_config[f"LatticeFiles"]
    system_config = general_config[f"{dic_key}System"]
    beam_config = general_config[f"BeamParameters"]
    tracking_config = general_config[f"TrackParameters"]

print(f"Read from LatticeFiles and {dic_key}System entry")

# -- General options: where are the lattices, what elements to get, etc...

# -- Paths to lattice files 
input_main_path = lattice_config["main_input_path"]
sequence_path = input_main_path + lattice_config["sequence_path"]   # Where the sequence is located
parameters_path = input_main_path + lattice_config["params_path"]   # Parameters such as the beam energy and so


# -- Definitionsquads strengths and lengths
sequence_name = lattice_config["sequence_name"]             # The name of the sequence
measure_class = lattice_config["measure_class"]             # Which elements will be used as measurement points
optics_class = lattice_config["optics_class"]               # Which elements will be used as measurement points
twiss_parameters = lattice_config["measure_parameters"]     # What parameters we want to save from the twiss file
quad_parameters = lattice_config["optics_parameters"]       # What parameters we need for the quads (optics + integrals)
debug = lattice_config["debug"]                             # Debug flag

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

        simulate_system(parameters_path, sequence_path, create_measurement_twiss, create_quads_data, twiss_path, quadrupoles_path, None, debug, make_integrals=True)

    # In general, the errors and corrections system should behave the same
    case _:
        # -- E/C twiss options
        create_measurement_twiss = system_config["measure_twiss"] 
        create_quads_data = system_config["optics_twiss"]
        track_flag = system_config["TBT_track"]

        # Where the E/C .madx file is 
        errors_path = system_config["modifications_path"]
    

        # -- Paths to E/C output files
        main_out = system_config["main_output_path"]
        twiss_path = main_out + system_config["measure_path"]                # Output file for the errors twiss
        quadrupoles_path = main_out + system_config["optics_path"]           # Output file for the quads strengths and lengths
        track_path = main_out + system_config["track_path"]                   # Output file for the trackone

        simulate_system(parameters_path, sequence_path, create_measurement_twiss, create_quads_data, twiss_path, quadrupoles_path, track_path, debug, errors_path, track_flag, tracking_config)


print(f"""
F i n i s h e d   s y s t e m   s i m u l a t i o n
written {dic_key} system in {main_out}
""")
