"""
This script is in charge of reading the .madx sequence file 
for a given accelerator and create the necessary twiss files 
"""

import numpy as np                  # Used for random number generation
import tfs                          # I'll read the measurement twiss to generate the PTC observers 
import csv
import pandas as pd

# Packages for command input
import argparse                 # Will be used mainly to select if we're creating nominal, errors or errors+corrections systems
import tomllib                  # This is to parse the config file

# Libraries for xsuite and pa calculations
from xobjects import Method
import xtrack as xt
import xobjects as xo
from xtrack.twiss import strengths
from modules.calculateBetaIntegrals import calculate_integrals
from modules.data_tools.simulateSystem_parser import create_parser_args, parse  # To parse the code's parameters
import time 

# I will import this one to get rid of a troublesome file 
# TODO: I should make it create the whole directory too
import os



# ---------------------------
# MAIN FUNCTIONS
# ---------------------------
def simulate_system(beam_params ,sequence_path, sequence_name,
                    create_measurement_twiss, create_quads_data,
                    main_path, twiss_path, quadrupoles_path, track_path,
                    debug, errors_path = None, 
                    track_flag = False, tracking_config = None,
                    make_integrals = False, _save_tfs = False, _integrals_path = False):

    # -- Create a mad-x connection
    # Read beam parameters
    
    print("\n----| Executing XTrack |------ \n")
    print("Creating system...")
    print(f"  -> Loading sequence in {sequence_path}")

    env = xt.load(sequence_path)        # Read the sequence
    line = None                         # Create the line
    if sequence_name != False:
        # If a name sequence name was given, use it
        line = env[sequence_name]
    else: 
        # If not, take the first one
        sequence_name = list(env.lines)[0]
        line = env[sequence_name]

    print(f"  \\___ ring selected: {sequence_name}")

    # If we have errors or corrections, apply them
    if errors_path != None:
        print(f"  \\__ adding modifications in {errors_path}")
        with open(errors_path, 'r') as f:
            rows = f.readlines()
            for row in rows:
                quad, err = row.split("\t")
                line[quad.lower()].k1 = line[quad.lower()].k1 + float(err)
   

    print(f"  -> Creating reference particle")
    print(f"""  \\__ Energy = {beam_params["energy"]:.3g} eV
  \\__ Mass = {beam_params["mass"]:.3g} eV
  \\__ Charge = {beam_params["charge"]} e 
  \\__ Direction = {beam_params["dir"]}
  \\__ Radiation = {beam_params["radiate"]}""")

    # TODO: gotta add it when the direction is like the other way around ykwim
    line.particle_ref = xt.Particles(
                energy0 = beam_params["energy"],
                mass0 = beam_params["mass"],
                q0 = beam_params["dir"]
            )

    twiss_method = "6d"     # We define the original twiss method, if don't use radiation, then change it

    # -- Here we'll add wether we're radiating or not
    if beam_params["radiate"]:
        print("  -> Adding radiation")
        line.configure_radiation(model = "mean")        # We'll add radiation
        line.compensate_radiation_energy_loss()         # We'll add tapering
    else:
        print("  --> Deactivating RF cavities")
        tab = line.get_table()
        tab_cav = tab.rows[tab.element_type == "Cavity"]
        for nn in tab_cav.name:
            line[nn].voltage = 0
        twiss_method = "4d"
        

    #   ---------------------------
    #       CREATE TWISS FILES 
    #   ---------------------------
    
    # First, we save all the needed names
    element_names = line.element_names

    if create_measurement_twiss:
        
        print(f"\nCreating twiss for measuring elements...")
        # We perform the twiss
        # HACK: We only want certain twiss parameters, not all the table
        
        # Tell which elements to keep
        selected_names = [name for name in element_names if line[name].__class__.__name__ in measure_class]

        # Perform the twiss
        print(f"  -> Performing twiss calculation.")
        twiss = line.twiss(method = twiss_method)

        # Filter the columns
        elem_mask = np.isin(twiss.name, selected_names)
        filtered_twiss = twiss.rows[elem_mask].cols[twiss_parameters]

        # Create pandas data for easy write and all that
        pandas_data = filtered_twiss.to_pandas()
        pandas_data.columns = pandas_data.columns.str.upper()
        
        print(f"  -> Saving twiss file to {main_path}/{twiss_path}.parquet")
        pandas_data.to_parquet(f"{main_path}/{twiss_path}.parquet", engine="pyarrow")       # Save dataframe in this format to read faster

        # Save to tfs if needed
        if _save_tfs:
            print(f"  \\__ Saving twiss file to {main_path}/tfs/{twiss_path}.tfs")
            filtered_twiss.to_tfs(f"{main_path}/tfs/{twiss_path}.tfs")


    if create_quads_data:
        print(f"\nCreating twiss for optics elements...")
        # -- CREATE QUADRUPOLE STRENGTHS AND LENGTHS PATHS
        # Tell which elements to keep
        selected_names = [name for name in element_names if line[name].__class__.__name__ in optics_class or name.lower().startswith("ip") ]
                                                                                                            # TODO: maybe would be good if I generalized this

        # Perform the twiss
        print(f"  -> Performing twiss calculation.")
        twiss = line.twiss(method = twiss_method)

        # Filter the columns
        elem_mask = np.isin(twiss.name, selected_names)
        filtered_twiss = twiss.rows[elem_mask].cols[quad_parameters]

        # Create pandas data for easy write and all that
        pandas_data = filtered_twiss.to_pandas()
        pandas_data.columns = pandas_data.columns.str.upper()
        

        # -- Calculate the integrals 
        if make_integrals: 
            print("  -> Calculating beta integrals on optics elements")
            print("  \\__ Getting elements' physical parameters")

            # Create a dictionary where we'll store the physical parameters data
            names = list()
            strengths = list()
            lengths = list()

            for name in selected_names:
                names.append(name)
                
                # We check atributes just because the ip's do not have these k1 or length
                strengths.append(getattr(line[name], "k1", 0.0))
                lengths.append(getattr(line[name], "length", 0.0))

            # Here I will convert this to a pandas dataframe
            phys_params = pd.DataFrame({'NAME': names, 'K1': strengths, 'L': lengths})

            # Here I add both columns to each the pandas data 
            pandas_data = pandas_data.merge(
                    phys_params[['NAME', "K1", "L"]],
                    on = 'NAME', how = 'left')
            
            print("  \\__ Calculating integrals")
            # With this, I have everything needed to calculate the integrals
            integrals_data = calculate_integrals(pandas_data, beam_parameters)
            
            # We save the dataframe
            print(f"  \\__ Saving integrals data to {main_path}/{integrals_path}.parquet")
            integrals_data.to_parquet(f"{main_path}/{integrals_path}.parquet", engine="pyarrow")
            
            # Also save the tfs
            if _save_tfs:
                print(f"  \\__ Saving integrals data to {main_path}/tfs/{integrals_path}.tfs")
                tfs.write(f"{main_path}/tfs/{integrals_path}.tfs", integrals_data)
            


        print(f"  -> Saving twiss file to {main_path}/{quadrupoles_path}.parquet")
        pandas_data.to_parquet(f"{main_path}/{quadrupoles_path}.parquet", engine="pyarrow")       # Save dataframe in this format to read faster

        
        # Save to tfs if needed
        if _save_tfs:
            print(f"  \\__ Saving twiss file to {main_path}/tfs/{quadrupoles_path}.tfs")
            filtered_twiss.to_tfs(f"{main_path}/tfs/{quadrupoles_path}"+".tfs")       



    #   -------------------------------
    #       CREATE TRACKING TABLE
    #   -------------------------------

    if track_flag:
         
        print(f"\nCreating track on measuring elements...")
        print(f"  -> Preparing tracking.") 
        
        # Time trackers just because I'm obsessed with them
        tracking_start = time.time()

        # -- Create the initial conditions of the tracking
        _x = tracking_config["x"]
        _y = tracking_config["y"]
        _px = tracking_config["px"]
        _py = tracking_config["py"]

        if tracking_config["random_positions"]:
            print("  -> Created random injection positions")
            x = np.random.normal(0.0, tracking_config["x_sigma"])
            y = np.random.normal(0.0, tracking_config["y_sigma"])
        
        if tracking_config["random_velocities"]:
            print("  -> Created random injection velocities")
            px = np.random.normal(0.0, tracking_config["px_sigma"])
            py = np.random.normal(0.0, tracking_config["py_sigma"])   
        

        # -- Create the monitors we'll use (basically all of them)
        # TODO: Maybe selecting it in the arcs would be better
        
        print("  -> Selecting monitors from the measure classes")
        # We'll read the twiss file with the measurement class as it is fast and can use pandas to sort
        measure_data = pd.read_parquet(f"{main_path}/{twiss_path}.parquet")
        
        # TODO: Maybe could be nice to have another criteria, but so far I think this is what we'll be used
        monitors_names = measure_data["NAME"].to_list()

        print(f"  -> Initilizing track with context {tracking_config["context"]}")
        
        # We will select the XTrack context
        con = tracking_config["context"]
        context = None 

        if con == "CUPY":
            context = xo.ContextCupy()
        elif con == "OPENCL":
            context = xo.ContextPyopencl()
        elif con == "CPU":
            context = xo.ContextCpu()
        elif con.startswith("CPU"):
            N = int(con.split("CPU")[1])
            context = xo.ContextCpu(omp_num_threads = N)       # TODO: check this param because I forgot
        
        print(f"  \\__ Building tracker")
        line.build_tracker(_context = context)
        
        print(f"  \\__ Creating initial particle")
        particle = line.build_particles(x = _x, y = _y, px = _px, py = _py)

        print(f"  -> Running track \n")
        line.track(particle, num_turns = int(tracking_config["N"]), multi_element_monitor_at = monitors_names, with_progress = 10)
         
        tracking_end = time.time()

        print(f"---| Tracking time: {tracking_end - tracking_start:.2f} seconds. \n")

        
        # -- Here we'll save up this data to something we can read later
        print(f"  -> Formatting tracking data for saving")
        mon = line.record_multi_element_last_track          # With this we retreive the monitor's data with shape (turns, particles [1], elements)
        X = mon.get("x")
        Y = mon.get("y")
        S = measure_data["S"].to_numpy()
        NAMES = mon.obs_names                   # TODO: are the names ordered?
        

        # Now we have to convert this into a pandas dataframe
        turns_range = X.shape[0]
        elements_range = X.shape[2]
        track_data = list()             # With this list we'll create the pandas dataframe

        for element in range(elements_range):
            # Create the first dicts
            rowx = {}                   # Reinitialize them just in case
            rowy = {}

            rowx = {
                    "NAME": NAMES[element],
                    "PLANE": "x"}
            rowy = {
                    "NAME": NAMES[element],
                    "PLANE": "y"}

            for turn in range(turns_range):
                # Then we add  the turns
                rowx[f"{turn}"] = X[turn, 0, element]
                rowy[f"{turn}"] = Y[turn, 0, element]
            
            track_data.append(rowx)
            track_data.append(rowy)

        # Now we convert it to a dataframe
        trackone = pd.DataFrame(track_data)
    
        # Here we save the data
        print(f"  \\__ Saving trackone file to {main_path}/{track_path}.parquet")
        trackone.to_parquet(f"{main_path}/{track_path}.parquet", engine="pyarrow")
        
        # Save it to tfs 
        if _save_tfs:
            print(f"  \\__ Saving trackone file to {main_path}/tfs/{track_path}.tfs")
            tfs.write(f"{main_path}/tfs/{track_path}.tfs", trackone)




# ---------------------------
# MAIN SCRIPT
# ---------------------------

# Time tracking
general_start = time.time()

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


# -- Definitions quads strengths and lengths
sequence_name = lattice_config["sequence_name"]             # The name of the sequence
measure_class = lattice_config["measure_class"]             # Which elements will be used as measurement points
optics_class = lattice_config["optics_class"]               # Which elements will be used as measurement points
twiss_parameters = lattice_config["measure_parameters"]     # What parameters we want to save from the twiss file
quad_parameters = lattice_config["optics_parameters"]       # What parameters we need for the quads (optics + integrals)
debug = lattice_config["debug"]                             # Debug flag


# -- Beam parameters
beam_parameters = dict()

beam_parameters["energy"] = beam_config["energy"] * 10**(9) # Remmember to convert all units to eV
beam_parameters["mass"] = beam_config["mass"] * 10**(6)     # Remmember to convert all units to eV
beam_parameters["charge"] = beam_config["charge"]           # Charge in units of e
beam_parameters["radiate"] = beam_config["radiate"]         # TODO: add this when it's time
beam_parameters["dir"] = beam_config["direction"]           # 1 if moving along the acc, -1 if moving antialong


# -- We will check the different systems now and create the corresponding twiss
match system:

    case "N":
        # -- Nominal twiss options
        create_measurement_twiss = system_config["measure_twiss"] 
        create_quads_data = system_config["optics_twiss"]

        # -- Paths to nominal output files
        main_out = system_config["main_output_path"]
        twiss_path = system_config["measure_path"]                # Output file for the nominal twiss
        quadrupoles_path = system_config["optics_path"]           # Output file for the quads strengths and lengths
        integrals_path = system_config["integrals_path"]

        # -- Save .tfs files 
        save_tfs = system_config["save_tfs"]
        
        # TODO: I may have to add the make integrals option
        simulate_system(beam_parameters, sequence_path, sequence_name, create_measurement_twiss, create_quads_data, main_out, twiss_path, quadrupoles_path, None, debug, make_integrals=True, _save_tfs = save_tfs, _integrals_path = integrals_path)

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
        twiss_path = system_config["measure_path"]                # Output file for the errors twiss
        quadrupoles_path = system_config["optics_path"]           # Output file for the quads strengths and lengths
        track_path = system_config["track_path"]                   # Output file for the trackone

        # -- Save .tfs files 
        save_tfs = system_config["save_tfs"]

        simulate_system(beam_parameters, sequence_path, sequence_name, create_measurement_twiss, create_quads_data, main_out, twiss_path, quadrupoles_path, track_path, debug, errors_path, track_flag, tracking_config, _save_tfs = save_tfs)


# Time tracking
general_end = time.time()

print(f"""
F i n i s h e d   s y s t e m   s i m u l a t i o n
written {dic_key} system in {main_out}
Took {general_end - general_start:.2f} seconds :D
""")



