"""
This script will call the package's pipeline for easy setup
"""

import argparse
import sys, subprocess

# -- Create the general parser 
parser = argparse.ArgumentParser()


# Create all the parsing
parser.add_argument(
        "-m", "--mode",
        help = """Which mode to select:
        -> all: recalls everything, including new errors generation.
        -> from_errs: performs the simulation from the errors simulation, including errors generation
        -> from_errs_we: save as above but without errors generation
        -> from_corrs: performs the simulation from the corrections simulation, including corrections generation
        -> from_corrs_wc: save as above but without corrections generation""",
        choices = ["all", "from_errs", "from_errs_we", "from_corrs", "from_corrs_wc"],
        default="from_errs",
        dest = "mode"
    )


# TODO: maybe add stuff to select which config file to use

args = parser.parse_args()


# It's not needed to use match because I can reuse a lot of things
mode = args.mode
if mode == "all":

    # I'd have to run the nominal system sim
    subprocess.run([sys.executable, "simulateSystem.py", "-s", "nom"],
                   check = True, 
                   text = True)


if mode in ["all", "from_errs"]:

    # I'll run the errors generator
    # TODO: make it possible to change the errors from here or something

    # I'd have to run the nominal system sim
    subprocess.run([sys.executable, 
                    "modules/CreateErrors.py", 
                    "-op", "inputs",
                    "-if", "outputs/nominal/nominal_quads_twiss.parquet",
                    "-ip", "ipd",
                    "-w", "300",
                    "-of", "XTrack",
                    "-re", "-rs", "5e-6"],
                   check = True, 
                   text = True)
       

# Run the simulation of the errors
if mode in ["all", "from_errs", "from_errs_we"]:

    # Run the simulation
    subprocess.run([sys.executable, "simulateSystem.py", "-s", "err"],
                   check = True, 
                   text = True)
    # Run the APJ calculation
    subprocess.run([sys.executable, "calculateAPJ.py", "-s", "err"],
                   check = True, 
                   text = True)


# Run the corrections finder
if mode in ["all", "from_errs", "from_errs_we", "from_corrs"]:

    # Run the simulation
    subprocess.run([sys.executable, "calculateCorrections.py"],
                   check = True, 
                   text = True)


# Run the corrections finder
if mode in ["all", "from_errs", "from_errs_we", "from_corrs", "from_corrs_wc"]:
    # Run the simulation
    subprocess.run([sys.executable, "simulateSystem.py", "-s", "corr"],
                   check = True, 
                   text = True)
    # Run the APJ calculation
    subprocess.run([sys.executable, "calculateAPJ.py", "-s", "corr"],
                   check = True, 
                   text = True)

subprocess.run([sys.executable, "modules/calculateBetaBeating.py"],
               check = True, 
               text = True)

# TODO: add beta-beating calculation
