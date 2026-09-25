"""
This script contain the parsing utilities for the main simulateSystem.py script
"""

import argparse                 # Will be used mainly to select if we're creating nominal, errors or errors+corrections systems


def create_parser(parser):
    """ In this function we create the needed command line parsing args """

    # --------------------------------
    #  General command line parsing
    # --------------------------------

    # -- Choose which system to simulate
    parser.add_argument(
        "-s", "--system",
        help="""Which system to handle: nominal, errors or corrections. 
        1. Nominal (nom) will only perform the respective twiss for the nominal lattice elements. 
        2. Errors (err) and corrections (corr) will perform the twiss* and the particle tracking""",
        dest="system",
        default = "Invalid"
    )

    # -- Choose the configuration file
    parser.add_argument(
        "-cf", "--configuration_file",
        help="Path to the configuration file to use",
        dest = "config_file",
        default = "configuration.toml"
    )

    # -- Choose if to ignore the parsing
    parser.add_argument(
            "-C", "--use_config_file",
            help = "Ignore the command line entries and use only what's on then config files",
            required = False,
            dest = "use_config",
            action = 'store_true'
    )





    # --------------------------------
    #  Specific command line parsing
    # --------------------------------
    
    # ----------------------
    # For system simulation

    # -- Choose the input sequence file
    parser.add_argument(
            "-id", "--input_dir",
            help = "Path to where the input sequence file is. Can be in any format supported by XTrack",
            dest = "sequence_path",
    )

    # -- Choose the sequence to use. 
    parser.add_argument(
            "-sn", "--sequence_name",
            help = "Name of the sequence. If not set, it will use the first sequence found in the input file.",
            dest = "sequence_name",
    )

    # -- Choose the output directory
    parser.add_argument(
            "-od", "--output_dir",
            help = "Output directory",
            dest = "main_output_path",
    )

    # -- Choose if to save the .tfs files or not
    parser.add_argument(
            "-tfs", "--save_tfs",
            help = "Save tfs files if flag is set",
            dest = "save_tfs",
            action = 'store_true'
    )

    # -- Choose the beam to work with (it's direction)
    parser.add_argument(
            "-b", "--beam_direction",
            help = "Direction of the beam",
            dest = "dir",
    )

    # -- Errors modifications path
    parser.add_argument(
            "-ef", "--errors_file",
            help = "Path to the errors/corrections files",
            dest = "modifications_path",
    )
    
    # -- Number of turns to simulate
    parser.add_argument(
        "-n", "--turns",
        help = "Number of turns to simulate",
        dest = "turns"
    )


    # ---------------------
    # For APJ calculation

    # -- Reference for avermax BPM
    parser.add_argument(
        "-ab", "--avermax_bpm",
        help = "Reference BPM for avermax calculation",
        dest = "ref_bpm"
    )

    # -- Track file
    parser.add_argument(
        "-sf", "--tracking_file",
        help = "Path to the tracking file to use (preferably .parquet)",
        dest = "trackone_path"
    )

    # -- Model dir
    parser.add_argument(
        "-md", "--model_dir",
        help = "Path of the model twiss files to use",
        dest = "twiss_path"
    )

    # -- Arcs definitions
    parser.add_argument(
        "-las", "--left_arc_start",
        help = "Left arc start position s for APJ calculation",
        dest = "left_arc_start"
    )
    parser.add_argument(
        "-lae", "--left_arc_end",
        help = "Left arc end position s for APJ calculation",
        dest = "left_arc_end"
    )
    parser.add_argument(
        "-ras", "--right_arc_start",
        help = "Right arc start position s for APJ calculation",
        dest = "right_arc_start"
    )
    parser.add_argument(
        "-rae", "--right_arc_end",
        help = "Right arc end position s for APJ calculation",
        dest = "right_arc_end"
    )


def parse_system(parsed_args, system = "Invalid"):
    """ This function takes the parsed args and checks that the entries are allright """

    # -- We first parse the system
    # We check the entriues are correct
    parse_system = parsed_args.system.lower()

    possible_systems = ["nominal", "errors", "corrections", "nom", "err", "corr"]

    if parse_system not in possible_systems:
        print("Errors parsing the system. Possible entries are: \n1. Nominal (nom) \n2. Errors (err) \n3. Corrections (corr)")
    else: 
        if parse_system in ["nominal", "nom"]:
            system = "N"
        elif parse_system in ["errors", "err"]:
            system = "E"
        elif parse_system in ["corrections", "corr"]:
            system = "C"

    return system


def parse_commands(parsed_args):
    """ Here we'll parse the command line arguments to overwrite the config file for quick modifications """

