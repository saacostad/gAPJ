"""
This script contain the parsing utilities for the main simulateSystem.py script
"""

import argparse                 # Will be used mainly to select if we're creating nominal, errors or errors+corrections systems


def create_parser_args(parser):
    """ In this function we create the needed command line parsing args """
    parser.add_argument(
        "-s", "--system",
        help="""Which system to handle: nominal, errors or corrections. 
        1. Nominal (nom) will only perform the respective twiss for the nominal lattice elements. 
        2. Errors (err) and corrections (corr) will perform the twiss* and the particle tracking""",
        required=True,
        dest="system"
    )


def parse(parsed_args):
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

        print("System selected...")
    
    return system
