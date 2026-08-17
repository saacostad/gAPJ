"""
This python script is made in order to create the TBT_from_PTC3.madx file for an arbitrary accelerator
"""

import argparse
import tfs as tfs 
import pandas as pd
import numpy as np
import os


# We start by parsing the inputs of the command like
parser = argparse.ArgumentParser()
parser.add_argument(
        "-op", "--output_path",
        help="Output path where the errors and corrections files will be written to",
        required=True,
        dest="out_path"
    )
parser.add_argument(
        "-if", "--input_file",
        help="Input file of twiss file of the given accelerator. NOTE: it has to have a keyword column",
        required=True,
        dest="input_file"
    )
parser.add_argument(
        "-ip", "--interaction_point",
        help="The interaction point IP we'll work around",
        required=True,
        dest="IP"
    )
parser.add_argument(
        "-cf", "--corrections_file",
        help="If there exists corrections, this is the file to them",
        dest="corr_file"
    )
parser.add_argument(
        "-w", "--window",
        help="Over which region [m] around the IP we'll select the quadrupoles",
        required=True,
        dest="window"
    )
parser.add_argument(
        "-re", "--random_errors",
        help="Flag that tells if we want to add random errors to the selected quadrupoles",
        action="store_true",
        dest="random_flag"
    )
parser.add_argument(
        "-rs", "--random_sigma",
        help="The average dispersion of the errors (sigma of a normal)",
        dest="random_sigma"
        )
parser.add_argument(
        "-bt", "--beta_test",
        help="Beta-Beating testing: resets the corrections but not the errors",
        dest="beta_test",
        action="store_true"
        )
parser.add_argument(
        "-of", "--output_format",
        help = "Select wheter we want this output to be in MAD-X or XTrack format",
        dest = "format",
        choices=["MADX","XTrack"]
        )

args = parser.parse_args()


def get_elements_around_ip(df, ip_name, window):
    """ This function will get the quadrupoles around a given IP """

    circumference = df["S"].iloc[-1]    # First, check the position of the last element 
    ip_s = df.loc[df["NAME"] == ip_name, "S"].values[0]     # Check the s position of the given IP 
    
    # We create a new column on our dataframe that accounts for the lattice centered at our IP
    # Basically maps all s around [-C/2, C/2]
    df["S_SHIFTED"] = (df["S"] - ip_s + circumference / 2) % circumference - circumference / 2
    
    # We select only those quadrupoles around the IP
    mask =(df["S_SHIFTED"].abs() <= window)

    return df[mask].sort_values("S_SHIFTED")

def get_qp(df):
    """ This function will get the quadrupoles around a given IP """

    return df.sort_values("S")


# We will read the lines of our .seq file
input_file_path = args.input_file

# We read the twiss file using tfs 
# We check which file we're passing, if a .parquet or a .tfs
if input_file_path.endswith(".parquet"):
    df = pd.read_parquet(input_file_path)
elif input_file_path.endswith(".tfs"):
    df = tfs.read(input_file_path)
else:
    print("Invalid input format. It has to be either a `.parquet` or a `.tfs`")


# We'll select only the quadrupoles and we'll check for the IPs 
filtered_df = get_elements_around_ip(df, args.IP, float(args.window))


errors = args.out_path + "/IR_errors.madx"
corrs = args.out_path + "/IR_errors+corrections.madx"


if args.beta_test: 
    # If we'll check the beta-beating, then we will just copy the contents or the errors file to the errors+corrections 
    e_f = open(errors,'r')
    c_f = open(corrs, 'w')

    c_f.write(e_f.read())

else:
    # Now, we'll write the errors and corrections files 
    e_f = open(errors, 'w')
    c_f = open(corrs, 'w')

    if args.random_flag:
        for _, QP in filtered_df.iterrows():

            name = QP["NAME"] # We get the name of the QP  
            
            if "IP" in name.upper():
                continue

            # We'll write the errors for each of the quadrupoles 
            err = np.random.normal(0.0, float(args.random_sigma))
            sign = "+" if err > 0.0 else "-"

            if args.format == "MADX":
                print(f"{name}->K1 = {name}->K1{sign}{abs(err)};", file = e_f)
                print(f"{name}->K1 = {name}->K1{sign}{abs(err)};", file = c_f)
            elif args.format == "XTrack":
                print(f"{name}\t{err}", file = e_f)
                print(f"{name}\t{err}", file = c_f)
            else: 
                print("No valid format. Try `--output_format MADX/XTrack`")
    else:
        for _, QP in filtered_df.iterrows():
            
            name = QP["NAME"] # We get the name of the QP  

            if "IP" in name.upper():
                continue

            if args.format == "MADX":
                # We'll write the errors for each of the quadrupoles 
                print(f"{name}->K1 = {name}->K1+0.0;", file = e_f)
                print(f"{name}->K1 = {name}->K1+0.0;", file = c_f)
            elif args.format == "XTrack":
                print(f"{name}\t0.0", file = e_f)
                print(f"{name}\t0.0", file = c_f)
            else: 
                print("No valid format. Try `--output_format MADX/XTrack`")

e_f.close()
c_f.close()


