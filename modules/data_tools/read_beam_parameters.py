"""
This script has the function to read the lattice parameters and all that's needed
in order to be able to take the .madx files (or .seq) and perform the TWISS
and track calculations 
"""

import numpy as np 
import tfs 



def read_parameters(path, 
                    initial_dict = {"PARTICLE": "ELECTRON", "RADIATE": "FALSE", "BV": "+1"}, 
                    params = ["NPART", "KBUNCH", "ENERGY", "EX", "EY"],
                    file_params = ["bunch population", "number of bunches", "beam energy", "horizontal emittance", "vertical emittance"],
                    ):
    """ This function will read the `parameters.txt` file and return a dictionary with the
    ready to use parameters to create the beam to use """

    with open(path, 'r') as file:
        
        # We read the file
        lines = file.readlines()

    # We'll parse the file
    for line in lines: 

        splitted = line.split(" ")

        # -- Parse each line 
        # We'll read the param name with a complex (because this is not easily readable) script
        if splitted[0] == "":
            # Skip the first line
            continue 
        
        # Create the param name and the value holders
        param_name = splitted[0]
        value = None 
        i = 0               # An iterable to keep track of the list
        value_flag = False  # This will keep track of wether we want the name or the value]

        # Iterate
        while True: 
            
            # update the iterable
            i += 1

            # Get the new word 
            new_word = splitted[i]
            
            # This is the loop to get the param name
            if new_word != "" and not value_flag:
                param_name += " " + new_word
                continue
            elif not value_flag: 
                value_flag = True 
                continue
            
            # This is the loop to get the value
            if value_flag and new_word != "":
                try: 
                    value = float(new_word)
                except:
                    value = new_word
                break
        
        # -- Units handling
        # Now, we'll handle units

        unit = splitted[-1]
        
        # TODO: I only write those that I will actually use

        # for emittances
        if "nm" in unit:
            value *= 1e-9
        elif "pm" in unit:
            value *= 1e-12
        # For NPART
        elif "1E" in unit:
            scale = float(unit)
            value *= scale
        # For energy
        elif "c^2" in unit: 
            metric = unit[0]
            
            if metric != "G":
                letters = ["e", "K", "M", "G", "T", "P"] # eV, KeV, MeV, GeV, TeV, PeV
                index = letters.index(metric)
                
                # Scale respectively
                scale = 10**(-(3 - index) * 3)

                value *= scale
        
        # Now, we'll check if the current read parameter is in the selected params list 
        if param_name in file_params:
            index = file_params.index(param_name)
            key = params[index]
            initial_dict[key] = value 

    return initial_dict

