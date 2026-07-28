# gAPJ
New code for Action and Phase Jump (APJ) method for magnetic error corrections using newer software, modularization and the general APJ algorithm.

## What it has so far
* `configuration.toml` is the general config file, from where one is able to tweak the parameters of the simulation. Documentation is in the comments.
* `simulateSystem.py` is in charge of creating the twiss and track files of a system.
* `calculateAPJ.py` will calculate the APJ variables $J$ and $\delta$.
* `modules.py/` contain auxiliar scripts. A useful one here that's meant to be called is `createErrors.py`, which will create `IR_errors.madx` and `IR_errors+corrections.madx` files.

##  Consider 
* The general APJ algorithm can so far deal with gradient quadrupole errors, but the avermax algorithms can deal with dipole and skew errors components. 
* So far, the project needs the `.madx` file defining the __sequence__ of the accelerator and the `parameters.txt` from where one gets the beam parameters.

## TODO 
1. Merge the corrections script so it works well with the project. 
2. Add measurement noise to the BPMs.
3. Add filters, specially on the APJ plots (and calculations).
4. Create an overall orchestrator to merge all scripts. 
5. Generalize the beam parameters entries.
6. Check if using Xtrack would improve runtimes.
7. Make iteratibility possible (whatever it is spelled).

# Run the project

Every main script will take a command like parameter `system`, which will tell the script which part of the pipeline it is supossed to execute. There are 3 main steps:

1. __Nominal system:__ [`simulateSystem.py -s nom/nominal`] to create nominal twiss files without having added errors
2. __Errors system:__ [`simulateSystem.py -s err/errors`] to create the twiss and tracking files for the system having added errors. [`calculateAPJ.py -s err/errors`] to calculate the APJ variables using the tracking file. [`createCorrections.py -s err/errors`] to create the corrections using the APJ general method.
3. __Corrections system:__ the same pipeline as *errors* but using `-s corr/corrections` instead. 

## Before running the project

* Make sure you have the sequence of the accelerator in `.madx` format. 
* So far, the code will read a `parameters.txt` file to get the beam it will use: particle, energy, emmisitivity....
* Create a custom `errors.madx` file. You can use the `modules/createErrors.py` script for this.


# Extras

## Use of errorsCreator.py

This script creates two files, `IR_errors.madx` and `IR_errors+corrections.madx`, which will have errors on the quadrupoles around an interaction region. The general use of this script is

`python createErrors.py -if [path to a twiss file with the quadrupoles] -ip [name of the IP] -w [how many meters around the IP to place errors] -re [if activated, will create errors] -rs [errors sigma (normal distribution)] -op [output path]`
