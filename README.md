# gAPJ
New code for Action and Phase Jump (APJ) method for magnetic error corrections using newer software, modularization and the general APJ algorithm.


## What it has so far
1. `NominalFiles.py` processes the `.madx` sequence of the given accelerator to get 2 outputs: `nominal_measurement_twiss.tfs`, the twiss file of the measurement points; and `quads_twiss.tfs`, the nominal twiss files of the quadrupoles and IPs.
2. `modules` contain different scripts to perform different jobs. The most important are the `ActionPhaseJump.py` and `orbit_tools.py`, which perform the avermax calculation from TBT data and the APJ variables calculation. 
3. `CalculateAPJ.py` takes the `trackone` file and performs the APJ variables calculation.


##  Consider 

1. The general APJ algorithm can so far deal with gradient quadrupole errors.
2. So far, the project needs the `.madx` file defining the sequence of the accelerator and the `parameters.txt` from where one gets the beam parameters.


## TODO 
1. Create the TBT simulator. 
