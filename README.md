# Simple WPA Simulator

This is a simple weak-phase approximation simulator for cryo-EM particles from an atomic structure. The code takes an atomic structure and a set of imaging parameters, and outputs an image stack in STAR file format.




## Installation 

1. I recommend creating a virtual environment before installing this package. Make sure the python version is at least 3.9 and that you have pip installed. 

2. Install package `pip install .`. Add the flag `-e` is you want to install it in editable mode.

3. Packages that will be installed

    * Numpy
    * Matplotlib
    * Pandas
    * MDAnalysis
    * Natsort
    * Starfile
    * Mrcfile
    * Jax
    * Torch


## Usage

Please check the tutorial provided in the repo for specific instructions on how to run the simulator.

After installing the package you will be able to generate stacks by running
```bash
generate_stack --config config_generator --overwrite 1
```
The default value for `overwrite` is 0, which will make the command raise an error if the directory specified by `output_path` already exists.

Note: this will only work if you have the environment where you installed the package is activated.

You can run `generate_stack --help` to get details on the entries of the config files, and the command line arguments for this script.