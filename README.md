# Simple WPA Simulator

This is a simple weak-phase approximation simulator for cryo-EM particles from an atomic structure. The code takes an atomic structure and a set of imaging parameters, and outputs an image stack in STAR file format.




## Installation 

1. I HEAVILY recommend creating a virtual environment before installing this package. This package requires poetry, which could create issues if not isollated from the rest of your system through a virtual environment. Make sure the python version is at least 3.9 and that you have pip installed. 

2. Install [Poetry](https://python-poetry.org/docs/). Poetry is also available in [PyPi](https://pypi.org/project/poetry/).

3. This package can run with both GPU and CPU.

    * To Run with CPU install with `poetry install --with cpu`
    * To Run with GPU install with `poetry install --with gpu`

4. Packages that will be installed

    * Numpy
    * Matplotlib
    * Pandas
    * MDAnalysis
    * Natsort
    * Starfile
    * Mrcfile
    * Jax (cpu or gpu)
    * Torch (cpu or gpu)


## Usage

Please check the tutorial provided in the repo for specific instructions on how to run the simulator.

After installing the package you will be able to generate stacks by running
```bash
generate_stack --config config_generator --overwrite 1
```
The default value for `overwrite` is 0, which will make the command raise an error if the directory specified by `output_path` already exists.

Note: this will only work if you have the environment where you installed the package is activated.

You can run 