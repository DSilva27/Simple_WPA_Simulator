# write a json config loader

import json
import os
import sys
from numbers import Number
import numpy as np
import textwrap

def load_config(config_file):
    """
    Load a json config file

    Parameters
    ----------
    config_file : str
        Path to the config file

    Returns
    -------
    config : dict
    """
    if not os.path.isfile(config_file):
        print("Config file not found: {}".format(config_file))
        sys.exit(1)

    with open(config_file, "r") as f:
        config = json.load(f)

    config = parse_config(config)

    return config


def parse_config(config):

    req_keys = {
        "experiment_name": str,
        "mode": str,
        "particles_per_model": (Number, list),
        "box_size": Number,
        "resolution": Number,
        "pixel_size": Number,
        "defocus_u": (Number, list),
        "working_dir": str,
        "output_path": str,
        "models_fname": str,
        "starfile_fname": str,
        "batch_size": Number,
    }

    for key in req_keys:
        if key not in config:
            raise ValueError("{} not found in config file".format(key))

        if not isinstance(config[key], req_keys[key]):
            raise ValueError("{} must be of type {}".format(key, req_keys[key]))

    optional_keys = {
        "defocus_v": [(Number, list), config["defocus_u"]],
        "defocus_ang": [(Number, list), 0],
        "bfactor": [(Number, list), 0.0],
        "scalefactor": [Number, 1.0],
        "phaseshift": [Number, 0.0],
        "amp_contrast": [Number, 0.01],
        "volt": [Number, 300],
        "spherical_aberr": [Number, 2.7],
        "seed": [Number, 0],
    }

    for key in optional_keys:
        if key not in config:
            config[key] = optional_keys[key][1]
        elif not isinstance(config[key], optional_keys[key][0]):
            raise ValueError("{} must be of type {}".format(key, optional_keys[key][0]))

    return config


def validate_config_generator(config):

    if config["mode"] not in ["all-atom", "resid", "cg"]:
        raise ValueError("Invalid mode, must be 'all-atom', 'resid' or 'cg'")

    if config["box_size"] <= 0:
        raise ValueError("Box size must be greater than 0")

    if config["resolution"] <= 0:
        raise ValueError("Resolution must be greater than 0")

    if config["pixel_size"] <= 0:
        raise ValueError("Pixel size must be greater than 0")

    if np.all(np.array(config["particles_per_model"]) <= 0):
        raise ValueError("Particles per model must be greater than 0")

    if np.all(np.array(config["defocus_u"]) <= 0):
        raise ValueError("Defocus u must be greater than 0")

    if np.all(np.array(config["defocus_v"]) <= 0):
        raise ValueError("Defocus v must be greater than 0")

    if config["batch_size"] <= 0:
        raise ValueError("Batch size must be greater than 0")

    if config["amp_contrast"] < 0 or config["amp_contrast"] > 1:
        raise ValueError("Amplitude contrast must be between 0 and 1")

    if config["volt"] <= 0:
        raise ValueError("Voltage must be greater than 0")

    if config["spherical_aberr"] <= 0:
        raise ValueError("Spherical aberration must be greater than 0")

    return


def help_config_generator():

    string = textwrap.dedent("""\
        Required keys in config file:
            experiment_name: name of the experiment, used for logging
            mode: detail level of the atomic models (all-atom, resid, cg)
            particles_per_model: number of particles per model as a list [n_particles_1, n_particles_2, ...]
            box_size: box size of the particles
            resolution: resolution of the particles
            pixel_size: pixel size of the particles
            defocus_u: defocus_u of the particles (Angstrom)
            working_dir: path where the atomic models (pdb files) are located
            output_path: output path for the generated data
            models_fname: path to the atomic models (pdb files). The path should be relative to working_dir, you can use * to indicate multiple files
            starfile_fname: name of the starfile for output
            batch_size: batch size of the output dataset

        Optional keys in config file:
            defocus_v: defocus_v of the particles (Angstrom) (default: defocus_u)
            defocus_ang: defocus_ang of the particles (degrees) (default: 0°)
            bfactor: bfactor of the particles (Angstrom^2) (default: 0.0 Angstrom^2)
            scalefactor: scalefactor of the particles (default: 1.0)
            phase_shift: phase_shift of the particles (default: 0.0)
            amp_contrast: amp_contrast of the particles (default: 0.01)
            volt: volt of the particles (kv) (default: 300 kV)
            spherical_aberr: spherical_aberr of the particles (mm) (default: 2.7 mm)
            seed: seed for parameter generation (default: 0)
 
    """)

    return string
