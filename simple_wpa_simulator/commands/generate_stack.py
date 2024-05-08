import os
import sys
import glob
import natsort
import argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
import warnings
import logging
import json

from simple_wpa_simulator.simulator import simulate_stack
from simple_wpa_simulator.parser import pdb_parser
from simple_wpa_simulator.utils import load_config, help_config_generator

warnings.filterwarnings("ignore", module="MDAnalysis")


def add_args(parser):
    parser.add_argument("--config", type=str, help="Config file name", required=True)
    parser.add_argument(
        "--overwrite",
        type=int,
        help="Overwrite previous results (1 for yes). Default is 0.",
        default=0,
        required=False,
    )

    return parser


def load_configfile(config_fname, overwrite):

    config = load_config(config_fname)

    if not os.path.exists(config["output_path"]):
        os.makedirs(config["output_path"])

    else:
        if overwrite == 0:
            raise FileExistsError(
                f"Output path {config['output_path']} already exists. Set overwrite to 1 to overwrite."
            )

    with open(os.path.join(config["output_path"], config_fname), "w") as f:
        json.dump(config, f, indent=4)

    return config


def load_models(config):

    logging.info("Checking that the models exist...")
    logging.info(f"Looking for models {config['models_fname']}")

    if "*" in config["models_fname"]:
        models_fname = natsort.natsorted(glob.glob(config["models_fname"]))
        if len(models_fname) == 0:
            raise FileNotFoundError(
                f"No files found with pattern {config['models_fname']}"
            )
    else:
        models_fname = [config["models_fname"]]

    logging.info("Models found.")

    logging.info("Using models: {}".format(config["models_fname"]))

    config["defocus_ang"] = list(np.radians(config["defocus_ang"]))
    n_models = len(models_fname)
    # unit_cell = [34., 34., 34., 90., 90., 90.]

    models = []
    print(models_fname[0])
    model_0 = mda.Universe(models_fname[0])
    model_0.atoms.translate(-model_0.atoms.center_of_mass())

    logging.info(f"Using model {models_fname[0]} as reference.")
    path_ref_model = (
        os.path.join(config["output_path"], "ref_model.")
        + models_fname[0].split(".")[-1]
    )
    model_0.atoms.write(path_ref_model)
    logging.info(f"Reference model written to {path_ref_model}")

    if config["mode"] == "resid":
        atom_list_filter = "protein and name CA"

    elif config["mode"] == "all-atom":
        atom_list_filter = "protein and not name H*"

    elif config["mode"] == "cg":
        atom_list_filter = "protein and name CA"

    for i in range(0, n_models):

        model_path = os.path.join(config["working_dir"], models_fname[i])
        uni = mda.Universe(model_path)
        align.alignto(uni, model_0, select=atom_list_filter, weights="mass")
        models.append(uni.select_atoms(atom_list_filter).positions.T)

    models = np.array(models)
    struct_info = pdb_parser(model_path, mode=config["mode"])

    return models, struct_info


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, epilog=help_config_generator())
    parser = add_args(parser)
    args = parser.parse_args()

    config = load_configfile(args.config, args.overwrite)

    logger = logging.getLogger()
    logger_fname = os.path.join(
        config["output_path"], config["experiment_name"] + ".log"
    )
    fhandler = logging.FileHandler(filename=logger_fname, mode="a")
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    fhandler.setFormatter(formatter)
    logger.addHandler(fhandler)
    logger.setLevel(logging.INFO)

    logging.info(
        "A copy of the used config file has been written to {}".format(
            os.path.join(config["output_path"], "config.json")
        )
    )

    models, struct_info = load_models(config)

    logging.info("Simulating particle stack...")
    simulate_stack(
        working_dir=config["output_path"],
        starfile_fname=config["starfile_fname"],
        models=models,
        struct_info=struct_info,
        images_per_model=config["particles_per_model"],
        config=config,
        batch_size=config["batch_size"],
        dtype=float,
        seed=config["seed"],
    )
    logging.info("Simulation complete.")
    logging.info("Output written to {}".format(config["output_path"]))

    return


if __name__ == "__main__":
    main()
