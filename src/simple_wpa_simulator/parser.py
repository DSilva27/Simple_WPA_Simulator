import MDAnalysis as mda
import numpy as np
import jax.numpy as jnp
from jaxtyping import Array

from simple_wpa_simulator.gemmi_utils import read_gemmi_atoms, extract_atomic_parameter


def pdb_parser(input_file: str) -> Array:
    """
    Parses a pdb file and returns an atomic model of the protein. The atomic model is a 5xN array, where N is the number of atoms or residues in the protein. The first three rows are the x, y, z coordinates of the atoms or residues. The fourth row is the atomic number of the atoms or the density of the residues. The fifth row is the variance of the atoms or residues, which is the resolution of the cryo-EM map divided by pi squared.

    Parameters
    ----------
    input_file : str
        The path to the pdb file.

    Returns
    -------
    struct_info : Array
        The atomic model of the protein.
    """

    atoms = read_gemmi_atoms(input_file)
    ff_a = np.array(extract_atomic_parameter(atoms, "form_factor_a"))
    ff_b = np.array(extract_atomic_parameter(atoms, "form_factor_b"))
    b_factor = np.array(extract_atomic_parameter(atoms, "B_factor"))
    struct_info = {
        "ff_a": ff_a,
        "ff_b": ff_b,
        "b_factor": b_factor,
    }

    return struct_info
