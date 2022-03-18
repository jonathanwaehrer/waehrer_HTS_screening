"""
Purpose of this script is to randomly sample a specified amount of ligands from a directory containing .mol2 files. The
sampled files will be copied from the initial directory to a target directory. Another option would have been to return
a list. However, since Autodock Vina requires the conversion of .mol2 to .pdbt and Smina requires a single concatenated
multi-.mol2 file, the less messy option is to have a specific directory containing only the sampled ligands.
"""
# ====================================================== IMPORTS ===================================================== #
import os
import random
import shutil

from glob import glob
from tqdm import tqdm


# ===================================================== FUNCTIONS ==================================================== #
def get_mol2_sample(initial_dir, sample_size, seed, target_dir):
    """
    Function to randomly copy a specified amount of .mol2 files to a given directory.

    Parameters
    ----------
    initial_dir: Path to the original ligand library.
    sample_size: Integer specifying sample size.
    seed: (int) Random seed.
    target_dir: Directory to copy the sample in.
    """
    # ---- set up directory ---- #
    if not os.path.exists(target_dir):
        os.mkdir(target_dir)

    # ---- Remove existing files in target directory (to avoid potential problems in the future) ----#
    print("...cleaning sample target directory...")
    if len(glob(target_dir+"*")) != 0:
        os.system("rm %s*" % target_dir)

    # ---- Get file paths of all ligands and copy sample ----
    full_ligand_library = glob(initial_dir + "*.mol2")
    print(initial_dir)
    print("...sampling %d out of %d ligands (random seed: %d)..." % (sample_size, len(full_ligand_library), seed))
    random.seed(seed)
    ligand_sample = random.sample(full_ligand_library, sample_size)
    for file in tqdm(ligand_sample, desc="...Sampling...                                 "):
        shutil.copy(file, target_dir)


# ==================================================================================================================== #
def run(ligand_library, sample_directory, sample_size, seed):
    get_mol2_sample(initial_dir=ligand_library, target_dir=sample_directory, sample_size=sample_size, seed=seed)
