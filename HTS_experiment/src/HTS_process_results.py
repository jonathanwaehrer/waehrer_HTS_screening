"""
Purpose is to pass the results of the HTS experiment in order to extract the scores (affinities in kcal/mol) and docking
times resulting from the experiment. Also, the ligand efficiency of every ligand is calculated.
"""

# ====================================================== IMPORTS ===================================================== #
import os
import pandas as pd

from glob import glob
from pathlib import Path
from pymol import cmd
from tqdm import tqdm


# ===================================================== FUNCTIONS ==================================================== #
def count_non_hydrogen_atoms(ligand_dir):
    """
    For all .mol2 files in provided directory, counts non-hydrogen atoms.

    Parameters
    ----------
    ligand_dir: Path to .mol2 files

    Returns
    ----------
    List of atom counts (non-hydrogen) and a list of ligand names in the order they were processed.
    """
    ligand_files = sorted(glob(ligand_dir + "[!concatenated_ligands]*.mol2"))  # single .mol2 files only!
    atom_counts = []
    ligand_names = []
    for ligand_mol2 in tqdm(ligand_files, desc="...counting non-hydrogen atoms in ligands...   "):
        # -- load file -- #
        cmd.reinitialize()
        cmd.load(ligand_mol2)
        # -- count non H -- #
        cmd.select("non_hydrogen", "not e. H")
        atom_counts.append(cmd.count_atoms("non_hydrogen"))
        # -- get ligand_name -- #
        ligand_names.append(Path(ligand_mol2).stem.split('.')[0])  # removes filepath and extension

    return atom_counts, ligand_names


def get_smina_scores(input_sdf):
    """
    Extracts the score of all poses in a multi-.sdf resulting from Smina.

    Smina .sd-files contain all predicted poses in a single file. The score of a predicted ligand/pose follows right
    after the line '> <minimizedAffinity>'. This means that this pattern can be used to append the desired line
    containing the score to a list.

    Parameters
    ----------
    input_sdf: Path to .sdf.

    Returns
    ----------
    List containing the score for each pose.
    """
    scores = []
    next_line = False

    predicted_poses = open(input_sdf, 'r')
    for line in predicted_poses.readlines():
        if next_line:
            scores.append(float(line))
        if line.startswith("> <minimizedAffinity>"):
            # This means that the next iteration will contain the score => set next_line to true to get the score
            next_line = True
        else:
            next_line = False
    return scores


def parse_results(results_dir, ligand_dir):
    """
    Iterates over results directories and calls corresponding functions to extract score and calculate ligand
    efficiencies

    Parameters
    ----------
    results_dir: Directory containing docking results.
    ligand_dir: Directory containing .mol2 files of the prepared ligands (for atom count)
    """
    # ---- Get ligand atom counts and names ---- #
    counts, names = count_non_hydrogen_atoms(ligand_dir)

    # ---- Get smina scores from multi .sd-files ---- #
    smina_results = sorted(glob(results_dir + "smina_output/*/*.sdf"))
    smina_scores = []
    protein_names = []
    for smina_sdf in tqdm(smina_results, desc="...parsing Smina output (multi .sd-files)...   "):
        smina_scores.extend(get_smina_scores(smina_sdf))
        protein_names.extend([smina_sdf.split('/')[-2]] * len(names))  # for long dataframe format

    # ---- Create dataframe of scores and calculate ligand efficiency ---- #
    smina_scores_df = pd.DataFrame(
        {'name': names * len(smina_results), 'docked_protein': protein_names, 'score': smina_scores,
         'no_atoms': counts * len(smina_results)})
    smina_scores_df['ligand_efficiency'] = smina_scores_df['score'] / smina_scores_df['no_atoms']
    # Save as .tsv
    smina_data_file = os.path.join(results_dir, "smina_output", "Smina_HTS_results.tsv")
    smina_scores_df.to_csv(smina_data_file, index=False)


# ==================================================================================================================== #
def run(results_dir, ligand_dir):
    parse_results(results_dir, ligand_dir)
