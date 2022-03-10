"""
Purpose of this script is to calculate the RMSD between the best predicted ligand positions and the original position
of the ligand in the .pdb file.
"""
# ====================================================== IMPORTS ===================================================== #
import math
import numpy as np
import pandas as pd
import os

from glob import glob
from pymol import cmd

from rdkit import Chem
from rdkit.Chem import rdFMCS, PandasTools
from rdkit import RDLogger


# ===================================================== FUNCTIONS ==================================================== #
def get_first_sdf_entry(sdf_file):
    """
    Extracts the first entry of an .sdf file to .mol2 using pymol.

    Parameters
    ----------
    sdf_file: Path to .sdf-file containing multiple structures
    """
    mol2_file = sdf_file[:-3] + 'mol2'
    # ---- Load .sdf file, save first pose (happens automatically when saving to .mol2) ---- #
    cmd.load(sdf_file, object='predictions')
    cmd.save(mol2_file, 'predictions')
    cmd.reinitialize()

    return mol2_file


def get_inplace_rmsd(ref, target):
    """
    Function to calculate the RMSD between two structures using RDKit based on their maximum common substructure (simple
    calculation of the RMSD was not possible, because the docking tools tend to change the order of atom entries in the
    output without any hint w.r.t. the original position in the file).

    Function taken and adapted from: https://github.com/AngelRuizMoreno/Jupyter_Dock/tree/main/utilities
    The original function

    Parameters
    ----------
    ref: Filepath to the reference structure as .mol2 file
    target: Filepath to target poses as .sdf file

    Returns
    ----------
    List containing RMSD of all poses in target.sdf
    """
    # ---- Read reference and calculate RMSD ---- #
    ref = Chem.MolFromMol2File(ref)
    ref.SetProp('_Name', 'Ref')

    r = rdFMCS.FindMCS([ref, target])

    a = ref.GetSubstructMatch(Chem.MolFromSmarts(r.smartsString))
    b = target.GetSubstructMatch(Chem.MolFromSmarts(r.smartsString))
    amap = list(zip(a, b))

    distances = []
    for atomA, atomB in amap:
        pos_A = ref.GetConformer().GetAtomPosition(atomA)
        pos_B = target.GetConformer().GetAtomPosition(atomB)
        coord_A = np.array((pos_A.x, pos_A.y, pos_A.z))
        coord_B = np.array((pos_B.x, pos_B.y, pos_B.z))
        dist_numpy = np.linalg.norm(coord_A - coord_B)
        distances.append(dist_numpy)

    rmsd = math.sqrt(1 / len(distances) * sum([i * i for i in distances]))

    return rmsd


def get_RMSD(true_ligands, predicted_positions):
    """
    Iterates over .sdf files in a directory and calls a function to calculate the RMSD between the predicted poses and
    the original poses of the ligand using RDKit.

    Parameters
    ----------
    true_ligands: Directory containing the original ligand as .mol2 file
    predicted_positions: Directory containing .sdf files of the re-docked ligands. Only the first entry is considered.
    """
    tool = os.path.basename(predicted_positions).split('_')[0]  # name of current docking tool
    # ---- Get file paths of all ligands ---- #
    true_ligands = sorted(glob(true_ligands + "*.mol2"))
    predicted_ligands = sorted(glob(predicted_positions + "/*.sdf"))
    # ---- Initialize lists for pandas dataframe (will be used for plotting in the analysis) ---- #
    """
    Seaborn and ggplot require long format. Desired format is as follows:
    
    RMSD    Score        Protein    Tool
    ----    -----        -------    ----
    0.4      -8.4        5yabx      vina
    0.3      -7.3        5yab       ledock
    ...     ...          ...        ...   
    
    => 3 separate lists to keep track of RMSD, current protein and current tool
    """
    rmsd = []
    scores = []
    protein_indicator = []
    tool_indicator = []

    # ---- Iterate over ligands ---- #
    print("...calculating RMSD for %s poses..." % tool)
    for truth, predicted in zip(true_ligands, predicted_ligands):
        protein_name = os.path.basename(truth).split("_")[0]
        print("   Protein: %s" % protein_name)

        # ---- Read target poses (might fail. In this case, this ligand will be skipped) ---- #
        RDLogger.DisableLog('rdApp.*')
        predicted_poses = PandasTools.LoadSDF(predicted)

        if len(predicted_poses) == 0:
            print("      Note: Skipped calculation for %s because of sanitation error." % protein_name)
            continue
        else:
            # ---- Calculate RMSD and update lists for pandas DF ---- #
            rmsd_current = [get_inplace_rmsd(ref=truth, target=predicted_poses.loc[i, 'ROMol']) for i in
                            predicted_poses.index]

            rmsd.extend(rmsd_current)
            protein_indicator.extend([protein_name] * len(rmsd_current))
            tool_indicator.extend([tool] * len(rmsd_current))

            if 'Score' in predicted_poses.columns:
                # Name of calculated affinity differs in tools. Vina and LeDock have 'Score' entry...
                scores.extend(predicted_poses['Score'])
            elif 'minimizedAffinity' in predicted_poses.columns:
                # ...Smina has 'minimizedAffinity' entry
                scores.extend(predicted_poses['minimizedAffinity'])
            else:
                # ...SeeSar has no comparable scoring output
                scores.extend(['NaN'] * len(rmsd_current))

    rmsd_df = pd.DataFrame(zip(rmsd, scores, protein_indicator, tool_indicator),
                           columns=['RMSD', 'Score', 'Protein', 'Tool'])
    return rmsd_df


# ==================================================================================================================== #
def run(original_lig_dir, predicted_lig_dir, rmsd_dir):
    """
    Run RMSD calculation pipeline.

    Parameters
    ----------
    original_lig_dir: Path to .mol2 of original ligand positions
    predicted_lig_dir: Path to .sdf files containing docked poses
    rmsd_dir: Output directory for the .tsv file containing results
    """
    if not os.path.exists(rmsd_dir):
        os.mkdir(rmsd_dir)

    # ---- iterate over subdirectories corresponding to different docking tools ---- #
    sub_dirs = sorted(glob(predicted_lig_dir + "*"))
    rmsd_data = pd.DataFrame()
    for directory in sub_dirs:
        rmsd_data = pd.concat([rmsd_data, get_RMSD(true_ligands=original_lig_dir, predicted_positions=directory)])
        print("------------------\n")
    rmsd_data.to_csv(rmsd_dir + "RMSD_experiment_results.tsv", index=False, sep='\t')
