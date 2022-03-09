"""
Purpose of this script is to calculate the RMSD between the best predicted ligand positions and the original position
of the ligand in the .pdb file.
"""
# ====================================================== IMPORTS ===================================================== #
from biopandas.mol2 import PandasMol2
from glob import glob
from pymol import cmd
from tqdm import tqdm
import os


# ===================================================== FUNCTIONS ==================================================== #
def get_first_sdf_entry(sdf_file):
    """
    Extracts the first entry of an .sdf file to .mol2 using pymol

    Parameters
    ----------
    sdf_file: Path to .sdf-file containing multiple structures
    """
    mol2_file = sdf_file[:-3] + 'pdb'
    # ---- Load .sdf file, save first pose (happens automatically when saving to .mol2) ---- #
    cmd.load(sdf_file, object='predictions')
    cmd.save(mol2_file, 'predictions')
    cmd.reinitialize()

    return mol2_file


def get_RMSD(true_ligands, predicted_positions):
    """
    Iterates over .sdf files in a directory and calculates the RMSD to the original ligand position. The first entry of
    the .sdf file containing the predicted positions is exported to .mol2 using PyMOL. Coordinates of both .mol2 files
    are parsed using BioPandas.

    Parameters
    ----------
    true_ligands: Directory containing the original ligand as .mol2 file
    predicted_positions: Directory containing .sdf files of the re-docked ligands. Only the first entry is considered.
    """
    tool = os.path.basename(predicted_positions).split('_')[0]  # name of current docking tool
    # ---- Get file paths of all ligands ---- #
    true_ligands = sorted(glob(true_ligands+"*.mol2"))
    predicted_ligands = sorted(glob(predicted_positions+"/*.sdf"))
    # ---- Iterate over ligands ---- #
    "...preparing ligands for vina...         "
    "...calculating RMSD for %s...            "
    for truth, predicted in tqdm(zip(true_ligands, predicted_ligands), total=(len(predicted_ligands)),
                                 desc="...calculating RMSD for %s..." % tool, position=0, leave=True):
        # Export first position of predictions.sdf, parse coordinates
        predicted = get_first_sdf_entry(sdf_file=predicted)
        #truth_coords = PandasMol2().read_mol2(truth).df[['x', 'y', 'z']]
        #predicted_coords = PandasMol2().read_mol2(predicted).df[['x', 'y', 'z']]


    '''truth = true_ligands[0]
    print(truth)
    predicted = predicted_ligands[0]

    
    print(truth_coords)'''

    '''cmd.load(truth, object='truth')
    cmd.load(predicted, object='predicted', format='sdf')
    cmd.split_states('predicted')

    cmd.align((predicted+'_0001'), truth)'''


# ==================================================================================================================== #
def run(original_lig_dir, predicted_lig_dir, rmsd_dir):
    if not os.path.exists(rmsd_dir):
        os.mkdir(rmsd_dir)
    #print(original_lig_dir)
    #print(predicted_lig_dir)
    #print(rmsd_dir)
    # ---- iterate over subdirectories corresponding to different docking tools ---- #
    sub_dirs = sorted(glob(predicted_lig_dir+"*"))
    get_RMSD(true_ligands=original_lig_dir, predicted_positions=sub_dirs[3])
    '''for directory in sub_dirs:
        get_RMSD(true_ligands=original_lig_dir, predicted_positions=directory)
        #print(os.path.basename(directory).split('_')[0])'''
