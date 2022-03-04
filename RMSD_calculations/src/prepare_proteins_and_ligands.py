"""
Purpose of this script is the preparation of the proteins using LePro and the preparation of ligands using open babel
 LePro to prepare proteins and openbabel to prepare ligands
"""
# ====================================================== IMPORTS ===================================================== #
import os

from openbabel import pybel
from pymol import cmd
from tqdm import tqdm


# ===================================================== FUNCTIONS ==================================================== #
def prepare_proteins(pdb_dir: str, bin_dir: str, system='mac'):
    """
    LePro automatically adds hydrogen by considering protonation state of histidine. Also removes all crystal water,
    ions, small ligands and cofactors. While LePro was especially designed especially for subsequent use of LeDock,
    the prepared proteins can also be used by Autodock vina and smina.

    Parameters
    ----------
    pdb_dir: Directory containing an arbitrary amount of .pdb files for preparation.
    bin_dir: Binary executables.
    system: Operating system. Choose mac (default) or linux.
    """
    # ---- Check input system (LePro is only available for mac and linux) ---- #
    systems = ['mac', 'linux']
    if system not in systems:
        raise ValueError("Invalid OS. Expected one of: %s" % systems)

    # ---- Get file paths, create output directory ---- #
    pdb_dir = os.path.dirname(pdb_dir)
    bin_dir = os.path.dirname(bin_dir)

    current_dir = os.path.dirname(os.path.realpath(__file__))
    output_dir = current_dir + "/../out/prepared_proteins"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # ---- Iteratively use LePro (either for mac or linux) on pdb files---- #
    if system == 'mac':
        for pdb in tqdm(os.listdir(pdb_dir), desc="...preparing receptors...                      "):
            file = os.path.join(pdb_dir, pdb)
            if pdb.endswith(".pdb"):
                os.system("%s/lepro_mac %s" % (bin_dir, file))
                os.system("mv pro.pdb %s/%s" % (output_dir, pdb))  # LePro always saves output as 'pro.pdb'
    else:
        for pdb in tqdm(os.listdir(pdb_dir)):
            file = os.path.join(pdb_dir, pdb)
            if pdb.endswith(".pdb"):
                os.system("%s/lepro_mac %s" % (bin_dir, file))
                os.system("mv pro.pdb %s/%s" % (output_dir, pdb))  # LePro always saves output as 'pro.pdb'
    os.system("rm dock.in")

    # ---- Vina requires .pdbqt files for docking (AutoDock provides a tool for that) ---- #
    for pdb in tqdm(os.listdir(output_dir), desc="...adding charges and atom types for Vina...   "):
        if pdb.endswith(".pdb"):
            file = os.path.join(output_dir, pdb)
            os.system("%s/prepare_receptor -r %s -o %s/%sqt >/dev/null 2>&1" % (bin_dir, file, output_dir, pdb))

    '''t1 = "/Users/jonathanwahrer/Desktop/Msc/3_Semester/Drug_Design_Praktikum/waehrer_HTS_screening/RMSD_calculations/out/prepared_proteins/2jam.pdb"
    t2 = "/Users/jonathanwahrer/Desktop/Msc/3_Semester/Drug_Design_Praktikum/waehrer_HTS_screening/RMSD_calculations/out/prepared_proteins/2jam.pdbqt"
    os.system("%s/prepare_receptor -r %s -o %s" % (bin_dir, t1, t2))'''


def prepare_ligands(pdb_dir):
    """
    Uses Open Babel's Pybel to add hydrogens to ligands.

    Parameters
    ----------
    pdb_dir: Directory containing an arbitrary amount of .pdb files for ligand preparation.
    """
    # ---- Get file paths, create output directory ---- #
    pdb_dir = os.path.dirname(pdb_dir)

    current_dir = os.path.dirname(os.path.realpath(__file__))
    output_dir = current_dir + "/../out/prepared_ligands"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # ---- get ligands from chain A of pdb files ---- #
    for pdb in tqdm(os.listdir(pdb_dir), desc="...extracting ligands from proteins...         "):
        '''cmd.load(pdb)
        cmd.select("sele_ligand", "organic and chain A")
        cmd.save("%s/%s_ligand.mol2" % (pdb, output_dir), format='mol2', selection='sele_ligand')
        cmd.reinitialize()'''


def run(input_proteins, path_to_bin):
    prepare_proteins(pdb_dir=input_proteins, bin_dir=path_to_bin, system="mac")
    prepare_ligands(pdb_dir=input_proteins)
