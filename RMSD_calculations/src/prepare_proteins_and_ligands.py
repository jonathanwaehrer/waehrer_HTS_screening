"""
Purpose of this script is the preparation of the proteins using LePro and the preparation of ligands using open babel
 LePro to prepare proteins and openbabel to prepare ligands
"""
# ====================================================== IMPORTS ===================================================== #
import os
from openbabel import pybel


# ===================================================== FUNCTIONS ==================================================== #
def prepare_proteins(pdb_dir: str, lepro_dir: str, system='mac'):
    """
    LePro automatically

    Parameters
    ----------
    pdb_dir: Directory containing an arbitrary amount of .pdb files for preparation
    lepro_dir: Path to lepro binary executable
    system: Operating system. Choose mac (default) or linux.
    """
    # ---- Check input system (LePro is only available for mac and linux) ---- #
    systems = ['mac', 'linux']
    if system not in systems:
        raise ValueError("Invalid OS. Expected one of: %s" % systems)

    # ---- Get file paths, create output directory ---- #
    pdb_dir = os.path.dirname(pdb_dir)
    lepro_dir = os.path.dirname(lepro_dir)

    current_dir = os.path.dirname(os.path.realpath(__file__))
    output_dir = current_dir + "/../out/prepared_proteins"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # ---- Iteratively use LePro (either for mac or linux) on pdb files---- #
    if system == 'mac':
        for file in os.listdir(pdb_dir):
            if file.endswith(".pdb"):
                os.system("%s/lepro_mac %s" % (lepro_dir, file))
                os.system("mv pro.pdb %s/%s" % (output_dir, file))  # LePro always saves output as 'pro.pdb'
    else:
        for file in os.listdir(pdb_dir):
            if file.endswith(".pdb"):
                os.system("%s/lepro_mac %s" % (lepro_dir, file))
                os.system("mv pro.pdb %s/%s" % (output_dir, file))  # LePro always saves output as 'pro.pdb'
    os.system("rm dock.in")


def prepare_ligands():
    print()


def run(input_proteins, path_to_bin):
    prepare_proteins(pdb_dir=input_proteins, lepro_dir=path_to_bin, system="mac")
    prepare_ligands()
