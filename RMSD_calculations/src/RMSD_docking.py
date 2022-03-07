"""
Purpose of this script is to perform docking using Vina, Smina and LeDock given a set of (preprocessed) .pdb files with
their corresponding ligand. This script assumes that the ligands are provided with the following naming scheme:
    - my_protein_ligand.mol2
    - my_protein_ligand.pdbqt
"""

# ====================================================== IMPORTS ===================================================== #
import os
import time

from openbabel import pybel
from pymol import cmd
from tqdm import tqdm
from vina import Vina


# ===================================================== FUNCTIONS ==================================================== #
def get_box(ligand: str, extending=6.0, software='vina'):
    """
    Gets the coordinates of the docking box in a format such that it can be used by Vina/Smina or LeDock. It uses the
    original position of the ligand and extends the coordinates by specified amount of units in Angstrom.

    Function taken and adapted from: https://github.com/AngelRuizMoreno/Jupyter_Dock/tree/main/utilities

    Parameters
    ----------
    ligand: .mol2 file containing the ligand
    extending: Specifies how much the box is extended around the original ligand coordinates
    software: ('vina', 'ledock', 'both') Sets the output-format of this function.
    """
    cmd.load(ligand, object='lig')
    ([minX, minY, minZ], [maxX, maxY, maxZ]) = cmd.get_extent('lig')

    minX = minX - float(extending)
    minY = minY - float(extending)
    minZ = minZ - float(extending)
    maxX = maxX + float(extending)
    maxY = maxY + float(extending)
    maxZ = maxZ + float(extending)

    SizeX = maxX - minX
    SizeY = maxY - minY
    SizeZ = maxZ - minZ
    CenterX = (maxX + minX) / 2
    CenterY = (maxY + minY) / 2
    CenterZ = (maxZ + minZ) / 2

    cmd.delete('all')

    if software == 'vina':
        return {'center_x': CenterX, 'center_y': CenterY, 'center_z': CenterZ}, {'size_x': SizeX, 'size_y': SizeY,
                                                                                 'size_z': SizeZ}
    elif software == 'ledock':
        return {'minX': minX, 'maxX': maxX}, {'minY': minY, 'maxY': maxY}, {'minZ': minZ, 'maxZ': maxZ}
    elif software == 'both':
        return ({'center_x': CenterX, 'center_y': CenterY, 'center_z': CenterZ},
                {'size_x': SizeX, 'size_y': SizeY, 'size_z': SizeZ}), (
                   {'minX': minX, 'maxX': maxX}, {'minY': minY, 'maxY': maxY}, {'minZ': minZ, 'maxZ': maxZ})

    else:
        print('software options must be "vina", "ledock" or "both"')


def pdbqt_to_sdf(pdbqt_file=None, output=None):
    """
    Converts pdbqt file to .sdf file using open babel.

    Function taken from: https://github.com/AngelRuizMoreno/Jupyter_Dock/tree/main/utilities

    Parameters
    ----------
    pdbqt_file:
    output: .sdf file
    """
    results = [m for m in pybel.readfile(filename=pdbqt_file, format='pdbqt')]
    out = pybel.Outputfile(filename=output, format='sdf', overwrite=True)
    for pose in results:
        pose.data.update({'Pose': pose.data['MODEL']})
        pose.data.update({'Score': pose.data['REMARK'].split()[2]})
        del pose.data['MODEL'], pose.data['REMARK'], pose.data['TORSDO']

        out.write(pose)
    out.close()


def vina_docking(ligand_dir, protein_dir, protein, docking_coordinates, results_dir):
    """
    Docks a single protein with its corresponding ligand using the python implementation of Autodock Vina.

    Parameters
    ----------
    ligand_dir: Filepath to ligands
    protein_dir: Filepath to proteins
    protein: Name of the protein
    docking_coordinates: Dictionary containing coordinates of box center and box size
    results_dir: Filepath to docked ligands
    """
    # ---- paths, directories, general variables ---- #
    ligand_file = ligand_dir + protein + "_ligand.pdbqt"
    docked_ligand_file = results_dir + "vina_output/" + protein + "_docked_ligand.pdbqt"
    docked_ligand_sdf = results_dir + "vina_output/" + protein + "_docked_ligand.sdf"
    protein_file = protein_dir + protein + ".pdbqt"
    box_center, box_size = docking_coordinates[0]

    # ---- setting up vina ---- #
    v = Vina(sf_name='Vina', seed=42)
    v.set_receptor(protein_file)
    v.set_ligand_from_file(ligand_file)
    v.compute_vina_maps(center=[box_center['center_x'], box_center['center_y'], box_center['center_z']],
                        box_size=[box_size['size_x'], box_size['size_y'], box_size['size_z']])

    # ---- docking ---- #
    v.dock(exhaustiveness=10, n_poses=10)
    v.write_poses(docked_ligand_file, n_poses=10, overwrite=True)
    pdbqt_to_sdf(pdbqt_file=docked_ligand_file, output=docked_ligand_sdf)


def dock_ligands(ligand_dir: str, protein_dir: str, results_dir: str, system='mac'):
    """
    Performs docking of all proteins in a directory with their respective ligand using Vina, Smina and LeDock.

    Parameters
    ----------
    ligand_dir: Filepath to ligands
    protein_dir: Filepath to proteins
    results_dir: Filepath to docked ligands
    system: Operating system. Choose mac (default) or linux.
    """
    # ---- Check input system (LeDock and Smina are only available for mac and linux) ---- #
    systems = ['mac', 'linux']
    if system not in systems:
        raise ValueError("Invalid OS. Expected one of: %s" % systems)

    # ---- calculate docking boxes ---- #
    protein_names = [i.split('.')[0] for i in os.listdir(protein_dir) if i.endswith(".pdb")]  # unique names
    dock_boxes = []
    for protein in tqdm(protein_names, desc="...calculating docking boxes...          "):
        ligand_file = ligand_dir + protein + "_ligand.mol2"
        dock_boxes.append(get_box(ligand=ligand_file, software='both'))  # index 0: vina format; index 1: LeDock format

    # ---- docking using Vina ---- #
    print("...docking %d proteins with Autodock Vina..." % len(protein_names))
    for i, protein_setup in enumerate(zip(protein_names, dock_boxes)):
        print("Protein: ", protein, " (%d/%d)" % (i+1, len(protein_names)))

        # paths, directories, general variables
        protein = protein_setup[0]
        dock_coords = protein_setup[1]
        vina_docking(ligand_dir, protein_dir, protein, dock_coords, results_dir)
        print("------------------------------------------------------\n")

    # ---- docking using Smina ---- #
    print("...docking %d proteins with Smina..." % len(protein_names))
    for i, protein_setup in enumerate(zip(protein_names, dock_boxes)):
        print("Protein: ", protein, " (%d/%d)" % (i+1, len(protein_names)))

        # paths, directories, general variables
        protein = protein_setup[0]
        dock_coords = protein_setup[1]
        vina_docking(ligand_dir, protein_dir, protein, dock_coords, results_dir)
        print("------------------------------------------------------\n")


# ==================================================================================================================== #
def run(ligand_dir, protein_dir, results_dir, system='mac'):
    # ---- Check input system (LeDock and Smina are only available for mac and linux) ---- #
    systems = ['mac', 'linux']
    if system not in systems:
        raise ValueError("Invalid OS. Expected one of: %s" % systems)
    dock_ligands(ligand_dir=ligand_dir, protein_dir=protein_dir, results_dir=results_dir, system=system)

    # t1 = "/Users/jonathanwahrer/Desktop/Msc/3_Semester/Drug_Design_Praktikum/waehrer_HTS_screening/RMSD_calculations/out/prepared_ligands/2jam_ligand.mol2"
    # test = get_box(t1, software="both")
