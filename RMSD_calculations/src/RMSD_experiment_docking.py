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
    ligand: .mol2 file containing the ligand.
    extending: Specifies how much the box is extended around the original ligand coordinates.
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
    output: .sdf file.
    """
    results = [m for m in pybel.readfile(filename=pdbqt_file, format='pdbqt')]
    out = pybel.Outputfile(filename=output, format='sdf', overwrite=True)
    for pose in results:
        pose.data.update({'Pose': pose.data['MODEL']})
        pose.data.update({'Score': pose.data['REMARK'].split()[2]})
        del pose.data['MODEL'], pose.data['REMARK'], pose.data['TORSDO']

        out.write(pose)
    out.close()


def dok_to_sdf(dok_file=None, output=None):
    """
    LeDock outputs docked ligands as .dok file which is not supported by many 3D visualization software including PyMol.
    This function converts .dok fiels to .sdf files.

    Function taken from: https://github.com/AngelRuizMoreno/Jupyter_Dock/tree/main/utilities

    Parameters
    ----------
    dok_file: str or path-like ; dok file from ledock docking
    output: str or path-like ; outfile from ledock docking, extension must be sdf
    """
    out = pybel.Outputfile(filename=output, format='sdf', overwrite=True)

    with open(dok_file, 'r') as f:
        doc = [line for line in f.readlines()]

    doc = [line.replace(line.split()[2], line.split()[2].upper()) if 'ATOM' in line else line for line in doc]

    start = [index for (index, p) in enumerate(doc) if 'REMARK Cluster' in p]
    finish = [index - 1 for (index, p) in enumerate(doc) if 'REMARK Cluster' in p]
    finish.append(len(doc))

    interval = list(zip(start, finish[1:]))
    for num, i in enumerate(interval):
        block = ",".join(doc[i[0]:i[1]]).replace(',', '')

        m = pybel.readstring(format='pdb', string=block)

        m.data.update({'Pose': m.data['REMARK'].split()[4]})
        m.data.update({'Score': m.data['REMARK'].split()[6]})
        del m.data['REMARK']

        out.write(m)

    out.close()


def generate_ledock_file(receptor='pro.pdb', rmsd=1.0, x=[0, 0], y=[0, 0], z=[0, 0], n_poses=10, l_list=[],
                         l_list_outfile='', out='dock.in'):
    """
    LeDock requires a single dock.in file as input for the docking process. This file contains all the necessary
    instructions for the docking process including the receptor, docking box coordinates and the path to another
    ligand.list file which contains the file paths for all ligands to be docked. It is also possible to use LePro
    with the protein as input in order to automatically generate a dock.in file with an automatic docking box
    recognition. However, if the binding pocket is already known, one has to provide the coordinates manually.

    Function taken from: https://github.com/AngelRuizMoreno/Jupyter_Dock/tree/main/utilities

    Parameters
    ----------
    receptor: .pdb file.
    rmsd: Minimum RMSD between different poses.
    x, y, z: List containing the minimum and maximum coordinate for the docking box in the corresponding dimension.
    n_poses: Amount of best poses per ligand to be provided in the output.
    l_list: List containing paths to the ligands.
    l_list_outfile: Output path for the ligand.list.
    out: Output path for the dock.in file.
    """
    rmsd = str(rmsd)
    x = [str(x) for x in x]
    y = [str(y) for y in y]
    z = [str(z) for z in z]
    n_poses = str(n_poses)

    with open(l_list_outfile, 'w') as l_out:
        for element in l_list:
            l_out.write(element)
    l_out.close()

    file = [
        'Receptor\n',
        receptor + '\n\n',
        'RMSD\n',
        rmsd + '\n\n',
        'Binding pocket\n',
        x[0], ' ', x[1], '\n',
        y[0], ' ', y[1], '\n',
        z[0], ' ', z[1], '\n\n',
        'Number of binding poses\n',
        n_poses + '\n\n',
        'Ligands list\n',
        l_list_outfile + '\n\n',
        'END']

    with open(out, 'w') as output:
        for line in file:
            output.write(line)
    output.close()


def vina_docking(ligand_dir, protein_dir, protein, docking_coordinates, results_dir):
    """
    Docks a single protein with its corresponding ligand using the python implementation of Autodock Vina.

    Parameters
    ----------
    ligand_dir: Filepath to ligands.
    protein_dir: Filepath to proteins.
    protein: Name of the protein.
    docking_coordinates: Dictionary containing coordinates of box center and box size.
    results_dir: Filepath to docked ligands.
    """
    # ---- paths, directories, general variables ---- #
    # Input
    ligand_file = ligand_dir + protein + "_ligand.pdbqt"
    protein_file = protein_dir + protein + ".pdbqt"
    box_center, box_size = docking_coordinates[0]
    # Output
    docked_ligand_file = results_dir + "vina_output/" + protein + "_docked_ligand.pdbqt"
    docked_ligand_sdf = results_dir + "vina_output/" + protein + "_docked_ligand.sdf"

    # ---- Set up directory ---- #
    if not os.path.exists(results_dir + "vina_output"):
        os.mkdir(results_dir + "vina_output")

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


def smina_docking(ligand_dir, protein_dir, protein, docking_coordinates, results_dir, bin_dir, system='mac'):
    """
    Docks a single protein with its corresponding ligand using the command line tool of Smina.

    Parameters
    ----------
    ligand_dir: Filepath to ligands.
    protein_dir: Filepath to proteins.
    protein: Name of the protein.
    docking_coordinates: Dictionary containing coordinates of box center and box size.
    results_dir: Filepath to docked ligands.
    bin_dir: Binary executables.
    system: Operating system. Choose mac (default) or linux.
    """
    # ---- Check input system (LeDock and Smina are only available for mac and linux) ---- #
    systems = ['mac', 'linux']
    if system not in systems:
        raise ValueError("Invalid OS. Expected one of: %s" % systems)

    # ---- Set up directory ---- #
    if not os.path.exists(results_dir + "smina_output"):
        os.mkdir(results_dir + "smina_output")

    # ---- paths, directories, general variables ---- #
    # Input
    ligand_file = ligand_dir + protein + "_ligand.mol2"
    protein_file = protein_dir + protein + ".pdb"
    box_center, box_size = docking_coordinates[0]
    # Output
    docked_ligand_sdf = results_dir + "smina_output/" + protein + "_docked_ligand.sdf"
    # Executable
    if system == 'mac':
        smina_binary = bin_dir + "smina.osx"
    else:
        smina_binary = bin_dir + "smina"
    # ---- docking ---- #
    smina_input = " -r %s -l %s -o %s" % (protein_file, ligand_file, docked_ligand_sdf)
    smina_box_center = " --center_x %.3f --center_y %.3f --center_z %.3f" % (box_center['center_x'], box_center['center_y'], box_center['center_z'])
    smina_box_size = " --size_x %.3f --size_y %.3f --size_z %.3f" % (box_size['size_x'], box_size['size_y'], box_size['size_z'])
    smina_specs = " --exhaustiveness 10 --num_modes 10 --seed -42"
    os.system(smina_binary + smina_input + smina_box_center + smina_box_size + smina_specs)


def ledock_docking(ligand_dir, protein_dir, protein, docking_coordinates, results_dir, bin_dir, system='mac'):
    """
    Docks a single protein with its corresponding ligand using the command line tool of LeDock. Before running the tool,
    an instructions dock.in and ligand.list is created containing necessary parameters for the docking run and paths to
    the ligands.

    Parameters
    ----------
    ligand_dir: Filepath to ligands.
    protein_dir: Filepath to proteins.
    protein: Name of the protein.
    docking_coordinates: Dictionary containing coordinates of box center and box size.
    results_dir: Filepath to docked ligands.
    bin_dir: Binary executables.
    system: Operating system. Choose mac (default) or linux.
    """
    # ---- Check input system (LeDock and Smina are only available for mac and linux) ---- #
    systems = ['mac', 'linux']
    if system not in systems:
        raise ValueError("Invalid OS. Expected one of: %s" % systems)

    # ---- Set up directory ---- #
    if not os.path.exists(results_dir + "ledock_output"):
        os.mkdir(results_dir + "ledock_output")

    # ---- paths, directories, general variables ---- #
    # >>>> NOTE: Some paths must be relative to the binary executable of LeDock instead of the main script
    # Input
    ligand_file = ligand_dir + protein + "_ligand.mol2"
    protein_file = protein_dir + protein + ".pdb"
    X, Y, Z = docking_coordinates[1]
    # Output
    docking_instructions = results_dir + "ledock_output/" + protein + "_dock.in"
    ligand_list_out = results_dir + "ledock_output/" + protein + "_ligand.list"
    docked_ligand_sdf = results_dir + "ledock_output/" + protein + "_docked_ligand.sdf"
    # Executable
    if system == 'mac':
        ledock_binary = bin_dir + "ledock_mac "
    else:
        ledock_binary = bin_dir + "ledock_linux_x86 "

    # ---- generating dock.in file and ligand.list ---- #
    generate_ledock_file(receptor=protein_file, x=[X['minX'], X['maxX']], y=[Y['minY'], Y['maxY']],
                         z=[Z['minZ'], Z['maxZ']], l_list=ligand_file, l_list_outfile=ligand_list_out,
                         out=docking_instructions)

    # ---- docking ---- #
    os.system(ledock_binary + docking_instructions)  # LeDock saves .dok files in same directory as the ligands
    os.system("mv %s*.dok %s" % (ligand_dir, results_dir+"ledock_output/"))

    dok_to_sdf(dok_file=results_dir+"ledock_output/"+protein+"_ligand.dok",
               output=docked_ligand_sdf)


def dock_ligands(ligand_dir: str, protein_dir: str, results_dir: str, bin_dir: str, system='mac'):
    """
    Performs docking of all proteins in a directory with their respective ligand using Vina, Smina and LeDock.

    Parameters
    ----------
    ligand_dir: Filepath to ligands.
    protein_dir: Filepath to proteins.
    results_dir: Filepath to docked ligands.
    bin_dir: Binary executables.
    system: Operating system. Choose mac (default) or linux.
    """
    # ---- Set up directory ---- #
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
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
    '''print("...docking %d proteins with Autodock Vina..." % len(protein_names))
    for i, protein_setup in enumerate(zip(protein_names, dock_boxes)):
        print("Protein: ", protein, " (%d/%d)" % (i + 1, len(protein_names)))

        # paths, directories, general variables
        protein = protein_setup[0]
        dock_coords = protein_setup[1]
        vina_docking(ligand_dir, protein_dir, protein, dock_coords, results_dir)
        print("------------------------------------------------------\n")'''

    # ---- docking using Smina ---- #
    '''print("...docking %d proteins with Smina..." % len(protein_names))
    for i, protein_setup in enumerate(zip(protein_names, dock_boxes)):
        print("Protein: ", protein, " (%d/%d)" % (i + 1, len(protein_names)))

        # paths, directories, general variables
        protein = protein_setup[0]
        dock_coords = protein_setup[1]
        smina_docking(ligand_dir, protein_dir, protein, dock_coords, results_dir, bin_dir, system=system)
        print("------------------------------------------------------\n")'''

    # ---- docking using LeDock ---- #
    print("...docking %d proteins with LeDock..." % len(protein_names))
    for i, protein_setup in enumerate(zip(protein_names, dock_boxes)):
        print("Protein: ", protein, " (%d/%d)" % (i + 1, len(protein_names)))

        # paths, directories, general variables
        protein = protein_setup[0]
        dock_coords = protein_setup[1]
        ledock_docking(ligand_dir, protein_dir, protein, dock_coords, results_dir, bin_dir, system=system)
        print("------------------------------------------------------\n")


# ==================================================================================================================== #
def run(ligand_dir, protein_dir, results_dir, bin_dir, system='mac'):
    # ---- Check input system (LeDock and Smina are only available for mac and linux) ---- #
    systems = ['mac', 'linux']
    if system not in systems:
        raise ValueError("Invalid OS. Expected one of: %s" % systems)
    dock_ligands(ligand_dir=ligand_dir, protein_dir=protein_dir, results_dir=results_dir, bin_dir=bin_dir,
                 system=system)
