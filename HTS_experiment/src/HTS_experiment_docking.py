"""
Purpose of this script is to perform high throughput screening using Vina, Smina and LeDock given a set of
(preprocessed) .pdb files and a set of (preprocessed) ligands.

This script assumes that the ligands are provided with the following naming scheme:
    - my_ligand.pdbqt (for Vina)
    - my_ligand.mol2 (for LeDock)
    - all_ligands_multi.mol2 (for Smina, all ligands concatenated to single .mol2 file)

For docking box calculations, this script assumes that the original ligands of the proteins are contained in a directory
waehrer_HTS_screening/data/original_protein_ligands/*.mol2
"""

# ====================================================== IMPORTS ===================================================== #
import os
import pandas as pd
import time

from glob import glob
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

    Function taken and adapted from: https://github.com/AngelRuizMoreno/Jupyter_Dock/tree/main/utilities

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
            element = element[element.find("/out"):][1:]  # Modified, because of LeDock path problem
            l_out.write(element + "\n")
    l_out.close()

    # Modified, because of LeDock path problem:
    l_list_outfile = l_list_outfile[l_list_outfile.find("/out"):][1:]
    receptor = receptor[receptor.find("/out"):][1:]

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


def vina_HTS_docking(ligand_dir, protein_dir, protein, docking_coordinates, results_dir):
    """
    Docks a single protein with a library of prepared ligands using the python implementation of Autodock Vina. For this
    purpose, all .pdbqt files in the directory are chosen as docking candidates. HTS with vina can only be done
    iteratively over the ligands by calling the bash script (or python function) once for each ligand. An advantage of
    the python implementation is that the docking box must be calculated only once before iterative docking. When
    calling Vina from the command line, the same docking box would be calculated once per ligand which slows down the
    screening process.

    Parameters
    ----------
    ligand_dir: Filepath to ligands.
    protein_dir: Filepath to proteins.
    protein: Name of the protein.
    docking_coordinates: Dictionary containing coordinates of box center and box size.
    results_dir: Filepath to docked ligands.

    Returns
    ----------
    Lists containing the docking time for each ligand that was docked to the receptor.
    """
    # ---- paths, directories, general variables ---- #
    # Input
    ligand_files = sorted(
        glob(ligand_dir + "[!concatenated_ligands]*.pdbqt"))  # this excludes concatenated_ligands.mol2
    protein_file = protein_dir + protein + ".pdbqt"
    box_center, box_size = docking_coordinates[0]

    # ---- Set up directories ---- #
    protein_docked_dir = os.path.join(results_dir, "vina_output", protein)
    if not os.path.exists(protein_docked_dir):
        os.makedirs(protein_docked_dir, exist_ok=True)  # create one directory per protein

    # ---- setting up vina ---- #
    v = Vina(sf_name='Vina', seed=42)
    v.set_receptor(protein_file)
    v.compute_vina_maps(center=[box_center['center_x'], box_center['center_y'], box_center['center_z']],
                        box_size=[box_size['size_x'], box_size['size_y'], box_size['size_z']])

    # ---- iterative docking ---- #
    # tracking arrays for .tsv output of time consumption
    ligand_docking_times = []  # will contain docking times per ligand for the current protein
    docked_ligands = []
    docked_protein = []
    for ligand in tqdm(ligand_files, desc="...%s: docking ligands...                  " % protein):
        # update tracking arrays
        docked_ligands.append(os.path.basename(ligand)[:-6])
        docked_protein.append(protein)

        # Output variables
        docked_ligand_file = os.path.join(protein_docked_dir, os.path.basename(ligand)[:-6] + "_docked.pdbqt")
        docked_ligand_sdf = os.path.join(protein_docked_dir, os.path.basename(ligand)[:-6] + "_docked.sdf")

        # Run docking
        start = time.time()  # start timer
        v.set_ligand_from_file(ligand)
        v.dock(exhaustiveness=5, n_poses=1)
        v.write_poses(docked_ligand_file, n_poses=1, overwrite=True)
        ligand_docking_times.append(time.time() - start)  # end timer
        pdbqt_to_sdf(pdbqt_file=docked_ligand_file, output=docked_ligand_sdf)

    return ligand_docking_times, docked_ligands, docked_protein


def smina_HTS_docking(ligand_dir, protein_dir, protein, docking_coordinates, results_dir, bin_dir, system='mac'):
    """
    Docks a single protein with arbitrarily many ligands using Smina. The ligands must be contained in a single
    concatenated .mol2 file. Given the ligand-directory as input parameter, this functions assumes that the file is
    named 'concatenated_ligands.mol2'.

    Parameters
    ----------
    ligand_dir: Filepath to ligands.
    protein_dir: Filepath to proteins.
    protein: Name of the protein.
    docking_coordinates: Dictionary containing coordinates of box center and box size.
    results_dir: Filepath to docked ligands.
    bin_dir: Binary executables.
    system: Operating system. Choose mac (default) or linux.

    Returns
    ----------
    Time in seconds required to dock all ligands with the input receptor.
    """
    # ---- Check input system (LeDock and Smina are only available for mac and linux) ---- #
    systems = ['mac', 'linux']
    if system not in systems:
        raise ValueError("Invalid OS. Expected one of: %s" % systems)

    # ---- Set up directory ---- #
    protein_docked_dir = os.path.join(results_dir, "smina_output", protein)
    if not os.path.exists(protein_docked_dir):
        os.makedirs(protein_docked_dir, exist_ok=True)  # create one directory per protein

    # ---- paths, directories, general variables ---- #
    # Input
    multi_mol_name = "concatenated_ligands.mol2"  # Smina requires multi-ligand .mol2 file
    ligand_file = os.path.join(ligand_dir, multi_mol_name)
    if not os.path.exists(ligand_file):
        raise FileNotFoundError(ligand_file + " does not exists, but is required for Smina HTS-docking.")
    protein_file = os.path.join(protein_dir, protein + ".pdb")
    box_center, box_size = docking_coordinates[0]
    # Output
    docked_ligand_sdf = os.path.join(results_dir, "smina_output/", protein, protein + "_docked_ligands_smina.sdf")
    # Executable
    if system == 'mac':
        smina_binary = bin_dir + "smina.osx"
    else:
        smina_binary = bin_dir + "smina"
    # ---- docking ---- #
    # Smina is called via the command line (os.system). Partial strings for better readability:
    smina_input = " -r %s -l %s -o %s" % (protein_file, ligand_file, docked_ligand_sdf)
    smina_box_center = " --center_x %.3f --center_y %.3f --center_z %.3f" % (
        box_center['center_x'], box_center['center_y'], box_center['center_z'])
    smina_box_size = " --size_x %.3f --size_y %.3f --size_z %.3f" % (
        box_size['size_x'], box_size['size_y'], box_size['size_z'])
    smina_specs = " --exhaustiveness 5 --num_modes 1 --seed -42"
    smina_log_out = " --log %s" % (docked_ligand_sdf.rstrip(".sdf") + ".log")

    start = time.time()  # start timer
    os.system(smina_binary + smina_input + smina_box_center + smina_box_size + smina_specs + smina_log_out)
    return time.time() - start  # end timer


def ledock_HTS_docking(ligand_dir, protein_dir, protein, docking_coordinates, results_dir, bin_dir, system='mac'):
    """
    Docks a single protein with a library of ligands using the command line tool of LeDock. Before running the tool,
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

    Returns
    ----------
    Time in seconds required to dock all ligands with the input receptor.
    """
    # ---- Check input system (LeDock and Smina are only available for mac and linux) ---- #
    systems = ['mac', 'linux']
    if system not in systems:
        raise ValueError("Invalid OS. Expected one of: %s" % systems)

    # ---- Set up directory ---- #
    protein_docked_dir = os.path.join(results_dir, "ledock_output", protein)
    if not os.path.exists(protein_docked_dir):
        os.makedirs(protein_docked_dir, exist_ok=True)  # create one directory per protein

    # ---- paths, directories, general variables ---- #
    # >>>> NOTE: Some paths must be relative to the binary executable of LeDock or the main script in order to work
    # Input
    protein_file = os.path.join(protein_dir, protein + ".pdb")
    ligand_files = sorted(glob(ligand_dir + "[!concatenated_ligands]*.mol2"))  # excludes concatenated_ligands.mol2
    X, Y, Z = docking_coordinates[1]
    # Output
    docking_instructions = os.path.join(results_dir, "ledock_output/", protein, protein + "_dock.in")
    ligand_list_out = os.path.join(results_dir, "ledock_output/", protein, protein + "_ligand.list")
    # Executable
    if system == 'mac':
        ledock_binary = bin_dir + "ledock_mac "
    else:
        ledock_binary = bin_dir + "ledock_linux_x86 "

    # ---- generating dock.in file and ligand.list ---- #
    generate_ledock_file(receptor=protein_file, x=[X['minX'], X['maxX']], y=[Y['minY'], Y['maxY']],
                         z=[Z['minZ'], Z['maxZ']], l_list=ligand_files, l_list_outfile=ligand_list_out,
                         out=docking_instructions, n_poses=1)

    # ---- docking ---- #
    start = time.time()  # start timer
    os.system(ledock_binary + docking_instructions)  # LeDock saves .dok files in same directory as the ligands
    end = time.time() - start  # end timer
    current_protein_out = os.path.join(results_dir, "ledock_output/", protein)
    os.system("mv %s*.dok %s" % (ligand_dir, current_protein_out))  # move to out/docking_results/ledock_output/protein/

    # ---- iteratively, convert .dok files to .sdf ---- #
    for ligand_dok in glob(os.path.join(current_protein_out, "*.dok")):
        docked_ligand_sdf = ligand_dok.rstrip(".dok") + "_docked.sdf"
        dok_to_sdf(dok_file=ligand_dok, output=docked_ligand_sdf)

    return end


def dock_ligands(ligand_dir: str, protein_dir: str, results_dir: str, bin_dir: str, original_ligands_dir: str,
                 system='mac'):
    """
    Performs docking of all proteins in a directory with their respective ligand using Vina, Smina and LeDock.
    Additionally, this function also saves the required docking time per protein as .tsv file.

    Parameters
    ----------
    ligand_dir: Filepath to docking ligands.
    protein_dir: Filepath to proteins.
    results_dir: Filepath to docked ligands.
    bin_dir: Binary executables.
    original_ligands_dir: Path to original ligands of the proteins (used for docking box_calculation)
    system: Operating system. Choose mac (default) or linux.
    """
    # ---- Set up directories ---- #
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    current_path = os.path.dirname(os.path.realpath(__file__))
    original_ligands_dir = original_ligands_dir

    # ---- Check input system (LeDock and Smina are only available for mac and linux) ---- #
    systems = ['mac', 'linux']
    if system not in systems:
        raise ValueError("Invalid OS. Expected one of: %s" % systems)

    # ---- calculate docking boxes ---- #
    protein_names = sorted([i.split('.')[0] for i in os.listdir(protein_dir) if i.endswith(".pdb")])[10:12]  # unique
    original_ligand_names = sorted([i for i in os.listdir(original_ligands_dir) if i.endswith(".mol2")])[10:12]
    dock_boxes = []
    for protein, ligand in tqdm(zip(protein_names, original_ligand_names),
                                desc="...calculating docking boxes...          ",
                                total=len(protein_names)):
        ligand = original_ligands_dir + ligand
        dock_boxes.append(get_box(ligand=ligand, software='both'))  # index 0: vina format; index 1: LeDock format

    # ---- docking using Vina ---- #
    print("...docking %d protein(s) with Autodock Vina..." % len(protein_names))
    # Lists for time comparison
    vina_times = []  # keeps track of elapsed time in seconds for each ligand that was docked to a single receptor
    vina_ligands = []
    vina_proteins = []
    for i, protein_setup in enumerate(zip(protein_names, dock_boxes)):
        # paths, directories, general variables
        protein = protein_setup[0]
        dock_coords = protein_setup[1]
        # run Vina
        print("Protein: ", protein, " (%d/%d)" % (i + 1, len(protein_names)))
        hts_docking_data = vina_HTS_docking(ligand_dir, protein_dir, protein, dock_coords, results_dir)
        print("------------------------------------------------------\n")
        # Extend lists for time comparison
        vina_times.extend(hts_docking_data[0])
        vina_ligands.extend(hts_docking_data[1])
        vina_proteins.extend(hts_docking_data[2])

    # save measured time as .tsv:
    vina_times_df = pd.DataFrame({'dock_time': vina_times, 'ligand': vina_ligands, 'protein': vina_proteins})
    vina_times_df.to_csv(os.path.join(results_dir, "vina_output", "vina_docking_times.tsv"), sep="\t", index=False)

    # ---- docking using Smina ---- #
    print("...docking %d protein(s) with Smina..." % len(protein_names))
    # Lists for time comparison
    smina_times = []  # keeps track of elapsed time in seconds for each ligand that was docked to a single receptor
    smina_proteins = []
    for i, protein_setup in enumerate(zip(protein_names, dock_boxes)):
        # paths, directories, general variables
        protein = protein_setup[0]
        dock_coords = protein_setup[1]
        # run Smina
        print("Protein: ", protein, " (%d/%d)" % (i + 1, len(protein_names)))
        smina_times.append(smina_HTS_docking(ligand_dir, protein_dir, protein, dock_coords, results_dir,
                                             bin_dir, system=system))
        print("------------------------------------------------------\n")
        # Extend lists for time comparison
        smina_proteins.append(protein)

    # save measured time as .tsv:
    smina_times_df = pd.DataFrame({'dock_time': smina_times, 'protein': smina_proteins})
    smina_times_df.to_csv(os.path.join(results_dir, "smina_output", "smina_docking_times.tsv"), sep="\t", index=False)

    # ---- docking using LeDock ---- #
    print("...docking %d protein(s) with LeDock..." % len(protein_names))
    # Lists for time comparison
    ledock_times = []  # keeps track of elapsed time in seconds for each ligand that was docked to a single receptor
    ledock_proteins = []
    for i, protein_setup in enumerate(zip(protein_names, dock_boxes)):
        # paths, directories, general variables
        protein = protein_setup[0]
        dock_coords = protein_setup[1]
        # run LeDock
        print("Protein: ", protein, " (%d/%d)" % (i + 1, len(protein_names)))
        ledock_times.append(ledock_HTS_docking(ligand_dir, protein_dir, protein, dock_coords, results_dir, bin_dir,
                                               system=system))
        print("------------------------------------------------------\n")
        # Extend lists for time comparison
        ledock_proteins.append(protein)

    # save measured time as .tsv:
    ledock_times_df = pd.DataFrame({'dock_time': ledock_times, 'protein': ledock_proteins})
    ledock_times_df.to_csv(os.path.join(results_dir, "ledock_output", "ledock_docking_times.tsv"), sep="\t",
                           index=False)


# ==================================================================================================================== #
def run(ligand_dir, protein_dir, results_dir, bin_dir, original_ligands_dir, system='mac'):
    """
    Parameters
    ----------
    ligand_dir: Filepath to docking ligands.
    protein_dir: Filepath to proteins.
    results_dir: Filepath to docked ligands.
    bin_dir: Binary executables.
    original_ligands_dir: Path to original ligands of the proteins (used for docking box_calculation)
    system: Operating system. Choose mac (default) or linux.
    """
    # ---- Check input system (LeDock and Smina are only available for mac and linux) ---- #
    systems = ['mac', 'linux']
    if system not in systems:
        raise ValueError("Invalid OS. Expected one of: %s" % systems)
    # ---- Run docking ---- #
    dock_ligands(ligand_dir=ligand_dir, protein_dir=protein_dir, results_dir=results_dir, bin_dir=bin_dir,
                 original_ligands_dir=original_ligands_dir, system=system)
