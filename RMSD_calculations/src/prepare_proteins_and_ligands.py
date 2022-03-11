"""
Purpose of this script is the preparation of proteins and ligands using LePro and Autodock's provided tools for the
necessary preprocessing steps for Vina (conversion to .pdbqt files).
"""
# ====================================================== IMPORTS ===================================================== #
import os
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

    current_dir = os.path.dirname(os.path.realpath(__file__))  # leads to /src folder
    output_dir = current_dir + "/../out/prepared_proteins"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # ---- Iteratively use LePro (either for mac or linux) on pdb files---- #
    if system == 'mac':
        for pdb in tqdm(os.listdir(pdb_dir), desc="...preparing receptors...                      "):
            file = os.path.join(pdb_dir, pdb)
            if pdb.endswith(".pdb"):  # problem with mac: often, invisible .DS_Store files are saved causing errors.
                os.system("%s/lepro_mac %s" % (bin_dir, file))
                os.system("mv pro.pdb %s/%s" % (output_dir, pdb))  # LePro always saves output as 'pro.pdb'
    else:
        for pdb in tqdm(os.listdir(pdb_dir)):
            file = os.path.join(pdb_dir, pdb)
            if pdb.endswith(".pdb"):
                os.system("%s/lepro_linux_x86 %s" % (bin_dir, file))
                os.system("mv pro.pdb %s/%s" % (output_dir, pdb))  # LePro always saves output as 'pro.pdb'
    os.system("rm dock.in")

    # ---- Vina requires .pdbqt files for docking (AutoDock provides a tool for that) ---- #
    for pdb in tqdm(os.listdir(output_dir), desc="...adding charges and atom types for Vina...   "):
        if pdb.endswith(".pdb"):
            file = os.path.join(output_dir, pdb)
            os.system("%s/prepare_receptor -r %s -o %s/%sqt >/dev/null 2>&1" % (bin_dir, file, output_dir, pdb))


def prepare_ligands_for_vina(bin_dir: str):
    """
    Vina requires ligands in .pdbqt format as well. Can be acquire by simply applying Autodock's prepare_ligands tool.

    Parameters
    ----------
    bin_dir: Binary executables.
    """
    # ---- Get file path of output directory ---- #
    current_dir = os.path.dirname(os.path.realpath(__file__))
    ligand_dir = current_dir + "/../out/prepared_ligands"
    if not os.path.exists(ligand_dir):
        os.mkdir(ligand_dir)
    bin_dir = os.path.dirname(bin_dir)

    # ---- get ligands and apply tool for pdbqt conversion ---- #
    for ligand in tqdm(os.listdir(ligand_dir), desc="...preparing ligands for vina...         "):
        if ligand.endswith(".mol2"):  # ignores potential invisible .DS_Store files from mac
            """
            For some reason, the provided prepare_ligand binary only works if the ligand .mol2 is in the same directory
            as the executing script (which is waehrer_HTS_screening/RMSD_calculations/prediction_RMSD.py). 
            This means that the working directory must be changed before applying the tool.
            """
            file = os.path.join(ligand_dir, ligand)
            os.chdir("out/prepared_ligands/")
            os.system("%s/prepare_ligand -l %s -o %s/%spdbqt" % (bin_dir, file, ligand_dir, ligand[:-4]))
            os.chdir("../../")


# ==================================================================================================================== #
def run(input_proteins, path_to_bin, system="mac"):
    prepare_proteins(pdb_dir=input_proteins, bin_dir=path_to_bin, system=system)
    prepare_ligands_for_vina(bin_dir=path_to_bin)


'''def prepare_ligands(pdb_dir):
    """
    Extracts ligands from .pdb files and uses Open Babel's Pybel to add hydrogens to ligands. Additionally for Autodock
    Vina, the corresponding tool 'prepare_ligands' is applied to convert the .mol2 files to .pdbqt files.

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

    # ---- get ligands from chain A of pdb files and use PyBel ---- #
    for pdb in tqdm(os.listdir(pdb_dir), desc="...extracting ligands from proteins...         "):

        # Extract and save ligand:
        if pdb.endswith(".pdb"):
            extracted_ligand_path = "%s/%s_ligand.mol2" % (output_dir, pdb[:-4])
            file = os.path.join(pdb_dir, pdb)
            cmd.load(file)
            ligand_atoms = cmd.select("sele_ligand", "organic and chain A")
            if not ligand_atoms:
                # In one of the provided .pdb files, there was no ligand in chain A
                ligand_atoms = cmd.select("sele_ligand", "organic")
            cmd.save(extracted_ligand_path, format='mol2', selection='sele_ligand')
            cmd.reinitialize()

            # Add hydrogens:
            mol = [m for m in pybel.readfile(filename=extracted_ligand_path, format='mol2')][0]
            mol.addh()
            out = pybel.Outputfile(filename=extracted_ligand_path, format='mol2', overwrite=True)
            out.write(mol)
            out.close()'''