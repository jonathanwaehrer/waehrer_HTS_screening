"""
Purpose of this script if to perform virtual screening on a set of target proteins with a (defined random) sample of a
ligand library (directory of .mol2 files) using Autodock vina, smina and LeDock. This pipeline also includes a random
sampling step of the ligand library (user specifies sample size) as well as preprocessing of proteins and ligands.
"""
# ====================================================== IMPORTS ===================================================== #
import argparse
import os

import src.HTS_sample_ligands as sample
import src.HTS_prepare_proteins_and_ligands as preparation
import src.HTS_experiment_docking as docking
import src.HTS_process_results as postprocessing


# ==================================================================================================================== #
def main():
    # ---- Argument input ---- #
    parser = argparse.ArgumentParser(description="Frontiers in applied drug design - Jonathan Waehrer, HTS experiment")
    parser._action_groups.pop()
    required_args = parser.add_argument_group("Required parameters")
    required_args.add_argument("--proteins", required=True,
                               help="Path to .pdb files for docking. Ligands (if present) will be removed.")
    required_args.add_argument("--bin", required=True, help="Path to required binary executables.")
    required_args.add_argument("--initial_ligands", required=True, help="Path to the original ligands' .mol2 files for each protein. The file names must start with the corresponding protein name.")
    optional_args = parser.add_argument_group("Optional parameters")
    optional_args.add_argument("--ligands",
                        help="Path to ligand .mol2 files for docking. Will lead to sampling and preprocessing. Can be omitted if a prepared library (including .pdbqt files for Vina) already exists under ./out/prepared_ligands. NOTE: Adding this parameter will remove all ligands (.mol2 and .pdbqt) in HTS_experiment/out/prepared_ligands before resampling and start of preprocessing. This can NOT be undone.")
    optional_args.add_argument("--os", default='mac',
                        help="Affects choice of binary executables. Choose between mac (default) or linux.")
    optional_args.add_argument("--sample", help="Integer specifying sample size.", type=int, default=500)
    optional_args.add_argument("--seed", help="Seed for sampling during preprocessing (=42).", type=int, default=42)
    args = parser.parse_args()

    # ---- Assign variables from arguments ---- #
    proteins = args.proteins
    bin_path = args.bin
    original_ligand_path = args.initial_ligands

    system = args.os
    seed = args.seed
    lig_library = args.ligands
    sample_size = args.sample

    # ---- Check input system (LePro is only available for mac and linux) ---- #
    systems = ['mac', 'linux']
    if system not in systems:
        raise ValueError("Invalid OS. Expected one of: %s" % systems)

    # ---- Change and set directories ---- #
    current_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(current_path)
    if not os.path.exists(current_path + "/out"):
        os.mkdir(current_path + "/out")
    prepared_ligand_dir = current_path + "/out/prepared_ligands/"
    prepared_protein_dir = current_path + "/out/prepared_proteins/"
    docking_results_dir = current_path + "/out/docking_results/"

    # ---- Start run ---- #
    print("# ================== VIRTUAL SCREENING WITH DOCKING TOOLS ================= #")
    # Preprocessing:
    if args.ligands is not None:
        print("# ------------------------- Sampling ligands -------------------------- #")
        sample.run(ligand_library=lig_library, sample_directory=prepared_ligand_dir, sample_size=sample_size, seed=seed)
        print("\n# ---- Receptor- and Ligand-preparation for Vina, Smina and LeDock ---- #")
        preparation.run(input_proteins=proteins, path_to_bin=bin_path, system=system)
    # Actual HTS:
    print("# ------------ Performing HTS using Vina, Smina and LeDock ------------ #")
    docking.run(ligand_dir=prepared_ligand_dir, protein_dir=prepared_protein_dir, original_ligands_dir=original_ligand_path, results_dir=docking_results_dir, bin_dir=bin_path, system=system)
    print("# ------------------- Parsing HTS docking results --------------------- #")
    # postprocessing.run(results_dir=docking_results_dir, ligand_dir=prepared_ligand_dir)


if __name__ == "__main__":
    main()

    '''proteins = "/Users/jonathanwahrer/Desktop/Msc/3_Semester/Drug_Design_Praktikum/waehrer_HTS_screening/data/proteins/"
    bin_path = "/Users/jonathanwahrer/Desktop/Msc/3_Semester/Drug_Design_Praktikum/waehrer_HTS_screening/bin/"
    system = "mac"
    lig_library = "/Users/jonathanwahrer/Desktop/Msc/3_Semester/Drug_Design_Praktikum/waehrer_HTS_screening/data/ligands/"
    sample_size = 3
    seed = 42
    pre = False
    '''
