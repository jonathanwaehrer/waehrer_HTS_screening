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


# ==================================================================================================================== #
def main():
    # ---- Argument input ---- #
    parser = argparse.ArgumentParser(description="Frontiers in applied drug design - Jonathan Waehrer, HTS experiment")
    parser.add_argument("--proteins", required=True,
                        help="Path to .pdb files for docking. Ligands (if present) will be removed.")
    parser.add_argument("--bin", required=True, help="Path to required binary executables.")
    parser.add_argument("--os", default='mac',
                        help="Affects choice of binary executables. Choose between mac (default) or linux.")
    parser.add_argument("--pre", action='store_true',
                        help="OPTIONAL: Add this parameter if preprocessing (including sampling) is desired. Can be omitted if a desired set of (preprocessed) ligands and proteins already exists in directory ./out/prepared_ligands.\n"
                             "NOTE: Adding this parameter will remove all ligands (.mol2 and .pdbqt) in HTS_experiment/out/prepared_ligands before resampling and start of preprocessing. This can NOT be undone.")
    parser.add_argument("--ligands", help="Required if sampling is to be performed (--pre). Path to ligand .mol2 file for sampling.")
    parser.add_argument("--sample", help="Integer specifying sample size.", type=int, default=500)
    parser.add_argument("--seed", help="Seed for sampling during preprocessing (=42).", type=int, default=42)
    args = parser.parse_args()

    # Check for missing inputs
    if args.pre:
        if args.ligands not in args:
            raise ValueError("Preprocessing requires '--ligands' as input directory.")
        if args.sample not in args:
            raise ValueError("Preprocessing requires '--sample' for sample size.")
        lig_library = args.ligands
        sample_size = args.sample

    proteins = args.proteins
    bin_path = args.bin
    system = args.os
    seed = args.seed

    # Check input system (LePro is only available for mac and linux)
    systems = ['mac', 'linux']
    if system not in systems:
        raise ValueError("Invalid OS. Expected one of: %s" % systems)

    # change and set directories:
    current_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(current_path)
    if not os.path.exists(current_path + "/out"):
        os.mkdir(current_path + "/out")
    prepared_ligand_dir = current_path + "/out/prepared_ligands/"
    prepared_protein_dir = current_path + "/out/prepared_proteins/"
    docking_results_dir = current_path + "/out/docking_results/"

    print("# ================== VIRTUAL SCREENING WITH DOCKING TOOLS ================= #")
    # Preprocessing:
    if args.pre:
        print("# ------------------------- Sampling ligands -------------------------- #")
        sample.run(ligand_library=lig_library, sample_directory=prepared_ligand_dir, sample_size=sample_size, seed=seed)
        print("\n# ---- Receptor- and Ligand-preparation for Vina, Smina and LeDock ---- #")
        preparation.run(input_proteins=proteins, path_to_bin=bin_path, system=system)
    # Actual HTS:
    print("# ------------ Performing HTS using Vina, Smina and LeDock ------------ #")
    docking.run(ligand_dir=prepared_ligand_dir, protein_dir=prepared_protein_dir, results_dir=docking_results_dir,
                bin_dir=bin_path, system=system)


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