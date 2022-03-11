"""
Main script for re-docking a set of proteins with their corresponding ligand with Autodock vina, smina and LeDock in
order to determine the RMSD of the prediction from the actual ligand position. This pipeline also includes preprocessing
of the ligands and proteins.
"""
# ====================================================== IMPORTS ===================================================== #
import argparse
import os

import src.prepare_proteins_and_ligands as preparation
import src.RMSD_experiment_docking as docking
import src.calculate_RMSD as calculate_RMSD


# ==================================================================================================================== #
def main():
    # ---- Argument input ---- #
    parser = argparse.ArgumentParser(description="Frontiers in applied drug design - Jonathan Waehrer, RMSD experiment")
    parser.add_argument("--proteins", required=True,
                        help="Path to .pdb files for docking. Ligands (if present) will be removed.")
    parser.add_argument("--bin", required=True, help="Path to required binary executables.")
    parser.add_argument("--os", default='mac',
                        help="Affects choice of binary executables. Choose between mac (default) or linux.")
    args = parser.parse_args()

    proteins = args.proteins
    bin_path = args.bin
    system = args.os
    # Check input system (LePro is only available for mac and linux)
    systems = ['mac', 'linux']
    if system not in systems:
        raise ValueError("Invalid OS. Expected one of: %s" % systems)

    # change and set directories:
    current_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(current_path)
    if not os.path.exists("out"):
        os.mkdir("out")

    prepared_ligand_dir = "out/prepared_ligands/"
    prepared_protein_dir = "out/prepared_proteins/"
    docking_results_dir = "out/docking_results/"
    rmsd_dir = "out/RMSD/"

    # ---- Run pipeline ---- #
    print("# ================== RMSD COMPARISON OF DOCKING TOOLS ================= #")
    print("# ---- Receptor- and Ligand-preparation for Vina, Smina and LeDock ---- #")
    preparation.run(input_proteins=proteins, path_to_bin=bin_path, system=system)
    print("# ------- Docking original ligands using Vina, Smina and LeDock ------- #")
    docking.run(ligand_dir=prepared_ligand_dir, protein_dir=prepared_protein_dir, results_dir=docking_results_dir,
                bin_dir=bin_path, system=system)
    print("# ------ Computing RMSD between predicted and original positions ------ #")
    calculate_RMSD.run(original_lig_dir=prepared_ligand_dir, predicted_lig_dir=docking_results_dir, rmsd_dir=rmsd_dir)


if __name__ == "__main__":
    main()
