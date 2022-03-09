"""
Main script for re-docking a set of proteins with their corresponding ligand with Autodock vina, smina and LeDock in
order to determine the RMSD of the prediction from the actual ligand position.
"""
# ====================================================== IMPORTS ===================================================== #
import src.prepare_proteins_and_ligands as preparation
import src.RMSD_experiment_docking as docking
import src.calculate_RMSD as calculate_RMSD
import os


def main():
    proteins = "/Users/jonathanwahrer/Desktop/Msc/3_Semester/Drug_Design_Praktikum/waehrer_HTS_screening/data/proteins/"
    bin_path = "/Users/jonathanwahrer/Desktop/Msc/3_Semester/Drug_Design_Praktikum/waehrer_HTS_screening/bin/"
    system = "mac"

    # change and set directories:
    current_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(current_path)
    if not os.path.exists("out"):
        os.mkdir("out")

    prepared_ligand_dir = "out/prepared_ligands/"
    prepared_protein_dir = "out/prepared_proteins/"
    docking_results_dir = "out/docking_results/"
    rmsd_dir = "out/RMSD/"

    print("# ================== RMSD COMPARISON OF DOCKING TOOLS ================= #")
    print("# ---- Receptor- and Ligand-preparation for Vina, Smina and LeDock ---- #")
    # preparation.run(input_proteins=proteins, path_to_bin=bin_path, system = system)
    print("# ------- Docking original ligands using Vina, Smina and LeDock ------- #")
    # docking.run(ligand_dir=prepared_ligand_dir, protein_dir=prepared_protein_dir, results_dir=docking_results_dir, bin_dir=bin_path, system=system)
    print("# ------ Computing RMSD between predicted and original positions ------ #")
    calculate_RMSD.run(original_lig_dir=prepared_ligand_dir, predicted_lig_dir=docking_results_dir, rmsd_dir=rmsd_dir)


if __name__ == "__main__":
    main()
