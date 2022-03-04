"""
Main script for re-docking a set of proteins with their corresponding ligand with Autodock vina, smina and LeDock in
order to determine the RMSD of the prediction from the actual ligand position.
"""
# ====================================================== IMPORTS ===================================================== #
import src.prepare_proteins_and_ligands as preparation
import os


def main():
    inp = "/Users/jonathanwahrer/Desktop/Msc/3_Semester/Drug_Design_Praktikum/waehrer_HTS_screening/data/proteins/"
    lepr = "/Users/jonathanwahrer/Desktop/Msc/3_Semester/Drug_Design_Praktikum/waehrer_HTS_screening/bin/"
    preparation.run(inp, lepr)


if __name__ == "__main__":
    main()
