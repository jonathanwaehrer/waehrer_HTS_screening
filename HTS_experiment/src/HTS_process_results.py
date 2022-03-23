"""
Purpose is to pass the results of the HTS experiment in order to extract the scores (affinities in kcal/mol) and docking
times resulting from the experiment. Also, the ligand efficiency of every ligand is calculated.
"""

# ====================================================== IMPORTS ===================================================== #
import numpy as np
import os
import pandas as pd

from glob import glob
from pathlib import Path
from pymol import cmd
from tqdm import tqdm


# ===================================================== FUNCTIONS ==================================================== #
def count_non_hydrogen_atoms(ligand_dir):
    """
    For all .mol2 files in provided directory, counts non-hydrogen atoms.

    Parameters
    ----------
    ligand_dir: Path to .mol2 files

    Returns
    ----------
    List of atom counts (non-hydrogen) and a list of ligand names in the order they were processed.
    """
    ligand_files = sorted(glob(ligand_dir + "[!concatenated_ligands]*.mol2"))  # single .mol2 files only!
    atom_counts = []
    ligand_names = []
    for ligand_mol2 in tqdm(ligand_files, desc="...counting non-hydrogen atoms in ligands...   "):
        # -- load file -- #
        cmd.reinitialize()
        cmd.load(ligand_mol2)
        # -- count non H -- #
        cmd.select("non_hydrogen", "not e. H")
        atom_counts.append(cmd.count_atoms("non_hydrogen"))
        # -- get ligand_name -- #
        ligand_names.append(Path(ligand_mol2).stem.split('.')[0])  # removes filepath and extension

    return atom_counts, ligand_names


def get_scores_from_sdf(input_sdf, score_entry_string):
    """
    Extracts the score of all poses from .sd-files. Depending on the docking tool, the score either follows the line
    ">  <Score>" (Vina and LeDock) or "> <minimizedAffinity>"

    Smina .sd-files contain all predicted poses in a single file. The score of a predicted ligand/pose follows right
    after the line '> <minimizedAffinity>'. This means that this pattern can be used to append the desired line
    containing the score to a list.

    Parameters
    ----------
    input_sdf: Path to .sdf.
    score_entry_string: String preceding the score. ">  <Score>" (Vina and LeDock) or "> <minimizedAffinity>" (Smina)

    Returns
    ----------
    List containing the scores for each pose.
    """
    scores = []
    next_line = False

    predicted_poses = open(input_sdf, 'r')
    for line in predicted_poses.readlines():
        if next_line:
            scores.append(float(line))
        if line.startswith(score_entry_string):
            # This means that the next iteration will contain the score => set next_line to true to get the score
            next_line = True
        else:
            next_line = False
    return scores


def get_scores(results_dir, ligand_dir):
    """
    Iterates over results directories and calls corresponding functions to extract score and calculate ligand
    efficiencies. Depending on the docking tool and its corresponding output, the results must be handled differently
    (for an example, Smina outputs a multi-.sd-file per protein, while Vina and LeDock have one .sdf per docked ligand
    as output).

    Parameters
    ----------
    results_dir: Directory containing docking results (out/docking_results/).
    ligand_dir: Directory containing .mol2 files of the prepared ligands (for atom count)
    """
    # ---- Get ligand atom counts and names ---- #
    counts, names = count_non_hydrogen_atoms(ligand_dir)

    # ---- Get vina scores from single .sd-files ---- #
    vina_results = sorted(glob(results_dir + "vina_output/*/*.sdf"))
    scores = []
    protein_names = []
    for vina_sdf in tqdm(vina_results, desc="...parsing Vina output (single .sd-files)...   "):
        scores.extend(get_scores_from_sdf(input_sdf=vina_sdf, score_entry_string=">  <Score>"))
        protein_names.append(vina_sdf.split('/')[-2])

    # Create dataframe of scores and calculate ligand efficiency
    multiplication_factor = int(len(vina_results) / len(names))  # ratio: #docked ligands versus ligand library size
    vina_scores_df = pd.DataFrame(
        {'name': names * multiplication_factor, 'docked_protein': protein_names, 'score': scores,
         'no_atoms': counts * multiplication_factor, 'Tool': 'Vina'})
    vina_scores_df['ligand_efficiency'] = vina_scores_df['score'] / vina_scores_df['no_atoms']
    # Save as .tsv
    vina_data_file = os.path.join(results_dir, "vina_output", "vina_HTS_results.tsv")
    vina_scores_df.to_csv(vina_data_file, index=False, sep='\t')

    # ---- Get smina scores from multi .sd-files ---- #
    smina_results = sorted(glob(results_dir + "smina_output/*/*.sdf"))
    scores = []
    protein_names = []
    for smina_sdf in tqdm(smina_results, desc="...parsing Smina output (multi .sd-files)...   "):
        scores.extend(get_scores_from_sdf(input_sdf=smina_sdf, score_entry_string="> <minimizedAffinity>"))
        protein_names.extend([smina_sdf.split('/')[-2]] * len(names))  # protein name for long dataframe format

    # Create dataframe of scores and calculate ligand efficiency 
    smina_scores_df = pd.DataFrame(
        {'name': names * len(smina_results), 'docked_protein': protein_names, 'score': scores,
         'no_atoms': counts * len(smina_results), 'Tool': 'Smina'})
    smina_scores_df['ligand_efficiency'] = smina_scores_df['score'] / smina_scores_df['no_atoms']
    # Save as .tsv
    smina_data_file = os.path.join(results_dir, "smina_output", "smina_HTS_results.tsv")
    smina_scores_df.to_csv(smina_data_file, index=False, sep='\t')

    # ---- Get LeDock scores from single .sd-files (similar to vina)---- #
    ledock_results = sorted(glob(results_dir + "ledock_output/*/*.sdf"))
    scores = []
    protein_names = []
    for ledock_sdf in tqdm(ledock_results, desc="...parsing LeDock output (single .sd-files)... "):
        scores.extend(get_scores_from_sdf(input_sdf=ledock_sdf, score_entry_string=">  <Score>"))
        protein_names.append(ledock_sdf.split('/')[-2])

    # Create dataframe of scores and calculate ligand efficiency
    multiplication_factor = int(len(vina_results) / len(names))  # ratio: #docked ligands versus ligand library size
    ledock_scores_df = pd.DataFrame(
        {'name': names * multiplication_factor, 'docked_protein': protein_names, 'score': scores,
         'no_atoms': counts * multiplication_factor, 'Tool': 'LeDock'})
    ledock_scores_df['ligand_efficiency'] = ledock_scores_df['score'] / ledock_scores_df['no_atoms']
    # Save as .tsv
    ledock_data_file = os.path.join(results_dir, "ledock_output", "ledock_HTS_results.tsv")
    ledock_scores_df.to_csv(ledock_data_file, index=False, sep='\t')


def top_ligands_per_protein(hts_results: str, scoring_criteria='score', amount=50, separator='\t'):
    """
    Extracts the top N (default: 50) ligands (based on the lowest score or ligand efficiency) from a single
    hts_results.tsv file.

    Parameters
    ----------
    hts_results: Path to data containing hts_results.
    scoring_criteria: Sorting criterum to base the top 50 on. Choose 'score' (default) or 'ligand_efficiency'.
    amount: Integer specifying the amount of top ligands to be extracted. Per default, the top 50 ligands are retrieved.
    separator: Separator for file read-in (default is tab separated).

    Returns
    ----------
    pd.DataFrame containing the specified amount of top ligands per protein.
    """
    # ---- Read data and sort scores protein-wise ---- #
    hts_results = pd.read_csv(hts_results, sep=separator)
    hts_results = hts_results.sort_values(by=[scoring_criteria, 'docked_protein'])
    # ---- Get top N ligands with lowest score per protein ---- #
    top_ligands = []
    ranked_ligands = []
    for protein in hts_results['docked_protein'].unique():
        top_ligands.append(hts_results[hts_results['docked_protein'] == protein].iloc[0:amount])
        # reset index of sorted array to get rank for each ligand
        ranked_ligands.append(hts_results[hts_results['docked_protein'] == protein].reset_index(drop=True))

    # return the best ligands per protein as single long data frame and ranks for each ligand
    return pd.concat(top_ligands), pd.concat(ranked_ligands)


def rank_consensus(dataframe, amount=50, keyword='Progress'):
    """
    Takes input dataframe with rank-sorted index per protein and calculates the mean rank of each ligand (protein-wise)
    and the corresponding standard deviation over all tools that were used.

    Parameters
    ----------
    dataframe: pd.Dataframe containing ranks as index.
    amount: Amount of top ligands to return per protein.
    keyword: String for progress bar.

    Returns
    ----------
    A full pd.DataFrame containing mean rank and variance for each ligand as well as the top 50 ligands per protein
    based on the lowest mean of ranks
    """
    # ---- Index (column 'Unnamed: 0') corresponds to rank ----
    dataframe['Index'] = dataframe.index + 1  # in order to have ranks from e.g. 1 to 500 instead of 0 to 499
    mean_rank = dataframe.groupby(by=['name', 'docked_protein'], as_index=False)['Index'].mean()
    mean_rank['SD'] = dataframe.groupby(by=['name', 'docked_protein'], as_index=False)['Index'].std()['Index']
    mean_rank.rename(columns={'Index': 'Mean rank'}, inplace=True)
    mean_rank.sort_values(by=['docked_protein', 'Mean rank', 'SD'], inplace=True)

    # ---- Get top N based on sum of ranks ---- #
    top_ligands = []
    for protein in tqdm(mean_rank['docked_protein'].unique(), desc=keyword):
        top_ligands.append(mean_rank[mean_rank['docked_protein'] == protein].sort_values(by='Mean rank').iloc[0:amount])
    return mean_rank, pd.concat(top_ligands)


def overlap_heatmap_data(mean_rank_data, scored_data):
    """
    Function that calculates the protein-wise fraction of overlap in the top N ligands between all 3 tools and
    their resulting rank-mean-consensus data.

    Parameters
    ----------
    mean_rank_data: Consensus data containing the top N protein-wise mean ranks over all tools based on either score or ligand efficiency.
    scored_data: Data containing the top N ligands per protein based on either based on either score or ligand efficiency.

    Returns
    ----------
    Heatmap data of fraction of overlapping ligands.
    """
    # ----- Input dataframes are in long format => reshape to wide format and concatenate horizontally ----- #
    # Make unique identifier for each consensus set of ligand names per protein
    mean_rank_data['identifier'] = mean_rank_data['docked_protein'] + "_consensus"
    mean_rank_data['idx'] = mean_rank_data.groupby('identifier').cumcount()

    # Make unique identifier for each tool-wise set of ligand names per protein
    scored_data['identifier'] = scored_data['docked_protein'] + '_' + scored_data['Tool']
    scored_data['idx'] = scored_data.groupby('identifier').cumcount()

    # reshape to wide format and concatenate horizontally
    mean_rank_data = mean_rank_data.pivot(columns='identifier', values='name', index='idx')
    scored_data = scored_data.pivot(columns='identifier', values='name', index='idx')
    wide_data = pd.concat([mean_rank_data, scored_data], axis=1)

    # ----- Calculate overlap between columns i and j ----- #
    heatmap = np.zeros((len(wide_data.columns), len(wide_data.columns)))  # quadratic heatmap
    normalization_factor = len(wide_data)  # used to calculate fraction of overlap

    for i in range(len(wide_data.columns)):
        for j in range(len(wide_data.columns)):
            overlap = len(set(wide_data.iloc[:, i]).intersection(set(wide_data.iloc[:, j]))) / normalization_factor
            if i == 10000:
                print(wide_data.iloc[:, i].name, i)
                print(wide_data.iloc[:, j].name, j)
                print(overlap)
                print("------------\n")
            heatmap[j, i] = overlap

    heatmap = pd.DataFrame(heatmap, columns=wide_data.columns, index=wide_data.columns)
    return heatmap


def top_ligands_per_tool(results_dir, top_amount=50):
    """
    Extracts the N (default = 50) most affine ligands per protein for each tool that is used.

    Parameters
    ----------
    results_dir: Directory containing docking results (out/docking_results/).
    top_amount: Amount of ligands retrieved per protein for each tool.

    Returns
    ----------
    pd.DataFrame containing the specified amount of top ligands per protein for all tools.
    pd.DataFrame containing only the OVERLAPPING ligands per protein in the top N for all tools.
    """
    # ---- Make directory for plot data ---- #
    data_output = os.path.join(results_dir, "HTS_plots")
    if not os.path.exists(data_output):
        os.mkdir(data_output)

    # ---- Get top N based on score ---- #
    results_data_files = glob(results_dir + "/*/*HTS_results.tsv")
    top_N_per_tool_score = []  # top N ligands for each tool per protein
    for tool_results in tqdm(results_data_files,
                             desc="...extracting top %d ligands per protein (based on score)...            " % top_amount):
        top_N_per_tool_score.append(top_ligands_per_protein(hts_results=tool_results, amount=top_amount))
    # Save top N .tsv
    top_ligands_frame = pd.concat([dataframe[0] for dataframe in top_N_per_tool_score])
    top_ligands_frame.to_csv(os.path.join(data_output, "top_%d_ligands.tsv" % top_amount), sep='\t', index=False)

    # Only overlapping ligands (occur in top N for at least 2 tools) to .tsv
    top_ligands_frame_ov = top_ligands_frame[top_ligands_frame.duplicated(subset=['name', 'docked_protein'], keep=False)]
    top_ligands_frame_ov.to_csv(os.path.join(data_output, "top_%d_ligands_overlap.tsv" % top_amount), sep='\t',
                             index=False)

    # ---- Get top N based on ligand efficiency ---- #
    top_N_per_tool_LE = []  # top N ligands for each tool per protein
    for tool_results in tqdm(results_data_files,
                             desc="...extracting top %d ligands per protein (based on ligand efficiency)..." % top_amount):
        top_N_per_tool_LE.append(
            top_ligands_per_protein(hts_results=tool_results, scoring_criteria='ligand_efficiency', amount=top_amount))
    # Save to .tsv
    top_ligands_frame_LE = pd.concat([dataframe[0] for dataframe in top_N_per_tool_LE])
    top_ligands_frame_LE.to_csv(os.path.join(data_output, "top_%d_ligands_ligand_efficiency.tsv" % top_amount),
                                sep='\t', index=False)

    # Only overlapping ligands (occur in top N for at least 2 tools) to .tsv
    top_ligands_frame_LE_ov = top_ligands_frame_LE[
        top_ligands_frame_LE.duplicated(subset=['name', 'docked_protein'], keep=False)]
    top_ligands_frame_LE_ov.to_csv(os.path.join(data_output, "top_%d_ligands_overlap_ligand_efficiency.tsv" % top_amount),
                                sep='\t', index=False)

    # ---- Mean of ranks per ligand to .csv ---- #
    # based on score
    score_ranks, score_ranks_top_N = rank_consensus(pd.concat([dataframe[1] for dataframe in top_N_per_tool_score]),
                                                    amount=top_amount,
                                                    keyword="...extracting mean ligand ranks based on affinity score...              ")

    score_ranks.to_csv(os.path.join(data_output, "ligand_ranks_score.tsv"), sep='\t', index=False)
    score_ranks_top_N.to_csv(os.path.join(data_output, "ligand_ranks_score_top_%d.tsv" % top_amount), sep='\t',
                             index=False)

    # based on ligand efficiency
    ligand_efficiency_ranks, ligand_efficiency_ranks_top_N = rank_consensus(
        pd.concat([dataframe[1] for dataframe in top_N_per_tool_LE]), amount=top_amount,
        keyword="...extracting mean ligand ranks based on ligand efficiency...           ")
    ligand_efficiency_ranks.to_csv(os.path.join(data_output, "ligand_ranks_ligand_efficiency.tsv"), sep='\t',
                                   index=False)
    ligand_efficiency_ranks_top_N.to_csv(
        os.path.join(data_output, "ligand_ranks_ligand_efficiency_top_%d.tsv" % top_amount), sep='\t', index=False)

    # ---- Heatmap of overlapping ligands ---- #
    # score
    heatmap_based_on_score = overlap_heatmap_data(mean_rank_data=score_ranks_top_N, scored_data=top_ligands_frame)
    heatmap_based_on_score.to_csv(os.path.join(data_output, "top_%d_score_ligands_overlap_heatmap.tsv" % top_amount),
                                  sep='\t', index=True)
    # ligand efficiency
    heatmap_based_on_LE = overlap_heatmap_data(mean_rank_data=ligand_efficiency_ranks_top_N, scored_data=top_ligands_frame_LE)
    heatmap_based_on_LE.to_csv(os.path.join(data_output, "top_%d_ligand_efficiency_ligands_overlap_heatmap.tsv" % top_amount),
                               sep='\t', index=True)


# ==================================================================================================================== #
def run(results_dir, ligand_dir, top_amount):
    get_scores(results_dir, ligand_dir)
    top_ligands_per_tool(results_dir=results_dir, top_amount=top_amount)
