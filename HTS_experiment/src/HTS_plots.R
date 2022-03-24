"
Purpose of this script is to plot the results of the HTS experiment using ggplot2. This script is for plotting only and
will not be automatically called by the main function of the HTS experiment (virtual_screening.py).
"

# ---- Packages, imports ---- #
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("tidytext")
library(ggplot2)
library(pheatmap)
library(tidytext)

# ---- Read scores ---- #
script_path = dirname(rstudioapi::getSourceEditorContext()$path)       # path of THIS script
data_path = paste0(script_path, '/../out/docking_results/HTS_plots/')  # directory containing data for plots


# ===== TOP 50 PLOTS BASED ON SCORE ===== #
# ---- 1) Lines: Top 50 ligand scores per protein (facets) per tool (color) ---- #
full_top_50 = read.csv(paste0(data_path, "top_50_ligands.tsv"), sep = '\t')
p1 = ggplot(full_top_50, aes(x=name, y=score, color=Tool)) + 
  # line
  geom_line() + geom_point(alpha=0.3) +
  # facets
  facet_wrap(~full_top_50[,'docked_protein'], ncol = 2) +
  # theme, labels
  labs(x = "Ligand", y= "Score [kcal/mol]", title = "Top 50 ligands per protein")

ggsave(paste0(data_path, "1_top_50_ligands.pdf"), plot = p1, width = 11, height = 9)

# ---- 2) Lines: Scores over overlapping top 50 ligands in all tools ---- #
overlap_top_50 = read.csv(paste0(data_path, "top_50_ligands_overlap.tsv"), sep = '\t')
p2 = ggplot(overlap_top_50, aes(x=as.factor(name), y=score, color=Tool, group=Tool)) + 
  # lines & points
  geom_line() + geom_point(alpha=0.5) +
  # facets
  facet_wrap(~overlap_top_50[,'docked_protein'], ncol = 2, scales = "free_x") +
  # theme, labss
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Ligand", y= "Score [kcal/mol]", title = "Overlapping ligands between tools in top 50")

ggsave(paste0(data_path, "2_top_50_ligands_overlap.pdf"), plot = p2, width = 11, height = 9)


# ===== TOP 50 PLOTS BASED ON LIGAND EFFICIENCY ===== #
# ---- 3) Lines: Top 50 ligand scores per protein (facets) per tool (color) ---- #
full_top_50_LE = read.csv(paste0(data_path, "top_50_ligands_ligand_efficiency.tsv"), sep = '\t')
p3 = ggplot(full_top_50_LE, aes(x=name, y=ligand_efficiency, color=Tool)) + 
  # line
  geom_line() + geom_point(alpha=0.3) +
  # facets
  facet_wrap(~full_top_50_LE[,'docked_protein'], ncol = 2) +
  # theme, labels
  labs(x = "Ligand", y= "Ligand Efficiency", title = "Top 50 ligands per protein (LE)")

ggsave(paste0(data_path, "3_top_50_ligands_ligand_efficiency.pdf"), plot = p3, width = 11, height = 9)

# ---- 4) Lines: Scores over overlapping top 50 ligands in all tools ---- #
overlap_top_50_LE = read.csv(paste0(data_path, "top_50_ligands_overlap_ligand_efficiency.tsv"), sep = '\t')
p4 = ggplot(overlap_top_50_LE, aes(x=as.factor(name), y=ligand_efficiency, color=Tool, group=Tool)) + 
  # lines & points
  geom_line() + geom_point(alpha=0.5) +
  # facets
  facet_wrap(~overlap_top_50_LE[,'docked_protein'], ncol = 2, scales = "free_x") +
  # theme, labss
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Ligand", y= "Ligand Efficiency", title = "Overlapping ligands between tools in top 50 (LE)")

ggsave(paste0(data_path, "4_top_50_ligands_overlap_ligand_efficiency.pdf"), plot = p4, width = 11, height = 9)


# ===== LINES: RANK CONSENSUS (TOP 50 LIGANDS) OF USED TOOLS ===== #
# ---- 5) Based on score ---- #
score_rank_consensus = read.csv(paste0(data_path, "ligand_ranks_score_top_50.tsv"), sep = '\t')
p5 = ggplot(data = score_rank_consensus, aes(y=Mean.rank, x=reorder_within(name, Mean.rank, docked_protein), group=1)) + 
  # lines & points
  geom_line() + geom_point(alpha=0.2) +
  facet_wrap(~docked_protein, ncol = 2, scales = "free_x") +
  scale_x_reordered() + 
  # error bars
  geom_errorbar(aes(ymin=Mean.rank-SD, ymax=Mean.rank+SD), width=1, alpha=0.4) + 
  # axis labels, title
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = "Ligand", y = "Mean rank", title = "Top 50 of Mean rank consensus based on score")
ggsave(paste0(data_path, "5_top_50_rank_consensus_score.pdf"), plot = p5, width = 11, height = 9)

# ---- 6) Based on LE ---- #
LE_rank_consensus = read.csv(paste0(data_path, "ligand_ranks_ligand_efficiency_top_50.tsv"), sep = '\t')
LE_rank_consensus$name = as.factor(LE_rank_consensus$name)
p6 = ggplot(data = LE_rank_consensus, aes(y=Mean.rank, x=reorder_within(name, Mean.rank, docked_protein), group=1)) + 
  # lines & points
  geom_line() + geom_point(alpha=0.2) +
  facet_wrap(~docked_protein, ncol = 2, scales = "free_x") +
  scale_x_reordered() + 
  # error bars
  geom_errorbar(aes(ymin=Mean.rank-SD, ymax=Mean.rank+SD), width=1, alpha=0.4) + 
  # axis labels, title
  #scale_x_discrete(labels=LE_rank_consensus$name)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = "Ligand", y = "Mean rank", title = "Top 50 of Mean rank consensus based on LE")
ggsave(paste0(data_path, "6_top_50_rank_consensus_LE.pdf"), plot = p6, width = 11, height = 9)


# ===== HEATMAP OF OVERLAP ===== #
# ---- 7) Based on score ---- #
score_heat = read.table(paste0(data_path, "top_50_score_ligands_overlap_heatmap.tsv"), sep = '\t', header = TRUE, row.names = 1)
colnames(score_heat) = row.names(score_heat)
p7 = pheatmap(score_heat, main = "Fraction of overlap in top 50 ligands (by score)", border_color = NA)
ggsave(paste0(data_path, "7_overlap_heatmap_score.pdf"), plot = p7, width = 8.5, height = 8.5)

# ---- 8) Based on ligand efficiency ---- #
le_heat = read.table(paste0(data_path, "top_50_ligand_efficiency_ligands_overlap_heatmap.tsv"), sep = '\t', header = TRUE, row.names = 1)
colnames(le_heat) = row.names(le_heat)
p8 = pheatmap(le_heat, main = "Fraction of overlap in top 50 ligands (by LE)", border_color = NA)
ggsave(paste0(data_path, "8_overlap_heatmap_LE.pdf"), plot = p8, width = 8.5, height = 8.5)


# ===== 9) LINES: PROCESSING TIME ===== #
times = read.csv(paste0(data_path, "docking_times.tsv"), sep = '\t')
times$dock_time = times$dock_time/60  # seconds => minutes

p9 = ggplot(data = times, aes(x=factor(protein), y=dock_time, group=Tool, color = Tool)) + 
  geom_line() + geom_point(alpha=0.3) +
  labs(x='Protein', y='Docking Time [min]', title = "Docking time per protein") +
  scale_y_continuous(breaks = seq(0, 180, 10))

ggsave(paste0(data_path, "9_docking_time_per_protein.pdf"), plot = p9, width = 6, height = 3)


