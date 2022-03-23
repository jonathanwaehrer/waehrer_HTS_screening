"
Purpose of this script is to plot the results of the HTS experiment using ggplot2. This script is for plotting only and
will not be automatically called by the main function of the HTS experiment (virtual_screening.py).
"

# ---- Packages, imports ---- #
#install.packages("ggplot2")
library(ggplot2)

# ---- Read scores ---- #
script_path = dirname(rstudioapi::getSourceEditorContext()$path)       # path of THIS script
data_path = paste0(script_path, '/../out/docking_results/HTS_plots/')  # directory containing data for plots


# ===== TOP 50 PLOTS BASED ON SCORE ===== #
# ---- Lines: Top 50 ligand scores per protein (facets) per tool (color) ---- #
full_top_50 = read.csv(paste0(data_path, "top_50_ligands.tsv"), sep = '\t')
p1 = ggplot(full_top_50, aes(x=name, y=score, color=Tool)) + 
  # line
  geom_line() + geom_point(alpha=0.3) +
  # facets
  facet_wrap(~full_top_50[,'docked_protein'], ncol = 2) +
  # theme, labels
  labs(x = "Ligand", y= "Score [kcal/mol]", title = "Top 50 ligands per protein")

ggsave(paste0(data_path, "1_top_50_ligands.pdf"), plot = p1, width = 11, height = 7)

# ---- Lines: Scores over overlapping top 50 ligands in all tools ---- #
overlap_top_50 = read.csv(paste0(data_path, "top_50_ligands_overlap.tsv"), sep = '\t')
p2 = ggplot(overlap_top_50, aes(x=as.factor(name), y=score, color=Tool, group=Tool)) + 
  # lines & points
  geom_line() + geom_point(alpha=0.5) +
  # facets
  facet_wrap(~overlap_top_50[,'docked_protein'], ncol = 2, scales = "free_x") +
  # theme, labss
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Ligand", y= "Score [kcal/mol]", title = "Overlapping ligands between tools in top 50")

ggsave(paste0(data_path, "2_top_50_ligands_overlap.pdf"), plot = p2, width = 11, height = 7)


# ===== TOP 50 PLOTS BASED ON LIGAND EFFICIENCY ===== #
# ---- Lines: Top 50 ligand scores per protein (facets) per tool (color) ---- #
full_top_50_LE = read.csv(paste0(data_path, "top_50_ligands_ligand_efficiency.tsv"), sep = '\t')
p3 = ggplot(full_top_50_LE, aes(x=name, y=ligand_efficiency, color=Tool)) + 
  # line
  geom_line() + geom_point(alpha=0.3) +
  # facets
  facet_wrap(~full_top_50_LE[,'docked_protein'], ncol = 2) +
  # theme, labels
  labs(x = "Ligand", y= "Ligand Efficiency", title = "Top 50 ligands per protein (LE)")

ggsave(paste0(data_path, "3_top_50_ligands_ligand_efficiency.pdf"), plot = p3, width = 11, height = 7)

# ---- Lines: Scores over overlapping top 50 ligands in all tools ---- #
overlap_top_50_LE = read.csv(paste0(data_path, "top_50_ligands_overlap_ligand_efficiency.tsv"), sep = '\t')
p4 = ggplot(overlap_top_50_LE, aes(x=as.factor(name), y=ligand_efficiency, color=Tool, group=Tool)) + 
  # lines & points
  geom_line() + geom_point(alpha=0.5) +
  # facets
  facet_wrap(~overlap_top_50_LE[,'docked_protein'], ncol = 2, scales = "free_x") +
  # theme, labss
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Ligand", y= "Ligand Efficiency", title = "Overlapping ligands between tools in top 50 (LE)")

ggsave(paste0(data_path, "4_top_50_ligands_overlap_ligand_efficiency.pdf"), plot = p4, width = 11, height = 7)
