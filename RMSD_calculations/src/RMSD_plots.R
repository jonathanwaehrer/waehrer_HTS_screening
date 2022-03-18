"
Purpose of this script is to plot the results of the RMSD experiment using ggplot2. This script is for plotting only and
will not be automatically called by the main function of the RMSD experiment (prediction_RMSD.py).
This script assumes that you run it from the same directory where it is located in. If the filepath can not be found,
either adjust the path or change the directory accordingly.
"

# ---- Packages, imports ---- #
install.packages("ggplot2")
library(ggplot2)

# ---- Read RMSD file ---- #
filepath = paste0(dirname(rstudioapi::getSourceEditorContext()$path), '/../out/RMSD/RMSD_experiment_results.tsv')
df = read.csv(filepath, sep = '\t', stringsAsFactors = TRUE)

# ---- Boxplots of RMSD per protein ---- #
p1 = ggplot(data=df, aes(y=RMSD, x=Tool, fill=Tool)) + 
  geom_boxplot() + facet_wrap(~df[,'Protein'], scales = 'free_x') +
  theme(legend.position = "None") +
  ggtitle("RMSD of predicted poses")

ggsave('../out/RMSD/1_RMSD_boxplots_per_protein.pdf', plot = p1, width = 9, height = 7)

# ---- Boxplots of best scoring position per docking tool ---- #
# keep only first entry of each protein-tool kombination
df$merged = paste(df$Protein, df$Tool)
df2 = df[!duplicated(df$merged),]

# Plot
p2 = ggplot(df2, aes(y=RMSD, x=Tool, fill=Tool)) + 
  geom_boxplot(alpha=0.7) +
  labs(title = "RMSD for best scoring positions") +
  theme(legend.position = "None")
ggsave('../out/RMSD/2_RMSD_top_positions.pdf', plot = p2, width = 5, height = 3)

# ---- Lines: Best RMSD per protein ---- #
p3 = ggplot(df2, aes(y=RMSD, x=Protein, color=Tool, group=Tool)) + 
  geom_point() + 
  geom_line(alpha=0.5) +
  labs(title = "RMSD for best scoring positions per protein")
ggsave('../out/RMSD/3_RMSD_top_vs_protein.pdf', plot = p3, width = 6, height = 3)

# ---- Lines: RMSD vs pose rank ---- #
df$Rank=NA
df$Rank = with(df,ave(Rank,df$Protein,df$Tool,FUN = seq_along))

p4 = ggplot(df, aes(y=RMSD, x=Rank, color=Tool)) +
  geom_point(alpha=0.5) +
  geom_line() + 
  facet_wrap(~df[,'Protein'], scales = 'free_y', ncol = 4) +
  ggtitle("RMSD per pose rank") +
  scale_x_continuous(breaks = 1:10)
ggsave('../out/RMSD/4_RMSD_vs_pose_rank.pdf', plot = p4, width = 9, height = 7)


# ---- Lines: Score versus RMSD ---- #
# SeeSar has no scoring output => omit 
p5 = ggplot(na.omit(df), aes(y=Score, x=RMSD, color=Tool)) +
  geom_point(alpha=0.5) +
  geom_line() + 
  facet_wrap(~na.omit(df)[,'Protein'], scales = 'free_y', ncol = 4) +
  ggtitle("Score vs RMSD for tested proteins") +
  ylab("Score [Kcal/mol]")

ggsave('../out/RMSD/5_Score_vs_RMSD.pdf', plot = p5, width = 9, height = 7)

