install.packages("ggplot2")
library(ggplot2)

filepath = '../out/RMSD/RMSD_experiment_results.tsv'
df = read.csv(filepath, sep = '\t', stringsAsFactors = TRUE)

# ---- Boxplot of RMSD per protein ---- #
p1 = ggplot(data=df, aes(y=RMSD, x=Tool, fill=Tool)) + 
  geom_boxplot() + facet_wrap(~df[,'Protein'], scales = 'free_x') +
  theme(legend.position = "None") + ylab("Score [Kcal/mol]") + 
  ggtitle("RMSD of predicted poses")

ggsave('../out/RMSD/RMSD_boxplots_per_protein.pdf', plot = p1, width = 9, height = 7)

# ---- Boxplot of best scoring position ---- #


# ---- Score versus RMSD ---- #
df = na.omit(df)  # SeeSar has no scoring output
p3 = ggplot(data = df, aes(y=Score, x=RMSD, color=Tool)) + 
  geom_line() + 
  facet_wrap(~df[,'Protein'], scales = 'free_y', ncol = 4) +
  ggtitle("Score vs RMSD for tested proteins") +
  ylab("Score [Kcal/mol]")

ggsave('../out/RMSD/Score_vs_RMSD.pdf', plot = p3, width = 9, height = 7)

