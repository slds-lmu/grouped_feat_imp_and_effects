# This script creates Figures 2, 3 and 4 of the paper (section 3.2)

library(ggplot2)
library(dplyr)
source("code/simulations/helper.R")
theme_set(theme_bw())


# number of simulations
vec = 1:20

# create Figure 2: comparison of RF and SVM with same correlations of features within groups

# RF plot without correlations
path = "results/simulation_results/sim_corr_rf_base/results/"
gfi = get_gfi(vec, path)

colors = c("GOPFI" = "cadetblue1", "GPFI" = "cadetblue3", "GSI" = "cadetblue4", "LOGI" = "coral", "LOGO" = "coral3")
p1 = ggplot(gfi[-which(gfi$features=="G1"),], aes(x = features, y = rel_imp, fill = method)) + 
  geom_boxplot(position=position_dodge(1)) + 
  ylim(-0.5,1.5) +
  scale_fill_manual(values = colors) +
  labs(x = "Group", y = "Relative Importance", fill = "Legend") +
  scale_x_discrete(labels = c(expression("G"[2]), expression("G"[3]), expression("G"[4]))) +
  ggtitle("Random Forest")


# SVM plot without correlations
path = "results/simulation_results/sim_corr_svm_base/results/"
gfi = get_gfi(vec, path)
  
colors = c("GOPFI" = "cadetblue1", "GPFI" = "cadetblue3", "GSI" = "cadetblue4", "LOGI" = "coral", "LOGO" = "coral3")
p2 = ggplot(gfi[-which(gfi$features=="G1"),], aes(x = features, y = rel_imp, fill = method)) + 
  geom_boxplot(position=position_dodge(1)) + 
  ylim(-0.5,1.5) +
  scale_fill_manual(values = colors) +
  labs(x = "Group", y = "Relative Importance", fill = "Legend") +
  scale_x_discrete(labels = c(expression("G"[2]), expression("G"[3]), expression("G"[4]))) +
  ggtitle("SVM")
  


p = gridExtra::grid.arrange(p1,p2, nrow = 1)
ggsave("results/figures/corr_comparison_same.png", p, width = 9, height = 5)



# create Figure 3: comparison of RF and SVM with differing correlations of features within groups

# RF plot without correlations
path = "results/simulation_results/sim_corr_rf_diff/results/"
gfi = get_gfi(vec, path)

colors = c("GOPFI" = "cadetblue1", "GPFI" = "cadetblue3", "GSI" = "cadetblue4", "LOGI" = "coral", "LOGO" = "coral3")
p3 = ggplot(gfi[-which(gfi$features=="G1"),], aes(x = features, y = rel_imp, fill = method)) + 
  geom_boxplot(position=position_dodge(1)) + 
  ylim(-0.5,1.5) +
  scale_fill_manual(values = colors) +
  labs(x = "Group", y = "Relative Importance", fill = "Legend") +
  scale_x_discrete(labels = c(expression("G"[2]), expression("G"[3]), expression("G"[4]))) +
  ggtitle("Random Forest")


# SVM plot without correlations
path = "results/simulation_results/sim_corr_svm_diff/results/"
gfi = get_gfi(vec, path)

colors = c("GOPFI" = "cadetblue1", "GPFI" = "cadetblue3", "GSI" = "cadetblue4", "LOGI" = "coral", "LOGO" = "coral3")
p4 = ggplot(gfi[-which(gfi$features=="G1"),], aes(x = features, y = rel_imp, fill = method)) + 
  geom_boxplot(position=position_dodge(1)) + 
  ylim(-0.5,1.5) +
  scale_fill_manual(values = colors) +
  labs(x = "Group", y = "Relative Importance", fill = "Legend") +
  scale_x_discrete(labels = c(expression("G"[2]), expression("G"[3]), expression("G"[4]))) +
  ggtitle("SVM")



p = gridExtra::grid.arrange(p3,p4, nrow = 1)
ggsave("results/figures/corr_comparison_diff.png", p, width = 9, height = 5)




# Figure 4: analyse choice of RF for differing correlations
path = "results/simulation_results/sim_corr_rf_diff/results/"

df = get_split_features(vec, path)

df %>%
  group_by(group, splitvarName, repetition) %>%
  summarize(percentage = sum(leftChild)/10000) %>%
  ggplot(aes(x = splitvarName, y = percentage)) + 
  geom_boxplot(aes(fill = group)) + 
  scale_x_discrete(breaks = paste0("V", 1:40), labels = paste0("X", 1:40)) +
  labs(x = "Features", y = "Percentag", fill = "Group")

ggsave("results/figures/splitting_distr.png", width = 10, height = 5)


