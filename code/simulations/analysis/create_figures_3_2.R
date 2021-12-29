# This script creates Figures 2, 3 and 4 of the paper (section 3.2)

library(ggplot2)
library(dplyr)
library(batchtools)
source("code/simulations/helper.r")
theme_set(theme_bw())


# get experimental design
reg = loadRegistry("results/simulation_results/sim_corr")
summarizeExperiments()

# number of simulations
vec = 1:2000

# create Figure 2: comparison of RF and SVM with same correlations of features within groups

# RF plot without correlations
path = "results/simulation_results/sim_corr/results/"
gfi = get_gfi(vec, path)
params = unwrap(getJobPars())
gfi = ijoin(gfi, params, by = c("run"= "job.id"))

colors = c("GOPFI" = "cadetblue1", "GPFI" = "cadetblue3", "GSI" = "cadetblue4", "LOGI" = "coral", "LOGO" = "coral3")
# RF same correlation
gfi_rf_base = gfi[which(gfi$algo=="rf" & gfi$corr =="base"),]
p1 = ggplot(gfi_rf_base[-which(gfi_rf_base$features=="G1"),], aes(x = features, y = rel_imp, fill = method)) + 
  geom_boxplot(position=position_dodge(1)) + 
  ylim(-2,3.5) +
  scale_fill_manual(values = colors) +
  labs(x = "Group", y = "Relative Importance", fill = "Legend") +
  scale_x_discrete(labels = c(expression("G"[2]), expression("G"[3]), expression("G"[4]))) +
  ggtitle("Random Forest")


# SVM plot with same correlation
#path = "results/simulation_results/sim_corr_svm_base/results/"
#gfi = get_gfi(vec, path)
  
colors = c("GOPFI" = "cadetblue1", "GPFI" = "cadetblue3", "GSI" = "cadetblue4", "LOGI" = "coral", "LOGO" = "coral3")
gfi_svm_base = gfi[which(gfi$algo=="svm" & gfi$corr =="base"),]
p2 = ggplot(gfi_svm_base[-which(gfi_svm_base$features=="G1"),], aes(x = features, y = rel_imp, fill = method)) + 
  geom_boxplot(position=position_dodge(1)) + 
  ylim(-2,3.5) +
  scale_fill_manual(values = colors) +
  labs(x = "Group", y = "Relative Importance", fill = "Legend") +
  scale_x_discrete(labels = c(expression("G"[2]), expression("G"[3]), expression("G"[4]))) +
  ggtitle("SVM")
  


p = gridExtra::grid.arrange(p1,p2, nrow = 1)
ggsave("results/figures/corr_comparison_same.png", p, width = 9, height = 5)



# create Figure 3: comparison of RF and SVM with differing correlations of features within groups

# RF plot 
#path = "results/simulation_results/sim_corr_rf_diff/results/"
#gfi = get_gfi(vec, path)

colors = c("GOPFI" = "cadetblue1", "GPFI" = "cadetblue3", "GSI" = "cadetblue4", "LOGI" = "coral", "LOGO" = "coral3")
gfi_rf_diff = gfi[which(gfi$algo=="rf" & gfi$corr =="diff"),]
p3 = ggplot(gfi_rf_diff[-which(gfi_rf_diff$features=="G1"),], aes(x = features, y = rel_imp, fill = method)) + 
  geom_boxplot(position=position_dodge(1)) + 
  ylim(-1,2) +
  scale_fill_manual(values = colors) +
  labs(x = "Group", y = "Relative Importance", fill = "Legend") +
  scale_x_discrete(labels = c(expression("G"[2]), expression("G"[3]), expression("G"[4]))) +
  ggtitle("Random Forest")


# SVM plot 
#path = "results/simulation_results/sim_corr_svm_diff/results/"
#gfi = get_gfi(vec, path)

colors = c("GOPFI" = "cadetblue1", "GPFI" = "cadetblue3", "GSI" = "cadetblue4", "LOGI" = "coral", "LOGO" = "coral3")
gfi_svm_diff = gfi[which(gfi$algo=="svm" & gfi$corr =="diff"),]
p4 = ggplot(gfi_svm_diff[-which(gfi_svm_diff$features=="G1"),], aes(x = features, y = rel_imp, fill = method)) + 
  geom_boxplot(position=position_dodge(1)) + 
  ylim(-1,2) +
  scale_fill_manual(values = colors) +
  labs(x = "Group", y = "Relative Importance", fill = "Legend") +
  scale_x_discrete(labels = c(expression("G"[2]), expression("G"[3]), expression("G"[4]))) +
  ggtitle("SVM")



p = gridExtra::grid.arrange(p3,p4, nrow = 1)
ggsave("results/figures/corr_comparison_diff.png", p, width = 9, height = 5)




# Figure 4: analyse choice of RF for differing correlations
path = "results/simulation_results/sim_corr_rf_diff/results/"

df = get_split_features( which(params$algo=="rf" & params$corr=="diff"), path)

df %>%
  group_by(group, splitvarName, repetition) %>%
  summarize(percentage = sum(leftChild)/10000) %>%
  ggplot(aes(x = splitvarName, y = percentage)) + 
  geom_boxplot(aes(fill = group)) + 
  scale_x_discrete(breaks = paste0("V", 1:40), labels = paste0("X", 1:40)) +
  labs(x = "Features", y = "Percentag", fill = "Group")

ggsave("results/figures/splitting_distr.png", width = 10, height = 5)


