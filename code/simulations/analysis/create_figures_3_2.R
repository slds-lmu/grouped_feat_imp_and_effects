# This script creates Figures 2, 3 and 4 of the paper (section 3.2)

library(ggplot2)
library(dplyr)
library(batchtools)
library(RColorBrewer)
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

display.brewer.all(colorblindFriendly = T)
colors = c("GOPFI" = brewer.pal(11, "RdYlBu")[8], "GPFI" = brewer.pal(11, "RdYlBu")[9], "GSI" = brewer.pal(11, "RdYlBu")[10],
           "LOGI" = brewer.pal(11, "RdYlBu")[4], "LOGO" = brewer.pal(11, "RdYlBu")[3])
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
  
gfi_svm_base = gfi[which(gfi$algo=="svm" & gfi$corr =="base"),]
p2 = ggplot(gfi_svm_base[-which(gfi_svm_base$features=="G1"),], aes(x = features, y = rel_imp, fill = method)) + 
  geom_boxplot(position=position_dodge(1)) + 
  ylim(-2,3.5) +
  scale_fill_manual(values = colors) +
  labs(x = "Group", y = "Relative Importance", fill = "Legend") +
  scale_x_discrete(labels = c(expression("G"[2]), expression("G"[3]), expression("G"[4]))) +
  ggtitle("SVM")
  

p = ggpubr::ggarrange(p1,p2, ncol = 2, common.legend = TRUE, legend = "right")
ggsave("results/figures/corr_comparison_same.png", p, width = 9, height = 5)



# create Figure 3: comparison of RF and SVM with differing correlations of features within groups

# RF plot 

gfi_rf_diff = gfi[which(gfi$algo=="rf" & gfi$corr =="diff"),]
p3 = ggplot(gfi_rf_diff[-which(gfi_rf_diff$features=="G1"),], aes(x = features, y = rel_imp, fill = method)) + 
  geom_boxplot(position=position_dodge(1)) + 
  ylim(-1,2) +
  scale_fill_manual(values = colors) +
  labs(x = "Group", y = "Relative Importance", fill = "Legend") +
  scale_x_discrete(labels = c(expression("G"[2]), expression("G"[3]), expression("G"[4]))) +
  ggtitle("Random Forest")


# SVM plot 

gfi_svm_diff = gfi[which(gfi$algo=="svm" & gfi$corr =="diff"),]
p4 = ggplot(gfi_svm_diff[-which(gfi_svm_diff$features=="G1"),], aes(x = features, y = rel_imp, fill = method)) + 
  geom_boxplot(position=position_dodge(1)) + 
  ylim(-1,2) +
  scale_fill_manual(values = colors) +
  labs(x = "Group", y = "Relative Importance", fill = "Legend") +
  scale_x_discrete(labels = c(expression("G"[2]), expression("G"[3]), expression("G"[4]))) +
  ggtitle("SVM")



p = ggpubr::ggarrange(p3,p4, ncol = 2, common.legend = TRUE, legend = "right")
ggsave("results/figures/corr_comparison_diff.png", p, width = 9, height = 5)




# Figure 4: analyse choice of RF for differing correlations
df = get_split_features( which(params$algo=="rf" & params$corr=="diff"), path)
colors_groups <- brewer.pal(4,"Set2")
names(colors_groups) <- levels(df$group)

df %>%
  group_by(group, splitvarName, repetition) %>%
  summarize(percentage = sum(leftChild)/10000) %>%
  ggplot(aes(x = splitvarName, y = percentage)) + 
  geom_boxplot(aes(fill = group)) + 
  scale_fill_manual(values = colors_groups) +
  scale_x_discrete(breaks = paste0("V", 1:40), labels = paste0("X", 1:40)) +
  labs(x = "Features", y = "Percentage", fill = "Group")

ggsave("results/figures/splitting_distr.png", width = 10, height = 5)


