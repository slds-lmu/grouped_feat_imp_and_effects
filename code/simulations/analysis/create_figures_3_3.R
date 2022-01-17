# This script creates Figure 5 of the paper (section 3.3)

library(ggplot2)
library(gridExtra)

# parameters
path = "results/simulation_results/sim_groupsize/results/"
vec = 1:500

# create Figure 5: Comparison of Shapley importance on group and feature level
df_shap = get_shapley_imp(vec, path)
shap_feat = df_shap[[1]]
shap_group = df_shap[[2]]

colors_groups <- brewer.pal(2,"Set2")
names(colors_groups) <- levels(df$group)
p1 = ggplot(shap_group, aes(x = group, y = mse, fill = group)) + geom_boxplot() +
  scale_x_discrete(breaks = c("G1", "G2"), labels = c(expression("G"[1], "G"[2]))) +
  scale_fill_manual(values = colors_groups, labels = c(expression("G"[1]), expression("G"[2]))) +
  labs(x = "Group", y = "MSE", fill = "Group") 

p2 = ggplot(shap_feat, aes(x = feature, y = mse, fill = group)) + geom_boxplot() +
  scale_fill_manual(values = colors_groups, labels = c(expression("G"[1]), expression("G"[2]))) +
  scale_x_discrete(breaks = paste0("V", 1:8), labels = c(expression("X"[1]), expression("X"[2]),
                                                         expression("X"[3]), expression("X"[4]),
                                                         expression("X"[5]), expression("X"[6]),
                                                         expression("X"[7]), expression("X"[8]))) +
  labs(x = "Feature", y = "MSE", fill = "Group") 


p = ggpubr::ggarrange(p1,p2, ncol = 2, common.legend = TRUE, legend = "right")
ggsave("results/figures/sim_varying_size_shap.png", p, width = 9, height = 5)

