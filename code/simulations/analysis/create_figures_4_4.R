# This script creates Figures 7, 8 and 9 of the paper (section 4)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(reshape2)
library(ggExtra)

source("code/simulations/helper.R")
theme_set(theme_bw())


#---------------------------------------------------------------------------------------------------
# Figure 7: one factor/PC with linear trend

# parameters
path = "results/simulation_results/sim1latent/results/"
vec = 1:100

# CFEP plots

# linear
formula = pd ~ z
formula_true = pd_true ~ z

# data
df_1 = get_data_cfep(vec, path)


# sspca

# prepare data for plotting of estimated CFEP and confidence interval
df_sspca = df_1[["df_detailed_sspca"]]

list_sspca_f1 = data_prep_cfep(df_sspca, "z", formula, formula_true) 
result_sspca_aggr = list_sspca_f1[[1]]
pd_true = list_sspca_f1[[2]]

# plot for sspca
colors = c("Estimated CFEP (Sparse SPCA)" = "red", "Ground Truth" = "blue")
p1 = ggplot(df_sspca, aes(x = z, y = pd)) + geom_point(alpha = 0) +
  geom_line(data = result_sspca_aggr, aes(x = quantile, y = mean_prediction, color = "Estimated CFEP (Sparse SPCA)"), lty = 1, lwd = 1.2) + 
  geom_line(data = pd_true, aes(x = quantile, y = mean_prediction, color = "Ground Truth"), lty= 2, lwd = 1.2) + 
  geom_ribbon(data = result_sspca_aggr,aes(x = quantile, y = mean_prediction, ymin = lower, ymax = upper), alpha = 0.3, fill = "grey") + 
  scale_color_manual(values = colors) + 
  ylim(-20,20) + xlim(-25,25) +
  labs(x = "Z", y = "Mean Prediction", color = "") +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
p1 = ggMarginal(p1, type = "histogram", margins = "x")



# spca
df_pca = df_1[["df_detailed_pca"]]
list_pca_f1 = data_prep_cfep(df_pca, "z", formula, formula_true) 
result_pca_aggr = list_pca_f1[[1]]
pd_true = list_pca_f1[[2]]

# plot for spca
colors = c("Estimated CFEP (Sparse PCA)" = "red", "Ground Truth" = "blue")
p2 = ggplot(df_pca, aes(x = z, y = pd)) + geom_point(alpha = 0) +
  geom_line(data = result_pca_aggr, aes(x = quantile, y = mean_prediction, color = "Estimated CFEP (Sparse PCA)"), lty = 1, lwd = 1.2) + 
  geom_line(data = pd_true, aes(x = quantile, y = mean_prediction, color = "Ground Truth"), lty = 2, lwd = 1.2) + 
  geom_ribbon(data = result_pca_aggr,aes(x = quantile, y = mean_prediction, ymin = lower, ymax = upper), alpha = 0.3, fill = "grey") + 
  scale_color_manual(values = colors) +
  ylim(-20,20) + xlim(-25,25) +
  labs(x = "Z", y = "Mean Prediction", color = "") +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
p2 = ggMarginal(p2, type = "histogram", margins = "x")


# totalvis plot
#path = "results/simulation_results/sim1latent.totalvis/results/"

list_df = data_prep_totalvis(vec, path)
df_1 = list_df[[1]]$df

# first factor/PC:
df = df_1

# linear - estimated Totalvis value
formula = avg_pred ~ x_vals
q = list_sspca_f1[[3]]
mod = lm(formula = formula, data = df)
avg_pred = pred = mod$coefficients[1] + mod$coefficients[2] %*% q
data = data.frame(z = q, pred = t(avg_pred))

colors = c("Estimated Totalvis (PCA)" = "red", "Ground Truth" = "blue")

p3 = ggplot(df, aes(x = x_vals, y = avg_pred)) + geom_point(alpha = 0) + geom_line(aes(group = rep))  +# xlim(min(df$x_vals), max(df$x_vals)) +
  geom_line(data = data , aes(x = z, y = pred, color = "Estimated Totalvis (PCA)"), lty = 1, lwd = 1.2) +
  geom_line(data = pd_true, aes(x = quantile, y = mean_prediction, color = "Ground Truth"), lty= 2, lwd = 1.2) +
  theme_bw() + scale_color_manual(values = colors) +
  ylim(-20,20) + xlim(-25,25) +
  labs(x = "Z", y = "Mean Prediction", color = "") +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
p3 = ggMarginal(p3, type = "histogram", margins = "x")


# Figure 7
p = grid.arrange(p1,p2,p3, nrow = 1)

ggsave("results/figures/sim_effect_1factor.png",plot = p,  width = 9, height = 5)



#---------------------------------------------------------------------------------------------------
# Figure 8: Feature loading weights distribution 

# data
df_features_pca = df_1$df_features_pca
df_features_sspca = df_1$df_features_sspca

# calculate mean of feature loadings for spca and sspca and define true features vector
true_features = factor(colnames(df_features_pca[,1:50]) %in% c("V5","V8","V25","V47","V49") )
true_features_adj = vector("character", length = length(true_features))
true_features_adj[which(true_features==TRUE)] = "Yes"
true_features_adj[which(true_features==FALSE)] = "No"
mean_pca = data.frame("variable" = colnames(df_features_pca[,1:50]), "mean_value" = colMeans(df_features_pca[,1:50]), "true_feat" = true_features_adj)
mean_sspca = data.frame("variable" = colnames(df_features_pca[,1:50]), "mean_value" = colMeans(df_features_sspca[,1:50]), "true_feat" = true_features_adj)

# melt dataframes depending on simulation number for plotting
colnames(df_features_pca) = paste0("X",1:50)
colnames(df_features_sspca) = paste0("X",1:50)
mean_pca$variable = paste0("X",1:50)
mean_sspca$variable = paste0("X",1:50)
df_feat_pca_plot = melt(df_features_pca)
df_feat_sspca_plot = melt(df_features_sspca)

# plot feature loading weights as boxplots for all simulations for spca and sspca
p1 = ggplot(df_feat_pca_plot, aes(x = variable, y = value)) + 
  geom_boxplot()  + 
  geom_point(data = mean_pca, shape = 23, size = 2.3, aes(x = variable, y = mean_value, fill = true_feat)) +
  theme_bw() +
  labs(x = "Features", y = "Loading Weights", fill = "Feature in Z") 
p2 = ggplot(df_feat_sspca_plot, aes(x = variable, y = value)) + 
  geom_boxplot() + 
  geom_point(data = mean_sspca, shape = 23, size = 2.3, aes(x = variable, y = mean_value, fill = true_feat)) +
  theme_bw() + #ylim(-1,1) +
  labs(x = "Features", y = "Loading Weights", fill = "Feature in Z") 


# Figure 8
p = grid.arrange(p2, p1, nrow = 2)
ggsave("results/figures/feat_linear.png",p, width = 13, height = 6.5)



#---------------------------------------------------------------------------------------------------
# Figure 9: two factors/PCs

# parameters
path = "results/simulation_results/sim2latent/results/"

# CFEP plots


# data
df_1 = get_data_cfep(vec, path)

# extract sspca and pca data
df_sspca = df_1[["df_detailed_sspca"]]
df_pca = df_1[["df_detailed_pca"]]

# first factor: linear trend
formula = pd1 ~ z1
formula_true = pd_true1 ~ z1

# sspca
list_sspca_f1 = data_prep_cfep(df_sspca, "z1", formula, formula_true) 
result_sspca_aggr = list_sspca_f1[[1]]
pd_true = list_sspca_f1[[2]]


# plot for sspca
colors = c("Estimated CFEP (Sparse SPCA)" = "red", "Ground Truth" = "blue")
p1 = ggplot(df_sspca, aes(x = z1, y = pd1)) + geom_point(alpha = 0) +
  geom_line(data = result_sspca_aggr, aes(x = quantile, y = mean_prediction, color = "Estimated CFEP (Sparse SPCA)"), lwd = 1.2) + 
  geom_line(data = pd_true, aes(x = quantile, y = mean_prediction, color = "Ground Truth"), lty = 2, lwd = 1.2) + 
  geom_ribbon(data = result_sspca_aggr,aes(x = quantile, y = mean_prediction, ymin = lower, ymax = upper), alpha = 0.3, fill = "grey") + 
  theme_bw() + #ylim(min(result_sspca_aggr$lower), max(result_sspca_aggr$upper)) + 
  scale_color_manual(values = colors) + 
  ylim(0,40) + xlim(-25,25) +
  labs(x = expression("Z"[1]), y = "Mean Prediction", color = "") +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
p1 = ggMarginal(p1, type = "histogram", margins = "x")



# spca
list_pca_f1 = data_prep_cfep(df_pca, "z1", formula, formula_true) 
result_pca_aggr = list_pca_f1[[1]]
pd_true = list_pca_f1[[2]]

# plot for spca
colors = c("Estimated CFEP (Sparse PCA)" = "red", "Ground Truth" = "blue")
p2 = ggplot(df_pca, aes(x = z1, y = pd1)) + geom_point(alpha = 0) +
  geom_line(data = result_pca_aggr, aes(x = quantile, y = mean_prediction, color = "Estimated CFEP (Sparse PCA)"), lwd = 1.2) + 
  geom_line(data = pd_true, aes(x = quantile, y = mean_prediction, color = "Ground Truth"), lty = 2, lwd = 1.2) + 
  geom_ribbon(data = result_pca_aggr,aes(x = quantile, y = mean_prediction, ymin = lower, ymax = upper), alpha = 0.3, fill = "grey") + 
  theme_bw() + #ylim(min(result_sspca_aggr$lower), max(result_sspca_aggr$upper)) + 
  scale_color_manual(values = colors) + 
  ylim(0,40) + xlim(-25,25) +
  labs(x = expression("Z"[1]), y = "Mean Prediction", color = "") +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
  #ylim(min(result_sspca_aggr$lower), max(result_sspca_aggr$upper))
p2 = ggMarginal(p2, type = "histogram", margins = "x")



# plot sspca and spca next to each other
grid.arrange(p1, p2, nrow = 1)



# second factor: quadratic trend

formula = pd2 ~ poly(z2, 2, raw = TRUE)
formula_true = pd_true2 ~ poly(z2, 2, raw = TRUE)

# sspca
list_sspca_f2 = data_prep_cfep(df_sspca, "z2", formula, formula_true) 
result_sspca_aggr = list_sspca_f2[[1]]
pd_true = list_sspca_f2[[2]]


# plot for sspca
colors = c("Estimated CFEP (Sparse SPCA)" = "red", "Ground Truth" = "blue")
p3 = ggplot(df_sspca, aes(x = z2, y = pd2)) + geom_point(alpha = 0) +
  geom_line(data = result_sspca_aggr, aes(x = quantile, y = mean_prediction, color = "Estimated CFEP (Sparse SPCA)"), lwd = 1.2) + 
  geom_line(data = pd_true, aes(x = quantile, y = mean_prediction, color = "Ground Truth"), lty = 2, lwd = 1.2) + 
  geom_ribbon(data = result_sspca_aggr,aes(x = quantile, y = mean_prediction, ymin = lower, ymax = upper), alpha = 0.3, fill = "grey") + 
  theme_bw() +
  scale_color_manual(values = colors) + 
  ylim(-25,180) + xlim(-20,20) +
  labs(x = expression("Z"[2]), y = "Mean Prediction", color = "") +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
p3 = ggMarginal(p3, type = "histogram", margins = "x")



# spca

list_pca_f2 = data_prep_cfep(df_pca, "z2", formula, formula_true) 
result_pca_aggr = list_pca_f2[[1]]
pd_true = list_pca_f2[[2]]

# plot for spca
colors = c("Estimated CFEP (Sparse PCA)" = "red", "Ground Truth" = "blue")
p4 = ggplot(df_pca, aes(x = z2, y = pd2)) + geom_point(alpha = 0) +
  geom_line(data = result_pca_aggr, aes(x = quantile, y = mean_prediction, color = "Estimated CFEP (Sparse PCA)"), lwd = 1.2) + 
  geom_line(data = pd_true, aes(x = quantile, y = mean_prediction, color = "Ground Truth"), lty = 2, lwd = 1.2) + 
  geom_ribbon(data = result_pca_aggr,aes(x = quantile, y = mean_prediction, ymin = lower, ymax = upper), alpha = 0.3, fill = "grey") + 
  theme_bw() +
  scale_color_manual(values = colors) + 
  ylim(-25,180) + xlim(-20,20) +
  labs(x = expression("Z"[2]), y = "Mean Prediction", color = "") +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
  #ylim(min(result_sspca_aggr$lower), max(result_sspca_aggr$upper))
p4 = ggMarginal(p4, type = "histogram", margins = "x")


# plot sspca and spca next to each other
grid.arrange(p3, p4, nrow = 1)




# TOTALVIS plots 
#path = "results/simulation_results/sim2latent.totalvis/results/"

list_df = data_prep_totalvis(vec, path)
df_1 = list_df[[1]]$df
df_2 = list_df[[2]]$df


# first factor/PC:

df = df_1
pd_true = list_pca_f1[[2]]

# linear - estimated Totalvis value
formula = avg_pred ~ x_vals
q = list_sspca_f1[[3]]
mod = lm(formula = formula, data = df)
avg_pred = pred = mod$coefficients[1] + mod$coefficients[2] %*% q
data = data.frame(z = q, pred = t(avg_pred))

colors = c("Estimated Totalvis (PCA)" = "red", "Ground Truth" = "blue")
p5 = ggplot(df, aes(x = x_vals, y = avg_pred)) + geom_point(alpha = 0) + geom_line(aes(group = rep))  + 
  geom_line(data = data , aes(x = z, y = pred, color = "Estimated Totalvis (PCA)"), lty = 1, lwd = 1.2) +
  geom_line(data = pd_true, aes(x = quantile, y = mean_prediction, color = "Ground Truth"), lty= 2, lwd = 1.2) +
  theme_bw() + scale_color_manual(values = colors) +
  ylim(0,40) + xlim(-25,25) +
  labs(x = expression("Z"[1]), y = "Mean Prediction", color = "") +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
p5 = ggMarginal(p5, type = "histogram", margins = "x")



# second factor/PC:
df = df_2
pd_true = list_pca_f2[[2]]

# linear - estimated Totalvis value
formula = avg_pred ~ poly(x_vals, 2, raw = TRUE)
q = list_sspca_f2[[3]]
mod = lm(formula = formula, data = df)
avg_pred = mod$coefficients[1] + mod$coefficients[2] %*% q + mod$coefficients[3] %*% q^2
data = data.frame(z = q, pred = t(avg_pred))


colors = c("Estimated Totalvis (PCA)" = "red", "Ground Truth" = "blue")

p6 = ggplot(df, aes(x = x_vals, y = avg_pred)) + geom_point(alpha = 0) + geom_line(aes(group = rep))  + #xlim(min(df$x_vals), max(df$x_vals)) +
  geom_line(data = data , aes(x = z, y = pred, color = "Estimated Totalvis (PCA)"), lty = 1, lwd = 1.2) +
  geom_line(data = pd_true, aes(x = quantile, y = mean_prediction, color = "Ground Truth"), lty= 2, lwd = 1.2) +
  theme_bw() + scale_color_manual(values = colors) + 
  ylim(-25,180) + xlim(-20,20) +
  labs(x = expression("Z"[2]), y = "Mean Prediction", color = "") +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
p6 = ggMarginal(p6, type = "histogram", margins = "x")


# create and save Figure 9
p = gridExtra::grid.arrange(p1, p2, p5, p3, p4, p6, nrow = 2)

ggsave("results/figures/sim_2factors.png",plot = p,  width = 9, height = 7)

