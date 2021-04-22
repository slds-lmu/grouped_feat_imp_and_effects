# This script creates Figure 1 and Table 1 (section 3.1)
source("code/simulations/helper.r")
library(dplyr)
theme_set(theme_bw())

gfs_sim = readRDS("results/simulation_results/gfs_sim.RDS")


df = helper_gfs(gfs_sim$gfs_refit)
df_alluv = data.frame(first = df[[1]]$features, second = df[[2]]$features, third = df[[3]]$features)

links = data.frame(
  source = df_alluv$first,
  target = df_alluv$second
)
links = data.frame(links %>% group_by(source, target) %>% dplyr::summarize(value = n()))
links$value[is.na(links$target)] = 0

links2 = data.frame(
  source = df_alluv$second,
  target = df_alluv$third
)
links2 = data.frame(links2 %>% group_by(source, target) %>% summarize(value = n()))
links2$value[is.na(links2$target)] = 0

links = bind_rows(links, links2)


links$source = as.character(links$source)
links$target = as.character(links$target)
links$source[links$source == "G1"] = paste0("G1, MSE = ", df[[1]] %>% dplyr::filter(features == "G1") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[1]] %>% dplyr::filter(features == "G1") %>% nrow())
links$source[links$source == "G2"] = paste0("G2, MSE = ", df[[1]] %>% dplyr::filter(features == "G2") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[1]] %>% dplyr::filter(features == "G2") %>% nrow())
links$target[links$target == "G1.G3"] = paste0("G1.G3, MSE = ", df[[2]] %>% dplyr::filter(features == "G1.G3") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[2]] %>% dplyr::filter(features == "G1.G3") %>% nrow())
links$target[links$target == "G3.G2"] = paste0("G2.G3, MSE = ", df[[2]] %>% dplyr::filter(features == "G3.G2") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[2]] %>% dplyr::filter(features == "G3.G2") %>% nrow())
links$source[links$source == "G3.G2"] = paste0("G2.G3, MSE = ", df[[2]] %>% dplyr::filter(features == "G3.G2") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[2]] %>% dplyr::filter(features == "G3.G2") %>% nrow())
links$source[links$source == "G1.G3"] = paste0("G1.G3, MSE = ", df[[2]] %>% dplyr::filter(features == "G1.G3") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[2]] %>% dplyr::filter(features == "G1.G3") %>% nrow())
links$target[links$target == "G1.G3.G2"] = paste0("G1.G2.G3, MSE = ", df[[3]] %>% dplyr::filter(features == "G1.G3.G2") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[3]] %>% dplyr::filter(features == "G1.G3.G2") %>% nrow())


nodes = data.frame(
  name = c(as.character(links$source), as.character(links$target)) %>% unique()
)
links$IDsource = match(links$source, nodes$name) - 1
links$IDtarget = match(links$target, nodes$name) - 1
my_color = 'd3.scaleOrdinal() .domain([""]) .range(["grey"])'

# Make the Network
library(networkD3)
p = sankeyNetwork(Links = links, Nodes = nodes,
  Source = "IDsource", Target = "IDtarget",
  Value = "value", NodeID = "name", 
  colourScale = my_color,
  # LinkGroup="group",
  sinksRight=FALSE, iterations = 1, fontSize = 24, nodePadding = 120, height = 720, width = 1200, fontFamily = "sans-serif")
p
library(htmlwidgets)
saveWidget(p, "sankey_refit.html")
library(webshot)
webshot("sankey_refit.html" , "sankey_sim_refit.pdf", delay = 0.1, zoom = 0.6)


t_test = function(res1, res2) {
  n2 = mean(unlist(lapply(res1$pred$instance$test.inds, length)))
  n1 = mean(unlist(lapply(res1$pred$instance$train.inds, length)))
  aij = res1$measures.test[, 2]
  bij = res2$measures.test[, 2]
  xij = aij - bij
  m = mean(xij)
  s2 = var(xij)
  n = length(res1$pred$instance$train.inds)
  t_val = m / sqrt((1 / n + (n2 / n1)) * s2)
  p_value = 2 * pt(-abs(t_val), df = n - 1)
  p_value
}

#data frame with p values
# alpha = 0.05
gopfi_df = data.frame(gfs_sim$gopfi)
# gopfi_df$p_val_adj = p.adjust(gopfi_df$p_val, method = "bonferroni", n = 3 * 4)
gopfi_df$GOPFI = round(gopfi_df$GOPFI, 2)
# gopfi_df$GOPFI = paste0(round(gopfi_df$GOPFI, 2), ifelse(gopfi_df$p_val_adj <= alpha, "*", ""))
gopfi_df = gopfi_df %>% dplyr::select(features, GOPFI) %>% dplyr::filter(features != "all")
colnames(gopfi_df) = c("Group", "GOPFI")

gpfi_df = data.frame(gfs_sim$gpfi)
# gpfi_df$p_val_adj = p.adjust(gpfi_df$p_val, method = "bonferroni", n = 3 * 4)
# gpfi_df$mse = paste0(round(gpfi_df$mse, 2), ifelse(gpfi_df$p_val_adj <= alpha, "*", ""))
gpfi_df$mse = round(gpfi_df$mse, 2)
gpfi_df = gpfi_df %>% dplyr::select(features, mse) %>% dplyr::filter(features != "all")
colnames(gpfi_df) = c("Group", "GPFI")

goi = gfs_sim$goi
goi_df = data.frame(Group = c("G1", "G2", "G3"), GOI = NA, p_val = NA)
goi_df$GOI[1] = -(goi$G1$aggr - goi$featureless$aggr)
goi_df$GOI[2] = -(goi$G2$aggr - goi$featureless$aggr)
goi_df$GOI[3] = -(goi$G3$aggr - goi$featureless$aggr)
# goi_df$p_val[1] = t_test(goi$G1, goi$featureless)
# goi_df$p_val[2] = t_test(goi$G2, goi$featureless)
# goi_df$p_val[3] = t_test(goi$G3, goi$featureless)
# goi_df$p_val_adj = p.adjust(goi_df$p_val, method = "bonferroni", n = 3 * 4)

# goi_df$GOI = paste0(round(goi_df$GOI, 2), ifelse(goi_df$p_val_adj <= alpha, "*", ""))
goi_df$GOI = round(goi_df$GOI, 2)
goi_df = goi_df %>% dplyr::select(Group, GOI)


dgi = gfs_sim$dgi
dgi_df = data.frame(Group = c("G1", "G2", "G3"), DGI = NA, p_val = NA)
dgi_df$DGI[1] = -(dgi$all$aggr - dgi$G1$aggr)
dgi_df$DGI[2] = -(dgi$all$aggr - dgi$G2$aggr)
dgi_df$DGI[3] = -(dgi$all$aggr - dgi$G3$aggr)
# dgi_df$p_val[1] = t_test(dgi$all, dgi$G1)
# dgi_df$p_val[2] = t_test(dgi$all, dgi$G2)
# dgi_df$p_val[3] = t_test(dgi$all, dgi$G3)
# dgi_df$p_val_adj = p.adjust(dgi_df$p_val, method = "bonferroni", n = 3 * 4)

# dgi_df$DGI = paste0(round(dgi_df$DGI, 2), ifelse(dgi_df$p_val_adj <= alpha, "*", ""))
dgi_df$DGI = round(dgi_df$DGI, 2)
dgi_df = dgi_df %>% dplyr::select(Group, DGI)

fi_table = gopfi_df %>% left_join(gpfi_df) %>% left_join(goi_df) %>% left_join(dgi_df) %>% left_join(gfs_sim$shap, by = c("Group" = "group"))

library(xtable)
xtable(fi_table[c(3,2,1), ])
