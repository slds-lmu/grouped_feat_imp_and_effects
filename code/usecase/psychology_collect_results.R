psi = readRDS("results/usecase_results/psychology_feat_imp.RDS")

#feat imp
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
alpha = 0.05
gopfi_df = data.frame(psi$gopfi)
gopfi_df$p_val_adj = p.adjust(gopfi_df$p_val, method = "bonferroni", n = 5 * 4)
gopfi_df$GOPFI = paste0(round(gopfi_df$GOPFI, 3), ifelse(gopfi_df$p_val_adj <= alpha, "*", ""))
gopfi_df = gopfi_df %>% dplyr::select(features, GOPFI) %>% dplyr::filter(features != "all")
colnames(gopfi_df) = c("Group", "GOPFI")

gpfi_df = data.frame(psi$gpfi)
gpfi_df$p_val_adj = p.adjust(gpfi_df$p_val, method = "bonferroni", n = 5 * 4)
gpfi_df$mse = paste0(round(gpfi_df$mse, 3), ifelse(gpfi_df$p_val_adj <= alpha, "*", ""))
gpfi_df = gpfi_df %>% dplyr::select(features, mse) %>% dplyr::filter(features != "all")
colnames(gpfi_df) = c("Group", "GPFI")

goi = psi$goi
goi_df = data.frame(Group = names(goi)[-6], GOI = NA, p_val = NA)
goi_df$GOI[1] = goi$app$aggr - goi$featureless$aggr
goi_df$GOI[2] = goi$communication_social$aggr - goi$featureless$aggr
goi_df$GOI[3] = goi$mobility$aggr - goi$featureless$aggr
goi_df$GOI[4] = goi$music$aggr - goi$featureless$aggr
goi_df$GOI[5] = goi$general$aggr - goi$featureless$aggr
goi_df$GOI = -1 * goi_df$GOI
goi_df$p_val[1] = t_test(goi$featureless, goi$app)
goi_df$p_val[2] = t_test(goi$featureless, goi$communication_social)
goi_df$p_val[3] = t_test(goi$featureless, goi$mobility)
goi_df$p_val[4] = t_test(goi$featureless, goi$music)
goi_df$p_val[5] = t_test(goi$featureless, goi$general)
goi_df$p_val_adj = p.adjust(goi_df$p_val, method = "bonferroni", n = 5 * 4)

goi_df$GOI = paste0(round(goi_df$GOI, 3), ifelse(goi_df$p_val_adj <= alpha, "*", ""))
goi_df = goi_df %>% dplyr::select(Group, GOI)


dgi = psi$dgi
dgi_df = data.frame(Group = names(dgi)[-6], DGI = NA, p_val = NA)
dgi_df$DGI[1] = dgi$all$aggr - dgi$app$aggr
dgi_df$DGI[2] = dgi$all$aggr - dgi$communication_social$aggr
dgi_df$DGI[3] = dgi$all$aggr - dgi$mobility$aggr
dgi_df$DGI[4] = dgi$all$aggr - dgi$music$aggr
dgi_df$DGI[5] = dgi$all$aggr - dgi$general$aggr
dgi_df$DGI = dgi_df$DGI * (-1)

dgi_df$p_val[1] = t_test(dgi$app, dgi$all)
dgi_df$p_val[2] = t_test(dgi$communication_social, dgi$all)
dgi_df$p_val[3] = t_test(dgi$mobility, dgi$all)
dgi_df$p_val[4] = t_test(dgi$music, dgi$all)
dgi_df$p_val[5] = t_test(dgi$general, dgi$all)
dgi_df$p_val_adj = p.adjust(dgi_df$p_val, method = "bonferroni", n = 5 * 4)

dgi_df$DGI = paste0(round(dgi_df$DGI, 3), ifelse(dgi_df$p_val_adj <= alpha, "*", ""))
dgi_df = dgi_df %>% dplyr::select(Group, DGI)

fi_table = gopfi_df %>% left_join(gpfi_df) %>% left_join(goi_df) %>% left_join(dgi_df)

fi_table = fi_table %>% left_join(shap, by = c("Group" = "group"))
library(xtable)
fi_table$meanShapImp = round(fi_table$meanShapImp, 3)
fi_table$meanShapImp = as.character(fi_table$meanShapImp)
xtable(fi_table)



#sankey
library(dplyr)

helper_gfs2 = function(gfs) {
  links = setNames(data.frame(matrix(ncol = 3)), c("source", "target", "value"))
  
  z = list()  
  y = max(unlist(lapply(gfs, function(x) length(x[[1]]$B_results))))
  for (j in 1:y) z[[j]] = data.frame()
  
  gfs2 = gfs
  for (i in 1:length(gfs2)) {
    if (length(gfs2[[i]][[1]]$B_results) == 0) next
    for (j in 1:y) {
      if (j > length(gfs2[[i]][[1]]$B_results)) {z[[j]] = bind_rows(z[[j]], data.frame(features = NA, mse = NA)); next}
      if (is.numeric(gfs2[[i]][[1]]$B_results[[j]])) gfs2[[i]][[1]]$B_results[[j]] = data.frame(features = names(gfs2[[i]][[1]]$B_results[[j]]), mse = as.numeric(gfs2[[i]][[1]]$B_results[[j]]))
      z[[j]] = bind_rows(z[[j]], gfs2[[i]][[1]]$B_results[[j]])
      rownames(z[[j]]) = NULL
    }
  }
  z
}
library(batchtools)
reg = loadRegistry("usecase_gfs_refit", work.dir = ".")
gfs = reduceResultsList()
#make plot for permutation method
df = helper_gfs2(gfs)
df[[3]][which(df[[3]]$features == "app.communication_social.general"), ] = NA
df[[3]][which(df[[3]]$features == "app.mobility.general"), ] = NA
df[[2]][which(df[[2]]$features == "app.communication_social"), ] = NA
df[[2]][which(df[[2]]$features == "app.music"), ] = NA
df[[2]][which(df[[2]]$features == "music.general"), ] = NA

df_alluv = data.frame(first = df[[1]]$features, second = df[[2]]$features, third = df[[3]]$features, fourth = df[[4]]$features)
df_alluv$first = gsub("communication_social", paste0("C"), df_alluv$first)
df_alluv$first = gsub("app", paste0("A"), df_alluv$first)
df_alluv$first = gsub("music", paste0("Mu"), df_alluv$first)

df_alluv$second = as.character(df_alluv$second)
df_alluv$second = gsub("app.music", paste0("A.Mu"), df_alluv$second)
df_alluv$second = gsub("app.general", paste0("A.G"), df_alluv$second)
df_alluv$second = gsub("app.communication_social", paste0("A.C"), df_alluv$second)
df_alluv$second = gsub("communication_social.music", paste0("C.Mu"), df_alluv$second)

df_alluv$third = as.character(df_alluv$third)
df_alluv$third = gsub("app.music.general", paste0("A.Mu.G"), df_alluv$third)
df_alluv$third = gsub("app.mobility.music", paste0("A.Mo.Mu"), df_alluv$third)
df_alluv$third = gsub("app.communication_social.music", paste0("A.C.Mu"), df_alluv$third)
df_alluv$third = gsub("communication_social.music.general", paste0("C.Mu.G"), df_alluv$third)
df_alluv$third = gsub("app.communication_social.mobility", paste0("A.C.Mo"), df_alluv$third)
df_alluv$third = gsub("app.communication_social.general", paste0("A.C.G"), df_alluv$third)


df_alluv$fourth = as.character(df_alluv$fourth)
df_alluv$fourth = gsub("app.communication_social.music.general", paste0("A.C.Mu.G"), df_alluv$fourth)
df_alluv$fourth = gsub("app.communication_social.mobility.music", paste0("A.C.Mo.Mu"), df_alluv$fourth)

links = data.frame(
  source = df_alluv$first,
  target = df_alluv$second
)
links = data.frame(links %>% group_by(source, target) %>% dplyr::summarize(value = n()))


links2 = data.frame(
  source = df_alluv$second,
  target = df_alluv$third
)
links2 = data.frame(links2 %>% group_by(source, target) %>% summarize(value = n()))
links2$value[is.na(links2$target)] = 0
links2 = links2[-which((is.na(links2$source) * is.na(links2$target)) == 1), ]
links = bind_rows(links, links2)

#reduce links to min size 5 at ending
# links = links %>% filter(target != "A.C.Mo.Mu.G")
# links = links %>% filter(target != "A.C.Mu.G")
# links = links %>% filter(target != "A.Mo.Mu")
# links = links %>% filter(target != "C.Mu.G")
# links = links %>% filter(target != "mobility.general")
links = links[-which(links$target == "communication_social.general"), ]
links = links[-which(links$target == "mobility.general"), ]
links = links[-which(links$target == "C.Mu.G"), ]

links = links[c(1, 2, 3, 4, 5), ]


links$source = as.character(links$source)
links$target = as.character(links$target)
links$source[links$source == "A"] = paste0("A, MSE = ", df[[1]] %>% dplyr::filter(features == "app") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[1]] %>% dplyr::filter(features == "app") %>% nrow())
links$source[links$source == "general"] = paste0("G, MSE = ", df[[1]] %>% dplyr::filter(features == "general") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[1]] %>% dplyr::filter(features == "general") %>% nrow())
links$source[links$source == "C"] = paste0("C, MSE = ", df[[1]] %>% dplyr::filter(features == "communication_social") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[1]] %>% dplyr::filter(features == "communication_social") %>% nrow())
links$source[links$source == "Mu"] = paste0("Mu, MSE = ", df[[1]] %>% dplyr::filter(features == "music") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[1]] %>% dplyr::filter(features == "music") %>% nrow())
links$target[links$target == "A.C"] = paste0("A.C, MSE = ", df[[2]] %>% dplyr::filter(features == "app.communication_social") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[2]] %>% dplyr::filter(features == "app.communication_social") %>% nrow())
links$target[links$target == "A.G"] = paste0("A.G, MSE = ", df[[2]] %>% dplyr::filter(features == "app.general") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[2]] %>% dplyr::filter(features == "app.general") %>% nrow())
links$source[links$source == "A.G"] = paste0("A.G, MSE = ", df[[2]] %>% dplyr::filter(features == "app.general") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[2]] %>% dplyr::filter(features == "app.general") %>% nrow())
links$source[links$source == "A.C"] = paste0("A.C, MSE = ", df[[2]] %>% dplyr::filter(features == "app.communication_social") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[2]] %>% dplyr::filter(features == "app.communication_social") %>% nrow())
links$source[links$source == "A.Mu"] = paste0("A.Mu, MSE = ", df[[2]] %>% dplyr::filter(features == "app.music") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[2]] %>% dplyr::filter(features == "app.music") %>% nrow())
links$target[links$target == "A.Mu"] = paste0("A.Mu, MSE = ", df[[2]] %>% dplyr::filter(features == "app.music") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[2]] %>% dplyr::filter(features == "app.music") %>% nrow())
links$source[links$source == "C.Mu"] = paste0("C.Mu, MSE = ", df[[2]] %>% dplyr::filter(features == "communication_social.music") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[2]] %>% dplyr::filter(features == "communication_social.music") %>% nrow())
links$target[links$target == "C.Mu"] = paste0("C.Mu, MSE = ", df[[2]] %>% dplyr::filter(features == "communication_social.music") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[2]] %>% dplyr::filter(features == "communication_social.music") %>% nrow())
links$source[links$source == "music.general"] = paste0("Mu.G, MSE = ", df[[2]] %>% dplyr::filter(features == "music.general") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[2]] %>% dplyr::filter(features == "music.general") %>% nrow())
links$target[links$target == "music.general"] = paste0("Mu.G, MSE = ", df[[2]] %>% dplyr::filter(features == "music.general") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[2]] %>% dplyr::filter(features == "music.general") %>% nrow())
links$target[links$target == "A.C.G"] = paste0("A.C.G, MSE = ", df[[3]] %>% dplyr::filter(features == "app.communication_social.general") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[3]] %>% dplyr::filter(features == "app.communication_social.general") %>% nrow())
links$source[links$source == "A.C.G"] = paste0("A.C.G, MSE = ", df[[3]] %>% dplyr::filter(features == "app.communication_social.general") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[3]] %>% dplyr::filter(features == "app.communication_social.general") %>% nrow())
links$source[links$source == "A.C.Mo"] = paste0("A.C.Mo, MSE = ", df[[3]] %>% dplyr::filter(features == "app.communication_social.mobility") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[3]] %>% dplyr::filter(features == "app.communication_social.mobility") %>% nrow())
links$target[links$target == "A.C.Mo"] = paste0("A.C.Mo, MSE = ", df[[3]] %>% dplyr::filter(features == "app.communication_social.mobility") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[3]] %>% dplyr::filter(features == "app.communication_social.mobility") %>% nrow())
links$target[links$target == "A.C.Mu"] = paste0("A.C.Mu, MSE = ", df[[3]] %>% dplyr::filter(features == "app.communication_social.music") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[3]] %>% dplyr::filter(features == "app.communication_social.music") %>% nrow())
links$source[links$source == "A.C.Mu"] = paste0("A.C.Mu, MSE = ", df[[3]] %>% dplyr::filter(features == "app.communication_social.music") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[3]] %>% dplyr::filter(features == "app.communication_social.music") %>% nrow())
links$source[links$source == "A.Mu.G"] = paste0("A.Mu.G, MSE = ", df[[3]] %>% dplyr::filter(features == "app.music.general") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[3]] %>% dplyr::filter(features == "app.music.general") %>% nrow())
links$target[links$target == "A.Mu.G"] = paste0("A.Mu.G, MSE = ", df[[3]] %>% dplyr::filter(features == "app.music.general") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[3]] %>% dplyr::filter(features == "app.music.general") %>% nrow())
links$target[links$target == "C.Mu.G"] = paste0("C.Mu.G, MSE = ", df[[3]] %>% dplyr::filter(features == "communication_social.music.general") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[3]] %>% dplyr::filter(features == "communication_social.music.general") %>% nrow())
links$source[links$source == "C.Mu.G"] = paste0("C.Mu.G, MSE = ", df[[3]] %>% dplyr::filter(features == "communication_social.music.general") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[3]] %>% dplyr::filter(features == "communication_social.music.general") %>% nrow())
links$target[links$target == "app.mobility.general"] = paste0("A.Mo.G, MSE = ", df[[3]] %>% dplyr::filter(features == "app.mobility.general") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[3]] %>% dplyr::filter(features == "app.mobility.general") %>% nrow())
links$source[links$source == "app.mobility.general"] = paste0("A.Mo.G, MSE = ", df[[3]] %>% dplyr::filter(features == "app.mobility.general") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[3]] %>% dplyr::filter(features == "app.mobility.general") %>% nrow())
links$target[links$target == "A.C.Mu.G"] = paste0("A.C.Mu.G, MSE = ", df[[4]] %>% dplyr::filter(features == "app.communication_social.music.general") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[4]] %>% dplyr::filter(features == "app.communication_social.music.general") %>% nrow())
links$target[links$target == "A.C.Mo.Mu"] = paste0("A.C.Mo.Mu, MSE = ", df[[4]] %>% dplyr::filter(features == "app.communication_social.mobility.music") %>% pull(mse) %>% mean() %>% round(3), ", n = ", df[[4]] %>% dplyr::filter(features == "app.communication_social.mobility.music") %>% nrow())

links = bind_rows(links, data.frame(source = NA, target = "A, MSE = 0.519, n = 82", value = 82))
links = bind_rows(links, data.frame(source = NA, target = "G, MSE = 0.526, n = 18", value = 18))
links = links[c(-2, -4), ]

links$source[2] = "O, MSE = 0.526, n = 18"
links$source[3] = "A.O, MSE = 0.513, n = 73"
links$target[1] = "A.O, MSE = 0.513, n = 73"
links$target[2] = "A.O, MSE = 0.513, n = 73"
links$target[5] = "O, MSE = 0.526, n = 18"
# mean(df[[3]][df[[3]]$features == "app.music.general" & df[[2]]$features == "app.general", ]$mse, na.rm = TRUE)
# 0.5078229
# sum(!is.na(df[[3]][df[[3]]$features == "app.music.general" & df[[2]]$features == "app.general", ]$mse))
links$target[3] = "A.O.Mu, MSE = 0.508, n = 9"

# links = bind_rows(links, data.frame(source = paste0("Empty, n = ", 100 - nrow(df[[1]])), value = 23))
nodes = data.frame(
  name = c(as.character(links$source), as.character(links$target)) %>% unique()
)
links$IDsource = match(links$source, nodes$name) - 1
links$IDtarget = match(links$target, nodes$name) - 1
links$group =  as.factor(c("type_c","type_c","type_c","type_b","type_b"))
nodes$group = as.factor(c("type_a","type_a","type_a","type_b","type_a"))
my_color <- 'd3.scaleOrdinal() .domain(["type_a", "type_b", "type_c"]) .range(["grey", "white", "#cccccc"])'

# my_color = 'd3.scaleOrdinal() .domain([""]) .range(["grey"])'

# Make the Network
library(networkD3)
p = sankeyNetwork(Links = links, Nodes = nodes,
                  Source = "IDsource", Target = "IDtarget",
                  Value = "value", NodeID = "name", 
                  colourScale = my_color,
                  LinkGroup="group",NodeGroup="group",
                  sinksRight=FALSE, iterations = 20, fontSize = 16, nodePadding = 250, height = 720, width = 1200, fontFamily = "sans-serif")
p
library(htmlwidgets)
saveWidget(p, "sankey_usecase.html")
library(webshot)
webshot("sankey_usecase.html" , "sankey_usecase.pdf", delay = 0.1)

