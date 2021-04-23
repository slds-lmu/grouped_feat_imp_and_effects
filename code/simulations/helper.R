library(ranger)
library(dplyr)
# helper functions for simulations on varying correlations within groups

helper_gfs = function(gfs) {
  links = setNames(data.frame(matrix(ncol = 3)), c("source", "target", "value"))
  
  z = list()  
  y = max(unlist(lapply(gfs, function(x) length(x$B_results))))
  for (j in 1:y) z[[j]] = data.frame()
  
  gfs2 = gfs
  for (i in 1:length(gfs2)) {
    if (length(gfs2[[i]]$B_results) == 0) next
    for (j in 1:y) {
      if (j > length(gfs2[[i]]$B_results)) {z[[j]] = bind_rows(z[[j]], data.frame(features = NA, mse = NA)); next}
      if (is.numeric(gfs2[[i]]$B_results[[j]])) gfs2[[i]]$B_results[[j]] = data.frame(features = names(gfs2[[i]]$B_results[[j]]), mse = as.numeric(gfs2[[i]]$B_results[[j]]))
      z[[j]] = bind_rows(z[[j]], gfs2[[i]]$B_results[[j]])
      rownames(z[[j]]) = NULL
    }
  }
  z
}


get_gfi = function(vec, path){
  for (i in vec) {
    if (i == 1) {
      list = readRDS(paste0(path,i,".rds"))
      gpfi = list$gpfi
      gopfi = list$gopfi
      dgi_list = list$dgi
      dgi = data.frame("features" = c("G1", "G2", "G3", "G4"), "rmse" = c(dgi_list$G1$aggr-dgi_list$all$aggr, dgi_list$G2$aggr-dgi_list$all$aggr, dgi_list$G3$aggr-dgi_list$all$aggr, dgi_list$G4$aggr-dgi_list$all$aggr))
      goi_list = list$goi
      goi = data.frame("features" = c("G1", "G2", "G3", "G4"), "rmse" = c(goi_list$G1$aggr-goi_list$featureless$aggr, goi_list$G2$aggr-goi_list$featureless$aggr, goi_list$G3$aggr-goi_list$featureless$aggr, goi_list$G4$aggr-goi_list$featureless$aggr))
      shap = list$shap_imp
      perf = list$perf
      
    }
    else {
      list = readRDS(paste0(path,i,".rds"))
      gpfi = rbind.data.frame(gpfi, list$gpfi)
      gopfi = rbind.data.frame(gopfi, list$gopfi)
      dgi_list = list$dgi
      dgi = rbind.data.frame(dgi, data.frame("features" = c("G1", "G2", "G3", "G4"), "rmse" = c(dgi_list$G1$aggr-dgi_list$all$aggr, dgi_list$G2$aggr-dgi_list$all$aggr, dgi_list$G3$aggr-dgi_list$all$aggr, dgi_list$G4$aggr-dgi_list$all$aggr)))
      goi_list = list$goi
      goi = rbind.data.frame(goi, data.frame("features" = c("G1", "G2", "G3", "G4"), "rmse" = c(goi_list$G1$aggr-goi_list$featureless$aggr, goi_list$G2$aggr-goi_list$featureless$aggr, goi_list$G3$aggr-goi_list$featureless$aggr, goi_list$G4$aggr-goi_list$featureless$aggr)))
      shap = rbind.data.frame(shap, list$shap_imp)
      perf = c(perf, list$perf)
    }
    
  }
  # add run and method columns
  gpfi$run = rep(vec, each = 4)
  gpfi$method = "GPFI"
  gopfi = gopfi[!is.na(gopfi$GOPFI),]
  gopfi = as.data.frame(gopfi)[,which(colnames(gopfi)!="rmse")]
  gopfi$run = rep(vec, each = 4)
  gopfi$method = "GOPFI"
  colnames(gopfi)[2] = "rmse"
  dgi$run = rep(vec, each = 4)
  dgi$method = "LOGO"
  goi$run = rep(vec, each = 4)
  goi$method = "LOGI"
  
  # GSI
  colnames(shap) = c("G1", "G2", "G3", "G4")
  library(reshape2)
  shap$run = vec
  shap = melt(shap, id.vars = c("run"), variable.name = "features", value.name = "rmse")
  shap$method = "GSI"
  
  
  gpfi$rel_imp = gpfi$rmse/rep(gpfi$rmse[which(gpfi$features=="G1")], each = 4)
  gopfi$rel_imp = gopfi$rmse/rep(gopfi$rmse[which(gopfi$features=="G1")], each = 4)
  dgi$rel_imp = dgi$rmse/rep(dgi$rmse[which(dgi$features=="G1")], each = 4)
  goi$rel_imp = goi$rmse/rep(goi$rmse[which(goi$features=="G1")], each = 4)
  shap$rel_imp = shap$rmse/rep(shap$rmse[which(shap$features=="G1")], each = 4)
  gfi = rbind.data.frame(gpfi, gopfi, dgi, goi, shap)
  
  return(gfi)
}


get_split_features = function(vec, path){
  for (i in vec) {
    
    list = readRDS(paste0(path, i,".rds"))
    res = list$res
    models = res$models
    df_split = lapply(1:length(models), function(mod){
      
      for(t in 1:models[[mod]]$learner.model$num.trees){
        
        tree_first_split = treeInfo(models[[mod]]$learner.model, tree = t)[1,]
        tree_first_split$tree = t
        tree_first_split$resample = mod
        tree_first_split$repetition = i
        if(t == 1) df_first_split = tree_first_split
        else df_first_split = rbind(df_first_split, tree_first_split)
      }
      df_first_split
    })
    if(i ==1) df = do.call("rbind", df_split)
    else df = rbind(df, do.call("rbind", df_split))
    
  }
  
  df$splitvarName = as.character(df$splitvarName)
  df$splitvarName <- factor(df$splitvarName, levels = paste0("V", 1:40))
  df$group = NA
  df$group[df$splitvarName %in% paste0("V", 1:10)] = "G1"
  df$group[df$splitvarName %in% paste0("V", 11:20)] = "G2"
  df$group[df$splitvarName %in% paste0("V", 21:30)] = "G3"
  df$group[df$splitvarName %in% paste0("V", 31:40)] = "G4"
  return(df)
}



get_shapley_imp = function(vec, path){
  for (i in vec) {
    list = readRDS(paste0(path, i,".rds"))
    if (i == 1) {
      shap_feat = list$shapleyValue
      shap_group = data.frame("group" = c("G1", "G2"), "mse" = c(as.numeric(list$impG1), as.numeric(list$impG2)))
      
    }
    else {
      shap_feat = rbind(shap_feat, list$shapleyValue)
      shap_group = rbind(shap_group, data.frame("group" = c("G1", "G2"), "mse" = c(as.numeric(list$impG1), as.numeric(list$impG2))))
    }
    
  }
  shap_feat$group = rep(c(rep("G1",6), rep("G2",2)), length(vec))
  
  return(list(shap_feat, shap_group))
}  
 


model_sim <- function(df, formula, q) {
  for (i in 1:length(unique(df$sim_no))) {
    data = df[which(df$sim_no == i),]
    model = lm(formula = formula, data = data)
    coef = model$coefficients
    if (length(coef) == 2) {
      pred = coef[1] + coef[2] %*% q # vector der länge q
    }
    else if (length(coef) == 3) {
      pred = coef[1] + coef[2] %*% q + coef[3] %*% q^2
    }
    
    if (i == 1) {
      result = data.frame(matrix(pred, nrow = 1))
    }
    else result = rbind.data.frame(result, data.frame(matrix(pred, nrow = 1)))
    
  }
  return(result)
}

# totalvis helpers
data_prep_totalvis = function(vec, path){
  for (i in vec) {
    list = readRDS(paste0(path, i,".rds"))
    if(names(list)[1] == "pred_df") list = list(list)
    if (i == 1) {
      
      for(j in 1:length(list)){
        df = as.data.frame(list[[j]]$pred_df)
        df$rep = i
        weights = as.data.frame(list[[j]]$pca_object$rotation)
        weights$rep = i
        assign(paste0("df_", j), df)
        assign(paste0("weights_", j), weights)
      }
    }
    else {
      for(j in 1:length(list)){
        
        df = as.data.frame(list[[j]]$pred_df)
        df$rep = i
        df = rbind(get(paste0("df_", j)), df)
        weights = as.data.frame(list[[j]]$pca_object$rotation)
        weights$rep = i
        weights = rbind(get(paste0("weights_", j)), weights)
        assign(paste0("df_", j), df)
        assign(paste0("weights_", j), weights)
      }
    }
    
  }
  list_df = list()
  for(j in 1:length(list)){
    list_df[[j]] = NA
    list_df[[j]]$df = get(paste0("df_", j))
    list_df[[j]]$weights = get(paste0("weights_", j))
  }
  
  return(list_df)
}


# CFEP

get_data_cfep = function(vec, path){
  for (i in vec) {
    if (i == 1) {
      list_1 = readRDS(paste0(path, i,".rds"))
      #df_sspca = list_1$df_sspca
      #df_pca = list_1$df_pca
      df_detailed_sspca = list_1$df_detailed_sspca
      df_detailed_pca = list_1$df_detailed_pca
      df_features_sspca = list_1$df_features_sspca
      df_features_pca = list_1$df_features_pca
    }
    else {
      list = readRDS(paste0(path, i,".rds"))
      #df_sspca = rbind.data.frame(df_sspca, list$df_sspca)
      #df_pca = rbind.data.frame(df_pca, list$df_pca)
      df_detailed_sspca = rbind.data.frame(df_detailed_sspca, list$df_detailed_sspca)
      df_detailed_pca = rbind.data.frame(df_detailed_pca, list$df_detailed_pca)
      df_features_sspca = rbind.data.frame(df_features_sspca, list$df_features_sspca)
      df_features_pca = rbind.data.frame(df_features_pca, list$df_features_pca)
    }
    
  }
  df_detailed_sspca$sim_no = as.factor(rep(1:length(vec), each = 500))
  df_detailed_pca$sim_no = as.factor(rep(1:length(vec), each = 500))
  
  return(list(df_detailed_sspca = df_detailed_sspca, df_detailed_pca = df_detailed_pca, 
              df_features_sspca = df_features_sspca, df_features_pca = df_features_pca))
}



data_prep_cfep = function(df, z, f, f_true){
  q = round(seq(from = min(df[,z]), to = max(df[,z]), length.out = 50),1)
  result = model_sim(df, f, q)
  colnames(result)[1:(ncol(result))] = q
  result = gather(result, quantile, mean_prediction, 1:(ncol(result)))
  result$quantile = as.factor(as.numeric(result$quantile))
  sd = c()
  for (i in result$quantile) {
    sd[i] = sd(result$mean_prediction[which(result$quantile == i)])
    
  }
  result_aggr = aggregate(mean_prediction ~ quantile, result, mean)
  result_aggr$upper = result_aggr$mean_prediction + 1.96*sd
  result_aggr$lower = result_aggr$mean_prediction - 1.96*sd
  result_aggr$quantile = as.numeric(as.character(result_aggr$quantile))
  
  # data prep for ground truth
  pd_true = model_sim(df, f_true, q)
  colnames(pd_true)[1:(ncol(pd_true))] = q
  pd_true = gather(pd_true, quantile, mean_prediction, 1:(ncol(pd_true)))
  pd_true$quantile = as.numeric(pd_true$quantile)
  pd_true = aggregate(mean_prediction ~ quantile, pd_true, mean)
  
  return(list(result_aggr, pd_true, q))
}

