greedy_fw_search_refit = function(gimp, group_df, delta = 0.001, 
  outer_resampling = cv10, inner_resampling = cv10,
  measure = mse) {
  outer_res = makeResampleInstance(outer_resampling, task = gimp$task)
  G = unique(group_df$group)
  B_all = list()
  for (i in 1:length(outer_res$train.inds)) {
    print(round(i / length(outer_res$train.inds), 2))
    B = list()
    B_results = list()
    outer_train = getTaskData(gimp$task, subset = outer_res$train.inds[[i]])
    outer_task = makeRegrTask(data = outer_train, target = getTaskTargetNames(gimp$task))
    outer_test = getTaskData(gimp$task, subset = outer_res$test.inds[[i]])
    for (k in 1:length(G)) {
      if (k == 1) {
        B_tilde = G
        R_null = numeric(0L)
        R_G = list()
        inner_res = makeResampleInstance(inner_resampling, task = outer_task)
        for (j in 1:length(inner_res$train.inds)) {
          mod_inner = train(makeLearner("regr.featureless"), outer_task, subset = inner_res$train.inds[[j]])
          pred_inner = predict(mod_inner, task = outer_task, subset = inner_res$test.inds[[j]])
          R_null[[j]] = measure$fun(pred = pred_inner)
          for (G_i in B_tilde) {
            td = getTaskData(outer_task)
            td = td %>% dplyr::select(group_df %>% dplyr::filter(group == G_i) %>% pull(feature), getTaskTargetNames(outer_task))
            task_G_i = makeRegrTask(data = td, target = getTaskTargetNames(outer_task))
            mod_G_i = train(gimp$lrn, task = task_G_i, subset = inner_res$train.inds[[j]])
            pred_G_i = predict(mod_G_i, task = task_G_i, subset = inner_res$test.inds[[j]])
            R_G[[G_i]] = c(R_G[[G_i]],  measure$fun(pred = pred_G_i))  
          }
        }#end inner resampling
        res_inner = unlist(lapply(R_G, function(x) mean(x)))
        G_star = names(res_inner[res_inner == min(res_inner)])
        if (mean(R_null) - res_inner[G_star] > delta) {
          B = G_star
          B_results[[k]] = res_inner[res_inner == min(res_inner)]
          R_k_minus_1 = res_inner[G_star]
        }
      }#end if k == 1
      if ((k > 1) & (length(B) != 0)) {
        comb = combn(G, k, simplify = FALSE)
        B_tilde = comb[unlist(lapply(comb, function(x) all(B %in% x)))]
        R_G = list()
        for (j in 1:length(inner_res$train.inds)) {
          for (G_i in B_tilde) {
            td = getTaskData(outer_task)
            td = td %>% dplyr::select(group_df %>% dplyr::filter(group %in% G_i) %>% pull(feature), getTaskTargetNames(outer_task))
            task_G_i = makeRegrTask(data = td, target = getTaskTargetNames(outer_task))
            mod_G_i = train(gimp$lrn, task = task_G_i, subset = inner_res$train.inds[[j]])
            pred_G_i = predict(mod_G_i, task = task_G_i, subset = inner_res$test.inds[[j]])
            R_G[[paste0(G_i, collapse = ".")]] = c(R_G[[paste0(G_i, collapse = ".")]],  measure$fun(pred = pred_G_i))  
          }
        }#end inner resampling
        res_inner = unlist(lapply(R_G, function(x) mean(x)))
        G_star = ifelse(measure$minimize, names(res_inner[res_inner == min(res_inner)]), names(res_inner[res_inner == max(res_inner)]))
        if (R_k_minus_1 - res_inner[G_star] > delta) {
          B = strsplit(G_star, split = "\\.")[[1]]
          B_results[[k]] = res_inner[res_inner == min(res_inner)]
          R_k_minus_1 = res_inner[G_star]
        } else {
          break
        }
      }#end if k > 1
    }#end k in 1:length(G)
    B_all[[i]] = list()
    B_all[[i]][["B"]] = B
    B_all[[i]][["B_results"]] = B_results
    
    # td = getTaskData(gimp$task)
    # td = td %>% dplyr::select(group_df %>% dplyr::filter(group %in% B) %>% pull(feature), getTaskTargetNames(outer_task))
    # task_B = makeRegrTask(data = td, target = getTaskTargetNames(outer_task))
    # mod_B = train(gimp$lrn, task = task_B, subset = outer_res$train.inds[[i]])
    # pred_B = predict(mod_B, task = task_B, subset = outer_res$test.inds[[i]])
    # B_all[[i]][["B_outer_result"]] =  measure$fun(pred = pred_G_i)
  }#end outer resampling
  B_all
}#end function


greedy_fw_search_permutation = function(gimp, group_df, delta = 0.001, 
  outer_resampling = cv10, inner_resampling = cv10,
  measure = mse) {
  outer_res = makeResampleInstance(outer_resampling, task = gimp$task)
  G = unique(group_df$group)
  B_all = list()
  for (i in 1:length(outer_res$train.inds)) {
    print(round(i / length(outer_res$train.inds), 2))
    B = list()
    B_results = list()
    outer_train = getTaskData(gimp$task, subset = outer_res$train.inds[[i]])
    outer_task = makeRegrTask(data = outer_train, target = getTaskTargetNames(gimp$task))
    outer_test = getTaskData(gimp$task, subset = outer_res$test.inds[[i]])
    for (k in 1:length(G)) {
      if (k == 1) {
        B_tilde = G
        R_null = numeric(0L)
        R_G = list()
        inner_res = makeResampleInstance(inner_resampling, task = outer_task)
        
        res_inner = resample(gimp$lrn, task = outer_task, resampling = inner_res, models = TRUE, measures = measure)
        gfeats = list()
        groups = unique(group_df$group)
        groups = groups[!is.na(groups)]
        for (j in 1:length(groups)) {
          gfeats[[j]] = group_df %>% dplyr::filter(group %in% setdiff(unique(group_df$group), unique(group_df$group)[i])) %>% pull(feature) %>% as.character()
        }
        names(gfeats) = groups
        gfeats$all = as.character(unique(group_df$feature))
        imp = featureImportance::featureImportance(res_inner, data = getTaskData(outer_task), n.feat.perm = 10L,
          features = gfeats, measures = switch(getTaskType(outer_task), classif = mmce, regr = measure), local = FALSE, 
          importance.fun = function(permuted, unpermuted) return(permuted))
        
        imp2 = summary(imp)
        imp2$GOPFI = NA
        for (group in groups) {
          imp2$GOPFI[which(group == imp2$features)] = (imp2 %>% dplyr::filter(features == "all") %>% pull(2)) - (imp2 %>% dplyr::filter(features == group) %>% pull(2))
        }
        
        G_star = (imp2 %>% arrange(desc(GOPFI)) %>% pull(features))[1]
        if (imp2 %>% dplyr::filter(features == G_star) %>% pull(GOPFI) > delta) {
          B = G_star
          B_results[[k]] = imp2 %>% dplyr::filter(features == G_star) 
          R_k_minus_1 = imp2 %>% dplyr::filter(features == G_star) %>% pull(2) 
        }
      }#end if k == 1
      if ((k > 1) & (length(B) != 0)) {
        comb = combn(G, k, simplify = FALSE)
        B_tilde = comb[unlist(lapply(comb, function(x) all(unlist(strsplit(B, split = "\\.")) %in% x)))]
        R_G = list()
        
        gfeats = list()
        for (j in 1:length(B_tilde)) {
          gfeats[[j]] = group_df %>% dplyr::filter(group %in% setdiff(unique(group_df$group), B_tilde[[j]])) %>% pull(feature) %>% as.character()
        }
        names(gfeats) = unlist(lapply(B_tilde, function(x) paste0(x, collapse = ".")))
        
        imp = featureImportance::featureImportance(res_inner, data = getTaskData(outer_task), n.feat.perm = 10L,
          features = gfeats, measures = switch(getTaskType(outer_task), classif = mmce, regr = measure), local = FALSE, 
          importance.fun = function(permuted, unpermuted) return(permuted))
        
        imp2 = summary(imp)  
      
        G_star = (imp2 %>% arrange((.[[2]])) %>% pull(features))[1]
        if (R_k_minus_1 - (imp2 %>% dplyr::filter(features == G_star) %>% pull(2)) > delta) {
          B = G_star
          B_results[[k]] = imp2 %>% dplyr::filter(features == G_star) 
          R_k_minus_1 = imp2 %>% dplyr::filter(features == G_star) %>% pull(mse) 
        } else {
          break
        }
      }#end if k > 1
    }#end k in 1:length(G)
    B_all[[i]] = list()
    B_all[[i]][["B"]] = B
    B_all[[i]][["B_results"]] = B_results
    # 
    # if (length(B) > 0) {
    #   feats = group_df %>% dplyr::filter(group %in% setdiff(unique(group_df$group),unlist(strsplit(B, split = "\\.")))) %>% pull(feature) %>% as.character()
    #   mod = train(gimp$lrn, task = outer_task)
    #   imp = featureImportance::featureImportance(mod, data = outer_test, n.feat.perm = 10L,
    #     features = gfeats, measures = switch(getTaskType(outer_task), classif = mmce, regr = measure), local = FALSE, 
    #     importance.fun = function(permuted, unpermuted) return(permuted))
    #   
    #   B_all[[i]][["outer_results"]] = summary(imp)
    # }
  }#end outer resampling
  B_all
}#end function
