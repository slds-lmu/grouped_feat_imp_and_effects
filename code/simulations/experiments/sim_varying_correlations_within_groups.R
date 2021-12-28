

library(kernlab)
library(batchtools)
library(stringi)
library(data.table)

# parameters
algorithm = c("rf", "svm")
correlation = list("base" = c(0.1,0.1,0.1,0.1), "diff" = c(0.1,0.3,0.6,0.1))

for(algo in algorithm){
  for(j in 1:length(correlation)){
    reg = makeExperimentRegistry(file.dir = paste0("results/simulation_results/sim_corr_",algo,"_",names(correlation)[j]), source = c("code/functions/gimp.r", "code/functions/gimp_sim.r","code/functions/helper.r"), seed = 123, packages = "kernlab")
    
    
    # add problem
    for (i in 1:20) addProblem(as.character(i), i) 
    
    # define algorithm
    addAlgorithm("sim", function(data, job, instance, correlation,...) {
      
      correlation = list("base" = c(0.1,0.1,0.1,0.1), "diff" = c(0.1,0.3,0.6,0.1))
      
      
      # generate data
      p = 10
      n = 1000
      
      rv1 = rnorm(n)
      rv2 = rnorm(n)
      rv3 = rnorm(n)
      rv4 = rnorm(n)
    
      
      X_G1 = make_correlated_feats_sim(10, rv1, name = "G1_feat", prop = correlation[[j]][1])
      X_G2 = make_correlated_feats_sim(10, rv2, name = "G2_feat", prop = correlation[[j]][2])
      X_G3 = make_correlated_feats_sim(10, rv3, name = "G3_feat", prop = correlation[[j]][3])
      X_G4 = make_correlated_feats_sim(10, rv4, name = "G4_feat", prop = correlation[[j]][4])
      epsilon = rnorm(n)
      
      Xsim = as.data.frame(cbind(X_G1, X_G2, X_G3, X_G4))
      colnames(Xsim) = paste0("V",1:(4*p))
      
      Z1 = 3*(X_G1[,3])^2 - 4*(X_G1[,5]) - 6*X_G1[,7] + ifelse(mean(X_G1[,8]) > 0, 5*X_G1[,9], 0)
      Z2 = 3*(X_G2[,3])^2 - 4*(X_G2[,5]) - 6*X_G2[,7] + ifelse(mean(X_G2[,8]) > 0, 5*X_G2[,9], 0)
      Z3 = 3*(X_G3[,3])^2 - 4*(X_G3[,5]) - 6*X_G3[,7] + ifelse(mean(X_G3[,8]) > 0, 5*X_G3[,9], 0)
      y = Z1 + Z2 + Z3 + epsilon
      
      data = data.frame(cbind(Xsim,y))
      X = as.matrix(Xsim)
      colnames(X) = colnames(Xsim)
      
      
      # Model
      if(algo == "rf") learner = makeLearner("regr.ranger", par.vals = list(num.trees = 2000))
      else if (algo == "svm") learner = makeLearner("regr.ksvm")
      
      task = makeRegrTask(data = data, target = "y")
      mod = train(learner = learner, task = task)
      res = resample(learner, task, cv5, measures = rmse, models = TRUE)
      gimp = Gimp$new(task = task, res = res, mod = mod, lrn = learner)
      
      
      
      group_df = data.frame(feature = colnames(Xsim), group = c(rep("G1", 10), rep("G2", 10), rep("G3", 10), rep("G4", 10)), stringsAsFactors = FALSE)
      
      gpfi = gimp$group_permutation_feat_imp(group_df, PIMP = FALSE, s = 5, regr.measure = rmse)
      dgi = gimp$drop_group_importance(group_df, measures = rmse)
      goi = gimp$group_only_importance(group_df, measures = rmse)
      gopfi = gimp$group_only_permutation_feat_imp(group_df, PIMP = FALSE, s = 5, regr.measure = rmse)
      
      
      
      # grouped importance with shapley
      shap_imp = data.frame(matrix(NA, nrow = length(res$models), ncol = length(unique(group_df$group))))
      colnames(shap_imp) = unique(group_df$group)
      for(j in 1:length(res$models)){
        train.ind = res$pred$instance$train.inds[[j]]
        test.ind = res$pred$instance$test.inds[[j]]
        mod = res$models[[j]]
        data = getTaskData(task, subset = test.ind)
        features = getTaskFeatureNames(task)
        target = getTaskTargetNames(task)
        
        # generate all permutations of groups
        group_df = data.frame(feature = colnames(Xsim), group = c(rep("G1", 10), rep("G2", 10), rep("G3", 10), rep("G4", 10)), stringsAsFactors = FALSE)
        
        groups = unique(group_df$group)
        perm = featureImportance:::generatePermutations(groups, n.shapley.perm = 120)
        
        # generate all marginal contribution sets for each feature
        mc = lapply(groups, function(x) featureImportance:::generateMarginalContribution(x, perm))
        mc = unlist(mc, recursive = FALSE)
        
        # get all unique sets
        values = unique(unname(unlist(mc, recursive = FALSE)))
        values_feat = c(list(group_df$feature[group_df$group=="G1"]), values[[2]], list(group_df$feature[group_df$group%in% c("G1","G2")]),
                        list(group_df$feature[group_df$group=="G2"]), list(group_df$feature[group_df$group%in% c("G1","G2","G3")]),
                        list(group_df$feature[group_df$group%in% c("G2","G3")]), list(group_df$feature[group_df$group%in% c("G1","G3")]),
                        list(group_df$feature[group_df$group=="G3"]), list(group_df$feature[group_df$group%in% c("G1","G2","G3", "G4")]),
                        list(group_df$feature[group_df$group%in% c("G2","G3","G4")]),list(group_df$feature[group_df$group%in% c("G1","G3","G4")]),
                        list(group_df$feature[group_df$group%in% c("G3","G4")]),list(group_df$feature[group_df$group%in% c("G1","G2","G4")]),
                        list(group_df$feature[group_df$group%in% c("G2","G4")]),list(group_df$feature[group_df$group%in% c("G1","G4")]),
                        list(group_df$feature[group_df$group%in% c("G4")]))
        
        
        # compute value function based on importance
        data_group = data
        names(data_group)[1:40] = group_df$group 
        value.function = lapply(values_feat, function(f) {
          featureImportance:::calculateValueFunctionImportance(object = mod, data = data, measures = rmse,
                                                               n.feat.perm = 5, features = f)
        })
        
        vf = rbindlist(value.function)
        vf$features = stri_paste_list(values_feat, ",")
        vf$features = stri_paste_list(values, ",")
        
        
        for(group in groups){
          mc = featureImportance:::generateMarginalContribution(group, perm)
          mc.val = featureImportance:::getMarginalContributionValues(mc, vf)
          imp = featureImportance:::getShapleyImportance(mc.val)
          shap_imp[j,group] = imp
          #assign(paste0("imp",group,j), imp)
        }
        
      }
      
      shap_imp = colMeans(shap_imp)
      
      
      return(list(gpfi = gpfi, gopfi = gopfi, dgi = dgi, goi = goi, perf = res$aggr, shap_imp = shap_imp, res = res))
      
    })
    
    
    # reg$cluster.functions = makeClusterFunctionsSocket(4) #for windows
    # reg$cluster.functions = makeClusterFunctionsMulticore(4) #for Linux
    
    # add experiments
    addExperiments(repls = 1)
    
    # submit jobs
    start_time <- Sys.time()
    submitJobs()
    end_time <- Sys.time()
    
    
  }
}


