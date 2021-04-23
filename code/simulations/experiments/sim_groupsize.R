
library(BBmisc)
library(mlr)
library(data.table)
library(stringi)
library(batchtools)
library(checkmate)
library(featureImportance)
# unlink("results/simulation_results/sim_groupsize", recursive = TRUE)
reg = makeExperimentRegistry(file.dir = "results/simulation_results/sim_groupsize", seed = 123)

# add problem
for (i in 1:20) addProblem(as.character(i), i) 

# define algorithm
addAlgorithm("sim", function(data, job, instance) {
  # generate data
  p1 = 6
  p2 = 2
  n = 2000
  
  X_G1 = cbind(runif(n), runif(n), runif(n), runif(n), runif(n), runif(n))
  X_G2 = cbind(runif(n), runif(n))
  epsilon = rnorm(n, sd = 1)
  
  Xsim = as.data.frame(cbind(X_G1, X_G2))
  colnames(Xsim) = paste0("V",1:8)
  
  Z1 = 2*(X_G1[,1]) + 2*(X_G1[,3]) 
  Z2 = 2*(X_G2[,1])
   
  y = Z1 + Z2 + epsilon
  
  data = data.frame(cbind(Xsim,y))
  X = as.matrix(Xsim)
  colnames(X) = colnames(Xsim)
  
  
  # Model
  learnerRF = makeLearner("regr.ranger", par.vals = list(num.trees = 2000))
  task = makeRegrTask(data = data, target = "y")
  train.ind = sample(1:getTaskSize(task), floor(getTaskSize(task)/2))
  test.ind = setdiff(1:getTaskSize(task), train.ind)
  mod = mlr::train(learnerRF, task, subset = train.ind)
  data = getTaskData(task, subset = test.ind)
  features = getTaskFeatureNames(task)
  target = getTaskTargetNames(task)
  #mod = train(learner = learnerRF, task = task)
  #res = resample(learnerRF, task, cv5, measures = rsq, models = TRUE)
  
  #imp = shapleyImportance(mod, data, features = c("V1", "V2","V3","V4","V5"), target = "y", n.feat.perm = 5, measures = mse)
  
  # generate all permutations of groups
  group_df = data.frame(feature = colnames(Xsim), group = c(rep("G1", p1), rep("G2", p2)), stringsAsFactors = FALSE)
  groups = unique(group_df$group)
  perm = featureImportance:::generatePermutations(groups, n.shapley.perm = 120)
  
  # generate all marginal contribution sets for each feature
  mc = lapply(groups, function(x) featureImportance:::generateMarginalContribution(x, perm))
  mc = unlist(mc, recursive = FALSE)
  
  # get all unique sets
  values = unique(unname(unlist(mc, recursive = FALSE)))
  values_feat = c(list(group_df$feature[group_df$group=="G1"]), values[[2]], list(group_df$feature[group_df$group%in% c("G1","G2")]),
                  list(group_df$feature[group_df$group=="G2"]))
  
  
  # compute value function based on importance
  data_group = data
  names(data_group)[1:8] = group_df$group 
  value.function = lapply(values_feat, function(f) {
    featureImportance:::calculateValueFunctionImportance(object = mod, data = data, measures = getDefaultMeasure(task),
                                     n.feat.perm = 20, features = f)
  })
  
  vf = rbindlist(value.function)
  vf$features = stri_paste_list(values_feat, ",")
  vf$features = stri_paste_list(values, ",")
  
  for(group in groups){
    mc = featureImportance:::generateMarginalContribution(group, perm)
    mc.val = featureImportance:::getMarginalContributionValues(mc, vf)
    imp = featureImportance:::getShapleyImportance(mc.val)
    assign(paste0("imp",group), imp)
  }
  
  
  # individual features
  imp = featureImportance::shapleyImportance(mod, data, features = features, n.feat.perm = 5,n.shapley.perm = 1000, measures = getDefaultMeasure(task))
  
  

  
  return(list(shapleyValue = imp$shapley.value, shapleyUncertainty = imp$shapley.uncertainty, marginalContributions = imp$marginal.contributions,
              impG1 = impG1, impG2 = impG2))
  
})


# reg$cluster.functions = makeClusterFunctionsSocket(4) #for windows
# reg$cluster.functions = makeClusterFunctionsMulticore(4) #for Linux

# add experiments
addExperiments(repls = 1)

# submit jobs
start_time <- Sys.time()
submitJobs()
end_time <- Sys.time()

