set.seed(1312314)
library(mlrCPO)
source("code/functions/gimp.R")
load("usecase_psychology.RData")

task = personalityTasks$C #choose Task (E, C, or O are ok)
target = getTaskTargetNames(task)
data = getTaskData(task)
data = mlr::impute(data, classes = list(integer = imputeMedian(), factor = imputeMode(), numeric = imputeMedian()))$data
data[, which(colnames(data) != target)] = scale(data[, which(colnames(data) != target)])
learnerRF = makeLearner("regr.ranger")
task = makeRegrTask(data = data, target = target)
mod = train(learner = learnerRF, task = task)
res = resample(learnerRF, task, cv10, measures = rsq, models = TRUE)
gimp = Gimp$new(task = task, res = res, mod = mod, lrn = learnerRF)

library(dplyr)
cat_table = feat_imp_total %>% dplyr::select(feature, category2) %>% distinct(feature, .keep_all = TRUE)
features = getTaskFeatureNames(task) %>% data.frame()
colnames(features) = "feature"
cat_table = features %>% left_join(cat_table)
colnames(cat_table)[2] = "group"
group_df = cat_table

# greedy
source("code/functions/greedy_forward_search.R")

library(batchtools) 
unlink("results/usecase_results/usecase_gfs_refit", recursive = TRUE)
reg = makeExperimentRegistry(file.dir = "results/usecase_results/usecase_gfs_refit", source = c("code/functions/greedy_forward_search.R", "code/functions/gimp.R"), seed = 123)

for (i in 1:100) addProblem(as.character(i), i)
addAlgorithm("gfs", function(data, job, instance) {
  greedy_fw_search_refit(gimp, cat_table, delta = 0.001, outer_resampling = makeResampleDesc("Subsample", iters = 1), inner_resampling = cv10, measure = mse)
})
addExperiments(repls = 1)
batchtools::batchExport(export = list(
  gimp = gimp, cat_table = cat_table
 )
)

reg$cluster.functions = makeClusterFunctionsMulticore(6)
submitJobs()
