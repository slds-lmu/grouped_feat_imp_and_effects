set.seed(1312314)
library(mlrCPO)
library(future.apply)
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
res = resample(learnerRF, task, cv10, measures = mse, models = TRUE)
gimp = Gimp$new(task = task, res = res, mod = mod, lrn = learnerRF)

library(dplyr)
cat_table = feat_imp_total %>% dplyr::select(feature, category2) %>% distinct(feature, .keep_all = TRUE)
features = getTaskFeatureNames(task) %>% data.frame()
colnames(features) = "feature"
cat_table = features %>% left_join(cat_table)
colnames(cat_table)[2] = "group"
group_df = cat_table

#gpfi
gpfi = gimp$group_permutation_feat_imp(cat_table, s = 50, n.feat.perm = 10, regr.measure = mse)

#gopfi
gopfi = gimp$group_only_permutation_feat_imp(cat_table, s = 50, n.feat.perm = 10, regr.measure = mse)

# dgi
dgi = gimp$drop_group_importance(cat_table, resampling = cv10, measures = mse)

# goi
goi = gimp$group_only_importance(cat_table, resampling = cv10, measures = mse)

# shapley
shap = gimp$shapley(group_df = cat_table, res = res, n.shapley.perm = 120)

saveRDS(list(gpfi = gpfi, gopfi = gopfi, dgi = dgi, goi = goi, shap = shap), file = "results/usecase_results/psychology_feat_imp.RDS")
