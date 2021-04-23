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


#sspca in group "app"
source("code/functions/sspca.R")
data = getTaskData(gimp$task)
feat = cat_table %>% dplyr::filter(group == "app") %>% pull(feature) %>% as.character()
X = data %>% dplyr::select(all_of(feat)) %>% as.matrix()
Y = data %>% dplyr::select(getTaskTargetNames(gimp$task)) %>% as.matrix()

library(kernlab)
kernel.y = "gaussian"
sigma.y = sqrt(sigest(Y~X)[2]/2)
pmd = sspca(X, Y, sigma.y = sigma.y, sumabsv = 1.8, K = 4)
i = 2
pmd$v[, i][pmd$v[, i] != 0] -> a; names(a) <- colnames(X)[pmd$v[, i] != 0]; a
gg = gimp$group_pdp2(features = a)

gg$p1

ggsave("results/figures/LPDP_usecase1.pdf", gg$p1, width = 10, height = 7)
a
