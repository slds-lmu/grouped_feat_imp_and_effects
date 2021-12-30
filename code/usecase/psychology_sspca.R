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

wid = 6.5
hei = 4

# ggsave("LPDP_usecase_app_1.pdf", gg1$p1 + theme_bw(), width = wid, height = hei)
# ggsave("LPDP_usecase_app_2.pdf", gg2$p1 + theme_bw(), width = wid, height = hei)
# ggsave("LPDP_usecase_general_1.pdf", gg3$p1 + theme_bw(), width = wid, height = hei)
# ggsave("LPDP_usecase_music_1.pdf", gg4$p1 + theme_bw(), width = wid, height = hei)

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
gg1 = gimp$group_pdp2(features = a, x = 0.9)
df_features = data.frame(feature = names(a), variable = paste0("x", 1:length(a)))
tb = paste0(df_features$variable, " = ",df_features$feature, collapse = "\n")
gg1$p1 + annotate(geom = "label", x = 0.7, y = -0.08, label = tb, hjust = 0, fill = "white")

library(ggplot2)
gg1$p1 + theme_bw()
ggsave("LPDP_usecase_app_1.pdf", gg1$p1 + theme_bw(), width = wid, height = hei)

#plot2
pmd = sspca(X, Y, sigma.y = sigma.y, sumabsv = 1.5, K = 4)
i = 1
pmd$v[, i][pmd$v[, i] != 0] -> a; names(a) <- colnames(X)[pmd$v[, i] != 0]; a
gg2 = gimp$group_pdp2(features = a, x = 1)
gg2$p1 + theme_bw()
ggsave("LPDP_usecase_app_2.pdf", gg2$p1 + theme_bw(), width = wid, height = hei)

#plot3
feat = cat_table %>% dplyr::filter(group == "general") %>% pull(feature) %>% as.character()
X = data %>% dplyr::select(all_of(feat)) %>% as.matrix()
Y = data %>% dplyr::select(getTaskTargetNames(gimp$task)) %>% as.matrix()

library(kernlab)
kernel.y = "gaussian"
sigma.y = sqrt(sigest(Y~X)[2]/2)
pmd = sspca(X, Y, sigma.y = sigma.y, sumabsv = 1.5, K = 4)
i = 1
pmd$v[, i][pmd$v[, i] != 0] -> a; names(a) <- colnames(X)[pmd$v[, i] != 0]; a
gg3 = gimp$group_pdp2(features = a, x = 1.5)
gg3$p1 + theme_bw()
ggsave("LPDP_usecase_general_1.pdf", gg3$p1 + theme_bw(), width = wid, height = hei)

#plot 4
# i = 2
# pmd$v[, i][pmd$v[, i] != 0] -> a; names(a) <- colnames(X)[pmd$v[, i] != 0]; a
# gg = gimp$group_pdp2(features = a)
# gg$p1 + theme_bw()
# ggsave("LPDP_usecase_general_2.pdf", gg$p1 + theme_bw(), width = 10, height = 7)

# plot 5
feat = cat_table %>% dplyr::filter(group == "music") %>% pull(feature) %>% as.character()
X = data %>% dplyr::select(all_of(feat)) %>% as.matrix()
Y = data %>% dplyr::select(getTaskTargetNames(gimp$task)) %>% as.matrix()

sigma.y = sqrt(sigest(Y~X)[2]/2)
pmd = sspca(X, Y, sigma.y = sigma.y, sumabsv = 1.5, K = 4)
i = 1
pmd$v[, i][pmd$v[, i] != 0] -> a; names(a) <- colnames(X)[pmd$v[, i] != 0]; a
gg4 = gimp$group_pdp2(features = a, x = 4.5, y = -0.082)
gg4$p1 + theme_bw()
ggsave("LPDP_usecase_music_1.pdf", gg4$p1 + theme_bw(), width = wid, height = hei)
