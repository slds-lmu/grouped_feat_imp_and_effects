#three groups G1, G2, G3
#G1 and G2 very similar, each have high importance
#G3 independent, lower importance
#question: If only 2 groups of features should be used, which ones?
source("code/functions/helper.R")
source("code/functions/gimp.R")
source("code/functions/greedy_forward_search.R")
library(future.apply)
set.seed(213456)
p = 20
n = 1000
Xsim = sim_data(N=n, n_groups= 2, y = c(2,1), size_groups = rep(p/2,2))
epsilon = rnorm(n)

df = data.frame(cbind(Xsim, Xsim[, 1:10] + rnorm(n, sd = 0.01)))

# Model
learner = makeLearner("regr.ksvm")
task = makeRegrTask(data = df, target = "target")

mod = train(learner = learner, task = task)
res = resample(learner, task, cv10, measures = mse, models = TRUE)
gimp = Gimp$new(task = task, res = res, mod = mod, lrn = learner)

group_df = data.frame(feature = colnames(df %>% dplyr::select(-target)), group = c(rep("G1", 10), rep("G3", 10), rep("G2", 10)), stringsAsFactors = FALSE)

gopfi = gimp$group_only_permutation_feat_imp(group_df, PIMP = FALSE, regr.measure = mse)
gpfi = gimp$group_permutation_feat_imp(group_df, PIMP = FALSE, regr.measure = mse)
goi = gimp$group_only_importance(group_df, measures = mse)
dgi = gimp$drop_group_importance(group_df, measures = mse)
shap = gimp$shapley(group_df = group_df, n.shapley.perm = 120)

gfs_refit = greedy_fw_search_refit(gimp, group_df, delta = 0.01, outer_resampling = makeResampleDesc("Subsample", iters = 100), measure = mse)

saveRDS(list(gopfi = gopfi, gpfi = gpfi, goi = goi, dgi = dgi, gfs_refit = gfs_refit, shap = shap), file = "results/simulation_results/gfs_sim.RDS")
