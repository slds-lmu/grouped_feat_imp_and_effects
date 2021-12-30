# install.packages("grpreg")
set.seed(123351)
library(grpreg)
data("Birthwt")
library("lars")
data(diabetes)
library(gglasso) #https://cran.r-project.org/web/packages/gglasso/gglasso.pdf
source("code/functions/gimp_classif.R")
library(mlr)

#bw task
bw.df = data.frame(Birthwt$X, bwt = Birthwt$bwt)
task_bw = makeRegrTask(id = "bw", data = bw.df, target = "bwt")

#colon task
data(colon)
group = rep(1:20,each=5)
x = data.frame(colon$x)
group_df_colon = data.frame(feature = as.character(names(x)), group = as.character(group), stringsAsFactors = FALSE)
task_colon = makeClassifTask(id = "colon", data = data.frame(x, y = as.factor(colon$y)), target = "y")

learners_regr = makeLearners(c("lm", "ranger", "h2o.deeplearning", "ksvm", "cvglmnet"), type = "regr")
learners_classif = makeLearners(c("logreg", "ranger", "h2o.deeplearning", "ksvm", "cvglmnet"), type = "classif")

bmr_bw = benchmark(learners_regr, task_bw, resamplings = makeResampleDesc("Subsample", iters = 100L), measures = rsq)
bmr_colon = benchmark(learners_classif, task_colon, resamplings = makeResampleDesc("Subsample", iters = 100L), measures = acc)


#bwt use case
learner_bwt = makeLearner("regr.ranger")
mod_bwt = train(learner_bwt, task1)
res_bwt = resample(learner_bwt, task1, cv10, measures = mse, models = TRUE)

learner_diabetes = makeLearner("regr.lm")
mod_diabetes = train(learner_diabetes, task2)
res_diabetes = resample(learner_diabetes, task2, cv10, measures = mse, models = TRUE)

#group lasso
group_df_bwt = data.frame(feature = as.character(getTaskFeatureNames(task1)), group = Birthwt$group)
group_df_bwt$feature = as.character(group_df_bwt$feature)
group_df_bwt$group = as.character(group_df_bwt$group)
group_df_diabetes = data.frame(feature = getTaskFeatureNames(task2), group = c(rep("body", 4), rep("lab", 6)))
group_df_diabetes$feature = as.character(group_df_diabetes$feature)
group_df_diabetes$group = as.character(group_df_diabetes$group)

x = as.data.frame(getTaskData(task1, target.extra = TRUE)$data)
x = as.matrix(x)
y = getTaskData(task1, target.extra = TRUE)$target
fit_bwt = grpreg(x,y,group_df_bwt$group, penalty="grLasso", nlambda = 400)
cv_bwt = cv.grpreg(x,y,group_df_bwt$group, penalty="grLasso", nlabmda = 400)
plot(fit_bwt, main = "Group Lasso")

#gimp analysis
library(future.apply)
source("code/functions/gimp.R")
gimp_bwt = Gimp$new(task = task1, res = res_bwt, mod = mod_bwt, lrn = learner_bwt)
gpfi_bwt = gimp_bwt$group_permutation_feat_imp(group_df = group_df_bwt, PIMP = FALSE)
gopfi_bwt = gimp_bwt$group_only_permutation_feat_imp(group_df = group_df_bwt, PIMP = FALSE)
LOGO_bwt = gimp_bwt$drop_group_importance(group_df = group_df_bwt, resampling = cv10, measures = mse)
  dgi = LOGO_bwt
  dgi_df = data.frame(Group = names(LOGO_bwt)[-9], DGI = NA)
  dgi_df$DGI[1] = dgi$all$aggr - dgi$age$aggr
  dgi_df$DGI[2] = dgi$all$aggr - dgi$lwt$aggr
  dgi_df$DGI[3] = dgi$all$aggr - dgi$race$aggr
  dgi_df$DGI[4] = dgi$all$aggr - dgi$smoke$aggr
  dgi_df$DGI[5] = dgi$all$aggr - dgi$ptl$aggr
  dgi_df$DGI[6] = dgi$all$aggr - dgi$ht$aggr
  dgi_df$DGI[7] = dgi$all$aggr - dgi$ui$aggr
  dgi_df$DGI[8] = dgi$all$aggr - dgi$ftv$aggr
  dgi_df$DGI = dgi_df$DGI * (-1)
LOGO_bwt = dgi_df

LOGI_bwt = gimp_bwt$group_only_importance(group_df = group_df_bwt, resampling = cv10, measures = mse)
  goi = LOGI_bwt
  goi_df = data.frame(Group = names(LOGI_bwt)[-9], goi = NA)
  goi_df$goi[1] = goi$age$aggr - goi$featureless$aggr
  goi_df$goi[2] = goi$lwt$aggr - goi$featureless$aggr
  goi_df$goi[3] = goi$race$aggr - goi$featureless$aggr
  goi_df$goi[4] = goi$smoke$aggr - goi$featureless$aggr
  goi_df$goi[5] = goi$ptl$aggr - goi$featureless$aggr
  goi_df$goi[6] = goi$ht$aggr - goi$featureless$aggr
  goi_df$goi[7] = goi$ui$aggr - goi$featureless$aggr
  goi_df$goi[8] = goi$ftv$aggr - goi$featureless$aggr
  goi_df$goi = goi_df$goi * (-1)
LOGI_bwt = goi_df

shap_bwt = gimp_bwt$shapley(group_df = group_df_bwt, res = res_bwt, n.shapley.perm = 120)

rank_df = data.frame(group = gpfi_bwt$features, stringsAsFactors = FALSE)
gpfi_bwt$rank_gpfi = rank(-gpfi_bwt$mse)
rank_df = rank_df %>% left_join(gpfi_bwt %>% dplyr::select(features, rank_gpfi), by = c("group" = "features"))
gopfi_bwt$rank_gopfi = rank(-gopfi_bwt$GOPFI)
rank_df = rank_df %>% left_join(gopfi_bwt %>% dplyr::select(features, rank_gopfi), by = c("group" = "features"))
LOGO_bwt$rank_LOGO = rank(-LOGO_bwt$DGI)
rank_df = rank_df %>% left_join(LOGO_bwt %>% dplyr::select(Group, rank_LOGO), by = c("group" = "Group"))
LOGI_bwt$rank_LOGI = rank(-LOGI_bwt$goi)
rank_df = rank_df %>% left_join(LOGI_bwt %>% dplyr::select(Group, rank_LOGI), by = c("group" = "Group"))
shap_bwt$rank_shap = rank(-shap_bwt$meanShapImp)
rank_df = rank_df %>% left_join(shap_bwt %>% dplyr::select(group, rank_shap), by = c("group" = "group"))

df_lasso = data.frame(group = c("ui", "smoke", "race", "ht", "ptl", "lwt", "age", "ftv"), rank_LASSO = 1:8)
rank_df = rank_df %>% left_join(df_lasso %>% dplyr::select(group, rank_LASSO), by = c("group" = "group"))
rank_df




#colon usecase
learner_colon = makeLearner("classif.ranger")
mod_colon = train(learner_colon, task_colon)
res_colon = resample(learner_colon, task_colon, cv10, measures = acc, models = TRUE)
gimp_colon = Gimp_classif$new(task = task_colon, res = res_colon, mod = mod_colon, lrn = learner_colon)
gpfi_colon = gimp_colon$group_permutation_feat_imp(group_df = group_df_colon, PIMP = FALSE)

LOGO_colon = gimp_colon$drop_group_importance(group_df = group_df_colon, resampling = cv10, measures = acc)
dgi = LOGO_colon
dgi_df = data.frame(Group = names(LOGO_colon)[-21], DGI = NA)
dgi_df$DGI[1] = dgi$all$aggr - dgi$`1`$aggr
dgi_df$DGI[2] = dgi$all$aggr - dgi$`2`$aggr
dgi_df$DGI[3] = dgi$all$aggr - dgi$`3`$aggr
dgi_df$DGI[4] = dgi$all$aggr - dgi$`4`$aggr
dgi_df$DGI[5] = dgi$all$aggr - dgi$`5`$aggr
dgi_df$DGI[6] = dgi$all$aggr - dgi$`6`$aggr
dgi_df$DGI[7] = dgi$all$aggr - dgi$`7`$aggr
dgi_df$DGI[8] = dgi$all$aggr - dgi$`8`$aggr
dgi_df$DGI[9] = dgi$all$aggr - dgi$`9`$aggr
dgi_df$DGI[10] = dgi$all$aggr - dgi$`10`$aggr
dgi_df$DGI[11] = dgi$all$aggr - dgi$`11`$aggr
dgi_df$DGI[12] = dgi$all$aggr - dgi$`12`$aggr
dgi_df$DGI[13] = dgi$all$aggr - dgi$`13`$aggr
dgi_df$DGI[14] = dgi$all$aggr - dgi$`14`$aggr
dgi_df$DGI[15] = dgi$all$aggr - dgi$`15`$aggr
dgi_df$DGI[16] = dgi$all$aggr - dgi$`16`$aggr
dgi_df$DGI[17] = dgi$all$aggr - dgi$`17`$aggr
dgi_df$DGI[18] = dgi$all$aggr - dgi$`18`$aggr
dgi_df$DGI[19] = dgi$all$aggr - dgi$`19`$aggr
dgi_df$DGI[20] = dgi$all$aggr - dgi$`20`$aggr
LOGO_colon = dgi_df
LOGO_colon %>% arrange(desc(DGI))



LOGI_colon = gimp_colon$group_only_importance(group_df = group_df_colon, resampling = cv10, measures = acc)
goi = LOGI_colon
goi_df = data.frame(Group = as.character(1:20), goi = NA)
goi_df$goi[1] = goi$featureless$aggr - goi[[1]]$aggr
goi_df$goi[2] = goi$featureless$aggr - goi[[2]]$aggr
goi_df$goi[3] = goi$featureless$aggr - goi[[3]]$aggr
goi_df$goi[4] = goi$featureless$aggr - goi[[4]]$aggr
goi_df$goi[5] = goi$featureless$aggr - goi[[5]]$aggr
goi_df$goi[6] = goi$featureless$aggr - goi[[6]]$aggr
goi_df$goi[7] = goi$featureless$aggr - goi[[7]]$aggr
goi_df$goi[8] = goi$featureless$aggr - goi[[8]]$aggr
goi_df$goi[9] = goi$featureless$aggr - goi[[9]]$aggr
goi_df$goi[10] = goi$featureless$aggr - goi[[10]]$aggr
goi_df$goi[11] = goi$featureless$aggr - goi[[11]]$aggr
goi_df$goi[12] = goi$featureless$aggr - goi[[12]]$aggr
goi_df$goi[13] = goi$featureless$aggr - goi[[13]]$aggr
goi_df$goi[14] = goi$featureless$aggr - goi[[14]]$aggr
goi_df$goi[15] = goi$featureless$aggr - goi[[15]]$aggr
goi_df$goi[16] = goi$featureless$aggr - goi[[16]]$aggr
goi_df$goi[17] = goi$featureless$aggr - goi[[17]]$aggr
goi_df$goi[18] = goi$featureless$aggr - goi[[18]]$aggr
goi_df$goi[19] = goi$featureless$aggr - goi[[19]]$aggr
goi_df$goi[20] = goi$featureless$aggr - goi[[20]]$aggr

goi_df$goi = goi_df$goi * (-1)
LOGI_colon = goi_df
LOGI_colon %>% arrange(desc(goi))

gopfi_colon = gimp_colon$group_only_permutation_feat_imp(group_df = group_df_colon, PIMP = FALSE)

shap_colon = gimp_colon$shapley(group_df = group_df_colon, res = res_colon, n.shapley.perm = 120)
