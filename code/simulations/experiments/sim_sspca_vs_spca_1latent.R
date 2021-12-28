library(kernlab)
library(batchtools)
library(totalvis)
# unlink("results/simulation_results/sim1latent", recursive = TRUE)
reg = makeExperimentRegistry(file.dir = "results/simulation_results/sim1latent", source = c("code/functions/gimp_sim.R", "code/functions/sspca.R", "code/functions/helper.R"), seed = 123, packages = "kernlab")

# add problem
for (i in 1:20) addProblem(as.character(i), i) 

# define algorithm
addAlgorithm("sim", function(data, job, instance) {
  # generate data
  p = 50
  n = 500
  
  Xsim = sim_data_sim(N = n, n_groups = 5, y = rep(0,5), size_groups = rep(10,5), prop = c(0.2,0.4,0.6,0.8,1))
  epsilon = rnorm(n)
  
  Xsim = as.data.frame(Xsim)
  colnames(Xsim) = paste0("V",1:p)
  
  Z1 = Xsim[,5] - 2*Xsim[,8] - 4*Xsim[,25] + 8*Xsim[,47] + 4*Xsim[,49]
  y = (Z1) + epsilon
 
  
  data = data.frame(cbind(Xsim,y))
  X = as.matrix(Xsim)
  colnames(X) = colnames(Xsim)
  
  
  # Model
  learnerRF = makeLearner("regr.ranger", par.vals = list(num.trees = 2000))
  task = makeRegrTask(data = data, target = "y")
  
  mod = train(learner = learnerRF, task = task)
  res = resample(learnerRF, task, cv10, measures = rsq, models = TRUE)
  gimp = Gimpsim$new(task = task, res = res, mod = mod, lrn = learnerRF)
  #Y = as.matrix(y)
  #Y = as.matrix(res$pred$data$response)[res$pred$data$id]
  data = getTaskData(task)
  data2 = data
  data2$id = as.numeric(rownames(data2))
  data2 = data2 %>% left_join(res$pred$data[which(colnames(res$pred$data) %in% c("id", "truth", "response"))])
  Y = as.matrix(data.frame(data2$response))
  
  
  #-----------------------------------------------------------------------------
  # SSPCA
  
  kernel.y = "gaussian"
  sigma.y = sigest(Y~X)[2]
  L <- kernelMatrix(kernel = rbfdot(sigma = sigma.y), Y)
  #SSPCA 1). Decompose L such that L = Δ^TΔ
  Delta = my_chol_psd(L)
  n.row <- nrow(X)
  # Delta = chol(L, pivot = TRUE)
  #SSPCA 2). H ← I−n^{−1}ee^T
  H = diag(1, n.row, n.row) - matrix(1, nrow = n.row, ncol = n.row) / n.row
  #SSPCA 3). Ψ ← ΔTHX
  Psi = t(Delta) %*% H %*% X
  # SSPCA 4). Compute the sparse basis based on the PMD method
  # Psi_1 = Psi
  
  library(PMA)
  Psi[1, 1] = NA
  ncomp = 1
  pmd = SPC(Psi, K = ncomp, sumabsv = 2, vpos = TRUE)
  #(pmd$u %*% diag(pmd$d) %*% t(pmd$v))[1, 4]
  pmd
  # PMD.cv(Psi)
  
  pmd$v[, 1][pmd$v[, 1] != 0] -> a; names(a) <- colnames(X)[pmd$v[, 1] != 0]
  a = list(a)
  b = as.numeric(c(1,-2,-4,8,4)); names(b) <- c("V5","V8","V25","V47","V49")
  #b = as.numeric(c(1,1)); names(b) <- c("V2","V4")
  gg = gimp$group_pdp(features = a, featuresTrue = b, parts = 4, x = 2, y = 0)
  
  df_sspca = data.frame("R2" = numeric(), "mse_z1" = numeric(), "mse_resp" = numeric(), "mse_resp_true" = numeric(),
    "mse_resp_true_true" = numeric(), "mse_resp_yhat" = numeric(),  "mse_resp_y" = numeric(),
    "mse_zhat_true" = numeric(), "feat_corr" = numeric(),
    "feat_mod" = numeric(), "feat_true" = numeric())
  
  df_sspca[1, "R2"] = res$aggr
  df_sspca[1, "mse_z1"] = measureMSE(scale(gg$pdp$dim_redTrue), scale(gg$pdp$dim_red))
  df_sspca[1, "mse_resp"] = measureMSE(gg$pdp$pdTrue, gg$pdp$pd)
  df_sspca[1, "mse_resp_true"] = measureMSE(gg$pdp$res_true, gg$pdp$pd)
  df_sspca[1, "mse_resp_true_true"] = measureMSE(gg$pdp$res_truetrue, gg$pdp$pd)
  df_sspca[1, "mse_resp_yhat"] = measureMSE(Y, gg$pdp$pd)
  df_sspca[1, "mse_resp_y"] = measureMSE(gg$pdp$res_truetrue, Y)
  df_sspca[1, "mse_zhat_true"] = measureMSE(gg$pdp$res_truetrue, gg$pdp$res_true)
  df_sspca[1, "feat_corr"] = gg$pdp$feat_numbBoth[1]
  df_sspca[1, "feat_mod"] = gg$pdp$feat_numb[1]
  df_sspca[1, "feat_true"] = gg$pdp$feat_numbTrue[1]
  
  df_detailed_sspca = data.frame("pd" = numeric(), "pd_true" = numeric(), "res_mod" = numeric(),
    "zhat_true" = numeric(), "z_true" = numeric(),"zhat" = numeric(), "z" = numeric())
  
  df_detailed_sspca[1:n, "pd"] = gg$pdp$pd
  df_detailed_sspca[1:n, "pd_true"] = gg$pdp$pdTrue
  df_detailed_sspca[1:n, "res_mod"] = Y
  df_detailed_sspca[1:n, "zhat_true"] = gg$pdp$res_true
  df_detailed_sspca[1:n, "z_true"] = gg$pdp$res_truetrue
  df_detailed_sspca[1:n, "z_true"] = gg$pdp$res_truetrue
  df_detailed_sspca[1:n, "zhat"] = gg$pdp$dim_red
  df_detailed_sspca[1:n, "z"] = gg$pdp$dim_redTrue
  
  
  # selected features with loadings
  df_features_sspca = data.frame(matrix(ncol = p, nrow = 0))
  colnames(df_features_sspca) = colnames(X)
  pmd$v[, 1][pmd$v[, 1] > 0.05] -> a; names(a) <- colnames(X)[pmd$v[, 1] > 0.05]
  df_features_sspca[1,] = as.numeric(colnames(X) %in% names(a))
  df_features_sspca[1,which(colnames(X) %in% names(a))] = a
  
  
  #-----------------------------------------------------------------------------
  # pca
  
  spca = elasticnet::spca(X, ncomp, sparse = "varnum", para = 5)
  
  spca$loadings[,1][spca$loadings[,1] != 0] -> a; names(a) <- colnames(X)[spca$loadings[,1] != 0]
  a = list(a)
  gg = gimp$group_pdp(features = a, featuresTrue = b, parts = 4, x = 2, y = 0)
  
  df_pca = data.frame("R2" = numeric(), "mse_z1" = numeric(), "mse_resp" = numeric(), "mse_resp_true" = numeric(),
    "mse_resp_true_true" = numeric(), "mse_resp_yhat" = numeric(),  "mse_resp_y" = numeric(),
    "mse_zhat_true" = numeric(), "feat_corr" = numeric(),
    "feat_mod" = numeric(), "feat_true" = numeric())
  df_pca[1,"R2"] = res$aggr
  df_pca[1,"mse_z1"] = measureMSE(scale(gg$pdp$dim_redTrue), scale(gg$pdp$dim_red))
  df_pca[1, "mse_resp"] = measureMSE(gg$pdp$pdTrue, gg$pdp$pd)
  df_pca[1, "mse_resp_true"] = measureMSE(gg$pdp$res_true, gg$pdp$pd)
  df_pca[1, "mse_resp_true_true"] = measureMSE(gg$pdp$res_truetrue, gg$pdp$pd)
  df_pca[1, "mse_resp_yhat"] = measureMSE(Y, gg$pdp$pd)
  df_pca[1, "mse_resp_y"] = measureMSE(gg$pdp$res_truetrue, Y)
  df_pca[1, "mse_zhat_true"] = measureMSE(gg$pdp$res_truetrue, gg$pdp$res_true)
  df_pca[1, "feat_corr"] = gg$pdp$feat_numbBoth[1]
  df_pca[1, "feat_mod"] = gg$pdp$feat_numb[1]
  df_pca[1, "feat_true"] = gg$pdp$feat_numbTrue[1]
  
  df_detailed_pca = data.frame("pd" = numeric(), "pd_true" = numeric(), "res_mod" = numeric(), 
    "zhat_true" = numeric(), "z_true" = numeric(),"zhat" = numeric(), "z" = numeric())
  
  df_detailed_pca[1:n, "pd"] = gg$pdp$pd
  df_detailed_pca[1:n, "pd_true"] = gg$pdp$pdTrue
  df_detailed_pca[1:n, "res_mod"] = Y
  df_detailed_pca[1:n, "zhat_true"] = gg$pdp$res_true
  df_detailed_pca[1:n, "z_true"] = gg$pdp$res_truetrue
  df_detailed_pca[1:n, "zhat"] = gg$pdp$dim_red
  df_detailed_pca[1:n, "z"] = gg$pdp$dim_redTrue 
  
  # selected features with loadings
  df_features_pca = data.frame(matrix(ncol = p, nrow = 0))
  colnames(df_features_pca) = colnames(X)
  
  spca$loadings[,1][abs(spca$loadings[,1]) > 0.05] -> a; names(a) <- colnames(X)[abs(spca$loadings[,1]) > 0.05]
  df_features_pca[1,] = as.numeric(colnames(X) %in% names(a))
  df_features_pca[1,which(colnames(X) %in% names(a))] = a
  
  
  #-------------------------------------------------------------------------------------------------
  # Totalvis
  df_totalvis = totalvis:::totalvis(mod, X)
  
  return(list(df_detailed_pca = df_detailed_pca, df_pca = df_pca, df_detailed_sspca = df_detailed_sspca, 
              df_sspca = df_sspca, df_features_pca = df_features_pca, df_features_sspca = df_features_sspca,
              df_totalvis = df_totalvis))
})


# reg$cluster.functions = makeClusterFunctionsSocket(4) #for windows
# reg$cluster.functions = makeClusterFunctionsMulticore(4) #for Linux

# add experiments
addExperiments(repls = 1)

# submit jobs
start_time <- Sys.time()
submitJobs()
end_time <- Sys.time()



