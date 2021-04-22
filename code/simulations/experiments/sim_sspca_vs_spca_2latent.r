
library(kernlab)
library(batchtools)
load_all("packages/totalvis")
# unlink("results/simulation_results/sim2latent", recursive = TRUE)
reg = makeExperimentRegistry(file.dir = "results/simulation_results/sim2latent", source = c("code/functions/gimp_sim.r", "code/functions/sspca.r", "code/functions/helper.r"), seed = 123, packages = "kernlab")

# add problem
for (i in 1:100) addProblem(as.character(i), i) 

# add algorithm
addAlgorithm("sim", function(data, job, instance) {
    p = 20
    n = 500
  
    Xsim = data.frame(sim_data_sim(N = n, n_groups = 2, y = rep(0,2), size_groups = rep(p/2,2), prop = c(0.15,0.35)),
                      sim_data_sim(N = n, n_groups = 2, y = rep(0,2), size_groups = rep(p/2,2), prop = c(0.55,0.8)))
    epsilon = rnorm(n)
    
    Xsim = as.data.frame(Xsim)
    colnames(Xsim) = paste0("V",1:(2*p))
    
    Z1 = 3*Xsim[,3] - 2*Xsim[,8] - 4*Xsim[,13] + 8*Xsim[,18] 
    Z2 =  2*Xsim[,25] + 4*Xsim[,35]
    y = (Z1) + (Z2)^2 + epsilon
    
    data = data.frame(cbind(Xsim,y))
    X = as.matrix(Xsim)
    colnames(X) = colnames(Xsim)
    
    
    # Model
    learnerRF = makeLearner("regr.ranger", par.vals = list(num.trees = 2000))
    task = makeRegrTask(data = data, target = "y")
    
    mod = train(learner = learnerRF, task = task)
    res = resample(learnerRF, task, cv10, measures = rsq, models = TRUE)
    gimp = Gimpsim$new(task = task, res = res, mod = mod, lrn = learnerRF)

    data = getTaskData(task)
    data2 = data
    data2$id = as.numeric(rownames(data2))
    data2 = data2 %>% left_join(res$pred$data[which(colnames(res$pred$data) %in% c("id", "truth", "response"))])
    Y = as.matrix(data.frame(data2$response))
    
    
    #---------------------------------------------------------------------------
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
    ncomp = 2
    pmd = SPC(Psi, K = ncomp, sumabsv = 2, vpos = TRUE)
    #(pmd$u %*% diag(pmd$d) %*% t(pmd$v))[1, 4]
    # PMD.cv(Psi)
    
    # firts latent Variable
    a = list()
    pmd$v[, 1][pmd$v[, 1] != 0] -> a[[1]]; names(a[[1]]) <- colnames(X)[pmd$v[, 1] != 0]
    pmd$v[, 2][pmd$v[, 2] != 0] -> a[[2]]; names(a[[2]]) <- colnames(X)[pmd$v[, 2] != 0]
    b = as.numeric(c(3,-2,-4,8)); names(b) <- c("V3","V8","V13","V18")
    gg = gimp$group_pdp(features = a, featuresTrue = b, parts = 4, x = 2, y = 0)
    
    df_sspca = data.frame("R2" = numeric(), "mse_resp_y" = numeric(), "mse_z1" = numeric(), "mse_resp_yhat1" = numeric(), 
      "feat_corr1" = numeric(),"feat_mod1" = numeric(), "feat_true1" = numeric(),
      "mse_z2" = numeric(), "mse_resp_yhat2" = numeric(), 
      "feat_corr2" = numeric(),"feat_mod2" = numeric(), "feat_true2" = numeric())
    
    df_sspca[1,"R2"] = res$aggr
    df_sspca[1, "mse_resp_y"] = measureMSE(y, Y)
    
    df_sspca[1,"mse_z1"] = measureMSE(scale(gg$pdp$dim_redTrue), scale(gg$pdp$dim_red))
    df_sspca[1, "mse_resp_yhat1"] = measureMSE(Y, gg$pdp$pd)
    df_sspca[1, "feat_corr1"] = gg$pdp$feat_numbBoth[1]
    df_sspca[1, "feat_mod1"] = gg$pdp$feat_numb[1]
    df_sspca[1, "feat_true1"] = gg$pdp$feat_numbTrue[1]
    
    df_detailed_sspca = data.frame("pd1" = numeric(), "pd_true1" = numeric(), "res_mod1" = numeric(), 
      "zhat1" = numeric(), "z1" = numeric(), "pd2" = numeric(), "pd_true2" = numeric(), 
      "res_mod2" = numeric(),  "zhat2" = numeric(), "z2" = numeric())
    
    df_detailed_sspca[1:n, "pd1"] = gg$pdp$pd
    df_detailed_sspca[1:n, "pd_true1"] = gg$pdp$pdTrue
    df_detailed_sspca[1:n, "res_mod1"] = Y
    df_detailed_sspca[1:n, "zhat1"] = gg$pdp$dim_red
    df_detailed_sspca[1:n, "z1"] = gg$pdp$dim_redTrue
    
    # second latent variable
    b = as.numeric(c(2,4)); names(b) <- c("V25","V35")
    gg = gimp$group_pdp(features = a, featuresTrue = b, parts = 4, x = 2, y = 0)
    
    
    df_sspca[1,"mse_z2"] = measureMSE(scale(gg$pdp$dim_redTrue), scale(gg$pdp$dim_red))
    df_sspca[1, "mse_resp_yhat2"] = measureMSE(Y, gg$pdp$pd)
    df_sspca[1, "feat_corr2"] = gg$pdp$feat_numbBoth[1]
    df_sspca[1, "feat_mod2"] = gg$pdp$feat_numb[1]
    df_sspca[1, "feat_true2"] = gg$pdp$feat_numbTrue[1]
    
    df_detailed_sspca[1:n, "pd2"] = gg$pdp$pd
    df_detailed_sspca[1:n, "pd_true2"] = gg$pdp$pdTrue
    df_detailed_sspca[1:n, "res_mod2"] = Y
    df_detailed_sspca[1:n, "zhat2"] = gg$pdp$dim_red
    df_detailed_sspca[1:n, "z2"] = gg$pdp$dim_redTrue
    
    
    
    # feature selection
    df_features_sspca = df <- data.frame(matrix(ncol = (2*p + 1), nrow = 0))
    colnames(df_features_sspca) = c(colnames(X), "latent")
    
    # latent 1
    a = list()
    pmd$v[, 1][pmd$v[, 1] > 0.05] -> a[[1]]; names(a[[1]]) <- colnames(X)[pmd$v[, 1] > 0.05]
    pmd$v[, 2][pmd$v[, 2] > 0.05] -> a[[2]]; names(a[[2]]) <- colnames(X)[pmd$v[, 2] > 0.05]
    b = as.numeric(c(3,-2,-4,8)); names(b) <- c("V3","V8","V13","V18")
    
    
    dim_redTrue = (data %>% data.frame()  %>% dplyr::select(names(b)) %>% as.matrix()) %*% b
    mse = c()
    df_feat = list()
    for (k in 1:length(a)) {
      df_feat[[k]] = (data %>% scale() %>% data.frame() %>% dplyr::select(names(a[[k]])) %>% as.matrix()) %*% a[[k]]
      mse[k] = measureMSE(truth = scale(dim_redTrue), response = scale(df_feat[[k]]))
    }
    best = which.min(mse)
    dim_red = df_feat[[best]]
    features = a[[best]]
    
    
    df_features_sspca[1,1:(2*p)] = as.numeric(colnames(X) %in% names(features))
    df_features_sspca[1,which(colnames(X) %in% names(features))] = features
    df_features_sspca$latent[1] = 1
    
    # second latent variable
    b = as.numeric(c(2,4)); names(b) <- c("V25","V35")
    
    dim_redTrue = (data %>% data.frame()  %>% dplyr::select(names(b)) %>% as.matrix()) %*% b
    mse = c()
    df_feat = list()
    for (k in 1:length(a)) {
      df_feat[[k]] = (data %>% scale() %>% data.frame() %>% dplyr::select(names(a[[k]])) %>% as.matrix()) %*% a[[k]]
      mse[k] = measureMSE(truth = scale(dim_redTrue), response = scale(df_feat[[k]]))
    }
    best = which.min(mse)
    dim_red = df_feat[[best]]
    features = a[[best]]
    
    
    df_features_sspca[2,1:(2*p)] = as.numeric(colnames(X) %in% names(features))
    df_features_sspca[2,which(colnames(X) %in% names(features))] = features
    df_features_sspca$latent[2] = 2
    
    
    
    #---------------------------------------------------------------------------
    # PCA
    
    spca = elasticnet::spca(X, ncomp, sparse = "varnum", para = c(7,7))
    
    
    # firts latent Variable
    a = list()
    spca$loadings[,1][spca$loadings[,1] != 0] -> a[[1]]; names(a[[1]]) <- colnames(X)[spca$loadings[,1] != 0]
    spca$loadings[,2][spca$loadings[,2] != 0] -> a[[2]]; names(a[[2]]) <- colnames(X)[spca$loadings[,2] != 0]
    b = as.numeric(c(1,-2,-4,8)); names(b) <- c("V3","V8","V13","V18")
    
    gg = gimp$group_pdp(features = a, featuresTrue = b, parts = 4, x = 2, y = 0)
    
    df_pca = data.frame("R2" = numeric(), "mse_resp_y" = numeric(), "mse_z1" = numeric(), "mse_resp_yhat1" = numeric(), 
      "feat_corr1" = numeric(),"feat_mod1" = numeric(), "feat_true1" = numeric(),
      "mse_z2" = numeric(), "mse_resp_yhat2" = numeric(), 
      "feat_corr2" = numeric(),"feat_mod2" = numeric(), "feat_true2" = numeric())
    
    df_pca[1,"R2"] = res$aggr
    df_pca[1, "mse_resp_y"] = measureMSE(y, Y)
    
    df_pca[1,"mse_z1"] = measureMSE(scale(gg$pdp$dim_redTrue), scale(gg$pdp$dim_red))
    df_pca[1, "mse_resp_yhat1"] = measureMSE(Y, gg$pdp$pd)
    df_pca[1, "feat_corr1"] = gg$pdp$feat_numbBoth[1]
    df_pca[1, "feat_mod1"] = gg$pdp$feat_numb[1]
    df_pca[1, "feat_true1"] = gg$pdp$feat_numbTrue[1]
    
    df_detailed_pca = data.frame("pd1" = numeric(), "pd_true1" = numeric(), "res_mod1" = numeric(), 
      "zhat1" = numeric(), "z1" = numeric(), "pd2" = numeric(), "pd_true2" = numeric(), 
      "res_mod2" = numeric(),  "zhat2" = numeric(), "z2" = numeric())
    
    df_detailed_pca[1:n, "pd1"] = gg$pdp$pd
    df_detailed_pca[1:n, "pd_true1"] = gg$pdp$pdTrue
    df_detailed_pca[1:n, "res_mod1"] = Y
    df_detailed_pca[1:n, "zhat1"] = gg$pdp$dim_red
    df_detailed_pca[1:n, "z1"] = gg$pdp$dim_redTrue
    
    # second latent variable
    b = as.numeric(c(2,4)); names(b) <- c("V25","V35")
    gg = gimp$group_pdp(features = a, featuresTrue = b, parts = 4, x = 2, y = 0)
    
    df_pca[1,"mse_z2"] = measureMSE(scale(gg$pdp$dim_redTrue), scale(gg$pdp$dim_red))
    df_pca[1, "mse_resp_yhat2"] = measureMSE(Y, gg$pdp$pd)
    df_pca[1, "feat_corr2"] = gg$pdp$feat_numbBoth[1]
    df_pca[1, "feat_mod2"] = gg$pdp$feat_numb[1]
    df_pca[1, "feat_true2"] = gg$pdp$feat_numbTrue[1]
    
    df_detailed_pca[1:n, "pd2"] = gg$pdp$pd
    df_detailed_pca[1:n, "pd_true2"] = gg$pdp$pdTrue
    df_detailed_pca[1:n, "res_mod2"] = Y
    df_detailed_pca[1:n, "zhat2"] = gg$pdp$dim_red
    df_detailed_pca[1:n, "z2"] = gg$pdp$dim_redTrue
    
    # feature selection
    df_features_pca = df <- data.frame(matrix(ncol = (2*p + 1), nrow = 0))
    colnames(df_features_pca) = c(colnames(X), "latent")
    
    # latent 1
    a = list()
    spca$loadings[,1][abs(spca$loadings[,1]) > 0.05] -> a[[1]]; names(a[[1]]) <- colnames(X)[abs(spca$loadings)[,1] > 0.05]
    spca$loadings[,2][abs(spca$loadings[,2]) > 0.05] -> a[[2]]; names(a[[2]]) <- colnames(X)[abs(spca$loadings)[,2] > 0.05]
    b = as.numeric(c(1,-2,-4,8)); names(b) <- c("V3","V8","V13","V18")
    
    
    dim_redTrue = (data %>% data.frame()  %>% dplyr::select(names(b)) %>% as.matrix()) %*% b
    mse = c()
    df_feat = list()
    for (k in 1:length(a)) {
      df_feat[[k]] = (data %>% scale() %>% data.frame() %>% dplyr::select(names(a[[k]])) %>% as.matrix()) %*% a[[k]]
      mse[k] = measureMSE(truth = scale(dim_redTrue), response = scale(df_feat[[k]]))
    }
    best = which.min(mse)
    dim_red = df_feat[[best]]
    features = a[[best]]
    
    
    df_features_pca[1,1:(2*p)] = as.numeric(colnames(X) %in% names(features))
    df_features_pca[1,which(colnames(X) %in% names(features))] = features
    df_features_pca$latent[1] = 1
    
    # second latent variable
    b = as.numeric(c(2,4)); names(b) <- c("V25","V35")
    
    dim_redTrue = (data %>% data.frame()  %>% dplyr::select(names(b)) %>% as.matrix()) %*% b
    mse = c()
    df_feat = list()
    for (k in 1:length(a)) {
      df_feat[[k]] = (data %>% scale() %>% data.frame() %>% dplyr::select(names(a[[k]])) %>% as.matrix()) %*% a[[k]]
      mse[k] = measureMSE(truth = scale(dim_redTrue), response = scale(df_feat[[k]]))
    }
    best = which.min(mse)
    dim_red = df_feat[[best]]
    features = a[[best]]
    
    
    df_features_pca[2,1:(2*p)] = as.numeric(colnames(X) %in% names(features))
    df_features_pca[2,which(colnames(X) %in% names(features))] = features
    df_features_pca$latent[2] = 2
    
    #-----------------------------------------------------------------------------------------------
    # Totalvis
    df_totalvis_1 = totalvis:::totalvis(mod, X, pc_num = 1)
    df_totalvis_2 = totalvis:::totalvis(mod, X, pc_num = 2)
    
    
    
    return(list(df_detailed_pca = df_detailed_pca, df_pca = df_pca, df_detailed_sspca = df_detailed_sspca, 
                df_sspca = df_sspca, df_features_pca = df_features_pca, df_features_sspca = df_features_sspca,
                df_totalvis_1 = df_totalvis_1, df_totalvis_2 = df_totalvis_2))
  
})


# add experiment
addExperiments(repls = 1)

# submit jobs
submitJobs()

