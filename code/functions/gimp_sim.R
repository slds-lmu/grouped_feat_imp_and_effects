# Gimp adjusted for Simulation Setting in Section 4 (only Effects)

library(dplyr)
library(mlr)
library(dplyr)
library(elasticnet)

library(R6)
Gimpsim = R6Class("grouped_imp_sim",
  public = list(
    initialize = function(task, res, mod, lrn) {
      private$.task = task
      private$.res = res
      private$.mod = mod
      private$.lrn = lrn
    },
    group_pdp = function(features, featuresTrue, parts = 4, x = 0, y = 0) {
      #df_features = data.frame(feature = names(features), variable = paste0("x", 1:length(features)))
      data = getTaskData(self$task)
      library(ggplot2)
      library(gridExtra)
      library(future)
      plan(multisession, workers = 2L) #or plan(sequential)
      mod = self$mod
      task = self$task
      mse = c()
      df = list()
      dim_redTrue = (data %>% data.frame()  %>% dplyr::select(names(featuresTrue)) %>% as.matrix()) %*% featuresTrue
      for(k in 1:length(features)){
        df[[k]] = (data %>% scale() %>% data.frame() %>% dplyr::select(names(features[[k]])) %>% as.matrix()) %*% features[[k]]
        mse[k] = measureMSE(truth = scale(dim_redTrue), response = scale(df[[k]]))
      }
      best = which.min(mse)
      dim_red = df[[best]]
      features = features[[best]]
      res = future.apply::future_lapply(1:nrow(data), 
        function(i) {
          
          # based on dim red weights
          data_pdp = data
          data_pdp[-i, as.character(names(features))] = data[i, as.character(names(features))]
          pred = predict(mod, newdata = data_pdp)
          resp = switch(getTaskType(task), classif = pred$data[[1]], regr = pred$data$response)
          
          # based on true weights
          data_pdp = data
          data_pdp[-i, as.character(names(featuresTrue))] = data[i, as.character(names(featuresTrue))]
          predTrue = predict(mod, newdata = data_pdp)
          respTrue = switch(getTaskType(task), classif = predTrue$data[[1]], regr = predTrue$data$response)
          
          
          pdp_i = data.frame(index = i, pd = mean(resp), pdTrue = mean(respTrue))
          pdp_i$dim_red = (data %>% scale() %>% data.frame() %>% slice(i) %>% dplyr::select(names(features)) %>% as.matrix()) %*% features
          pdp_i$res_true = (pdp_i$dim_red)
          pdp_i$dim_redTrue = (data %>% data.frame() %>% slice(i) %>% dplyr::select(names(featuresTrue)) %>% as.matrix()) %*% featuresTrue
          pdp_i$res_truetrue = (pdp_i$dim_redTrue)
          pdp_i$feat_numb = length(names(features)[which(abs(features)>0.05)])
          pdp_i$feat_numbTrue = length(names(featuresTrue))
          pdp_i$feat_numbBoth = length(which(names(features)[which(abs(features)>0.05)] %in% names(featuresTrue)))
          pdp_i$feat_v5 = length(which(names(features)[which(abs(features)>0.05)] %in% names(featuresTrue)[1]))
          pdp_i$feat_v25 = length(which(names(features)[which(abs(features)>0.05)] %in% names(featuresTrue)[2]))
          pdp_i$feat_v45 = length(which(names(features)[which(abs(features)>0.05)] %in% names(featuresTrue)[3]))
          pdp_i 
        }
      )
      plan(sequential)
      pdp = bind_rows(res)
      
      return(list(pdp = pdp))
    }
  ),
  private = list(
    .task = NULL,
    .res = NULL,
    .feat_imp = NULL,
    .mod = NULL,
    .lrn = NULL
  ),
  active = list(
    mod = function() private$.mod,
    task = function() private$.task,
    lrn = function() private$.lrn
  )
)


