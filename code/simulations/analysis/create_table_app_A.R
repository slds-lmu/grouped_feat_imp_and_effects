# Motivational example for grouped feature importance in the case of dummy-encoded features

library(reshape2)
library(mlr)
library(ranger)

# data generation
set.seed(123)

x1 = as.factor(sample(c(sample(1:2, size = 800, replace = TRUE), sample(3:4, size = 200, replace = TRUE)))) #round(rnorm(1000), 1)
x2 = sample(c(sample(1:2, size = 800, replace = TRUE), sample(3:4, size = 200, replace = TRUE)))
x2 = as.factor(x2)

data = data.frame(x1,model.matrix(~ x2 - 1))

data$y = 5*ifelse(data$x1=="2", 1, 0) + 
  5*ifelse(data$x1=="3", 1, 0) + 
  5*ifelse(data$x1=="4", 1, 0) +
  5*(data$x22 + data$x23 + data$x24) +
  rnorm(1000, sd = 1)

# train model
#set.seed(123)
train.ind = sample(1:nrow(data), 0.7*nrow(data))
data.test = data[setdiff(1:nrow(data),train.ind),]

mod = lm(y ~ x1*(x22 + x23 + x24), data = data[train.ind,])
predict.function = function(model, newdata) predict.lm(model, newdata)


summary(mod)


# function to calculate (grouped) permutation feature importance
gfi = function(gfeats, nperm = 50, mod, predict.function){
  gfi = lapply(gfeats, function(feat) {
    fi = c()
    set.seed(123)
    for(i in 1:nperm){
      perm = featureImportance:::permuteFeature(data.test, features = feat)
      pred.unperm = predict.function(mod, data.test)
      pred.perm = predict.function(mod, perm)
      perf.unperm = measureRMSE(truth = data.test$y, response = pred.unperm)
      perf.perm = measureRMSE(truth = data.test$y, response = pred.perm)
      
      fi[i] = perf.perm-perf.unperm
    }
    mean(fi)
  })
  names(gfi) = gfeats
  gfi
}

# encoded features as "one" feature
gfeats = list("x1",paste0("x2", 2:4))
gfi_grouped = gfi(gfeats, mod = mod, predict.function = predict.function)
unlist(gfi_grouped)

# encoded features as "single" features
gfeats = append(list("x1"), paste0("x2", 2:4), after = 1)
gfi_grouped = gfi(gfeats, mod = mod, predict.function = predict.function)
unlist(gfi_grouped)
mean(unlist(gfi_grouped)[2:4])