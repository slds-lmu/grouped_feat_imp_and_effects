#helper functions

# Greedy Forward Search
make_correlated_feats = function(n, rv, name) {
  for (i in 1:n) {
    if (i == 1) {
      U_i = rv
      N = length(rv)
      random_10 = sample(1:N, round(0.1*N))
      U_i[random_10] = rv[random_10] + rnorm(round(0.1*N), 0, 0.5)
      df = setNames(data.frame(U_i), paste0(name, i))
    } else {
      U_i = rv
      random_10 = sample(1:N, round(0.1*N))
      U_i[random_10] = rv[random_10] + rnorm(round(0.1*N), 0, 0.5)
      df = data.frame(df, setNames(data.frame(U_i), paste0(name, i)))
    }
  }
  df
}

sim_data = function(N = 100, n_groups, y = numeric(), size_groups = numeric()) {
  for (i in 1:length(size_groups)) {
    x = rnorm(N, 0, 1)
    if (i == 1) {
      df = make_correlated_feats(size_groups[i], rv = x, name = paste0("Group_", i, "_feat_"))
      rv = list(x)
      next
    }
    df = data.frame(df, make_correlated_feats(size_groups[i], rv = x, name = paste0("Group_", i, "_feat_")))
    rv = c(rv, list(x))
  }
  checkmate::assert_true(length(y) == n_groups)
  for (i in 1:length(y)) {
    if (i == 1) {
      a = list(y[i] * rv[[i]])
      z = a[[1]]
      next
    }
    a = c(a, list(y[i] * rv[[i]]))
    z = z + a[[i]]
  }
  target =  setNames(data.frame(z - mean(z) + rnorm(N, 0, 0.1)), "target")
  data.frame(df, target)
}


#---------------------------------------------------------------------------------------------------
# Other simulation examples and gimp_sim

get_df2 = function(imp) {
  output_imp = testthat::capture_output(imp, print = TRUE, width = 200)
  x = strsplit(output_imp, split = "\n")[[1]]
  x = x[7:length(x)]
  res = lapply(x, function(y) {
    while (grepl("  ", y)) y = gsub("  ", " ", y)
    z = strsplit(y, split = " ")[[1]]
    z2 = strsplit(y, split = ",")[[1]]
    z = z[c(length(z) - 1, length(z))]
    data.frame(feat = z[1], imp = as.numeric(z[2]))
  })
  bind_rows(res)
}

add_group2 = function(imp, gimp) {
  df_imp = get_df2(imp)
  output_gimp = testthat::capture_output_lines(imp, print = TRUE, width = 200)
  z = output_gimp[7:length(output_gimp)]
  grouped_imp = bind_rows(lapply(1:length(z), function(x) {
    print(x)
    x = z[x]
    while (grepl("  ", x)) x = gsub("  ", " ", x)
    y = strsplit(x, split = ",")[[1]]
    data.frame(feat = y[length(y) - 1])
  }
  )
  )
  grouped_imp$group = ""
  for (i in 1:nrow(grouped_imp)) {
    feat = grouped_imp$feat[i]
    grouped_imp$group[i] = gimp$group_df %>% filter(feature == feat) %>% pull(group) %>% as.character()
  }
  grouped_imp$imp = df_imp$imp[7:length(df_imp$imp)]
  grouped_imp$feat = NULL
  grouped_imp
}



library(mlr)
make_correlated_feats_sim = function(n, rv, name, prop) {
  for (i in 1:n) {
    if (i == 1) {
      U_i = rv
      random_20 = sample(1:length(U_i), prop*length(U_i))
      #U_i[random_20] = rnorm(10, 0, 1)
      U_i[random_20] = 0.2*rv[random_20] + 0.8*rnorm(prop*length(U_i), 0, 1)
      df = setNames(data.frame(U_i), paste0(name, i))
    } else {
      U_i = rv
      random_20 = sample(1:length(U_i), prop*length(U_i))
      U_i[random_20] = 0.2*rv[random_20] + 0.8*rnorm(prop*length(U_i), 0, 1)
      df = data.frame(df, setNames(data.frame(U_i), paste0(name, i)))
    }
  }
  df
}

make_correlated_feats_unif_sim = function(n, rv, name, prop) {
  for (i in 1:n) {
    if (i == 1) {
      U_i = rv
      random_20 = sample(1:length(U_i), prop*length(U_i))
      #U_i[random_20] = rnorm(10, 0, 1)
      U_i[random_20] = 0.2*rv[random_20] + 0.8*runif(prop*length(U_i), 0, 1)
      df = setNames(data.frame(U_i), paste0(name, i))
    } else {
      U_i = rv
      random_20 = sample(1:length(U_i), prop*length(U_i))
      U_i[random_20] = 0.2*rv[random_20] + 0.8*runif(prop*length(U_i), 0, 1)
      df = data.frame(df, setNames(data.frame(U_i), paste0(name, i)))
    }
  }
  df
}




sim_data_sim = function(N = 100, n_groups, y = numeric(), size_groups = numeric(), prop = numeric()) {
  x = rnorm(N, 0, 1)
  for (i in 1:length(size_groups)) {
    
    if (i == 1) {
      df = make_correlated_feats_sim(size_groups[i], rv = x, name = paste0("Group_", i, "_feat_"), prop = prop[i])
      rv = list(x)
      next
    }
    df = data.frame(df, make_correlated_feats_sim(size_groups[i], rv = x, name = paste0("Group_", i, "_feat_"), prop = prop[i]))
    rv = c(rv, list(x))
  }
  checkmate::assert_true(length(y) == n_groups)
  for (i in 1:length(y)) {
    if (i == 1) {
      a = list(y[i] * rv[[i]])
      z = a[[1]]
      next
    }
    a = c(a, list(y[i] * rv[[i]]))
    z = z + a[[i]]
  }
  #target =  setNames(data.frame(z - mean(z) + rnorm(N, 0, 0.1)), "target")
  data.frame(df)
}


