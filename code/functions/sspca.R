my_chol_psd = function(a){
  n = dim(a)[1];
  root = matrix(0,n,n);
  
  for (i in 1:n){
    a[i,i] = 3*a[i,i]
    sum = 0;
    if (i>1){
      sum = sum(root[i,1:(i-1)]^2);
    }
    
    x = a[i,i] - sum;
    
    if (x<0){
      x = 0;
    }
    
    root[i,i] = sqrt(x);
    
    if (i < n){
      for (j in (i+1):n){
        
        if (root[i,i] == 0){
          x=0;
        }
        else{
          sum = 0;
          if (i>1) {
            sum = root[i,1:(i-1)] %*% t(t(root[j,1:(i-1)]))
          }
          x = (a[i,j] - sum)/root[i,i];
        }
        
        root[j,i] = x;
      }
    }
  }
  return(root);
}

library(PMA)
sspca = function(X, Y, sigma.y, sumabsv = 2, K = 2) {
  if (nrow(X) != nrow(Y))
    stop("X and Y must be the same number of samples")
  

  L = kernelMatrix(kernel = rbfdot(sigma = sigma.y), Y)
  
  n.row = nrow(X)
  
  #SSPCA 1). Decompose L such that L = Δ^TΔ
  Delta = my_chol_psd(L)
  # Delta = chol(L, pivot = TRUE)
  #SSPCA 2). H ← I−n^{−1}ee^T
  H = diag(1, n.row, n.row) - matrix(1, nrow = n.row, ncol = n.row) / n.row
  #SSPCA 3). Ψ ← ΔTHX
  Psi = t(Delta) %*% H %*% X
  # SSPCA 4). Compute the sparse basis based on the PMD method
  
  pmd = SPC(Psi, K = K, orth = TRUE, sumabsv = sumabsv)  
  pmd
}