# Include RcppArmadillo dependency for C++ function
Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x); }", depends = "RcppArmadillo")

# Define Gaussian kernel function
gaussian_kernel <- function(x, y, sigma2) {
  if (length(x) > 1 || length(y) > 1) {
    return('Only scalars are accepted as inputs!')
  }
  return(exp(-((x - y)^2) / sigma2))
}

library('survival')

# Define the main function
count_me <- function(t, delta, z, x, lgrid, tmax, lambda, sigma2, tmin) {
  # Initialize variables
  t_grid = seq(tmin, tmax, length.out = lgrid)
  
  X_control = x[which(z == 0), ]
  X_treatment = x[which(z == 1), ]
  
  T_control = t[which(z == 0)]
  D_control = delta[which(z == 0)]
  
  n = length(T_control)
  m = nrow(X_treatment)
  
  # Reverse Kaplan-Meier estimator for control group
  km_control <- survfit(Surv(T_control, 1 - D_control) ~ 1)
  Ghat = km_control$surv[match(T_control, km_control$time)]
  
  weights = ifelse(D_control, 1 / Ghat, 0)
  W = diag(weights)
  
  # Kernel matrix for control group
  dists_squared = as.matrix(dist(X_control))^2
  sigma2x = median(dists_squared)
  K = exp(-(1 / sigma2x) * dists_squared)
  
  # Initialize Y matrix
  Y = matrix(0, nrow = n, ncol = lgrid)
  
  for (i in 1:n) {
    Y[i, ] = as.numeric(lapply(t_grid, function(t) gaussian_kernel(T_control[i], t, sigma2)))
  }
  
  dists_squared_treatment = as.matrix(dist(rbind(X_control, X_treatment)))[1:n, (n + 1):(n + m)]^2
  K_treatment = exp(-(1 / sigma2x) * dists_squared_treatment)
  
  # Calculate estimator
  inverse_matrix = armaInv(diag(rep(1, n)) * lambda + W %*% K)
  C = inverse_matrix %*% W %*% Y
  estimator = t(C) %*% K_treatment %*% rep(1 / m, m)
  
  # Calculate weights (actually, from which the estimator is computed)
  wf = rep(1 / m, m) %*% t(K_treatment) %*% inverse_matrix %*% W
  
  return(list(estimator = estimator, weights = as.numeric(wf)))
}
