simulate_tail <- function(n, m, beta0, delta0, shift, shift_time, delta1, beta1) {
  # Standard deviation of covariates
  sd_cov = 1
  
  # Generate initial data for control and treatment groups
  X_control = cbind(rnorm(n, sd = sd_cov), rnorm(n, sd = sd_cov))
  X_treatment = cbind(rnorm(m, sd = sd_cov) + shift, rnorm(m, sd = sd_cov))
  
  # Combine the control and treatment data
  X = rbind(X_control, X_treatment)
  
  # Compute the mu for control group
  mu_control = (beta0[1] * X[, 1] + beta0[2] * X[, 2])^2
  u_control = rexp(n + m, rate = 1 / (mu_control + 1))
  c_control = rchisq(n + m, df = mu_control + delta0)
  
  # Determine events for control group
  D_control = ifelse(u_control < c_control, 1, 0)
  T_control = ifelse(D_control, u_control, c_control)
  
  # Compute the mu for treatment group
  mu_treatment = (beta1[1] * X[, 1] + beta1[2] * X[, 2])^2
  u_treatment = shift_time + rexp(n + m, rate = 1 / (mu_treatment + 1))
  c_treatment = shift_time + rchisq(n + m, df = mu_treatment +delta1)
  
  # Determine events for treatment group
  D_treatment = ifelse(u_treatment < c_treatment, 1, 0)
  T_treatment = ifelse(D_treatment, u_treatment, c_treatment)
  
  # Print event percentages
  cat("Event percentage for control group (tunable with delta0): ", sum(D_control) / (n + m), "\n")
  cat("Event percentage for treatment group (tunable with delta1): ", sum(D_treatment) / (n + m), "\n")
  
  # Generate final dataset
  z = c(rep(0, n), rep(1, m))
  return(list(
    U0 = u_control,
    U1 = u_treatment,
    Tstar = (1 - z) * T_control + z * T_treatment,
    delta = (1 - z) * D_control + z * D_treatment,
    x = X,
    z = z,
    T0 = T_control,
    T1 = T_treatment,
    C0 = c_control,
    C1 = c_treatment,
    D0 = D_control
  ))
}
