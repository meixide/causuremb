# Causal Survival Embeddings

Code for the paper:

> **Causal survival embeddings: Non-parametric counterfactual inference under right-censoring**  
> Carlos García Meixide, Marcos Matabuena  
> *Statistical Methods in Medical Research*, 34(3):574–593, 2025  
> [10.1177/09622802241311455](https://doi.org/10.1177/09622802241311455) · [PMC11951469](https://pmc.ncbi.nlm.nih.gov/articles/PMC11951469/)

## Overview

Causal Survival Embeddings (CSE) is a non-parametric estimator of counterfactual survival functions under right-censoring. It leverages reproducing kernel Hilbert spaces (RKHS) and Kaplan–Meier inverse-probability weights to adjust for confounding in observational studies — without propensity score modeling or semi-parametric assumptions.

The method estimates the counterfactual parameter

$$P(\tilde{T}_1 > t) - P(\tilde{T}_0 > t), \quad t > 0$$

and decomposes observed survival differences into covariate composition effects and structural treatment effects (following Chernozhukov et al., 2013).

## Requirements

R (≥ 4.0) with the following packages:

```r
install.packages(c("Rcpp", "RcppArmadillo", "survival", "Matching"))
```

## Files

| File | Description |
|------|-------------|
| `count_me.R` | Core CSE estimator: takes observed times, event indicators, treatment assignments, and covariates; returns the estimated counterfactual survival embedding and its associated weights |
| `simulate_tail.R` | Data-generating process used in the paper's simulations |
| `main.R` | End-to-end example: simulate confounded data, tune λ via 5-fold cross-validation, compute the CSE estimate, and plot results |

## Usage

```r
source("main.R")
```

This script simulates a confounded observational survival dataset, selects the regularization parameter λ by cross-validation, and plots the result.

**Plot color code:** black = KM control arm · blue = KM treatment arm · red = causal survival embedding

To use the estimator on your own data:

```r
source("count_me.R")

result <- count_me(
  t      = observed_times,       # numeric vector of observed times T* = min(T, C)
  delta  = event_indicators,     # 1 = event, 0 = censored
  z      = treatment_assignment, # 1 = treatment, 0 = control
  x      = covariates,           # n x p covariate matrix
  lgrid  = 100,                  # number of time grid points
  tmax   = max(observed_times),
  lambda = 0.1,                  # regularization parameter (tune via CV)
  sigma2 = sigma2_kernel,        # bandwidth for the time kernel
  tmin   = 0
)

# result$estimator : (lgrid x 1) counterfactual mean embedding on the time grid
# result$weights   : (n_control x 1) weights for reconstructing the survival function
```

The survival function is then recovered as `1 - cumsum(weights[order(t_control)] ...)` — see `main.R` for a complete working example.

## Citation

```bibtex
@article{garciameixide2025causal,
  title   = {Causal survival embeddings: {N}on-parametric counterfactual inference under right-censoring},
  author  = {Garc{\'i}a Meixide, Carlos and Matabuena, Marcos},
  journal = {Statistical Methods in Medical Research},
  volume  = {34},
  number  = {3},
  pages   = {574--593},
  year    = {2025},
  doi     = {10.1177/09622802241311455}
}
```

## Contact

Open an issue for questions or collaboration inquiries.
