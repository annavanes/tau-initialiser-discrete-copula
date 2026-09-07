# tau-initialiser-discrete-copula
Implementation code for van Es, Anna and Cantoni, Eva. "Novel tau-informed initialization for maximum likelihood estimation of copulas with discrete margins" Dependence Modeling, vol. 14, no. 1, 2026, pp. 20250020. https://doi.org/10.1515/demo-2025-0020                     

## R
The folder "R" contains all simulation codes:

**d2_reparam_start_tau.R**

Functions needed for simulations with $d=2$. 

**d3plus.R**

Functions needed for simulations with $d>2$. 

**Option1d2.R**

Example of simulation for $d=2$, using the analytical gradient and Spherical-Cholesky reparametrization and starting from either the empirical correlation (set `params1$start = "rho.hat"`) or Option 1 (set `params1$start = "anal.tau"`). 

**Option1Alld.R**

Example of simulation for $d>2$, using the analytical gradient and Spherical-Cholesky reparametrization and starting from either the empirical correlation (set `params1$start = "rho.hat"`) or Option 1 (set `params1$start = "no.approx.tau"`). 


## Selected references
Perrone, E., Van Den Heuvel, E.R., and Zhan, Z. (2023). Kendall’s tau estimator for bivariate zero-inflated count data. Stat. Probab. Lett. 199:
109858,.

Pimentel, R.S. (2009). Kendall’s tau and Spearman’s rho for zero-inflated data, Ph.D. thesis. Kalamazoo, Western Michigan University.

Tong, Y.L. (2012). The multivariate normal distribution. Springer Science & Business Media, New York.

Lucchetti, R. and Pedini, L. (2024). The spherical parametrisation for correlation matrices and its computational advantages. Comput.
Econ. 64: 1023−1046,.