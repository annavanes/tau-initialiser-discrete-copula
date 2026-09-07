pacman::p_load(ggplot2,copula,Metrics,extraDistr,VineCopula,MASS,tictoc,mvtnorm,doRNG,reshape2,tidyverse,numDeriv,DescTools,zeallot,boot,Rmisc,utils,here,tinytex,weightedCL,pracma,doParallel,foreach,gtools,skellam,bench,pbivnorm,rlist) # load packages

rm(list=ls())   # clean environment

# Load all the necessary functions

source("d2_reparam_start_tau.r")
source("d3plus.r")

params1 = list(
  d = 2,                    # Number of dimensions (species)
  seed = 1
)

params1 = list.append(params1,
                      lambda0 = c(3,8),          # Poisson parameters
                      P = p2P(param = -0.5, d = params1$d), # Copula parameters
                      q = 0,                    # Number of observable covariates
                      NSim = 1,              # Number of simulations (toy example)
                      corr.str = "ex",         # corr structure
                      n = 1e2,
                      lam_setting = "big", # moderate, small, big
                      grad ="anal", # anal, num
                      optim_method = "optim",
                      start = "anal.tau" # rho.hat, anal.tau: rho.hat: start from empirical correlation; anal.tau: start from Option 1 -> see paper
)

# Reparametrized fit on theta(zeta) = (gamma1, gamma2, eta), example for d=2
fit_reparam <- function(Y, start = NULL) {
  obj <- function(t) negloglik_theta(t, Y)
  gr  <- function(t) grad_negloglik_theta(t, Y)
  t1  <- try(optim(start, obj, gr, method = "BFGS", control = list(maxit = 2000, reltol = 1e-8)), silent = TRUE)
  #t1.time = Sys.time()
  if (inherits(t1, "try-error")) return(list(converged = FALSE, par = rep(NA_real_, 3), par_nat = rep(NA_real_,3), obj = NA_real_, time = NA_real_))
  lam_hat <- exp(t1$par[1:2]); rho_hat <- rho_from_theta(t1$par[3])
  list(converged = (t1$convergence == 0), par = t1$par, par_nat = c(lam_hat, rho_hat), obj = t1$value, time = NA_real_)
}

# Simulation function:
mc_test <- function(par) {
  RTMB::getAll(par, warn = T)
  
  # Reconstruct the Correlation Matrix P
  
  rho0 = P2p(P)
  
  # Extract the true theta
  theta = par = c(lambda,rho0)
  
  # record start time
  t.start <- Sys.time()
  
  # Setup for parallel processing
  no_cores <- detectCores(logical = TRUE)
  cl <- makeCluster(min(no_cores, 4))  # Use 4 cores or fewer if not available
  registerDoParallel(cl)
  
  # Set up the progress bar
  pb = txtProgressBar(min = 0, max = NSim, style = 3)
  
  # Initialize
  na <- 0
  params.hat = matrix(NA, NSim, length(par))
  hessians = cov.matrices = list()
  Hess.mineig = Hess.kappa = numeric(NSim)
  standard_errors = matrix(NA, nrow = NSim, ncol = length(theta))
  
  Lambda_hat <- matrix(NA_real_, NSim, d)
  Rho_hat    <- rep(NA_real_, NSim)
  NLL        <- rep(NA_real_, NSim)
  Time       <- rep(NA_real_, NSim)
  Converged  <- rep(FALSE,   NSim)
  Starts <- matrix(NA_real_, NSim, d+1)
  
  # Run MC simulations
  for (b in 1:NSim) {
    
    set.seed(b)
    tryCatch({
      # Draw observations from a Gaussian copula with Poisson margins
      Y   <- rGausCopPois(n, lambda0, rho0)
      if (start == "rho.hat"){
        start0 = compute_start_corr(Y)
      } else {
        start0 <- compute_start_true_tau(Y)
      }
      
      t0 <- system.time(fr <- fit_reparam(Y, start = start0$start.psy))
      Time[b]      <- unname(t0["elapsed"])
      Converged[b] <- isTRUE(fr$converged)
      Lambda_hat[b,] <- fr$par_nat[1:2]
      Rho_hat[b]     <- fr$par_nat[3]
      NLL[b]         <- if (is.finite(fr$obj)) fr$obj else NA_real_
      Starts[b,] <- start0$start.nat
      
    }, error = function(e) {
      cat("Error in simulation", b, ": ", e$message, "\n")
    })
    setTxtProgressBar(pb, value = b)
  }
  
  # Stop cluster after use
  stopCluster(cl)
  close(pb)
  
  # Run time
  ptm.end = Sys.time()
  time.taken = ptm.end - t.start
  
  obj <- list(
    raw  = list(Lambda_hat = Lambda_hat, Rho_hat = Rho_hat, NLL = NLL, Time = Time, Converged = Converged, Starts = Starts),
    meta = list(B = NSim, n = n, lambda0 = lambda0, rho0 = rho0, method = "reparam", seed = seed),
    file = list(saved_as = NA_character_)
  )
  
  return(list(obj = obj,
              params1 = params1,
              theta = theta,
              true_lambda = params1$lambda,
              true_rhos = rho0,
              time.taken = time.taken
  ))
  as.list(environment())
}

raw_out = mc_test(par = params1)

# inspect results
raw_out$obj$raw
