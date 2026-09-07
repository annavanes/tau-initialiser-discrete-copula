pacman::p_load(ggplot2,copula,Metrics,extraDistr,VineCopula,MASS,tictoc,mvtnorm,doRNG,reshape2,tidyverse,numDeriv,DescTools,zeallot,boot,Rmisc,utils,here,tinytex,weightedCL,pracma,doParallel,foreach,gtools,skellam,bench,pbivnorm,rlist) # load packages

rm(list=ls())   # clean environment

# Load all the necessary functions

source("d2_reparam_start_tau.r")
source("d3plus.r")

# START FROM TRUE TAU WITH ANALYTICAL GRADIENT

params1 = list(
  d = 4,                    # Number of dimensions (species)
  seed = 1
)

params1 = list.append(params1,
                      lambda0 = generate_lambdas(params1$d, seed = 3),          # generate Poisson parameters
                      P = unstructured(dimension = params1$d, seed = params1$seed),              # generate Copula parameters
                      q = 0,                    # Number of observable covariates
                      NSim = 1,              # Number of simulations (toy example)
                      corr.str = "un", # unstructured corr
                      n = 1e2,
                      lam_setting = "small", # moderate, small
                      grad ="anal", # anal, num: analytical or numeric gradient
                      optim_method = "optim",
                      start = "no.approx.tau" # rho.hat, no.approx.tau: rho.hat: start from empirical correlation; no.approx.tau: start from Option 1 -> see paper
)

RTMB::getAll(params1, warn = T)

fit_reparam <- function(Y, start,
                        control = list(maxit = 300, reltol = 1e-8)) {
  d = ncol(Y)
  obj <- function(t) nll_tau_theta(t, Y)
  gr  <- function(t) grad_nll_tau_theta(t, Y)
  
  t1  <- try(optim(start$start.psy, 
                   obj, 
                   gr, 
                   method = "BFGS", 
                   control = control), 
             silent = TRUE)
  
  if (inherits(t1, "try-error")) return(list(converged = FALSE, par = rep(NA_real_, length(start)), par_nat = rep(NA_real_,length(start)), obj = NA_real_, time = NA_real_))
  lam_hat <- exp(t1$par[1:d]); rho_hat <- P2p(zeta_to_R(t1$par[(d+1):length(t1$par)], d))
  list(converged = (t1$convergence == 0), par = t1$par, par_nat = c(lam_hat, rho_hat), obj = t1$value, time = NA_real_, start_nat = start$start.nat)
}

# Simulation function:
mc_test <- function(par) {
  RTMB::getAll(par, warn = T)
  
  # Reconstruct the Correlation Matrix P
  
  rho0 = P2p(P)
  
  # Extract the true theta
  theta = par = c(lambda0,rho0)
  zeta0 = c(log(lambda0), R_to_zeta(P))
  
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
  Rho_hat    <- matrix(NA_real_, NSim, choose(d,2))
  NLL        <- rep(NA_real_, NSim)
  Time       <- rep(NA_real_, NSim)
  Converged  <- rep(FALSE,   NSim)
  Starts <- matrix(NA_real_, NSim, d+choose(d,2))
  
  # Run MC simulations
  for (b in 1:NSim) {
    
    set.seed(b)
    tryCatch({
      # Draw observations from a Gaussian copula with Poisson margins
      Y   <- rGausCopPoisAlld(n, lambda0, rho0)
      
      if (start == "no.approx.tau"){
        start0 <- compute_start_true_tau_perr_alld(Y)
      } else {
        start0 <- compute_start_corr_alld(Y)
      }
      
      
      t0 <- system.time(fr <- fit_reparam(Y, start0))
      Time[b]      <- unname(t0["elapsed"])
      Converged[b] <- isTRUE(fr$converged)
      Lambda_hat[b,] <- fr$par_nat[1:d]
      Rho_hat[b,]     <- fr$par_nat[(d+1):length(fr$par_nat)]
      NLL[b]         <- if (is.finite(fr$obj)) fr$obj else NA_real_
      Starts[b,] <- fr$start_nat
      
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
# inspect the results
raw_out$obj$raw
