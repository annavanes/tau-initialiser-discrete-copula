# generate a vector of lambdas for simulation
generate_lambdas <- function(d, seed){
  # Set the vector of lambdas we're interested in
  ls <- c(seq(0.2,0.9,by=0.1), 1:5)
  
  set.seed(seed)
  
  # Generate the lambdas vector
  lambdas <- sample(x = ls, size = d, replace = T)
  return(lambdas)
}

# generate an unstructured correlation matrix
unstructured <- function(dimension, seed) {
  set.seed(seed)
  temp <- matrix(NA,dimension,dimension)
  diag(temp) <- 1
  for (i in 2:dimension) {
    for (j in 1:(i-1)) {
      temp[i,j] <- runif(1,-.9,.9)
      temp[j,i] <- temp[i,j]
    }
  }
  return(round(pc2c(temp),2))
}

# from https://github.com/cran/corcounts
# helper function
pc2c <- function(Theta) {
  n <- dim(Theta)[2]
  Theta.alt <- Theta
  for (k in 2:(n-1)) {      # Anfangszeile
    for (i in (n-1):k) {    # Zeile (beginne unten)
      for (j in n:(i+1)) {  # Spalte (beginne rechts)
        Theta[i,j] <- Theta[i,j]*sqrt((1-Theta.alt[i-k+1,i]^2)*(1-Theta.alt[i-k+1,j]^2)) +
          Theta.alt[i-k+1,i]*Theta.alt[i-k+1,j]
      }
    }
  }
  
  C <- diag(rep(1,n))
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      C[i,j] <- Theta[j,i]
      C[j,i] <- Theta[j,i]
    }
  }
  return(C)
}

## ---------- utilities ----------

## zeta -> R (general d).  Angles w_{i,1:(i-1)} = pi * sigmoid(zeta_row)
zeta_to_R <- function(zeta, d) {
  stopifnot(length(zeta) == d*(d-1)/2)
  C <- diag(d)
  idx <- 1L
  for (i in 2:d) {
    wi <- pi * sigmoid(zeta[idx:(idx + i - 2)])
    idx <- idx + i - 1
    si <- sin(wi); ci <- cos(wi)
    
    ai <- numeric(i)
    prod_sin <- 1.0
    # j = 1..i-1
    for (j in 1:(i-1)) {
      ai[j] <- ci[j] * prod_sin
      prod_sin <- prod_sin * si[j]
    }
    # j = i
    ai[i] <- prod_sin
    
    C[i, 1:i] <- ai
  }
  R <- C %*% t(C)
  (R + t(R))/2
}

## R -> zeta (general d).  Compute L = lower Cholesky, read off a_{i,j} = L[i,j],
## then recover each row's angles via the recursive inverse, finally logit-map to R.
R_to_zeta <- function(R) {
  stopifnot(is.matrix(R), nrow(R) == ncol(R))
  d <- nrow(R)
  L <- t(chol(R))  # lower-triangular
  zeta <- numeric(d*(d-1)/2)
  idx <- 1L
  for (i in 2:d) {
    ai <- L[i, 1:i]
    wi <- numeric(i-1)
    denom <- 1.0
    for (j in 1:(i-1)) {
      x <- ai[j] / denom
      x <- pmin(pmax(x, -1), 1)            # clip for safety
      wi[j] <- acos(x)
      sj <- sin(wi[j])
      denom <- denom * if (sj < 1e-15) 1e-15 else sj
    }
    # map (0,pi) -> R via logit on w/pi
    wclip <- pmin(pmax(wi, 1e-9), pi - 1e-9)
    zeta[idx:(idx + i - 2)] <- log(wclip) - log(pi - wclip)
    idx <- idx + i - 1
  }
  zeta
}

pbvn <- function(z1, z2, rho) {
  if (requireNamespace("pbivnorm", quietly = TRUE)) {
    return(pbivnorm::pbivnorm(z1, z2, rho))
  } else {
    # Fallback: mvtnorm 2D
    return(mvtnorm::pmvnorm(
      lower = c(-Inf, -Inf), upper = c(z1, z2),
      mean  = c(0, 0),
      sigma = matrix(c(1, rho, rho, 1), 2, 2),
      algorithm = mvtnorm::Miwa()
    )[1])
  }
}
# theta (zeta)-variant spherical–Cholesky reparam for correlation matrices
# Stable analytic J_{rho,theta} and NLL/gradient on (tau, theta)
# Works for general d; d=3 path uses fast corner integrator.

# ------------------------ Utilities & indexing ------------------------------
sigmoid <- function(x) 1/(1+exp(-x))

pairs_P2p <- function(d) {
  # lower triangle in column order: (2,1),(3,1),...,(d,1),(3,2),...,(d,d-1)
  do.call(rbind, lapply(seq_len(d-1L), function(j) cbind(seq.int(j+1L, d), j)))
}

pair_to_index_P2p <- function(i, j, d) {
  if (i < j) { tmp <- i; i <- j; j <- tmp }
  offset <- sum(d - seq_len(j-1L))  # seq_len(0) -> integer(0) => sum 0
  offset + (i - j)
}

# -------------------- theta (zeta) <-> R (spherical–Cholesky), general d ---------------
# theta (zeta) maps unconstrained R^m to angles w = pi * sigmoid(theta (zeta)) in (0,pi)
# Row i has angles w_{i,1:(i-1)}, and C_{i,1:i} = a_i(w_i) 

theta_to_R <- function(theta, d, sin_floor = 1e-15) {
  stopifnot(length(theta) == d*(d-1)/2)
  C <- diag(d)
  idx <- 1L
  for (i in 2:d) {
    wi <- pi * sigmoid(theta[idx:(idx + i - 2L)])
    idx <- idx + i - 1L
    si <- pmax(sin(wi), sin_floor); ci <- cos(wi)
    ai <- numeric(i); ps <- 1.0
    for (j in 1:(i-1L)) { ai[j] <- ci[j] * ps; ps <- ps * si[j] }
    ai[i] <- ps
    C[i, 1:i] <- ai
  }
  R <- C %*% t(C)
  (R + t(R))/2
}

R_to_theta <- function(R, sin_floor = 1e-15) {
  d <- nrow(R); stopifnot(d == ncol(R))
  L <- t(chol(R))
  theta <- numeric(d*(d-1)/2); idx <- 1L
  for (i in 2:d) {
    ai <- L[i, 1:i]
    wi <- numeric(i-1L)
    denom <- 1.0
    for (j in 1:(i-1L)) {
      x <- ai[j] / denom; x <- pmin(pmax(x, -1), 1)
      wj <- acos(x); sj <- max(sin(wj), sin_floor)
      wi[j] <- wj; denom <- denom * sj
    }
    wclip <- pmin(pmax(wi, 1e-9), pi - 1e-9)
    theta[idx:(idx + i - 2L)] <- log(wclip) - log(pi - wclip)  # theta (zeta) = logit(w/pi)
    idx <- idx + i - 1L
  }
  theta
}

# ----------------------- Row a(w) and stable derivatives --------------------
# a_i components (length i) from row angles w (length i-1)
row_a_from_omega <- function(w) {
  i <- length(w) + 1L
  if (i == 1L) return(list(a = 1))
  S <- cumprod(c(1, sin(w)))              # S[j] = prod_{k<j} sin(w_k), length i
  a <- numeric(i)
  a[1:(i-1L)] <- cos(w) * S[1:(i-1L)]
  a[i]        <- S[i]
  list(a = a, S = S)
}

# STABLE derivative matrix J where J[k,s] = da_k / dw_s (i x (i-1)).
row_dalpha_domega_stable <- function(w) {
  i <- length(w) + 1L
  if (i == 1L) return(matrix(0, 1, 0))
  sw <- sin(w); cw <- cos(w)
  
  # prefix products of sin (for S_k = prod_{t<k} sin(w_t))
  S <- cumprod(c(1, sw))   # length i
  
  J <- matrix(0.0, nrow = i, ncol = i-1L)
  for (s in 1:(i-1L)) {
    # k == s < i: d = - sin(w_s) * prod_{t< s} sin(w_t)
    if (s <= (i-1L)) {
      J[s, s] <- - sw[s] * S[s]
    }
    # k in {s+1, ..., i-1}: d = cos(w_k) * [prod_{t<k, t!=s} sin(w_t)] * cos(w_s)
    if (s + 1L <= i - 1L) {
      for (k in (s+1L):(i-1L)) {
        prod_excl <- 1.0
        if (k > 2L) {
          for (t in 1:(k-1L)) if (t != s) prod_excl <- prod_excl * sw[t]
        }
        J[k, s] <- cw[k] * prod_excl * cw[s]
      }
    }
    # k == i: d = [prod_{t< i, t!= s} sin(w_t)] * cos(w_s)
    prod_excl_full <- 1.0
    if (i > 2L) {
      for (t in 1:(i-1L)) if (t != s) prod_excl_full <- prod_excl_full * sw[t]
    }
    J[i, s] <- prod_excl_full * cw[s]
  }
  J
}

# ----------------------- Analytic J_{rho,theta} (general d) -----------------
jacobian_rho_theta <- function(theta, d) {
  # Build per-row objects
  a_list <- vector("list", d)
  Jrow_list <- vector("list", d)
  idx <- 1L
  for (r in 1:d) {
    if (r == 1L) { a_list[[r]] <- 1; Jrow_list[[r]] <- matrix(0,1,0); next }
    wr <- pi * sigmoid(theta[idx:(idx + r - 2L)])
    idx <- idx + r - 1L
    ar <- row_a_from_omega(wr)
    a_list[[r]]   <- ar$a
    Jrow_list[[r]] <- row_dalpha_domega_stable(wr)
  }
  
  pairs <- pairs_P2p(d); m <- nrow(pairs)
  J_omega <- matrix(0.0, nrow = m, ncol = m)
  
  zcol <- 1L
  for (r in 2:d) {
    Jr <- Jrow_list[[r]]
    # angles in row r: s = 1..r-1
    for (s in 1:(r-1L)) {
      # (i=r, j<r)
      for (j in 1:(r-1L)) {
        idxp <- pair_to_index_P2p(r, j, d)
        J_omega[idxp, zcol] <- sum(Jr[1:j, s] * a_list[[j]][1:j])
      }
      # (i>r, j=r)
      if (r + 1L <= d) {
        for (i in (r+1L):d) {
          idxp <- pair_to_index_P2p(i, r, d)
          J_omega[idxp, zcol] <- sum(a_list[[i]][1:r] * Jr[1:r, s])
        }
      }
      zcol <- zcol + 1L
    }
  }
  
  # chain to theta (zeta): w = pi * sigmoid(theta (zeta))
  s <- sigmoid(theta)
  J_omega %*% diag(pi * s * (1 - s), nrow = length(theta))
}

## ---- h(y) in d=3 (Miwa + sigma=R) ----
h_gaus_pois_d3_fast <- function(corners, R) {
  n <- dim(corners$grid2)[1]; d <- nrow(R)
  probs <- numeric(n)
  for (i in 1:n) {
    ms <- corners$nonzero[[i]]
    if (length(ms) == 0L) { probs[i] <- 0; next }
    pr <- vapply(ms, function(m)
      mvtnorm::pmvnorm(lower = rep(-Inf, d),
                       upper = corners$grid2[i, , m],
                       mean  = rep(0, d),
                       sigma = R,
                       algorithm = ALG3D)[1],
      numeric(1L))
    probs[i] <- sum(corners$sign[ms] * pr)
  }
  pmax(pmin(probs, 1), 0)
}

# ---------------------- Likelihood & gradient on (tau (=log(lambda)), theta (zeta)) ---------------------
grad_loglik_fast <- function(Y, lambda, R, eps = 1e-12) {
  # stopifnot(ncol(Y) == 3L)
  d = ncol(Y)
  corners <- build_corner_data(Y, lambda, eps) # ok
  hvals   <- h_gaus_pois_d3_fast(corners, R) # ok
  wts     <- 1 / pmax(hvals, 1e-300)
  
  # d/d lambda
  lam_blocks <- precomp_lambda_blocks(R) # ok
  g_lambda <- numeric(d)
  for (k in 1:d) {
    dh <- dH_dlambda_k_fast(k, corners, lambda, lam_blocks[[k]], Y) # ok
    g_lambda[k] <- sum(wts * dh)
  }
  
  # d/d rho in EXACT copula::P2p order
  rho_pre <- precomp_rho_blocks(R)  # pairs = (1,2), (1,3), (2,3),...
  g_rho <- numeric(nrow(rho_pre$pairs))
  for (p in seq_len(nrow(rho_pre$pairs))) {
    dh <- dH_drho_pair_fast(corners, R, rho_pre$blocks[[p]]) #ok
    g_rho[p] <- sum(wts * dh)
  }
  
  list(grad_lambda = g_lambda, grad_rho = g_rho, pairs = rho_pre$pairs, hvals = hvals)
}

## ---- NLL wrapper that uses copula::p2P ----
negloglik_h_d3_fast <- function(par, Y) {
  d <- 3L
  lambda <- par[1:d]
  rhos   <- par[(d+1):length(par)]
  if (any(lambda <= 0)) return(Inf)
  R <- copula::p2P(rhos)                # <- keep the copula mapping
  # cheap checks
  if (!isTRUE(all.equal(diag(R), rep(1, d)))) return(Inf)
  if (any(abs(R - t(R)) > 1e-10)) return(Inf)
  if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 0) return(Inf)
  
  corners <- build_corner_data(Y, lambda, 1e-12)
  hvals   <- h_gaus_pois_d3_fast(corners, R)
  nll     <- -sum(log(pmax(hvals, 1e-300)))
  
  G <- grad_loglik_fast(Y, lambda, R)  # same P2p order
  attr(nll, "gradient") <- - c(G$grad_lambda, G$grad_rho)
  nll
}

loglik_and_grad_tau_theta <- function(Y, tau, theta) {
  d <- ncol(Y); m <- d*(d-1)/2
  lambda <- exp(tau)
  R <- theta_to_R(theta, d)
  
  G <- grad_loglik_fast(Y, lambda, R)
  gL <- G$grad_lambda; gR <- G$grad_rho
  
  g_tau_ll    <- gL * lambda
  J_rho_theta <- jacobian_rho_theta(theta, d)
  g_theta_ll  <- as.vector(t(J_rho_theta) %*% gR)
  
  list(grad_tau = g_tau_ll, grad_theta = g_theta_ll)
}

grad_nll_tau_theta <- function(par_psy, Y) {
  d <- ncol(Y)
  tau   <- par_psy[1:d]
  theta <- par_psy[(d+1):length(par_psy)]
  lg <- loglik_and_grad_tau_theta(Y, tau, theta)   # (+) log-lik grads
  c(-lg$grad_tau, -lg$grad_theta)                  # flip sign for NLL
}

# ----------------------- Gradient check helper ------------------------------
check_grad_theta <- function(Y, tau, theta) {
  f <- function(p) nll_tau_theta(p, Y)
  g_num <- numDeriv::grad(f, c(tau, theta))
  g_an  <- grad_nll_tau_theta(c(tau, theta), Y)
  list(max_abs_diff = max(abs(g_num - g_an)), g_num = g_num, g_an = g_an)
}

# ====================== Robust fallback for singular R in d=3 ======================
# pmvnorm Miwa fails on near-singular correlation matrices; provide a safe fallback
# to Genz–Bretz (higher tolerance, more robust) and wire it into the theta (zeta)-variant NLL.

safe_pmvnorm <- function(upper, R, abseps = 1e-9, maxpts = 1e6) {
  d=ncol(R)
  evmin <- try(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values), silent = TRUE)
  if (inherits(evmin, "try-error") || !is.finite(evmin) || evmin < 1e-10) {
    out = mvtnorm::pmvnorm(lower = rep(-Inf, d), upper = upper,
                           sigma = R,
                           algorithm = mvtnorm::GenzBretz(abseps = abseps, maxpts = maxpts))[1]
    #attr(out, "algo") = mvtnorm::GenzBretz(abseps = abseps, maxpts = maxpts)
    return(out)
  }
  out <- try(mvtnorm::pmvnorm(lower = rep(-Inf, d), upper = upper,
                              sigma = R,
                              algorithm = mvtnorm::Miwa())[1], silent = TRUE)
  # attr(out, "algo") = mvtnorm::Miwa()
  if (inherits(out, "try-error") || !is.finite(out)) {
    out <- mvtnorm::pmvnorm(lower = rep(-Inf, d), upper = upper,
                            sigma = R,
                            algorithm = mvtnorm::GenzBretz(abseps = abseps, maxpts = maxpts))[1]
    
  }
  as.numeric(out)
  # return(out)
}

h_gaus_pois_fast_safe <- function(corners, R) {
  d=ncol(R)
  
  n <- dim(corners$grid2)[1]
  probs <- numeric(n)
  for (i in 1:n) {
    ms <- corners$nonzero[[i]]
    if (length(ms) == 0L) { probs[i] <- 0; next }
    pr <- vapply(ms, function(m) safe_pmvnorm(upper = corners$grid2[i, , m], R = R), numeric(1L))
    probs[i] <- sum(corners$sign[ms] * pr)
  }
  pmax(pmin(probs, 1), 0)
}

## ---- corner scaffolding  ----
build_corner_data <- function(Y, lambda, eps = 1e-12) {
  n <- nrow(Y); d <- ncol(Y)
  
  grid1        <- expand.grid(rep(list(1:2), d)) |>
    dplyr::arrange(across(everything()))
  sign_factors <- (-1)^rowSums(grid1)
  n_opt        <- nrow(grid1)
  
  Fy  <- sapply(1:d, function(j) ppois(Y[, j],      lambda[j]))
  Fym <- sapply(1:d, function(j) ppois(Y[, j] - 1L, lambda[j]))
  ub  <- qnorm(pmin(pmax(Fy,  eps), 1 - eps))
  lb  <- qnorm(pmin(pmax(Fym, eps), 1 - eps))
  
  grid2 <- array(0.0, dim = c(n, d, n_opt))
  for (j in 1:d) {
    idx_lo <- which(grid1[, j] == 1L)
    idx_hi <- which(grid1[, j] == 2L)
    grid2[, j, idx_lo] <- lb[, j]
    grid2[, j, idx_hi] <- ub[, j]
  }
  nonzero_list <- lapply(1:n, function(i) which(colSums(is.finite(grid2[i, , ])) == d))
  
  # Poisson masses: y_k and y_k-1 (0 if y_k-1<0)
  dpm_y  <- sapply(1:d, function(k) dpois(Y[, k], lambda[k]))
  dpm_y1 <- sapply(1:d, function(k) ifelse(Y[, k] >= 1L, dpois(Y[, k]-1L, lambda[k]), 0.0))
  
  list(grid1 = grid1, sign = sign_factors, grid2 = grid2,
       nonzero = nonzero_list, ub = ub, lb = lb,
       dpm_y = dpm_y, dpm_y1 = dpm_y1)
}

grad_loglik <- function(par, data, abseps = 1e-9, maxpts = 1e6, eps = 1e-12, tiny = 1e-14) {
  lambda = par[1:d]
  rhos = par[(d+1):length(par)]
  R = p2P(rhos)
  Y <- as.matrix(data)
  n <- nrow(Y); d <- ncol(Y)
  hvals <- apply(Y, 1, function(y) h_one(y, lambda, R, abseps, maxpts, eps))
  wts <- 1 / pmax(hvals, tiny)
  
  ## grad wrt lambda (length d)
  g_lambda <- numeric(d)
  for (k in 1:d) {
    gk <- mapply(function(idx) dh_dlambda_k_one(Y[idx, ], lambda, R, k, abseps, maxpts, eps),
                 idx = seq_len(n))
    g_lambda[k] <- sum(wts * gk)
  }
  
  ## grad wrt rhos (P2p order)
  pairs <- t(utils::combn(d, 2))
  g_rho <- numeric(nrow(pairs))
  for (p in 1:nrow(pairs)) {
    i <- pairs[p, 1]; j <- pairs[p, 2]
    gij <- mapply(function(idx) dh_drho_ij_one(Y[idx, ], lambda, R, i, j, abseps, maxpts, eps),
                  idx = seq_len(n))
    g_rho[p] <- sum(wts * gij)
  }
  g_lambda = matrix(g_lambda, ncol = d)
  g_rho = matrix(g_rho, ncol = length(g_rho))
  grad = cbind(g_lambda, g_rho)
  out = colSums(grad)
  
  
  return(out)
}

# Replace the theta (zeta)-variant NLL to use the safe integrator
nll_tau_theta <- function(par_psy, Y) {
  d <- ncol(Y)
  tau   <- par_psy[1:d]
  theta <- par_psy[(d+1):length(par_psy)]
  lambda <- exp(tau)
  R <- theta_to_R(theta, d)
  
  corners <- build_corner_data(Y, lambda, 1e-12)
  hvals <- try(h_gaus_pois_fast_safe(corners, R), silent = TRUE)
  
  if (inherits(hvals, "try-error") || any(!is.finite(hvals))) {
    hvals <- apply(Y, 1, function(y) h_one(y, lambda, R))
  }
  return(-sum(log(pmax(hvals, 1e-300))))
}

h_one <- function(y, lambda, R, abseps = 1e-9, maxpts = 1e6, eps = 1e-12) {
  d <- length(y)
  a_hi <- qpois_norm_vec(y,     lambda, eps)
  a_lo <- qpois_norm_vec(y - 1, lambda, eps)
  total <- 0
  algo <- GB_algo(abseps, maxpts)
  for (m in 0:(2^d - 1)) {
    s <- as.integer(intToBits(m))[1:d]
    a <- ifelse(s == 1L, a_lo, a_hi)
    H_ms <- mvtnorm::pmvnorm(lower = rep(-Inf, d), upper = a, corr = R, algorithm = algo, seed = 1)[1]
    total <- total + ((-1) ^ sum(s)) * H_ms
  }
  as.numeric(total)
}

dh_drho_ij_one <- function(y, lambda, R, i, j, abseps = 1e-9, maxpts = 1e6, eps = 1e-12) {
  d <- length(y)
  a_hi <- qpois_norm_vec(y,     lambda, eps)
  a_lo <- qpois_norm_vec(y - 1, lambda, eps)
  total <- 0
  for (m in 0:(2^d - 1)) {
    s <- as.integer(intToBits(m))[1:d]
    a <- ifelse(s == 1L, a_lo, a_hi)
    term <- .dH_corner_drho_ij(a, R, i, j, abseps, maxpts)
    total <- total + ((-1) ^ sum(s)) * term
  }
  as.numeric(total)
}

dh_dlambda_k_one <- function(y, lambda, R, k, abseps = 1e-9, maxpts = 1e6, eps = 1e-12) {
  d <- length(y)
  a_hi <- qpois_norm_vec(y,     lambda, eps)
  a_lo <- qpois_norm_vec(y - 1, lambda, eps)
  total <- 0
  for (m in 0:(2^d - 1)) {
    s <- as.integer(intToBits(m))[1:d]
    y_corner <- y - s
    a <- ifelse(s == 1L, a_lo, a_hi)
    term <- .dH_corner_dlambda_k(a, y_corner, lambda, R, k, abseps, maxpts)
    total <- total + ((-1) ^ sum(s)) * term
  }
  as.numeric(total)
}

qpois_norm_vec <- function(y, lambda, eps = 1e-12) {
  u <- ppois(y, lambda)
  u <- pmin(pmax(u, eps), 1 - eps)
  qnorm(u)
}

## ---- dH/dlambda_k cornerwise (2D CDF via pbvn) ----
precomp_lambda_blocks <- function(R) {
  d = nrow(R)
  out <- vector("list", d)
  for (k in 1:d) {
    S    <- setdiff(1:d, k)          # d-1 indices
    RS   <- R[S, S, drop = FALSE]
    RSk  <- R[S, k, drop = FALSE]    # (d-1)×1
    Sig  <- RS - RSk %*% t(RSk)      # (d-1)x(d-1)
    sds   <- sqrt(diag(Sig))
    Dinv  <- diag(1 / sds, nrow = d-1)
    Rcond <- Dinv %*% Sig %*% Dinv        # Corr(Z_S | Z_k)  (d-1)x(d-1)
    out[[k]] <- list(k=k, S=S, RSk=RSk, sds = sds, rhoS=Rcond)
  }
  out
}

dH_dlambda_k_fast <- function(k, corners, lambda, blk, Y) {
  n <- nrow(Y); S <- blk$S
  out <- numeric(n)
  pmf_hi <- corners$dpm_y[,  k]   # uses y_k
  pmf_lo <- corners$dpm_y1[, k]   # uses y_k - 1
  
  for (idx in 1:n) {
    ms <- corners$nonzero[[idx]]
    if (length(ms) == 0L) { out[idx] <- 0; next }
    vals <- vapply(ms, function(m) {
      aS <- corners$grid2[idx, S, m]
      ak <- corners$grid2[idx, k, m]
      s_k <- (corners$grid1[m, k] == 1L)  # 1 -> lower -> y_k-1
      pmf <- if (s_k) pmf_lo[idx] else pmf_hi[idx]
      if (!is.finite(ak) || pmf == 0) return(0.0)
      mu  <- as.numeric(blk$RSk * ak)     # length d-1
      z = (aS - mu)/blk$sds
      - pmf * mvtnorm::pmvnorm(lower = rep(-Inf, d-1), upper = z, mean = rep(0,d-1), corr = blk$rhoS, algorithm = Miwa())[1]
    }, numeric(1L))
    out[idx] <- sum(corners$sign[ms] * vals)
  }
  out
}

dH_drho_pair_fast <- function(corners, R, block) {
  i <- block$i; j <- block$j; S <- block$S; rho <- block$rho
  n <- dim(corners$grid2)[1]
  out <- numeric(n)
  for (idx in 1:n) {
    ms <- corners$nonzero[[idx]]
    if (length(ms) == 0L) { out[idx] <- 0; next }
    vals <- vapply(ms, function(m) {
      a  <- corners$grid2[idx, , m]
      mu <- as.numeric(block$RSij %*% block$S22i %*% c(a[i], a[j]))  # (d-2)
      dbvnorm_rho(a[i], a[j], rho) * mvtnorm::pmvnorm(lower = rep(-Inf, d-2), 
                                                      upper = c(a[S]), mean = mu, sigma = block$Sigma_S, algorithm = Miwa())[1]
    }, numeric(1L))
    out[idx] <- sum(corners$sign[ms] * vals)
  }
  out
}

precomp_rho_blocks <- function(R) {
  d <- nrow(R)
  # P2p order: (2,1), (3,1), ..., (d,1), (3,2), ..., (d,d-1)
  pairs <- pairs_P2p_order(d)
  out <- vector("list", nrow(pairs))
  
  for (p in seq_len(nrow(pairs))) {
    i <- pairs[p, 1]; j <- pairs[p, 2]
    S <- setdiff(seq_len(d), c(i, j))
    
    rho  <- R[i, j]
    # 2x2 inverse of [[1, rho], [rho, 1]] (faster than solve)
    S22i <- (1 / (1 - rho^2)) * matrix(c(1, -rho, -rho, 1), 2, 2)
    
    RSij <- R[S, c(i, j), drop = FALSE]       # |S| x 2
    RSS  <- R[S, S, drop = FALSE]             # |S| x |S|
    Sigma_S <- RSS - RSij %*% S22i %*% t(RSij) # |S| x |S|  (Sigma_{S|ij})
    
    
    out[[p]] <- list(
      i = i, j = j, S = S,
      rho = rho,
      S22i = S22i,          # = R_{ij,ij}^{-1}
      RSij = RSij,          # = R_{S,ij}
      RSS  = RSS,           # = R_{SS}
      Sigma_S = Sigma_S    # = Sigma_{S|ij}
      #chol_Sigma_S = if (!inherits(L_S, "try-error")) L_S else NULL
    )
  }
  list(pairs = pairs, blocks = out)
}


pairs_P2p_order <- function(d) {
  # copula::P2p takes the lower triangle in column order:
  rbind(c(1L,2L), c(1L,3L), c(2L,3L))
  expand.grid(i=1:d,j=1:d) %>% filter(i<j) %>% arrange(i) %>% as.matrix()
}

dbvnorm_rho <- function(z1, z2, rho) {
  denom <- 2 * pi * sqrt(1 - rho^2)
  expo  <- - (z1^2 - 2*rho*z1*z2 + z2^2) / (2*(1 - rho^2))
  exp(expo) / denom
}

safe_cor <- function(Y) {
  R <- stats::cor(Y)
  if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) < 1e-8)
    R <- cov2cor(stats::cov(Y) + 1e-6 * diag(ncol(Y)))
  R
}

GB_algo <- function(abseps = 1e-9, maxpts = 1e6) mvtnorm::GenzBretz(abseps = abseps, maxpts = maxpts)
ALG3D <- mvtnorm::Miwa()   # deterministic & fast in d=3

simulate_Y_d <- function(n, lambda0, R0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  rhos0 <- copula::P2p(R0)
  nc <- normalCopula(param = rhos0, dim = length(lambda0), dispstr = "un")
  pm <- lapply(lambda0, function(x) list(lambda = x))
  dist <- mvdc(copula = nc, margins = rep("pois", length(lambda0)), paramMargins = pm)
  rMvdc(n, mvdc = dist)
}



compute_start_corr_alld <- function(data){
  # compute lambda start
  lam0 = pmax(colMeans(data), 0.1)
  r0 = pmax(-0.8, pmin(0.8, P2p(cor(data))))
  z1 <- R_to_zeta(p2P(r0))
  start <- c(log(lam0), z1)
  return(list(start.psy = start, start.nat = c(lam0,r0)))
}

rGausCopPoisAlld <- function(n, lambda, rho) {
  d = length(lambda)
  nc <- normalCopula(rho, dim = d, dispstr = "un")
  paramMargins.bp <- lapply(lambda, function(x) list(lambda = x))
  m  <- mvdc(copula = nc, margins = rep("pois", d), paramMargins = paramMargins.bp)
  rMvdc(n, m)
}

rho_hat_perr <- function(data, lam0){
  x = data[,1]; y = data[,2]
  L1 = lam0[1]; L2 = lam0[2]
  # compute rho start
  B1 = tie_prob_pois(L1); B2 = tie_prob_pois(L2)
  tau_b <- calculate_phz_2022_fast(x,y)
  A.hat = tau_b - B1 - B2 +1
  paramMargins <- list(list(lambda = L1),
                       list(lambda = L2))
  
  ymax = max(qpois(0.99999999999, L1),qpois(0.99999999999, L2))
  
  grid2 = expand_grid(x = 0:ymax, y = 0:ymax) %>% as.data.frame()
  
  
  f2 = function(r1) f(grid2, r1, paramMargins, A.hat)
  res = optimize(f2, lower = -0.95, upper = 0.95)

  res$minimum
}

compute_start_true_tau_perr_alld <- function(data){
  d = ncol(data)
  grid.ind = expand.grid(1:d,1:d) %>% filter(Var1 < Var2) %>% arrange(Var1)
  ng = nrow(grid.ind)
  
  starts.lam = pmax(colMeans(data), 0.1)
  starts.tau = log(starts.lam)
  
  starts.rho = matrix(NA, d, d)
  diag(starts.rho) = 1
  
  for (k in 1:ng){
    i = grid.ind$Var1[k]; j = grid.ind$Var2[k]
    
    Y1 = data[,c(i,j)]
    r0 = rho_hat_perr(Y1, starts.lam[c(i,j)])
    # starts1 = compute_start_true_tau_perr(Y1)
    starts.rho[i,j] = starts.rho[j,i] = r0
  }
  
  if (isposdef(starts.rho) == T){
    R = starts.rho
  } else {
    R = as.matrix(Matrix::nearPD(x = starts.rho, corr = T)$mat)
  }
  
  start.nat = c(starts.lam, P2p(R))
  start = c(starts.tau, R_to_zeta(R))
  
  return(list(start.psy = start, start.nat = start.nat))
}

calculate_phz_2022_robust <- function(A, B, na.rm = TRUE) {
  stopifnot(length(A) == length(B))
  
  if (na.rm) {
    keep <- is.finite(A) & is.finite(B)
    A <- A[keep]; B <- B[keep]
  }
  N <- length(A)
  if (N == 0L) return(NA_real_)
  
  # index sets
  ix01 <- which((A == 0) & (B > 0))
  ix11 <- which((A > 0) & (B > 0))
  ix10 <- which((B == 0) & (A > 0))
  
  # helpers to safely sum logicals with NA
  .sum_true <- function(x) sum(x, na.rm = TRUE)
  
  # counts for A: compare each A[ix10] to all A[ix11]
  x1_greater <- 0L; x1_tied <- 0L; x1_lower <- 0L
  if (length(ix10) > 0L && length(ix11) > 0L) {
    A11 <- A[ix11]
    for (a0 in A[ix10]) {
      x1_greater <- x1_greater + .sum_true(a0 >  A11)
      x1_tied    <- x1_tied    + .sum_true(a0 == A11)
      x1_lower   <- x1_lower   + .sum_true(a0 <  A11)
    }
  }
  
  # counts for B: compare each B[ix01] to all B[ix11]
  x2_greater <- 0L; x2_tied <- 0L; x2_lower <- 0L
  if (length(ix01) > 0L && length(ix11) > 0L) {
    B11 <- B[ix11]
    for (b0 in B[ix01]) {
      x2_greater <- x2_greater + .sum_true(b0 >  B11)
      x2_tied    <- x2_tied    + .sum_true(b0 == B11)
      x2_lower   <- x2_lower   + .sum_true(b0 <  B11)
    }
  }
  
  # convert to relative frequencies; define 0 when denominator is 0 to avoid NA*0 later
  den1 <- length(ix10) * length(ix11)
  den2 <- length(ix01) * length(ix11)
  .safe_div0 <- function(x, d) if (d > 0) x / d else 0
  
  p1       <- .safe_div0(x1_greater, den1)
  p1_ties  <- .safe_div0(x1_tied,    den1)
  p1_lower <- .safe_div0(x1_lower,   den1)
  
  p2       <- .safe_div0(x2_greater, den2)
  p2_ties  <- .safe_div0(x2_tied,    den2)
  p2_lower <- .safe_div0(x2_lower,   den2)
  
  # optional sanity checks (only when denominators > 0)
  if (den1 > 0) stopifnot(all.equal(p1 + p1_lower + p1_ties, 1, tol = 1e-12))
  if (den2 > 0) stopifnot(all.equal(p2 + p2_lower + p2_ties, 1, tol = 1e-12))
  
  # 2x2 for zeros vs positives (drop NAs already)
  x1pos <- A > 0
  x2pos <- B > 0
  p00 <- mean(!x1pos & !x2pos)
  p01 <- mean(!x1pos &  x2pos)
  p10 <- mean( x1pos & !x2pos)
  p11 <- mean( x1pos &  x2pos)
  
  # Kendall's tau on strictly positive pairs; use 0 if undefined to avoid NA * p11^2
  tau11 <- NA_real_
  if (p11 * N > 2) {
    adj <- x1pos & x2pos
    if (sum(adj) >= 2) {
      tau11 <- suppressWarnings(cor(A[adj], B[adj], method = "kendall"))
    }
  }
  tau11_safe <- ifelse(is.na(tau11), 0, tau11)
  
  tstar <- (p11^2) * tau11_safe +
    2 * (p00 * p11 - p01 * p10) +
    2 * p11 * ( p10 * (1 - 2 * p1 - p1_ties) + 
                  p01 * (1 - 2 * p2 - p2_ties) )
  return(tstar)
}

calculate_phz_2022_fast <- function(A, B, na.rm = TRUE) {
  stopifnot(length(A) == length(B))
  if (na.rm) {
    keep <- is.finite(A) & is.finite(B)
    A <- A[keep]; B <- B[keep]
  }
  N <- length(A); if (N == 0L) return(NA_real_)
  
  ix01 <- which((A == 0) & (B > 0))
  ix11 <- which((A > 0) & (B > 0))
  ix10 <- which((B == 0) & (A > 0))
  
  build_ref <- function(x) {
    if (length(x) == 0L) return(list(u = numeric(), freq = integer(), csum = integer(), tot = 0L))
    ux   <- sort(unique(x))
    idx  <- match(x, ux)
    freq <- tabulate(idx, nbins = length(ux))
    csum <- cumsum(freq)
    list(u = ux, freq = freq, csum = csum, tot = length(x))
  }
  # Return (greater, equal, lower) **for q vs ref**, i.e.
  # greater := #(q > ref) = #(ref < q) = lt
  # equal   := #(q == ref)
  # lower   := #(q < ref) = #(ref > q) = gt
  sum_rel <- function(q, ref) {
    if (length(q) == 0L || ref$tot == 0L) return(c(0L, 0L, 0L))
    idx  <- findInterval(q, ref$u)                 # #unique <= q
    mi   <- match(q, ref$u, nomatch = 0L)          # position for equality
    eq   <- ifelse(mi > 0L, ref$freq[mi], 0L)      # #(ref == q)
    leq  <- ifelse(idx > 0L, ref$csum[idx], 0L)    # #(ref <= q)
    lt   <- leq - eq                                # #(ref < q)
    gt   <- ref$tot - leq                           # #(ref > q)
    c(sum(lt), sum(eq), sum(gt))                    # (q>ref, q==ref, q<ref)
  }
  
  x1_greater <- x1_tied <- x1_lower <- 0L
  if (length(ix10) > 0L && length(ix11) > 0L) {
    refA <- build_ref(A[ix11])
    aggA <- sum_rel(A[ix10], refA)
    x1_greater <- aggA[1]; x1_tied <- aggA[2]; x1_lower <- aggA[3]
  }
  x2_greater <- x2_tied <- x2_lower <- 0L
  if (length(ix01) > 0L && length(ix11) > 0L) {
    refB <- build_ref(B[ix11])
    aggB <- sum_rel(B[ix01], refB)
    x2_greater <- aggB[1]; x2_tied <- aggB[2]; x2_lower <- aggB[3]
  }
  
  den1 <- length(ix10) * length(ix11)
  den2 <- length(ix01) * length(ix11)
  .safe_div0 <- function(x, d) if (d > 0) x / d else 0
  
  p1       <- .safe_div0(x1_greater, den1)
  p1_ties  <- .safe_div0(x1_tied,    den1)
  p1_lower <- .safe_div0(x1_lower,   den1)
  
  p2       <- .safe_div0(x2_greater, den2)
  p2_ties  <- .safe_div0(x2_tied,    den2)
  p2_lower <- .safe_div0(x2_lower,   den2)
  
  if (den1 > 0) stopifnot(all.equal(p1 + p1_lower + p1_ties, 1, tol = 1e-12))
  if (den2 > 0) stopifnot(all.equal(p2 + p2_lower + p2_ties, 1, tol = 1e-12))
  
  x1pos <- A > 0
  x2pos <- B > 0
  p00 <- mean(!x1pos & !x2pos)
  p01 <- mean(!x1pos &  x2pos)
  p10 <- mean( x1pos & !x2pos)
  p11 <- mean( x1pos &  x2pos)
  
  tau11 <- NA_real_
  if (p11 * N > 2) {
    adj <- x1pos & x2pos
    if (sum(adj) >= 2) {
      tau11 <- suppressWarnings(cor(A[adj], B[adj], method = "kendall"))
    }
  }
  tau11_safe <- ifelse(is.na(tau11), 0, tau11)
  
  tstar <- (p11^2) * tau11_safe +
    2 * (p00 * p11 - p01 * p10) +
    2 * p11 * ( p10 * (1 - 2 * p1 - p1_ties) +
                  p01 * (1 - 2 * p2 - p2_ties) )
  return(tstar)
}

sim_quick <- function(d, n = 10, seed = 1){
  lam = generate_lambdas(d, seed)
  P1 = unstructured(d, seed)
  rh = P2p(P1)
  Y   <- rGausCopPoisAlld(n, lam, rh)
  return(list(Y=Y, lambda = lam, P = P1, rho0 = rh))
}


