### REPARAM FNS

# Quick fn that only works for d=2 as we have w = acos(rho); zeta = log(w) - log(pi - w)

zeta_from_rho <- function(rho1) {
  w = acos(rho1)
  return(log(w) - log(pi-w))
}

# For all d
R_to_zeta <- function(R){
  L=chol(R) %>% t()
  Omega = compute_omega_mat(L)
  omega = P2p(Omega)
  zeta = omega_to_zeta(omega)
  return(zeta)
}

compute_omega_mat <- function(L){
  m = nrow(L)
  omega.mat = matrix(0,m,m)
  if (m == 2){
    omega.mat[2,1] = acos(L[2,1])
  } else {
    omega.mat[2:m,1] = acos(L[2:m,1])
    omega.mat[3:m,2] = acos(L[3:m,2]/sin(omega.mat[3:m,1]))
  }
  
  
  if( m > 4){
    for (i in 4:(m-1)){
      omega.mat[i:m,(i-1)] = acos(L[i:m,(i-1)]/apply(sin(omega.mat[i:m,1:(i-2)]), 1, prod))
      
    }
    omega.mat[m,(m-1)] = acos(L[m,(m-1)]/prod(sin(omega.mat[m,1:(m-2)])))
  }
  
  if( m == 4){
    omega.mat[m,(m-1)] = acos(L[m,(m-1)]/prod(sin(omega.mat[m,1:(m-2)])))
  }
  
  return(omega.mat)
}

omega_to_zeta <- function(omega){
  log(omega) - log(pi-omega)
}

zeta_to_R <- function(zeta, d){
  #indices <- expand.grid(i = 4:d, j = 3:d) %>% filter(i!=j, i > j) %>% arrange(i) # double check in d>3
  omega <- zeta_to_omega(zeta)
  omega.mat = p2P(omega); diag(omega.mat) = 0; omega.mat[upper.tri(omega.mat)] = 0
  L = build_L_from_omega(omega.mat)
  R = L %*% t(L)
  return(R)
}

zeta_to_omega <- function(zeta){
  e.zeta = exp(zeta)
  (e.zeta/(1+e.zeta))*pi
}

build_L_from_omega <- function(Omega) {
  d <- nrow(Omega)
  L <- matrix(0, d, d)
  
  for (i in 1:d) {
    # Use only the first (i-1) angles for row i
    
    if (i == 1) {
      w_row <- Omega[i,1] # 1 angle
      sines <- sin(w_row)
      cosines <- cos(w_row)
      L[i,1] = cosines[1]
    } else {
      w_row <- Omega[i, 1:(i-1)]
      sines <- sin(w_row)
      cosines <- cos(w_row)
      
      # j=1:
      L[i,1] = cosines[1]
      
      # For 2 <= j <= i-1:
      if (i > 2) {
        cp_sin <- cumprod(sines)
        for (j in 2:(i-1)) {
          # product of sin(w_{i,k}) for k=1 to j-1 = cp_sin[j-1]
          L[i,j] = cosines[j] * cp_sin[j-1]
        }
        # j = i:
        L[i,i] = cp_sin[i-1]
      } else {
        # i=2 case: we have only one angle w_{2,1}
        # a_{2,2} = sin(w_{2,1})
        L[i,i] = sines[1]
      }
    }
  }
  
  return(L)
}

# Reparametrized objective: theta(zeta) = (gamma1, gamma2, eta)
negloglik_theta <- function(theta, Y) {
  lam <- lambda_from_gamma(theta[1:2])
  rho <- rho_from_theta(theta[3])
  h   <- dGausCopPois_pbiv(Y, lam, rho)
  - sum(log(h))
}

# Analytic gradient for theta(zeta) via chain rule
grad_negloglik_theta <- function(theta, Y) {
  lam <- lambda_from_gamma(theta[1:2])
  rho <- rho_from_theta(theta[3])
  h   <- dGausCopPois_pbiv(Y, lam, rho)
  # base partials
  g1  <- dh_dlambda1(Y[,1], Y[,2], lam[1], lam[2], rho) / h
  g2  <- dh_dlambda2(Y[,1], Y[,2], lam[1], lam[2], rho) / h
  gr  <-      dh_drho(Y[,1], Y[,2], lam[1], lam[2], rho) / h
  # chain rule: dlam/dgamma = lam ; drho/deta = 1 - rho^2
  d_gamma1 <- - sum(g1) * lam[1]
  d_gamma2 <- - sum(g2) * lam[2]
  d_theta <- - sum(gr) * drho_dtheta(theta[3])
  c(d_gamma1, d_gamma2, d_theta)
}


# ------------------------------------------------------------
# Derivatives for analytic gradient
# ------------------------------------------------------------
# dPhi_2(a,b;rho)/drho = phi_2(a,b;rho)
phi2 <- function(a, b, rho) {
  denom <- 2*pi*sqrt(1 - rho^2)
  Q <- (a^2 - 2*rho*a*b + b^2) / (2*(1 - rho^2))
  exp(-Q) / denom
}


## GRADIENT
### gradient helpers
lambda_from_gamma <- function(gamma) exp(gamma)          # R -> (0,Inf)

rho_from_theta   <- function(theta) {
  w <- omega_from_theta(theta)
  cos(w)
}

# build Omega strictly lower-triangular from zeta
Omega_from_zeta <- function(zeta, d){
  omega <- zeta_to_omega(zeta)          # pi * logistic(zeta)
  Om <- copula::p2P(omega)              # fill symmetric matrix in copula order
  diag(Om) <- 0
  Om[upper.tri(Om)] <- 0                # keep strictly lower triangle as angles
  Om
}

omega_from_theta <- function(theta) pi * logistic(theta)

dlambda_dgamma <- function(gamma) exp(gamma)             # = lambda

logistic <- function(x) 1/(1 + exp(-x))

drho_dtheta <- function(theta) {
  s <- logistic(theta)
  w <- pi * s
  - sin(w) * (pi * s * (1 - s))
}

# dH/dlam for one margin (closed form):
# dH/dlam1 = - P(Y1 = y1; lam1) * Phi((b - rho a)/sqrt(1-rho^2))
# and symmetrically for lam2
.dH_dlambda <- function(y_fix, other_y, lambda_fix, lambda_other, rho, which = 1L, eps = 1e-12) {
  u_fix   <- clip01(ppois(y_fix,   lambda_fix), eps)
  u_other <- clip01(ppois(other_y, lambda_other), eps)
  # u_fix   <- ppois(y_fix,   lambda_fix)
  # u_other <- ppois(other_y, lambda_other)
  a <- qnorm(u_fix); b <- qnorm(u_other)
  pmf <- dpois(y_fix, lambda_fix)
  if (which == 1L)  t <- (b - rho * a) / sqrt(1 - rho^2)
  else               t <- (a - rho * b) / sqrt(1 - rho^2)
  - pmf * pnorm(t)
}

dH_dlambda1 <- function(y1, y2, lambda1, lambda2, rho, eps = 1e-12) .dH_dlambda(y1, y2, lambda1, lambda2, rho, which = 1L, eps = eps)

dH_dlambda2 <- function(y1, y2, lambda1, lambda2, rho, eps = 1e-12) .dH_dlambda(y2, y1, lambda2, lambda1, rho, which = 2L, eps = eps)

#Inclusion–exclusion for dh/dlamj (four corners of H)
dh_dlambda1 <- function(y1, y2, lambda1, lambda2, rho, eps = 1e-12) {
  dH_dlambda1(y1, y2, lambda1, lambda2, rho, eps) -
    dH_dlambda1(y1-1, y2, lambda1, lambda2, rho, eps) -
    dH_dlambda1(y1, y2-1, lambda1, lambda2, rho, eps) +
    dH_dlambda1(y1-1, y2-1, lambda1, lambda2, rho, eps)
}

dh_dlambda2 <- function(y1, y2, lambda1, lambda2, rho, eps = 1e-12) {
  dH_dlambda2(y1, y2, lambda1, lambda2, rho, eps) -
    dH_dlambda2(y1-1, y2, lambda1, lambda2, rho, eps) -
    dH_dlambda2(y1, y2-1, lambda1, lambda2, rho, eps) +
    dH_dlambda2(y1-1, y2-1, lambda1, lambda2, rho, eps)
}

# d/drho of h2 (rectangle mass)
dh_drho <- function(y1, y2, lambda1, lambda2, rho) {
  a1 <- qpois_norm(y1,     lambda1)
  a0 <- qpois_norm(y1 - 1, lambda1)
  b1 <- qpois_norm(y2,     lambda2)
  b0 <- qpois_norm(y2 - 1, lambda2)
  phi2(a1,b1,rho) - phi2(a0,b1,rho) - phi2(a1,b0,rho) + phi2(a0,b0,rho)
}

# Analytic gradient for direct parameterization (lam1, lam2, rho)
grad_negloglik_direct <- function(par, Y) {
  lam <- par[1:2]; rho <- par[3]
  h   <- dGausCopPois_pbiv(Y, lam, rho)
  # components for each observation
  g1  <- dh_dlambda1(Y[,1], Y[,2], lam[1], lam[2], rho) / h
  g2  <- dh_dlambda2(Y[,1], Y[,2], lam[1], lam[2], rho) / h
  gr  <- dh_drho(Y[,1], Y[,2], lam[1], lam[2], rho) / h
  - c(sum(g1), sum(g2), sum(gr))                               # gradient of negative log-lik
}

### helpers
# Helper: probit of CDFs at the two levels y and y-1
qpois_norm <- function(y, lambda, eps = 1e-12){
  qnorm(clip01(ppois(y, lambda), eps)) # NEEDED otherwise produces NaN
  # qnorm(ppois(y, lambda), eps) # produces NaN
} 

# Guard for correlation matrix; in d=2 PD <=> |rho|<1
R2 <- function(rho) { rho <- pmin(pmax(rho, -0.999999), 0.999999); matrix(c(1, rho, rho, 1), 2, 2) }

clip01 <- function(u, eps = 1e-12) pmin(pmax(u, eps), 1 - eps)

# ------------------------------------------------------------
# Simulation helpers
# ------------------------------------------------------------
# Simulate from a bivariate Gaussian copula with Poisson margins
rGausCopPois <- function(n, lambda, rho) {
  stopifnot(length(lambda) == 2)
  nc <- normalCopula(rho, dim = 2)
  m  <- mvdc(nc, margins = c("pois","pois"), paramMargins = list(list(lambda=lambda[1]), list(lambda=lambda[2])))
  rMvdc(n, m)
}

### fns

# Rectangle probability (pmf cell) via inclusion–exclusion
h2 <- function(y1, y2, lambda1, lambda2, rho, eps = 1e-12) {
  a1 <- qnorm(clip01(ppois(y1,     lambda1), eps)); b1 <- qnorm(clip01(ppois(y2,     lambda2), eps))
  a0 <- qnorm(clip01(ppois(y1 - 1, lambda1), eps)); b0 <- qnorm(clip01(ppois(y2 - 1, lambda2), eps))
  
  # a1 <- qnorm(ppois(y1,     lambda1), eps); b1 <- qnorm(ppois(y2,     lambda2), eps)
  # a0 <- qnorm(ppois(y1 - 1, lambda1), eps); b0 <- qnorm(ppois(y2 - 1, lambda2), eps)
  
  # Four corners
  F11 <- pbivnorm::pbivnorm(a1, b1, rho)
  F01 <- pbivnorm::pbivnorm(a0, b1, rho)
  F10 <- pbivnorm::pbivnorm(a1, b0, rho)
  F00 <- pbivnorm::pbivnorm(a0, b0, rho)
  s <- F11 - F01 - F10 + F00
  s[!is.finite(s)] <- .Machine$double.xmin
  pmax(s, .Machine$double.xmin)
}


# Vectorized pmf over a matrix Y (n x 2)
dGausCopPois_pbiv <- function(Y, lambda, rho, eps = 1e-12) {
  Y <- as.matrix(Y); stopifnot(ncol(Y) == 2)
  h2(Y[,1], Y[,2], lambda[1], lambda[2], rho, eps)
}

# Robust fallback using mvtnorm::pmvnorm on each rectangle (slower)
dGausCopPois_mvt <- function(Y, lambda, rho) {
  Y <- as.matrix(Y); n <- nrow(Y)
  a1 <- qnorm(ppois(Y[,1],     lambda[1])); a0 <- qnorm(ppois(Y[,1]-1, lambda[1]))
  b1 <- qnorm(ppois(Y[,2],     lambda[2])); b0 <- qnorm(ppois(Y[,2]-1, lambda[2]))
  Sig <- R2(rho)
  out <- numeric(n)
  for(i in seq_len(n)) {
    val <- mvtnorm::pmvnorm(lower = c(a0[i], b0[i]), upper = c(a1[i], b1[i]), mean = c(0,0), corr = Sig, seed = 1)[1]
    if(!is.finite(val) || val <= 0) val <- .Machine$double.xmin
    out[i] <- val
  }
  out
}


# Master pmf (choose fast pbivnorm, with optional fallback if needed)
dGausCopPois <- function(Y, lambda, rho, use_mvt_fallback = FALSE) {
  if (!use_mvt_fallback) return(dGausCopPois_pbiv(Y, lambda, rho))
  dGausCopPois_mvt(Y, lambda, rho)
}

### nll

# ------------------------------------------------------------
# Log-likelihoods and gradients
# ------------------------------------------------------------
negloglik_direct <- function(par, Y, use_mvt_fallback = FALSE) {
  lam <- par[1:2]
  rho <- par[3]
  if (any(!is.finite(lam)) || any(lam <= 0) || !is.finite(rho) || abs(rho) >= 1) return(Inf)
  h <- dGausCopPois(Y, lam, rho, use_mvt_fallback)
  - sum(log(h))
}

### Starting values and TAU
### Empirical tau
# Following two functions are adapted from from Perrone et.al. https://arxiv.org/abs/2208.03155 
calculate_pimtau_2015_debug <- function(A, B) {
  stopifnot(length(A) == length(B))
  # drop pairs with NA/Inf
  ok <- is.finite(A) & is.finite(B)
  A <- A[ok]; B <- B[ok]
  N <- length(A)
  if (N == 0L) return(NA_real_)
  
  # Indicators for zero vs >0
  x1new <- as.integer(A > 0)
  x2new <- as.integer(B > 0)
  
  # Robust 2x2 table with fixed levels
  tab <- table(factor(x1new, levels = 0:1),
               factor(x2new, levels = 0:1))
  p00 <- tab[1, 1] / N
  p01 <- tab[1, 2] / N
  p10 <- tab[2, 1] / N
  p11 <- tab[2, 2] / N
  
  # tau11 on the positive-positive sub-sample (only if at least 2 obs)
  tau11 <- NA_real_
  if (tab[2, 2] >= 2L) {
    idx11 <- (x1new == 1L) & (x2new == 1L)
    adjA <- A[idx11]; adjB <- B[idx11]
    # If adjA or adjB has zero variance, cor(..., "kendall") returns NA — that’s OK;
    # the p11^2 weight will zero out this contribution if p11==0.
    tau11 <- suppressWarnings(stats::cor(adjA, adjB, method = "kendall"))
  }
  
  # Helper to compute pairwise proportion in a memory-safe way
  pair_prop <- function(x, y, op = c("gt", "eq", "lt"), max_outer = 5e7) {
    op <- match.arg(op)
    nx <- length(x); ny <- length(y)
    if (nx == 0L || ny == 0L) return(NA_real_)
    # Use outer only if the matrix won't be too large
    if ((as.double(nx) * as.double(ny)) <= max_outer) {
      M <- switch(op,
                  gt = outer(x, y, ">"),
                  eq = outer(x, y, "=="),
                  lt = outer(x, y, "<"))
      return(mean(M))
    } else {
      cmp_fun <- switch(op,
                        gt = function(a) sum(a >  y),
                        eq = function(a) sum(a == y),
                        lt = function(a) sum(a <  y))
      s <- vapply(x, cmp_fun, integer(1L))
      return(sum(s) / (nx * ny))
    }
  }
  
  # Indices for the cells needed in p1 and p2
  ix10 <- which((x1new == 1L) & (x2new == 0L))  # A>0, B=0
  ix01 <- which((x1new == 0L) & (x2new == 1L))  # A=0, B>0
  ix11 <- which((x1new == 1L) & (x2new == 1L))  # A>0, B>0
  
  # p1 = P( A_{10} > A_{11} ), p2 = P( B_{01} > B_{11} )
  p1 <- if (length(ix10) > 0L && length(ix11) > 0L) {
    pair_prop(A[ix10], A[ix11], "gt")
  } else NA_real_
  p2 <- if (length(ix01) > 0L && length(ix11) > 0L) {
    pair_prop(B[ix01], B[ix11], "gt")
  } else NA_real_
  
  # Build the terms carefully so missing p1/p2 never propagate:
  term_tau11 <- (p11 ^ 2) * ifelse(is.na(tau11), 0, tau11)
  term_mass  <- 2 * (p00 * p11 - p01 * p10)
  term_p1    <- if (p11 > 0 && p10 > 0 && !is.na(p1)) 2 * p11 * p10 * (1 - 2 * p1) else 0
  term_p2    <- if (p11 > 0 && p01 > 0 && !is.na(p2)) 2 * p11 * p01 * (1 - 2 * p2) else 0
  
  tstar <- term_tau11 + term_mass + term_p1 + term_p2
  return(tstar)
}

calculate_phz_2022 <- function(A, B) {
  N <- length(A)
  # calculate p1 and p2 for Pimentel 2015, equation (8)
  ix01 <- which((A == 0) & (B > 0))
  ix11 <- which((A > 0) & (B > 0))
  ix10 <- which((B == 0) & (A > 0))
  x1_greater <- 0 # initial value for the frequency of X_{10}-X_{11}>0
  x2_greater <- 0 # initial value for the frequency of Y_{10}-Y_{11}>0
  x1_tied <- 0 # initial value for the frequency of X_{10}-X_{11}=0
  x2_tied <- 0 # initial value for the frequency of Y_{01}-Y_{11}=0
  x1_lower <- 0 # initial value for the frequency of X_{10}-X_{11}<0
  x2_lower <- 0 # initial value for the frequency of Y_{10}-Y_{11}<0
  # calculate x1_greater, x1_ties, and x1_lower
  for (i in 1:length(A[ix10])) {
    if (count(A[ix10[i]] > A[ix11])[1, 2] < length(ix11)) {
      x1_greater <- x1_greater + count(A[ix10[i]] > A[ix11])[2, 2]
    } else if (count(A[ix10[i]] > A[ix11])[1, 1] == TRUE){
      x1_greater <- x1_greater + length(ix11)
    }
    if (count(A[ix10[i]] == A[ix11])[1, 2] < length(ix11)) {
      x1_tied <- x1_tied + count(A[ix10[i]] == A[ix11])[2, 2]
    } else if (count(A[ix10[i]] == A[ix11])[1, 1] == TRUE){
      x1_tied <- x1_tied + length(ix11)
    } 
    if (count(A[ix10[i]] < A[ix11])[1, 2] < length(ix11)) {
      x1_lower <- x1_lower + count(A[ix10[i]] < A[ix11])[2, 2]
    } else if (count(A[ix10[i]] < A[ix11])[1, 1] == TRUE){
      x1_lower <- x1_lower + length(ix11)
    } 
  }
  # calculate x2_greater, x2_ties, and x2_lower
  for (i in 1:length(B[ix01])) {
    if (count(B[ix01[i]] > B[ix11])[1, 2] < length(ix11)) {
      x2_greater <- x2_greater + count(B[ix01[i]] > B[ix11])[2, 2]
    } else if(count(B[ix01[i]] > B[ix11])[1, 1]==TRUE){
      x2_greater <- x2_greater + length(ix11)
    }
    if (count(B[ix01[i]] == B[ix11])[1, 2] < length(ix11)) {
      x2_tied <- x2_tied + count(B[ix01[i]] == B[ix11])[2, 2]
    } else if(count(B[ix01[i]] == B[ix11])[1, 1]==TRUE){
      x2_tied <- x2_tied + length(ix11)
    }  
    if (count(B[ix01[i]] < B[ix11])[1, 2] < length(ix11)) {
      x2_lower <- x2_lower + count(B[ix01[i]] < B[ix11])[2, 2]
    } else if(count(B[ix01[i]] < B[ix11])[1, 1]==TRUE){
      x2_lower <- x2_lower + length(ix11)
    } 
  }
  # relative frequency of X_{10}-X_{11}>0
  p1 <- x1_greater / (length(ix10) * length(ix11))
  # relative frequency of Y_{10}-Y_{11}>0
  p2 <- x2_greater / (length(ix01) * length(ix11))
  # relative frequency of X_{10}-X_{11}==0
  p1_ties <- x1_tied / (length(ix10) * length(ix11))
  # relative frequency of Y_{10}-Y_{11}==0
  p2_ties <- x2_tied / (length(ix01) * length(ix11))
  # relative frequency of X_{10}-X_{11}<0
  p1_lower <- x1_lower / (length(ix10) * length(ix11))
  # relative frequency of Y_{10}-Y_{11}<0
  p2_lower <- x2_lower / (length(ix01) * length(ix11))
  # check that they sum up to 1
  stopifnot(all.equal((p1 + p1_lower + p1_ties),1) ,
            all.equal((p2 + p2_lower + p2_ties),1)) 
  # calculate the probabilities p00, p11, p01, p10
  x1new <- ifelse(A > 0, 1, 0)
  x2new <- ifelse(B > 0, 1, 0)
  if(dim(table(x1new, x2new))[1] == 2 & dim(table(x1new, x2new))[2] == 2){
    cont_table <- table(x1new, x2new)
  } else if (dim(table(x1new, x2new))[1] == 1){
    cont_table <- rbind(table(x1new, x2new), c(0, 0))
  } else if (dim(table(x1new, x2new))[2] == 1) {
    cont_table <- cbind(table(x1new, x2new), c(0, 0))
  } else{
    print(table(x1new,x2new))
  }
  p00 <- cont_table[1, 1] / N
  p01 <- cont_table[1, 2] / N
  p10 <- cont_table[2, 1] / N
  p11 <- cont_table[2, 2] / N
  if (p11 <= (2 / N)) { 
    tau11 <- NA
  } else {
    adjustedA <- A[A > 0 & B > 0]
    adjustedB <- B[B > 0 & A > 0]
    tau11 <- cor(adjustedA, adjustedB, method = "kendall")
  }
  tstar <- p11^2 * tau11 + 2 * (p00 * p11 - p01 * p10) + 2 * p11 * (p10 * (1 - 2 * p1 - p1_ties) + p01 * (1 - 2 * p2 - p2_ties)) # Zhan's derivations
  return(tstar)
}

# approx tau

tie_prob_pois <- function(lambda) exp(-2*lambda) * besselI(2*lambda, 0)

kappa_geom <- function(L1,L2){
  B1 <- tie_prob_pois(L1); B2 <- tie_prob_pois(L2)
  sqrt((1 - B1)*(1 - B2))
}
tau_cheapest <- function(rho,L1,L2){
  (2/pi) * kappa_geom(L1,L2) * asin(rho)
}

compute_rho_from_tau_pois <- function(tau, B1, B2){
  x = (pi/2) * (1/sqrt((1-B1)*(1-B2))) * tau
  sin(x)
}

rho_hat_cheapest <- function(x,y){
  stopifnot(length(x)==length(y))
  L1 <- mean(x); L2 <- mean(y)
  B1 = tie_prob_pois(L1); B2 = tie_prob_pois(L2)
  tau_b <- calculate_pimtau_2015_debug(x,y)
  rho  <- compute_rho_from_tau_pois(tau_b, B1, B2)
  pmax(-0.8, pmin(0.8, rho))
}

compute_start_tau <- function(data){
  x = data[,1]; y = data[,2]
  # compute lambda start
  lam0 = pmax(colMeans(data), 0.1)
  L1 = lam0[1]; L2 = lam0[2]
  
  # compute rho start
  B1 = tie_prob_pois(L1); B2 = tie_prob_pois(L2)
  tau_b <- calculate_pimtau_2015_debug(x,y)
  rho  <- compute_rho_from_tau_pois(tau_b, B1, B2)
  r0 = pmax(-0.8, pmin(0.8, rho))
  z1 <- zeta_from_rho(r0)
  start <- c(log(lam0), z1)
  return(start)
}

compute_start_corr <- function(data){
  x = data[,1]; y = data[,2]
  # compute lambda start
  lam0 = pmax(colMeans(data), 0.1)
  r0 = pmax(-0.8, pmin(0.8, cor(data)[1,2]))
  z1 <- zeta_from_rho(r0)
  start <- c(log(lam0), z1)
  return(list(start.psy = start, start.nat = c(lam0,r0)))
}

f = function(grid2, r, paramMargins, A.hat){
  nc <- normalCopula(r, dim = 2)
  m <- mvdc(copula = nc,
                  margins = c("pois", "pois"),
                  paramMargins = paramMargins)
  H = pMvdc(as.matrix(grid2-1),mvdc = m)
  # h = h_YY_pbvt(grid2, rho.1 = r, l1 = paramMargins[[1]]$lambda, l2 = paramMargins[[2]]$lambda)
  h = h2(y1 = grid2$x, y2 = grid2$y, lambda1 = paramMargins[[1]]$lambda, lambda2 = paramMargins[[2]]$lambda, rho = r)
  grid2 <- grid2 %>% mutate(H=H, h=h, Hh = 4*H*h, h2 = h^2, A=Hh-h2)
  A = sum(grid2$A)
  sum((A-A.hat)^2)
}

# function to compute MuZIC density for d=2
h_YY_pbvt <- function(counts, rho.1, l1, l2, df=100){
  
  ### INPUTS 
  # counts  [ n x 2 ] :                       matrix of counts with as columns Y1, Y2
  # l1, l2  scalars :                         λ1, λ2
  # rho.1   scalar :                          rho 
  
  ### OUTPUTS
  # hNN     [ n x 1 ] :                        values of pdf evaluated at (counts)
  
  n = nrow(counts)
  # compute lower and upper bounds
  lb = cbind(qnorm(ppois(q = counts[,1] - 1,lambda=l1)),qnorm(ppois(q = counts[,2] - 1,lambda=l2)))
  ub = cbind(qnorm(ppois(q = counts[,1],lambda=l1)),qnorm(ppois(q = counts[,2],lambda=l2)))
  
  # for numerical purposes
  lb[lb==-Inf]=-10
  lb[lb==Inf]=10
  ub[ub==-Inf]=-10
  ub[ub==Inf]=10
  
  # compute densities using rectangle probabilities and Student t approximation to Normal with df = 100
  hNN=pbvt(ub[,1], ub[,2], c(rho.1,100)) -
    pbvt(ub[,1], lb[,2], c(rho.1,100)) - pbvt(lb[,1], ub[,2], c(rho.1,100)) +
    pbvt(lb[,1], lb[,2], c(rho.1,100))
  return(hNN)
}

# compute starting values for zeta using reparametrization and argmin for rho for the expression of the true Kendall's tau (no approximation)
compute_start_true_tau <- function(data){
  x = data[,1]; y = data[,2]
  # compute lambda start
  lam0 = pmax(colMeans(data), 0.1)
  L1 = lam0[1]; L2 = lam0[2]
  # compute rho start
  B1 = tie_prob_pois(L1); B2 = tie_prob_pois(L2)
  tau_b <- calculate_pimtau_2015_debug(x,y)
  A.hat = tau_b - B1 - B2 +1
  paramMargins <- list(list(lambda = L1),
                       list(lambda = L2))
  
  ymax = max(qpois(0.99999999999, L1),qpois(0.99999999999, L2))
  
  grid2 = expand_grid(x = 0:ymax, y = 0:ymax) %>% as.data.frame()
  
  f2 = function(r1) f(grid2, r1, paramMargins, A.hat)
  res = optimize(f2, lower = -0.95, upper = 0.95)
  #r0 = pmax(-0.8, pmin(0.8, res$minimum))
  r0 = res$minimum
  z1 <- zeta_from_rho(r0)
  start <- c(log(lam0), z1)
  return(list(start.psy = start, start.nat = c(lam0,r0)))
}

# compute starting values for zeta using reparametrization and argmin for rho for the expression of the true Kendall's tau (no approximation)
compute_start_true_tau_pim <- function(data){
  x = data[,1]; y = data[,2]
  # compute lambda start
  lam0 = pmax(colMeans(data), 0.1)
  L1 = lam0[1]; L2 = lam0[2]
  # compute rho start
  B1 = tie_prob_pois(L1); B2 = tie_prob_pois(L2)
  tau_b <- calculate_pimtau_2015_debug(x,y)
  A.hat = tau_b - B1 - B2 +1
  paramMargins <- list(list(lambda = L1),
                       list(lambda = L2))
  
  ymax = max(qpois(0.99999999999, L1),qpois(0.99999999999, L2))
  
  grid2 = expand_grid(x = 0:ymax, y = 0:ymax) %>% as.data.frame()
  
  f <- function(r) f3(r = r, paramMargins = paramMargins, grid2 = grid2) - A.hat
  r0 = uniroot(f, lower = -.999, upper = 0.999)$root
  z1 <- zeta_from_rho(r0)
  start <- c(log(lam0), z1)
  return(list(start.psy = start, start.nat = c(lam0,r0)))
}

f3 = function(r,paramMargins,grid2){
  nc <- normalCopula(r, dim = 2)
  m <- mvdc(copula = nc,
            margins = c("pois", "pois"),
            paramMargins = paramMargins)
  H = pMvdc(as.matrix(grid2-1),mvdc = m)
  h = h2(y1 = grid2$x, y2 = grid2$y, lambda1 = paramMargins[[1]]$lambda, lambda2 = paramMargins[[2]]$lambda, rho = r)
  grid2 <- grid2 %>% mutate(H=H, h=h, Hh = 4*H*h, h2 = h^2, A=Hh-h2)
  A = sum(grid2$A)
  return(A)
}

compute_start_true_tau_perr <- function(data){
  x = data[,1]; y = data[,2]
  # compute lambda start
  lam0 = pmax(colMeans(data), 0.1)
  L1 = lam0[1]; L2 = lam0[2]
  # compute rho start
  B1 = tie_prob_pois(L1); B2 = tie_prob_pois(L2)
  tau_b <- calculate_phz_2022(x,y) 
  A.hat = tau_b - B1 - B2 +1
  paramMargins <- list(list(lambda = L1),
                       list(lambda = L2))
  
  ymax = max(qpois(0.99999999999, L1),qpois(0.99999999999, L2))
  
  grid2 = expand_grid(x = 0:ymax, y = 0:ymax) %>% as.data.frame()
  
  
  f2 = function(r1) f(grid2, r1, paramMargins, A.hat)
  res = optimize(f2, lower = -0.95, upper = 0.95)
  #r0 = pmax(-0.8, pmin(0.8, res$minimum))
  r0 = res$minimum
  
  z1 <- zeta_from_rho(r0)
  start <- c(log(lam0), z1)
  return(list(start.psy = start, start.nat = c(lam0,r0)))
}
