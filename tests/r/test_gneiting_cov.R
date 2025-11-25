# CorGneiting class: test of the covariance functions

knitr::opts_chunk$set(echo = TRUE)
library(gstlearn)
rm(list = ls())
# source(file = "./function_gneiting_model.R")

# auxiliary functions
kernel_power_exponential <- function(t, alpha = 2.0) {
  stopifnot((alpha > 0)&(alpha <= 2.0))
  exp(-abs(t)^alpha)
}
kernel_matern <- function(h, nu){
  val = rep(1.0, length(h))
  sel = (abs(h) > 0.0)
  val[sel] = abs(h[sel])^nu*besselK(abs(h[sel]), nu = nu)/gamma(nu)/2^(nu-1)
  val
}
kernel_gauss  <- function(h){exp(-h^2)}
kernel_generalized_cauchy <- function(u, alpha, beta) {(1+abs(u)^alpha)^{-beta}}
kernel_cauchy <- function(u, nu) {kernel_generalized_cauchy(u, alpha = 2, beta = nu)}

# computing a spectral model on R^d 
#' @param sp a list containing spectral waves
#'  - sp$omega a [np, d+1] matrix of space-time random frequencies
#'  - sp$phi a [np] vector of the random phases
#'  - sp$gamma a [np,nvar] matrix of the weights for the random waves
#'  @param st a [n, d+1] matrix of the d+1 coordinates of the n target space-time points
#'  @verbose a boolean
#'  @value a [n, nvar] matrix of the simulated values
compute_spectral <- function(sp, st, verbose = FALSE) {
  stopifnot(is.matrix(sp$omega))
  p <- nrow(sp$omega)
  d <- ncol(sp$omega)
  stopifnot(is.matrix(st)&(ncol(st) == d))
  n <- nrow(st)
  stopifnot(is.matrix(sp$lambda)&(nrow(sp$lambda) == d))
  nvar <- ncol(sp$lambda)
  stopifnot((length(sp$phi) == p)&(nrow(sp$gamma) == p))
  if(verbose) {
    print(paste0(">>> Computing using the spectral method:"))
    print(paste0(">>>     dimension of the space-time    : ", d))
    print(paste0(">>>     number of variables            : ", nvar))
    print(paste0(">>>     number of spectral components  : ", p))
    print(paste0(">>>     number of target points        : ", n))
  }
  t(cos(sp$omega %*% t(st) + sp$phi %*% t(rep(1.0, n)))) %*% sp$gamma
}

# rotations
rot1D <- function() {
  matrix(1,1,1)
}
rot2D <- function(angles = 0.0) {
  theta = angles[1] * pi / 180
  cs = cos(theta)
  sn = sin(theta)
  matrix(c(cs, sn, -sn, cs), 2, 2)
}
rot3D <- function(angles = rep(0.0, 3), flag.euler = TRUE) {
  # yaw
  alpha = angles[1] * pi / 180 
  cs = cos(alpha); sn = sin(alpha)
  Rz = matrix(c(cs, sn, 0, -sn, cs, 0, 0, 0, 1), 3, 3)
  # pitch
  beta = angles[2] * pi / 180 
  cs = cos(beta); sn = sin(beta)
  Ry = matrix(c(cs, 0, -sn, 0, 1, 0, sn, 0, cs), 3, 3)
  # roll
  gamma = angles[3] * pi / 180
  cs = cos(gamma); sn = sin(gamma)
  Rx = matrix(c(1, 0, 0, 0, cs, sn, 0, -sn, cs), 3, 3)
  if(flag.euler) {
    # extrinsic rotation whose (improper) Euler angles are (alpha, beta, gamma)
    Rx %*% Ry %*% Rz
  } else {
    # intrinsic rotation whose Tait-Bryan angles are (alpha, beta, gamma)
    Rz %*% Ry %*% Rx
  }
}

# Gneiting model for $d %\in% c(1,2,3)$
GneitingModel <- function(alpha = 1.0, beta = 1.0, nu = 1/2, type = 1,
                          ranges = c(1,1,1), rot = diag(length(ranges)-1), 
                          verbose = FALSE) {
  # checking parameters
  d = length(ranges) - 1
  stopifnot(d > 0)
  stopifnot(all(ranges > 0.0))
  stopifnot((alpha > 0.0)&(alpha <= 2.0))
  stopifnot((beta > 0.0)&(beta <= 1.0))
  stopifnot(is.null(nu)|(nu > 0))
  # type and trace covariance functions
  if (type == 0) {
    fn_type  <- function(){paste0("Gneiting-Gauss(", d, ")")}
    fn_cov_s <- function(h){kernel_gauss(h)}
    nu = NULL
  }
  if (type == 1) {
    fn_type  <- function(){paste0("Gneiting-Matern(", d, ")")}
    fn_cov_s <- function(h){kernel_matern(h, nu)}
  }
  if (type == 2) {
    fn_type  <- function(){paste0("Gneiting-Cauchy(", d, ")")}
    fn_cov_s <- function(h){kernel_generalized_cauchy(h, alpha = 2, beta = nu)}
  }
  fn_cov_t <- function(u) {kernel_generalized_cauchy(u, alpha = alpha, beta = beta*d/2)}
  
  if (d == 1) {
    ani_factor = matrix(1/ranges[1], nrow = 1, ncol = 1)  
  } else {
    ani_factor = diag(1/ranges[1:d])
  }
  
  dist_aniso = function(h){
    h = as.matrix(h)
    apply(X = as.matrix(h) %*% rot %*% ani_factor, 
          MARGIN = 1, 
          FUN = function(u){sqrt(sum(u^2))})
  }
  
  # space-time covariance function
  fn_cov <- function(H) {
    stopifnot(is.matrix(H)&(ncol(H) == d+1))
    ct = fn_cov_t(H[,d+1]/ranges[d+1])
    ct * fn_cov_s(dist_aniso(H[,1:d]) * ct^(1/d))
  }
  
  # simulation
  fn_sim <- function(n, seed = NULL, flag.double = TRUE) {
    if(!is.null(seed)) {set.seed(seed)}
    # simulation of the spectral components
    omega = matrix(NaN, nrow = n, ncol = d+1)
    # space spectral components with geometrical anisotropy
    if (type == 0) {
      R   = rep(1, n) 
    } else if (type == 1) {
      R   = 1/(4*rgamma(n = n, shape = nu))
    } else if (type == 2) {
      R   = rgamma(n = n, shape = nu)
    }
    omega[, 1:d] = sqrt(2*R) * matrix(rnorm(n*d), nrow = n, ncol = d) 
    # time spectral components
    lambda = (apply(X = as.matrix(omega[,1:d]), 
                   MARGIN = 1, 
                   FUN = function(u){sum(u^2)})/(4*R))^(1/beta)
    
    sim_S = rep(NaN, length(lambda))
    for (i in seq_along(lambda)) {
      sim_S[i] = simulate_stable_unilateral_exptilted(n = 1, alpha = beta, lambda = lambda[i], 
                                                   flag.double = flag.double)$X
    }
    sim_T = simulate_stable_bilateral(n = n, alpha = alpha)
    omega[, d+1] = (sim_S*lambda)^(1/alpha) * sim_T / ranges[d+1]
    omega[, 1:d] = omega[, 1:d] %*%
      ani_factor %*% t(rot)
    sp = list(
      omega = omega, 
      phi   = runif(n, min = 0.0, max = 2*pi),
      gamma = matrix(sqrt(-2*log(runif(n))/n), nrow = n, ncol = 1)
    )
    # computing the spectrum on the targets
    fn_compute = function(st, verbose = FALSE) {
      compute_spectral(sp = sp, st = st, verbose = verbose) 
    }
    
    list(
      spectrum = function(){sp},
      compute  = fn_compute
    )
  }
  if (verbose) {
    print(paste0("Creating a space-time ", fn_type(), " model:"))
    print(paste0(" - dimension of the space :", d))
    print(paste0(" - alpha                  :", alpha))
    print(paste0(" - beta                   :", beta))
    print(paste0(" - nu                     :", nu))
  }
  
  # return values
  list(
    type = fn_type,
    dimension = function(){d},
    alpha = function(){alpha},
    beta  = function(){beta},
    nu    = function(){nu},
    scales = function(){ranges},
    rotation = function(){rot},
    covariance = fn_cov, 
    simulate = fn_sim)
}

# comparison of the covariance values evaluated by R reference functions and by gstlearn methods
h   = seq(from = 0, to = 10, length.out = 100)
val_in_R = rep(NaN, length(h))
val_in_gstlearn = rep(NaN, length(h))
EPS_10 = 1.0e-10

# --------------------------------------------------------------
# Test of the generalized cauchy kernel using CorAniso on R^1
# --------------------------------------------------------------
setDefaultSpace(SpaceRN_create(ndim = 1))
corT = CorAniso(type = ECov_CAUCHY_GEN(), ctxt = CovContext(nvar = 1, ndim = 1))
for (a_t in c(0.5, 1.0, 1.5, 2.0) ){
for (alpha in c(0.25, 1.0, 1.5, 2.0)) {
for (beta in c(1.0, 3/2, 1.0, 2)) {
print(paste0(">>> Testing generalized Cauchy kernel with a_t = ", a_t, " alpha = ", alpha, " beta = ", beta))  
err = corT$setScaleDim(idim = 0, scale = a_t)
err = corT$setParam(param = alpha, ipar = 0)
err = corT$setParam(param = beta, ipar = 1)
val_in_gstlearn = corT$sample(h = h, codir = c(1.0), ivar = 0, jvar = 0)
val_in_R = kernel_generalized_cauchy(h/a_t, alpha, beta)
stopifnot(all(abs(val_in_gstlearn - val_in_R) <= EPS_10))
if (all(abs(val_in_gstlearn - val_in_R) <= EPS_10)) {
print(paste0(">>> Test OK"))  
} else {
print(paste0(">>> Test KO"))
break
}
}    
}
}

# --------------------------------------------------------------
# Test of the Gneiting space-time covariance using CorGneiting
# --------------------------------------------------------------
set.seed(seed = 122582)
# parameters
list_alpha = c(0.25, 1.0, 1.5, 2.0)
list_beta = c(0.5, 3/4, 1.0)
list_nu   = c(1/2, 1.0, 3/4, 2.0)
sep = 1.0
# anisotropies
scales_all = c(1, 2, 2.5, 3)
angles = list()
angles[[1+length(angles)]] = 0
angles[[1+length(angles)]] = c(30,0)
angles[[1+length(angles)]] = c(0,0,0)
type_label = c("Gneiting-Gauss", "Gneiting-Matern", "Gneiting-Cauchy")
for (type in 0:2 ){
print(paste0(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"))
print(paste0(">>> Testing the ", type_label[1+type]))
print(paste0(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"))
if (type == 0) {
  nus = 1.0
} else {
  nus = list_nu
}

for (d in 1:3 ){
print(paste0(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"))
print(paste0(">>> Testing the model in dimension d = ", d))
print(paste0(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"))

defineDefaultSpaceRnT(ndim = d)
displayDefaultSpace()
# rotation for the spatial anisotropy
scales = scales_all[c(1:d,4)]
rot = switch (d, rot1D(), rot2D(angles[[d]]), rot3D(angles[[d]]))
# directions of the spatial traces
if (d == 1) {
dir_u   = matrix(c(1,0), 1, 2)
}
if (d == 2) {
dir_u   = matrix(c(1,0,0, 1,1,0, 0,1,0, -1,1,0), nrow = 4, ncol = 3, byrow = TRUE)
}
if (d == 3) {
dir_u   = matrix(c(1,0,0,0, 0,1,0,0, 0,0,1,0), nrow = 3, ncol = 4, byrow = TRUE)
}
dir_u = dir_u / apply(X = dir_u, MARGIN = 1, FUN = function(u){sqrt(sum(u^2))})

# directions of the space-time traces
ndir = 4
dir_st = matrix(rnorm(ndir*(d+1)), nrow = ndir, ncol = d+1)
dir_st = dir_st / apply(X = dir_st, MARGIN = 1, FUN = function(u){sqrt(sum(u^2))})
for (alpha in list_alpha) {
for (beta in list_beta) {
for (nu in nus)
  if (type == 0) {
  print(paste0(">>> Testing Gneiting model with alpha = ", alpha, " beta = ", beta))  
  params = c(alpha, beta*d/2)
  } else {
  print(paste0(">>> Testing Gneiting model with alpha = ", alpha, " beta = ", beta, " nu = ", nu))  
  params = c(alpha, beta*d/2, nu)
  }
  # do the tests
  # Gneiting model defined in R
  mg_in_R = GneitingModel(alpha, beta, nu, ranges = scales, rot = rot, type = type, verbose = TRUE)
  # Gneiting model defined in gstlearn
  mg_in_gstlearn  = CorGneiting_create(
  ctxt = CovContext(nvar = 1, ndim = d+1), 
  type = switch(paste(type), "0" = ECov_GNEITING_G(), "1" = ECov_GNEITING_M(), "2" = ECov_GNEITING_C()),
  params = params, 
  ranges = scales,      # vector of length d+1
  angles = angles[[d]], # vector of length d (no angles in time)
  separability = 1.0,   # scalar in [0,1]
  flagRange = FALSE
  )
  
  res_test = rep(TRUE, 3)
  # temporal trace
  u0  = c(rep(0.0, d), 1.0)
  val_in_R = mg_in_R$covariance(H = h %*% t(u0))
  val_in_gstlearn = mg_in_gstlearn$sample(h = h, codir = u0, ivar = 0, jvar = 0)
  res_test[1] = all(abs(val_in_R - val_in_gstlearn) < EPS_10)
  # space trace
  for (idir in 1:nrow(dir_u)) {
  u0 = dir_u[idir,]
  val_in_gstlearn = mg_in_gstlearn$sample(h = h, codir = u0, ivar = 0, jvar = 0)
  val_in_R = mg_in_R$covariance(H = h %*% t(u0)) 
  res_test[2] = all(res_test[2], (all(abs(val_in_gstlearn - val_in_R) < EPS_10)))
  }

  # space-time trace
  for (idir in 1:nrow(dir_st)) {
  u0 = dir_st[idir,]
  val_in_gstlearn = mg_in_gstlearn$sample(h = h, codir = u0, ivar = 0, jvar = 0)
  val_in_R = mg_in_R$covariance(H = h %*% t(u0)) 
  res_test[3] = all(res_test[3], (all(abs(val_in_gstlearn - val_in_R) < EPS_10)))
  }
  # end of the test
  if (all(res_test)) {
    print(paste0(">>> Tests OK"))
  } else {
    print(paste0(">>> Tests KO"))
    break
}}}}}
