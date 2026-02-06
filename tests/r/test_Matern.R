suppressWarnings(suppressMessages(library(gstlearn)))
rm(list = ls())

#cat("PID =", Sys.getpid(), "\n")
#invisible(readLines(file("stdin"), n = 1))

flag.range = FALSE
ranges  = c(1,2)
angles  = c(30, 0) # en degrés
rescale = c(1, 2); rr = c(1.0, rescale)  #on va multiplier les ranges par 1. et 2.
params  = c(1/2, 1., 2.)

#' definition of the multi-variate parameters
nvar = length(rr)
r_ij  = outer(X = 1:nvar, Y = 1:nvar, FUN = function(i,j){sqrt((rr[i]^2+rr[j]^2)/2)})
nu_ij = outer(X = 1:nvar, Y = 1:nvar, FUN = function(i,j){(params[i]+params[j])/2})
compute_tau <- function(nu, r) {
  stopifnot(length(nu) == length(r))
  stopifnot(all(nu > 0)&all(r > 0))
  nv = length(nu)
  m = matrix(NaN, nv, nv)
  outer(X = 1:nv, Y = 1:nv, FUN = function(i,j) {
    nu_ij = (nu[i]+nu[j])/2 
    r_ij = sqrt((r[i]^2+r[j]^2)/2)
    gamma(nu_ij)/sqrt(gamma(nu[i])*gamma(nu[j])) *
      r[i]^nu[i] * r[j]^nu[j] / r_ij^(2*nu_ij)
  })
}
tau = compute_tau(nu = params, r = rr)
stopifnot(all(eigen(tau)$values >= 0)) # It should be a positive definite matrix

# definition of the test points
p0 = SpacePoint(c(0,0))
p1 = SpacePoint(c(1,0))
p2 = SpacePoint(sqrt(2)*c(1,1))
p3 = SpacePoint(c(0,1))
p4 = SpacePoint(sqrt(2)*c(-1,1))
pts = list(p1,p2,p3,p4)

#' Test mono-variable (check the anisotropy)
ivar = 1
sigma = matrix(1.0, nrow = 1, ncol = 1)
nvar = nrow(sigma)
gsSigma = MatrixSymmetric(nvar)
for (i in 1:nvar) {
  for (j in 1:nvar) {
    gsSigma$setValue(i-1, j-1, sigma[i, j])
  }
}
gsSigma$isDefinitePositive()

# XF: Old implementation
# cor_mono = CorMatern(
#   ranges = ranges, 
#   angle = angles, 
#   coeffScales = rr[ivar], 
#   params = params[ivar], 
#   sigma = gsSigma, 
#   flagRange = flag.range)

cor_mono = CorMatern_create(
  type = ECov_MATERN(),
  ctxt = CovContext(ndim = 2, nvar = 1), 
  params = params[ivar], 
  kappas =  rr[ivar], 
  ranges = ranges, 
  angle = angles, 
  flagRange = flag.range)

nvar = cor_mono$getNVar()
stopifnot(nvar == 1)

# compute the correlations
for (ivar in 1:nvar) {
  for (jvar in ivar:nvar) {
    mod_mat = Model_createFromParam(type = ECov_MATERN(), 
                                    ranges = ranges*r_ij[ivar, jvar], angles = angles, 
                                    sill = tau[ivar,jvar],
                                    param = nu_ij[ivar, jvar],
                                    flagRange = flag.range)
    for (i in 1:4)
    {
      pcur = pts[[i]]
      res = cor_mono$evalCov(p0, pcur, ivar = ivar-1, jvar = jvar-1)
      cat("ivar ",ivar," jvar ", jvar," pt: ", i," res ", round(res,4),"\n")
      stopifnot(abs(res - mod_mat$evalCov(p0, pcur)) < 1.e-12)
    }
  }
}
print("Mono-variable test is ok.")

#' Test tri-variable
nvar = nrow(tau)
gsSigma = MatrixSymmetric(nvar)
for (i in 1:nvar) {
  for (j in 1:nvar) {
    gsSigma$setValue(i-1, j-1, tau[i, j])
  }
}
gsSigma$isDefinitePositive()

# XF: Old implementation
# cor_tri = CorMatern(
#   ranges = ranges, 
#   angle = angles, 
#   coeffScales = rr, 
#   params = params, 
#   sigma = gsSigma, 
#   flagRange = flag.range)

cor_tri = CorMatern_create(
  type = ECov_MATERN(),
  ctxt = CovContext(ndim = 2, nvar = 3), 
  params = params, 
  kappas =  rr, 
  ranges = ranges, 
  angle = angles, 
  flagRange = flag.range)

nvar = cor_tri$getNVar()
stopifnot(nvar == 3)

C0 = matrix(NaN, nvar, nvar)
for (i in 1:nvar) {
  for (j in i:nvar) {
    C0[i,j] = cor_tri$getC0(ivar = i-1, jvar = j-1)
    C0[j,i] = C0[i,j]
  }
}
print(round(C0,3))
stopifnot(all(abs(C0 - tau) < 1.e-12))

# compute the correlations
for (ivar in 1:nvar) {
  for (jvar in ivar:nvar) {
    mod_mat = Model_createFromParam(type = ECov_MATERN(), 
                                    ranges = ranges/r_ij[ivar, jvar], angles = angles, 
                                    sill = tau[ivar,jvar],
                                    param = nu_ij[ivar, jvar],
                                    flagRange = flag.range)
    for (i in 1:4)
    {
 	    pcur = pts[[i]]
 	    res = cor_tri$evalCov(p0, pcur, ivar = ivar-1, jvar = jvar-1)
 	    cat("ivar ",ivar," jvar ", jvar," pt: ", i," res ", round(res,4),"\n")
 	    stopifnot(abs(res - mod_mat$evalCov(p0, pcur)) < 1.e-12)
    }
  }
}

ctxt = CovContext(nvar = 3)
model = ModelGeneric(ctxt)
a = model$setCov(cor_tri)
grid = DbGrid_create(c(100,100))
nmc = NamingConvention("Simu")
simuSpectral(dbin = NULL, dbout = grid, model = model, nbsimu = 1, seed = 43431, ns = 100, nd = 100, 
verbose = T, namconv = nmc)

print("All tests are ok.")
