suppressWarnings(suppressMessages(library(gstlearn)))
rm(list = ls())

flag.range = FALSE
ranges  = c(1,2)
angles  = c(30, 0) # en degrés
rescale = c(1/2, 2); rr = c(1.0, rescale)  #on va multiplier les ranges par 3. et 2.
params  = c(1/2, 1, 2)

p0 = SpacePoint(c(0,0))
p1 = SpacePoint(c(1,0))
p2 = SpacePoint(sqrt(2)*c(1,1))
p3 = SpacePoint(c(0,1))
p4 = SpacePoint(sqrt(2)*c(-1,1))

#' Test mono-variable (check the anisotropy)
ivar = 1
cor_mono = CorMatern(
  ranges = ranges, angle = angles, coeffScales = NULL, params = params[ivar], flagRange = flag.range)

stopifnot(cor_mono$getNVar() == 1)

mod_mono = Model_createFromParam(type = ECov_MATERN(), 
                             ranges = ranges*rr[ivar], angles = angles, 
                             sill = 1.0, param = c(params[ivar]),
                             flagRange = flag.range)

stopifnot(cor_mono$eval(p0, p1, ivar = ivar-1, jvar = ivar-1) == mod_mono$getCovAnisoList()$eval(p0, p1))
stopifnot(cor_mono$eval(p0, p2, ivar = ivar-1, jvar = ivar-1) == mod_mono$getCovAnisoList()$eval(p0, p2))
stopifnot(cor_mono$eval(p0, p3, ivar = ivar-1, jvar = ivar-1) == mod_mono$getCovAnisoList()$eval(p0, p3))
stopifnot(cor_mono$eval(p0, p4, ivar = ivar-1, jvar = ivar-1) == mod_mono$getCovAnisoList()$eval(p0, p4))

#' Test tri-variable
cor_tri = CorMatern(ranges = ranges, angle = angles, 
                    coeffScales = rescale, params = params, flagRange = flag.range)

nvar = cor_tri$getNVar()
# nvar = 3 # TODO: it does not work!
stopifnot(nvar == 3) 

#' Test de la correlation maximale
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
r_ij  = outer(X = 1:nvar, Y = 1:nvar, FUN = function(i,j){sqrt((rr[i]^2+rr[j]^2)/2)})
nu_ij = outer(X = 1:nvar, Y = 1:nvar, FUN = function(i,j){(params[i]+params[j])/2})
corMax = matrix(NaN, nvar, nvar)
for (i in 1:nvar) {
  for (j in i:nvar) {
    corMax[i,j] = cor_tri$getCorMax(ivar = i-1, jvar = j-1)
    corMax[j,i] = corMax[i,j]
  }
}
stopifnot(all(abs(corMax - tau) < 1.e-12))

# compute the correlations
for (ivar in 1:nvar) {
  for (jvar in ivar:nvar) {
    mod_mat = Model_createFromParam(type = ECov_MATERN(), 
                                    ranges = ranges*r_ij[ivar, jvar], angles = angles, 
                                    sill = tau[ivar,jvar],
                                    param = nu_ij[ivar, jvar],
                                    flagRange = flag.range)
    stopifnot(abs(cor_tri$eval(p0, p1, ivar = ivar-1, jvar = jvar-1) - mod_mat$getCovAnisoList()$eval(p0, p1)) < 1.e-12)
    stopifnot(abs(cor_tri$eval(p0, p2, ivar = ivar-1, jvar = jvar-1) - mod_mat$getCovAnisoList()$eval(p0, p2)) < 1.e-12)
    stopifnot(abs(cor_tri$eval(p0, p3, ivar = ivar-1, jvar = jvar-1) - mod_mat$getCovAnisoList()$eval(p0, p3)) < 1.e-12)
    stopifnot(abs(cor_tri$eval(p0, p4, ivar = ivar-1, jvar = jvar-1) - mod_mat$getCovAnisoList()$eval(p0, p4)) < 1.e-12)   
  }
}

print("All tests are ok.")
