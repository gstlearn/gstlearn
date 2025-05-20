#title: "Regression with spatially varying coefficients in gstlearn"
#author: "Xavier Freulon, Didier Renard, and the gstlearn Team"
#date: "2025-05-15"

knitr::opts_chunk$set(echo = FALSE)
rm(list = ls())
suppressWarnings(suppressMessages(library(gstlearn)))

set.seed(43243)
# global parameters
pal   = rev(RColorBrewer::brewer.pal(6, "RdYlBu"))
opers = EStatOption_fromKeys(c("NUM", "MINI", "MAXI", "MEAN", "STDV"))
err   = OptCst_defineByKey("ASP",0)
err = OptCst_define(ECst_NTROW(), -1)
err = OptCst_define(ECst_NTCOL(), -1)
ndim = 2
err = defineDefaultSpace(ESpaceType_RN(), ndim = ndim)

p_sel = 0.2

verbose  = TRUE
flag_est = TRUE
flag_std = TRUE
flag_varz = FALSE # TODO: it does not work for the estimation of the coefficient 

# Function SVC
ESTIM_KEY = c("obs", "res", "drift", "coeff")

#  Generate the input data (in the 1x1 square) with 'nvar' variables and 'nfex' external drifts
prefix = "SVC"
ndat   = 10
seedin = 13227
nout    = 10
seedout = 134484
nvar = 1
nfac = 2
dbin   = Db_createFillRandom(ndat = ndat, ndim = ndim, 
                             nvar = 1+nfac, nfex = nfac, ncode = 0, 
                             varmax = 0., selRatio = p_sel, heteroRatio = VectorDouble(), 
                             coormin = VectorDouble(), coormax = VectorDouble(), seed = seedin)
dbout   = Db_createFillRandom(ndat = nout, ndim = ndim, 
                              nvar = 0, nfex = nfac, ncode = 0, 
                              varmax = 0., selRatio = p_sel, heteroRatio = VectorDouble(), 
                              coormin = VectorDouble(), coormax = VectorDouble(), seed = seedout)
nm_fac = dbin$getNamesByLocator(locatorType = ELoc_F())
nm_var = dbin$getNamesByLocator(locatorType = ELoc_Z())
if (nfac > 0) {
  Z = dbin[nm_var[1]]
  for (i in seq_along(nm_fac)) {
    Z = dbin[nm_var[1+i]] * dbin[nm_fac[i]]
  }
  err = dbin$setColumn(tab = Z, name = "z")
  err = dbin$setLocators(names = "z", locatorType = ELoc_Z(), cleanSameLocator = TRUE)
  err = dbin$deleteColumns(names = paste("z", 1:(1+nfac), sep = "-"))
}
nm_var = dbin$getNamesByLocator(locatorType = ELoc_Z())

# defining the model of the effects
C0 = matrix(runif(nvar*nvar, min = -1, max = 1), ncol = nvar, nrow = nvar)
C0 = C0 %*% t(C0)
model = Model_createFromEnvironment(nvar = nvar, ndim = ndim)
err = model$addCovFromParam(type = ECov_MATERN(), sills = C0, param = 1.0, range = 2.0)
err = model$setDriftIRF(order = 0, nfex = nfac)
if (verbose) {
  err = model$display()
}

# neighborhood
# neigh = NeighUnique_create()
neigh = NeighMoving_create(nmini = 7, nmaxi = 7)

# list of estimation to be performed
estim_list = list()
estim_list[[1+length(estim_list)]] <- list(type = 1, idx = -1)
estim_list[[1+length(estim_list)]] <- list(type = 2, idx = -1)
estim_list[[1+length(estim_list)]] <- list(type = 3, idx = -1)
for (i in 0:nfac){
  estim_list[[1+length(estim_list)]] <- list(type = 4, idx = i)
}
nms = NULL
for (i in seq_along(estim_list)) {
  if (estim_list[[i]]$type %in% 1:3) {
    nms = c(nms, paste(sep = ".", ESTIM_KEY[estim_list[[i]]$type], nm_var))
  } else if (estim_list[[i]]$type %in% 4) {
    nms = c(nms, paste(sep = ".", ESTIM_KEY[estim_list[[i]]$type], estim_list[[i]]$idx, nm_var))
  }
}
# ----------------------------------------------------------------------------
# Using the test function
# ----------------------------------------------------------------------------
err = kriging_SVC(dbin = dbin, dbout = dbout, model = model, neigh = neigh,
                flag_est = flag_est, flag_std = flag_std, flag_varz = flag_varz, 
                estim_list = estim_list,
                prefix = "SVC", verbose = verbose
)
# ----------------------------------------------------------------------------
# Using the *kriging* function
# ----------------------------------------------------------------------------
OptCustom_define("NotOptimSimpleCase", 1)
pr = NULL
if (nfac == 0) {pr = "OK"}
if (nfac >  0) {pr = "KED"}
for (i in seq_along(estim_list)) {
  
  # estimation of the observed variable
  if (estim_list[[i]]$type == 1) {
    # Kriging the values
    prefix = paste(pr, ESTIM_KEY[estim_list[[i]]$type], sep = ".")
    err = kriging(dbin = dbin, dbout = dbout, model = model, neigh = neigh, 
                  flag_est = flag_est, flag_std = flag_std, flag_varz = flag_varz, 
                  namconv = NamingConvention(prefix)
    )
  }
  
  # estimation of the residual
  if (estim_list[[i]]$type == 2) {
    prefix = paste(pr, ESTIM_KEY[estim_list[[i]]$type], sep = ".")
    err = model$setCovFiltered(icov = 0, filter = FALSE)
    for (ifac in 0:nfac) {
      err = model$getDriftList()$setFiltered(i = ifac, filter = TRUE)
    }
    if(verbose) {
      err = model$display()
    }
    err = kriging(dbin = dbin, dbout = dbout, model = model, neigh = neigh, 
                  flag_est = flag_est, flag_std = flag_std, flag_varz = flag_varz, 
                  namconv = NamingConvention(prefix)
    )
  }
  
  # estimation of the drift
  if (estim_list[[i]]$type == 3) {
    prefix = paste(pr, ESTIM_KEY[estim_list[[i]]$type], sep = ".")
    err = model$setCovFiltered(icov = 0, filter = TRUE)
    for (ifac in 0:nfac) {
      err = model$getDriftList()$setFiltered(i = ifac, filter = FALSE)
    }
    if(verbose) {
      err = model$display()
    }
    err = kriging(dbin = dbin, dbout = dbout, model = model, neigh = neigh, 
                  flag_est = flag_est, flag_std = flag_std, flag_varz = flag_varz, 
                  namconv = NamingConvention(prefix)
    )
  }
  
  # Kriging the coefficients
  if (estim_list[[i]]$type == 4) {
    ifac = estim_list[[i]]$idx
    prefix = paste(pr, ESTIM_KEY[estim_list[[i]]$type], ifac, sep = ".")
    # selection of the drift to compute...
    err = model$setCovFiltered(icov = 0, filter = (estim_list[[i]]$idx > 0))
    for (jfac in 0:nfac) {
      flag = !(ifac == jfac)
      err = model$getDriftList()$setFiltered(i = jfac, filter = flag)
    }
    if(verbose) {
      err = model$display()
    }
    prefix = paste(pr, ESTIM_KEY[4], ifac, sep = ".")
    
    err = kriging(dbin = dbin, dbout = dbout, model = model, neigh = neigh, 
                  flag_est = flag_est, flag_std = flag_std, flag_varz = flag_varz, 
                  namconv = NamingConvention(prefix)
    )
    if (ifac > 0) {
      if(flag_est){
        nm_est = paste(prefix, nm_var, "estim", sep = ".")
        dbout[nm_est] = dbout[nm_est]/dbout[nm_fac[ifac]]
      }
      if (flag_std){
        nm_std = paste(prefix, nm_var, "stdev", sep = ".")
        dbout[nm_std] = dbout[nm_std]/abs(dbout[nm_fac[ifac]])
      }
      if (flag_varz){
        nm_varz = paste(prefix, nm_var, "varz", sep = ".")
        dbout[nm_varz] = dbout[nm_varz]/abs(dbout[nm_fac[ifac]])
      }
    }
  }
  
}
# QC of the results
if (length(nms)) {
  for (nm in nms) {
    for (cs in c("estim", "stdev", "varz")[c(flag_est, flag_std, flag_varz)]) {
      nm_1 = paste("SVC", nm, cs, sep = ".")
      nm_2 = paste(pr, nm, cs, sep = ".")
      if(all(abs(dbout$getColumn(name = nm_1, useSel = TRUE) - dbout$getColumn(name = nm_2, useSel = TRUE)) < 1e-6)) {
        res = paste0("Test passed: ", "SVC == ", pr)
      } else {
        res = paste0("Test failed: ", "SVC != ", pr)
      }
      print(paste0(res, "(nvar = ", nvar, " nfac  = ", nfac,"):", nm, "-", cs))
    } # loop over estim/stdev
  } # loop over the estimation
}

opers = EStatOption_fromKeys(c("NUM", "MINI", "MAXI", "MEAN", "STDV"))
knitr::kable(dbStatisticsMono(db = dbout, names = paste0(c(pr, "SVC"),".*"), 
                              opers = opers)$toTL(), digits = 4, caption = "Statistics")


# Test SVC with centred factors

# defining the model of the effects
# with constant drifts...
nfac   = 3
fac    = c(1,2,3)[1:nfac]
param  = c(1/2, 1, 3/2)[1:nfac]
ranges = c(0.25, 0.5, 1.0)[1:nfac]
nvar   = length(fac)
model = Model_createFromEnvironment(nvar = nvar, ndim = ndim)
model_mono = Model_createFromEnvironment(nvar = 1, ndim = ndim)
for (i in 1:nvar) {
  err = model$addCovFromParam(type = ECov_MATERN(), 
                              sills = diag(as.numeric(1:nvar == i)), 
                              param = param[i], range = ranges[i])
  err = model_mono$addCovFromParam(type = ECov_MATERN(), 
                                   sill = fac[i]^2, 
                                   param = param[i], range = ranges[i])
}
err = model$setDriftIRF(order = 0, nfex = nvar-1)
err = model_mono$setDriftIRF(order = 0, nfex = 0)

if (verbose) {
  err = model$display()
  err = model_mono$display()
}

# input/output data bases
ndat   = 13
seedin = 13227
nout    = 17
seedout = 134484
dbin   = Db_createFillRandom(ndat = ndat, ndim = ndim, 
                             nvar = nvar, nfex = 0, ncode = 0, 
                             varmax = 0., selRatio = p_sel, heteroRatio = VectorDouble(), 
                             coormin = VectorDouble(), coormax = VectorDouble(), seed = seedin)
dbout   = Db_createFillRandom(ndat = nout, ndim = ndim, 
                              nvar = 0, nfex = 0, ncode = 0, 
                              varmax = 0., selRatio = p_sel, heteroRatio = VectorDouble(), 
                              coormin = VectorDouble(), coormax = VectorDouble(), seed = seedout)

nm_var = dbin$getNamesByLocator(locatorType = ELoc_Z())
dbin["z"] = as.matrix(dbin[nm_var]) %*% fac
err = dbin$setLocators(names = "z", locatorType = ELoc_Z(), cleanSameLocator = TRUE)
nm_var = dbin$getNamesByLocator(locatorType = ELoc_Z())


for (ifac in seq_along(fac[-1])) {
  dbin[paste("f", ifac, sep = "-")] = rep(fac[-1][ifac], dbin$getNSample())
  dbout[paste("f", ifac, sep = "-")] = rep(fac[-1][ifac], dbout$getNSample())
}
err = dbin$setLocators(
  names = paste0("f-", 1:(nvar-1)), 
  locatorType = ELoc_F(), 
  cleanSameLocator = TRUE)
err = dbout$setLocators(
  names = paste0("f-", 1:(nvar-1)), 
  locatorType = ELoc_F(), 
  cleanSameLocator = TRUE)

nm_fac = dbin$getNamesByLocator(locatorType = ELoc_F())
nfac = length(nm_fac)

# neighborhood
neigh = NeighUnique_create()
# neigh = NeighMoving_create(nmini = 7, nmaxi = 7)

# list of estimation to be performed
estim_list = list()
estim_list[[1+length(estim_list)]] <- list(type = 1, idx = -1)
estim_list[[1+length(estim_list)]] <- list(type = 2, idx = -1)
estim_list[[1+length(estim_list)]] <- list(type = 3, idx = -1)
for (i in 0:nfac){
  estim_list[[1+length(estim_list)]] <- list(type = 4, idx = i)
}
nms = NULL
for (i in seq_along(estim_list)) {
  if (estim_list[[i]]$type %in% 1:3) {
    nms = c(nms, paste(sep = ".", ESTIM_KEY[estim_list[[i]]$type], nm_var))
  } else if (estim_list[[i]]$type %in% 4) {
    nms = c(nms, paste(sep = ".", ESTIM_KEY[estim_list[[i]]$type], estim_list[[i]]$idx, nm_var))
  }
}

# ----------------------------------------------------------------------------
# Using the test function
# ----------------------------------------------------------------------------
flag_centeredFactors = TRUE
err = kriging_SVC(dbin = dbin, dbout = dbout, model = model, neigh = neigh,
                flag_est = flag_est, flag_std = flag_std, flag_varz = flag_varz,
                estim_list = estim_list, flag_centeredFactors = flag_centeredFactors,
                prefix = "SVC", verbose = verbose
)
# ----------------------------------------------------------------------------
# Using the *kriging* function
# ----------------------------------------------------------------------------
pr = "OK"
for (i in seq_along(estim_list)) {
  # estimation of the observed variable
  if (estim_list[[i]]$type == 1) {
    # Kriging the values
    for (icov in 0:2) {
      err = model_mono$setCovFiltered(icov = icov, filter = FALSE)
    }
    err = model_mono$getDriftList()$setFiltered(i = 0, filter = FALSE)
    prefix = paste(pr, ESTIM_KEY[estim_list[[i]]$type], sep = ".")
    err = kriging(dbin = dbin, dbout = dbout, model = model_mono, neigh = neigh, 
                  flag_est = flag_est, flag_std = flag_std, flag_varz = flag_varz, 
                  namconv = NamingConvention(prefix)
    )
  }
  
  # estimation of the residual
  if (estim_list[[i]]$type == 2) {
    for (icov in 0:2) {
      err = model_mono$setCovFiltered(icov = icov, filter = FALSE)
    }
    err = model_mono$getDriftList()$setFiltered(i = 0, filter = TRUE)
    prefix = paste(pr, ESTIM_KEY[estim_list[[i]]$type], sep = ".")
    err = kriging(dbin = dbin, dbout = dbout, model = model_mono, neigh = neigh, 
                  flag_est = flag_est, flag_std = flag_std, flag_varz = flag_varz, 
                  namconv = NamingConvention(prefix)
    )  }
  
  # estimation of the drift
  if (estim_list[[i]]$type == 3) {
    for (icov in 0:2) {
      err = model_mono$setCovFiltered(icov = icov, filter = TRUE)
    }
    err = model_mono$getDriftList()$setFiltered(i = 0, filter = FALSE)
    prefix = paste(pr, ESTIM_KEY[estim_list[[i]]$type], sep = ".")
    err = kriging(dbin = dbin, dbout = dbout, model = model_mono, neigh = neigh, 
                  flag_est = flag_est, flag_std = flag_std, flag_varz = flag_varz, 
                  namconv = NamingConvention(prefix)
    )  
  }
  
  # Kriging the coefficients
  if (estim_list[[i]]$type == 4) {
    ifac = estim_list[[i]]$idx
    for (icov in 0:2) {
      err = model_mono$setCovFiltered(icov = icov, filter = !(icov == ifac))
    }
    err = model_mono$getDriftList()$setFiltered(i = 0, filter = (ifac != 0))
    prefix = paste(pr, ESTIM_KEY[estim_list[[i]]$type], ifac, sep = ".")
    err = kriging(dbin = dbin, dbout = dbout, model = model_mono, neigh = neigh, 
                  flag_est = flag_est, flag_std = flag_std, flag_varz = flag_varz, 
                  namconv = NamingConvention(prefix)
    )
    
    if (ifac > 0) {
      if(flag_est){
        nm_est = paste(prefix, nm_var, "estim", sep = ".")
        dbout[nm_est] = dbout[nm_est]/dbout[nm_fac[ifac]]
      }
      if (flag_std){
        nm_std = paste(prefix, nm_var, "stdev", sep = ".")
        dbout[nm_std] = dbout[nm_std]/abs(dbout[nm_fac[ifac]])
      }
      if (flag_varz){
        nm_varz = paste(prefix, nm_var, "varz", sep = ".")
        dbout[nm_varz] = dbout[nm_varz]/abs(dbout[nm_fac[ifac]])
      }
    }
  }
}

# QC of the results
epsilon = 1.e-3
if (length(nms)) {
  for (nm in nms) {
    for (cs in c("estim", "stdev", "varz")[c(flag_est, flag_std, flag_varz)]) {
      nm_1 = paste("SVC", nm, cs, sep = ".")
      nm_2 = paste(pr, nm, cs, sep = ".")
      if(all(abs(dbout$getColumn(name = nm_1, useSel = TRUE) - dbout$getColumn(name = nm_2, useSel = TRUE)) < epsilon)) {
        res = paste0("Test passed: ", "SVC == ", pr)
      } else {
        res = paste0("Test failed: ", "SVC != ", pr)
        print(dbout$getColumn(name = nm_1, useSel = TRUE))
        print(dbout$getColumn(name = nm_2, useSel = TRUE))
      }
      print(paste0(res, "(nvar = ", nvar, " nfac  = ", nfac,"):", nm, "-", cs))
    } # loop over estim/stdev
  } # loop over the estimation
}

opers = EStatOption_fromKeys(c("NUM", "MINI", "MAXI", "MEAN", "STDV"))
knitr::kable(dbStatisticsMono(db = dbout, names = paste0(c(pr, "SVC"),".*"), 
                              opers = opers)$toTL(), digits = 4, caption = "Statistics")
# SVC in general

# defining the model of the effects
P   = 1
L   = 2

sigma  = matrix(runif(L^2), nrow = L, ncol = L)
sigma  = sigma %*% t(sigma)
eigen(sigma)$values

param  = c(1/2, 1, 3/2)[1:L]
ranges = c(0.25, 0.5, 1.0)[1:L]

model = Model_createFromEnvironment(nvar = L, ndim = ndim)
err   = model$addCovFromParam(type = ECov_MATERN(), 
                              sills = sigma,
                              param = param[1], range = ranges[1])
err   = model$setDriftIRF(order = 0, nfex = P)

if (verbose) {
  err = model$display()
}

# input/output data bases
ndat   = 10
seedin = 13227
nout    = 3
seedout = 134484
dbin   = Db_createFillRandom(ndat = ndat, ndim = ndim, 
                             nvar = 1, nfex = P, ncode = 0, 
                             varmax = 0., selRatio = p_sel, heteroRatio = VectorDouble(), 
                             coormin = VectorDouble(), coormax = VectorDouble(), seed = seedin)
dbout   = Db_createFillRandom(ndat = nout, ndim = ndim, 
                              nvar = 0, nfex = P, ncode = 0, 
                              varmax = 0., selRatio = p_sel, heteroRatio = VectorDouble(), 
                              coormin = VectorDouble(), coormax = VectorDouble(), seed = seedout)

nm_var = dbin$getNamesByLocator(locatorType = ELoc_Z())
nm_fac = dbin$getNamesByLocator(locatorType = ELoc_F())
nfac = length(nm_fac)

# neighborhood
# neigh = NeighUnique_create()
neigh = NeighMoving_create(nmini = 7, nmaxi = 7)

# list of estimation to be performed
estim_list = list()
estim_list[[1+length(estim_list)]] <- list(type = 1, idx = -1)
nms = NULL

for (i in seq_along(estim_list)) {
  if (estim_list[[i]]$type %in% 1:3) {
    nms = c(nms, paste(sep = ".", ESTIM_KEY[estim_list[[i]]$type], nm_var))
  } else if (estim_list[[i]]$type %in% 4) {
    nms = c(nms, paste(sep = ".", ESTIM_KEY[estim_list[[i]]$type], estim_list[[i]]$idx, nm_var))
  }
}
# ----------------------------------------------------------------------------
# Using the test function
# ----------------------------------------------------------------------------
verbose = TRUE
flag_centeredFactors = FALSE
err = kriging_SVC(dbin = dbin, dbout = dbout, model = model, neigh = neigh,
                flag_est = flag_est, flag_std = flag_std, flag_varz = flag_varz,
                estim_list = estim_list, flag_centeredFactors = flag_centeredFactors,
                prefix = "SVC", verbose = verbose
)
# ----------------------------------------------------------------------------
# statistics
# ----------------------------------------------------------------------------
knitr::kable(dbStatisticsMono(db = dbout, names = paste0(c("SVC"),".*"), 
                              opers = opers)$toTL(), digits = 4, caption = "Statistics")


