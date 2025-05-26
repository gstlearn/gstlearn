################################################################################
#                                                                              #
#                            gstlearn R package                                #
#                                                                              #
# Copyright (c) (2023) MINES Paris / ARMINES                                   #
# Authors: gstlearn Team                                                       #
# Website: https://gstlearn.org                                                #
# License: GPL v3                                                              #
#                                                                              #
################################################################################
#
# This is a set of functions used by SVC implementation (from X. Freulon)
#

#title: "Regression with spatially varying coefficients in gstlearn"
#author: "Xavier Freulon, Didier Renard, and the gstlearn Team"
#date: "2025-05-15"
# -------------------------------------------------
#' Function to compute the kriging using Spatially Varying Coefficient process (SVC)
#' -------------------------------------------------
#' @param dbin  A Db containing the observations and the external drifts
#' @param dbout A Db containing the external drifts at the target locations
#' @param model A Model object defining the covariance function of the residual (mono-variate case)
#' @param neigh A Neigh object defining the neighbourhood to be used (UNIQUE or MOVING)
#' Computation of the estimation and the standard deviation of the estimation error
#' @param flag_est  A Boolean scalar. If TRUE, the estimate is computed
#' @param flag_std  A Boolean scalar. If TRUE, the standard deviation of the kriging error is computed
#' @param flag_varz A Boolean scalar. If TRUE, the variance of the kriging value is computed
#' Computation of the residual, the drift or the drift coefficients
#' @param estim_list  A list of variables to be estimated with
#' type %in% 1:4 with 
#'  1 - "observed value"
#'  2 - "residual"
#'  3 - "drift"
#'  4 - "coefficient"
#' idx  %in% c(NULL, 0:(nfac+1))
#' @param flag_centeredFactors A Boolean scalar
#'  If TRUE, the spatially varying factors are centered
#' @param prefix A string to build the name of the variables storing the results
#' @param verbose
#' @value a error code (not implemented) 
#' Notes: 
#' - error variance (ELoc_V()) not implemented
kriging_SVC <- function(dbin, dbout, model, neigh, 
                      flag_est = TRUE,  flag_std = FALSE, flag_varz = FALSE,
                      estim_list = list(list(type = 1, idx = NULL), 
                                        list(type = 2, idx = NULL),
                                        list(type = 3, idx = NULL)),
                      flag_centeredFactors = FALSE,
                      prefix = "SVC", verbose = FALSE) {
  #-----------------------------------------------------------------------
  # the observed variable
  #-----------------------------------------------------------------------
  nm_obs = dbin$getNamesByLocator(ELoc_Z())
  stopifnot(length(nm_obs) == 1)
  
  # the factors
  nm_fac = dbin$getNamesByLocator(ELoc_F())
  nfac = 1+length(nm_fac)       # the constant is added (ifac in 0:nfac)
  nvar = model$getNVar()
  stopifnot(nvar <= nfac)
  
  #-----------------------------------------------------------------------
  # Targets  
  #-----------------------------------------------------------------------
  iech_out = dbout$getColumns(names = "rank", useSel = TRUE)
  nout     = length(iech_out)
  
  #-----------------------------------------------------------------------
  # Estimations
  #-----------------------------------------------------------------------
  ESTIM_KEY = c("obs", "res", "drift", "coeff")
  n_estim   = length(estim_list)
  estim_nm  = NULL
  for (i in 1:n_estim) {
    stopifnot(estim_list[[i]]$type %in% 1:4)
    if (estim_list[[i]]$type %in% 1:3) {
      estim_nm = c(estim_nm, paste(sep = ".", ESTIM_KEY[estim_list[[i]]$type], nm_obs))
    } else if (estim_list[[i]]$type %in% 4) {
      estim_nm = c(estim_nm, paste(sep = ".", ESTIM_KEY[estim_list[[i]]$type], estim_list[[i]]$idx, nm_obs))
    }
  }
  if (flag_est) {
    res_est = matrix(NaN, nrow = dbout$getNSample(useSel = FALSE), ncol = n_estim)
  }
  if (flag_std) {
    res_std = matrix(NaN, nrow = dbout$getNSample(useSel = FALSE), ncol = n_estim)
  }
  if (flag_varz) {
    res_varz = matrix(NaN, nrow = dbout$getNSample(useSel = FALSE), ncol = n_estim)
  }
  
  if(verbose) {
    print(paste0("SVC: number of spatial effects     = ", nvar))
    print(paste0("SVC: number of constant effects    = ", nfac - nvar))
    if (flag_centeredFactors) {
      print(paste0("SVC: the spatially varying factors are centered (one universal condition)"))
    } else {
      print(paste0("SVC: the spatially varying factors are not centered (", nfac, " universal conditions)"))
    }
    print(paste0("SVC: number of estimated variables = ", n_estim * sum(c(flag_est, flag_std, flag_varz))))
    for (i_estim in 1:n_estim) {
      print(paste0("SVC: variable #", i_estim, " = ", estim_nm[i_estim]))
    }
  }
  
  # --------------------------
  # Initialization 
  # --------------------------
  model_mono = model$createReduce(0) # The SVC model is mono-variate (used to compute the drift matrix)
  err      = neigh$attach(dbin, dbout)
  nbghObs  = VectorInt()
  rankObs  = VectorVectorInt()
  
  Kcalc    = KrigingAlgebra()
  krigopt  = KrigOpt()
  
  # working matrices
  Sigma    = MatrixSymmetric() # covariance matrix of the observations
  matcov   = MatrixSymmetric() # working matrix for the covariance matrix evaluation
  Sigma0   = MatrixDense()     # covariance matrix Data x Target
  X        = MatrixDense()     # Drift Matrix of the data
  X0       = MatrixDense()     # Drift matrix of the target
  Sigma00  = model$eval0Mat()  # Covariance of the target
  
  # --------------------------
  # Loop on the target sites
  # --------------------------
  for (iout in 1:nout) {
    idxTarget = iech_out[iout]-1
    # ----------------------
    # computing LHS (iout == 1 or moving neighborhood)
    # ----------------------
    if ((iout == 1)|(neigh$getType()$getDescr() != "Unique Neighborhood")) {
      # data selection (ns observed values)
      err     = neigh$getNeigh(iech_out = idxTarget, ranks = nbghObs)
      rankObs = dbin$getSampleRanks(nbgh = nbghObs)
      Z       = dbin$getValuesByRanks(rankObs)
      ns      = length(Z)
      # computing LHS
      # Sigma is a matrix [ns x ns] and X is a matrix [ns x nfac]
      
      # The SVC model is mono-variate to compute the drift matrix
      err = model_mono$evalDriftMatInPlace(mat = X, db = dbin, nbgh = nbghObs, member = ECalcMember_fromKey("LHS"))
      
      # The SVC model is defined by a multi-variate to compute the covariance matrix
      err = dbin$clearLocators(locatorType = ELoc_Z())
      err = model$evalCovMatSymInPlace(matcov, db1 = dbin, nbgh1 = nbghObs, cleanOptim = FALSE)
      err = dbin$setLocators(nm_obs, locatorType = ELoc_Z(), cleanSameLocator = TRUE)
      
      err = computeCovMatSVCLHSInPlace(cov = Sigma, Sigma = matcov, F1 = X)
      if(flag_centeredFactors) {
        X_drift = MatrixDense(nrow = ns, ncol = 1)
        err = X_drift$setColumn(icol = 0, tab = rep(1.0, ns))
      } else {
        X_drift = X
      }
      err     = Kcalc$resetNewData()
      err     = Kcalc$setData(Z = Z, indices = rankObs)
      err     = Kcalc$setLHS(Sigma, X_drift)
    }
    
    # ----------------------
    # computing RHS (start)
    # ----------------------
    # RHS variables
    RHS_Sigma00 = MatrixSymmetric()
    RHS_Sigma0  = MatrixDense()
    RHS_X0 = MatrixDense()
    
    # drifts at target and covariance matrix Data-Target
    
    # The SVC model is mono-variate to compute the drift matrix
    err = model_mono$evalDriftMatInPlace(mat = X0, db = dbout, nbgh = idxTarget, member = ECalcMember_fromKey("LHS"))
    
    # The SVC model is defined by a multi-variate to compute the covariance matrix
    err = dbin$clearLocators(locatorType = ELoc_Z())
    err = model$evalCovMatInPlace(mat = Sigma0, 
                                  db1 = dbin,  nbgh1 = nbghObs,
                                  db2 = dbout, nbgh2 = idxTarget, 
                                  cleanOptim = FALSE)
    err = dbin$setLocators(nm_obs, locatorType = ELoc_Z(), cleanSameLocator = TRUE)
    
    for (i in 1:n_estim) {
      # variance of the target variable
      err = computeCovMatSVCRHSInPlace(cov = RHS_Sigma00, Sigma = Sigma00, F1 = X0, F2 = X0,
                                        type1 = estim_list[[i]]$type, idx1 = estim_list[[i]]$idx, 
                                        type2 = estim_list[[i]]$type, idx2 = estim_list[[i]]$idx
      )
      
      # covariance between data and target variable
      err = computeCovMatSVCRHSInPlace(cov = RHS_Sigma0, Sigma = Sigma0, F1 = X, F2 = X0,
                                        type1 = 1, idx1 = -1, 
                                        type2 = estim_list[[i]]$type, idx2 = estim_list[[i]]$idx
      )
      
      # drift between data and target variable
      # TODO: resizing the output matrix to be inserted in the function
      if (flag_centeredFactors) {
        err = RHS_X0$resize(nrows = 1, ncols = 1)
      } else {
        err = RHS_X0$resize(nrows = 1, ncols = X0$getNCols())
      }
      err = computeDriftMatSVCRHSInPlace(mat = RHS_X0, X0, 
                                          type = estim_list[[i]]$type, idx = estim_list[[i]]$idx,
                                          flagCenteredFactors = flag_centeredFactors)
      
      err = Kcalc$setVariance(RHS_Sigma00)
      err = Kcalc$setRHS(Sigma0 = RHS_Sigma0, X0 = RHS_X0)
      if(flag_est){
        res_est[iech_out[iout],i] = Kcalc$getEstimation()
      }
      if (flag_std) {
        res_std[iech_out[iout],i] = Kcalc$getStdv()
      }
      if (flag_varz) {
        res_varz[iech_out[iout],i] = Kcalc$getVarianceZstar()
      }
    } # loop over the estimations
  } # loop over the targets...
  
  # --------------------------
  # Storing the computed variables
  # --------------------------
  for (ivar in seq_along(estim_nm)) {
    if (flag_est) {
      nm  = paste(prefix, estim_nm[ivar], "estim", sep = ".")
      err = dbout$setColumn(tab = res_est[,ivar], name = nm, useSel = FALSE)
    }
    if (flag_std) {
      nm  = paste(prefix, estim_nm[ivar], "stdev", sep = ".")
      err = dbout$setColumn(tab = res_std[,ivar], name = nm, useSel = FALSE)
    }
    if (flag_varz) {
      nm  = paste(prefix, estim_nm[ivar], "varz", sep = ".")
      err = dbout$setColumn(tab = res_varz[,ivar], name = nm, useSel = FALSE)
    }
  }
  err
}

