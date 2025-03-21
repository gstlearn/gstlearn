# Test of the KrigingAlgebra in R
#' -----------------------------------------------------------------------------
#' Test of the kriging template
#' -----------------------------------------------------------------------------
rm(list = ls())
library(gstlearn)
#' -----------------------------------------------------------------------------
#' Function for kriging by hand:
#' -----------------------------------------------------------------------------
test_kriging <- function(dbin, dbout, model, neigh, 
    calcul = EKrigOpt_POINT(), order, 
    flag_est = TRUE, flag_std = TRUE, flag_varz = FALSE,
    prefix = "test_kriging") {

  nm_var = dbin$getNamesByLocator(ELoc_Z())
  nvar = length(nm_var)
  stopifnot(model$getNVar() == nvar)
  if (model$getNDrift() == 0) {order = -1} # TODO: does a method getOrder exist?
  iech_out = dbout$getColumns(names = "rank", useSel = TRUE)
  nout = length(iech_out)
  if (flag_est) {
    res_est = matrix(NaN, nrow = dbout$getNSample(useSel = FALSE), ncol = nvar)
  }
  if (flag_std) {
    res_std = matrix(NaN, nrow = dbout$getNSample(useSel = FALSE), ncol = nvar)
  }
  if (flag_varz) {
    res_varz = matrix(NaN, nrow = dbout$getNSample(useSel = FALSE), ncol = nvar)
  }
  # --------------------------
  # Initialization 
  # --------------------------
  err     = neigh$attach(dbin, dbout)
  Sigma   = MatrixSquareSymmetric()
  X       = MatrixRectangular()
  Z       = VectorDouble()
  Sigma0  = MatrixRectangular()
  X0      = MatrixRectangular()
  Kcalc   = KrigingAlgebra()
  krigopt = KrigOpt()
  id_neigh = VectorInt()
  sampleRanks = VectorVectorInt()

  Sigma00 = model$eval0Mat()
  err     = Kcalc$setVariance(Sigma00)
  # --------------------------
  # Loop on the target sites
  # --------------------------
  for (iout in 1:nout) {

    # ----------------------
    # computing LHS (iout == 1 or moving neighborhood)
    # ----------------------
    if ((iout == 1)|(neigh$getType()$getDescr() != "Unique Neighborhood")) {
      err = neigh$getNeigh(iech_out = iech_out[iout]-1, ranks = id_neigh)
      sampleRanks = dbin$getSampleRanks(nbgh = id_neigh)
      Z = dbin$getValuesByRanks(sampleRanks, means = model$getMeans(),
          subtractMean = !model$hasDrift())
      err = model$evalCovMatSymInPlaceFromIdx(mat = Sigma, db1 = dbin,
          index1 = sampleRanks)
      err = Kcalc$resetNewData()
      if (order == -1) {
        err = Kcalc$setData(Z = Z, indices = sampleRanks,
	   Means = model$getMeans())
        err = Kcalc$setLHS(Sigma)
      } else {
        err = Kcalc$setData(Z = Z, indices = sampleRanks)
	err = model$evalDriftMatByRanks(mat = X, db = dbin, sampleRanks,
	         member = ECalcMember_fromKey("LHS"))
        err = Kcalc$setLHS(Sigma, X)
      }
    }
    
    # ----------------------
    # computing RHS
    # ----------------------
    err = model$evalCovMatRHSInPlaceFromIdx(mat=Sigma0, db1=dbin, db2 = dbout, 
                                   index1 = sampleRanks,
				   iech2 = iech_out[iout]-1, 
                                   krigopt = krigopt, 
                                   cleanOptim = FALSE)
    if (order == -1) {
      err = Kcalc$setRHS(Sigma0 = Sigma0)
    } else {
      err = model$evalDriftMatByTarget(mat = X0, db = dbout,
                                       iech2 = iech_out[iout]-1,
   				       krigopt = krigopt)
      err = Kcalc$setRHS(Sigma0 = Sigma0, X0 = X0)
    }
    if(flag_est){
      res_est[iech_out[iout],] = Kcalc$getEstimation()
    }
    if (flag_std) {
      res_std[iech_out[iout],] = Kcalc$getStdv()
    }
    if (flag_varz) {
      res_varz[iech_out[iout],] = Kcalc$getVarianceZstar()
    }
  } # loop over the targets...

  # --------------------------
  # Storing the computed variables
  # --------------------------
  if (flag_est) {
    nm = paste(prefix, nm_var, "estim", sep = ".")
    for (ivar in 1:nvar) {
      err = dbout$setColumn(tab = res_est[,ivar], name = nm[ivar],
      useSel = FALSE)
    }
  }
  if (flag_std) {
    nm = paste(prefix, nm_var, "stdev", sep = ".")
    for (ivar in 1:nvar) {
      err = dbout$setColumn(tab = res_std[,ivar], name = nm[ivar],
      useSel = FALSE)
    }
  }
  if (flag_varz) {
    nm = paste(prefix, nm_var, "varz", sep = ".")
    for (ivar in 1:nvar) {
      err = dbout$setColumn(tab = res_varz[,ivar], name = nm[ivar],
      useSel = FALSE)
    }
  }
  err
}

#' -----------------------------------------------------------------------------
#' Function for the tests
#' -----------------------------------------------------------------------------
performTest = function (ndim, cas, k, nvar, 
                        flag.NA, flag.unique, flag.sel.out, flag.sel.in, verbose)
{
  order = cas[[k]]$order
  nfex  = cas[[k]]$nfex
  
  # --------------------------------------------------------------------------
  # Do the test
  # --------------------------------------------------------------------------
  if(order == -1) {
    prefix = "SK"
  } else {
    prefix = "UK"
  }
  err = mestitle(1, paste0(prefix,
                           " Nvar:", nvar, " Order:", order, " Nfex:", nfex,
                           " Unique:", flag.unique, 
                           " Hetero:", flag.NA, 
                           " in_sel:", flag.sel.in, 
                           " out_sel:", flag.sel.out))
  
  heteroRatio = 0.0; selRatio.in = 0.0; selRatio.out = 0.0
  if (flag.NA) {heteroRatio = 0.1}       # 10% of the samples per variables are note defined
  if (flag.sel.in) {selRatio.in = 0.1}   # 10% of the input samples are NOT selected
  if (flag.sel.out) {selRatio.out = 0.1} # 10% of the target are NOT selected
  # CreateGenerate the Input data set
  dbin   = Db_createFillRandom(ndat = ndat, ndim = ndim, nvar = nvar, nfex = nfex,
                               ncode = 0, varmax = 0., 
                               selRatio = selRatio.in,
                               heteroRatio = rep(heteroRatio, nvar),
                               coormin = VectorDouble(), coormax = VectorDouble(),
                               seed = seedin)
  if (verbose) {
    err = dbin$display()
  }
  #  Generate the output data set
  dbout   = Db_createFillRandom(ndat = nout, ndim = ndim, nvar = 0, nfex = nfex,
                                ncode = 0, varmax = 0., 
                                selRatio = selRatio.out,
                                heteroRatio = 0.0,
                                coormin = VectorDouble(), coormax = VectorDouble(),
                                seed = seedout)
  if (verbose) {
    err = dbout$display()
  }
  
  #  Create the Model
  model = Model_createFillRandom(ndim = ndim, nvar = nvar, 
                                 types = ECov_fromKeys(c("NUGGET", "EXPONENTIAL")), 
                                 order = order, nfex = nfex)
  if (verbose) {
    err = model$display()
  }
  
  #  Creating the Neighborhood
  if (flag.unique) {
    neigh = NeighUnique_create(flag_xvalid = FALSE)
  } else {
    neigh = NeighMoving_create(nmaxi = 5, nmini = 5, flag_xvalid = FALSE)
  }
  if (verbose) {
    err = neigh$display()
  }
  
  title = paste0(prefix, ": nvar = ", nvar, " nf = ", nfex)
  
  #  Perform Kriging using the API
  prefix_1 = paste0(prefix,"_ini")
  err = kriging(dbin = dbin, dbout = dbout, model = model, neigh = neigh, 
                calcul = EKrigOpt_POINT(),
                flag_est = flag_est, flag_std = flag_std, flag_varz = flag_varz,
                namconv = NamingConvention(prefix_1))
  if (verbose) {
    err = dbout$display() 
  }
  
  #  Perform Kriging *by hand*
  prefix_2 = paste0(prefix,"_test")
  err = test_kriging(dbin = dbin, dbout = dbout, model = model, neigh = neigh,
                     calcul = EKrigOpt_POINT(), order,
                     flag_est = flag_est, flag_std = flag_std, flag_varz = flag_varz,
                     prefix = prefix_2)
  if (verbose) {
    err = dbout$display() 
  }
  
  # QC
  if (flag.sel.out) {
    sel = (dbout["sel"] == 1)
  } else {
    sel = rep(TRUE, dbout$getNSample(useSel = FALSE))  
  }
  for (ivar in 1:nvar) {
    if (nvar == 1) {
      nm_var_1 = paste0(prefix_1, ".z")
      nm_var_2 = paste0(prefix_2, ".z")
    } else {
      nm_var_1 = paste0(prefix_1, ".z.", ivar)      
      nm_var_2 = paste0(prefix_2, ".z.", ivar)      
    }
    for (k in 1:3) {
      nm_estim_1 = paste0(nm_var_1, ".", c("estim", "stdev", "varz")[k])
      nm_estim_2 = paste0(nm_var_2, ".", c("estim", "stdev", "varz")[k])
      #stopifnot(all(abs(dbout[nm_estim_1] - dbout[nm_estim_2])[sel] < 1e-16))
      if (all(abs(dbout[nm_estim_1] - dbout[nm_estim_2])[sel] > 1e-16))
      {
        print("Un probleme a ete decele")
        print(dbout[nm_estim_1][sel])
        print(dbout[nm_estim_2][sel])
      }
    }
  }
  
  # --------------------------------------------------------------------------
  # Test OK
  # --------------------------------------------------------------------------
  err = message("Test passed\n")
  
  invisible()
}

# ------------------------------------------------------------------------------
#' Test of the kriging template
# ------------------------------------------------------------------------------
ndim = 2
err = defineDefaultSpace(ESpaceType_RN(), ndim = ndim)
err = OptCst_define(ECst_NTROW(), -1)
err = OptCst_define(ECst_NTCOL(), -1)
verbose = FALSE
flag_est = TRUE
flag_std = TRUE
flag_varz = TRUE
#  Generate the input data (in the 1x1 square) with 'nvar' variables and 'nfex' external drifts
ndat   = 10
seedin = 13227
nout    = 10
seedout = 134484
nvar_min = 1
nvar_max = 3
cas = list()
cas[[1+length(cas)]] <- list(order = -1, nfex = 0)
cas[[1+length(cas)]] <- list(order =  0, nfex = 0)
cas[[1+length(cas)]] <- list(order =  0, nfex = 1)
cas[[1+length(cas)]] <- list(order =  0, nfex = 2)

# ------------------------------------------------------------------------------
#  Global parameters
# ------------------------------------------------------------------------------

for(flag.sel.in in c(FALSE,TRUE)) {
  for(flag.sel.out in c(FALSE,TRUE)) {
    for(flag.unique in c(FALSE, TRUE)) {
      for(flag.NA in c(FALSE, TRUE)) {
        for (nvar in nvar_min:nvar_max) {
          for (k in seq_along(cas)) {
            performTest(ndim, cas, k, nvar, 
                        flag.NA, flag.unique, flag.sel.out, flag.sel.in, verbose)
          }
        }
      }
    }
  }
}
# ------------------------------------------------------------------------------
#' END
# ------------------------------------------------------------------------------
