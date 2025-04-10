/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Estimation/KrigingAlgebraSimpleCase.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"
#include "geoslib_old_f.h"

#include "Estimation/KrigingSystemSimpleCase.hpp"

#include "Enum/EKrigOpt.hpp"
#include "Enum/ECalcMember.hpp"
#include "Enum/ELoc.hpp"
#include "Enum/ENeigh.hpp"

#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Model/ModelGeneric.hpp"
#include "Neigh/ANeigh.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighImage.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/VectorHelper.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Space/SpaceRN.hpp"
#include "Estimation/KrigOpt.hpp"

#include <math.h>

KrigingSystemSimpleCase::KrigingSystemSimpleCase(Db* dbin,
                                                 Db* dbout,
                                                 const ModelGeneric* model,
                                                 ANeigh* neigh)
  : _dbin(dbin)
  , _dbout(dbout)
  , _model(nullptr)
  , _neigh(neigh)
  , _isReady(false)
  , _algebra()
  , _krigopt()
  , _means()
  , _meansTarget()
  , _iptrEst(-1)
  , _iptrStd(-1)
  , _iptrVarZ(-1)
  , _flagEst(false)
  , _flagStd(false)
  , _flagVarZ(false)
  , _flagDataChanged(false)
  , _iptrWeights(-1)
  , _flagWeights(false)
  , _flagSet(true)
  , _flagLTerm(false)
  , _neighUnique(false)
  , _iptrNeigh(-1)
  , _ndim(0)
  , _nvar(0)
  , _nech(0)
  , _nfeq(0)
  , _neq(0)
  , _dbinUidToBeDeleted()
  , _dboutUidToBeDeleted()
  , _flagVerr(false)
  , _flagNoStat(false)
{
  // _model is a copy of input model to allow modification (still used???)
  if (model != nullptr)
  {
    _model = (ModelGeneric*)model->clone();
  }
  if (model != nullptr)
    _flagNoStat = _model->isNoStat();

  // Reset the neighborhood
  if (neigh != nullptr)
    neigh->reset();

  _neighUnique = dynamic_cast<NeighUnique*>(_neigh) != nullptr;
  _algebra.setNeighUnique(_neighUnique);

  // Define local constants
  _flagVerr = _dbin->hasLocVariable(ELoc::V);

  _resetMemoryGeneral();
}

KrigingSystemSimpleCase::~KrigingSystemSimpleCase()
{
  // Turn OFF this option for future task
  OptDbg::setCurrentIndex(-1);

  // Clean elements from _dbin

  if (_dbin != nullptr)
  {
    if (!_dbinUidToBeDeleted.empty())
    {
      (void)_dbin->deleteColumnsByUID(_dbinUidToBeDeleted);
    }
  }

  // Clean elements from _dbout

  if (_dbout != nullptr)
  {
    if (!_dboutUidToBeDeleted.empty())
    {
      (void)_dbout->deleteColumnsByUID(_dboutUidToBeDeleted);
    }
  }
}

int KrigingSystemSimpleCase::_getNVar() const
{
  int nvar = 0;
  if (_model != nullptr)
  {
    if (nvar > 0 && nvar != _model->getNVar())
    {
      messerr("Inconsistent number of Variables - Value is returned as 0");
      return 0;
    }
    nvar = _model->getNVar();
  }

  // In the case of factor kriging, the number of Z-variables in the Data file
  // does not give the number of variables. Check should be avoided

  return nvar;
}

int KrigingSystemSimpleCase::_getNbfl() const
{
  if (_model == nullptr) return 0;
  return _model->getNDrift();
}

int KrigingSystemSimpleCase::_getNFeq() const
{
  if (_model == nullptr) return 0;
  return _model->getNDriftEquation();
}

int KrigingSystemSimpleCase::_getNeq() const
{
  int neq = _nvar * _nech + _nfeq;
  return neq;
}

void KrigingSystemSimpleCase::_resetMemoryGeneral()
{
  _setInternalShortCutVariablesGeneral();
}

void KrigingSystemSimpleCase::_dumpOptions()
{
  /* Kriging option */
  message("Punctual Estimation\n");
  message("\n");
}

void KrigingSystemSimpleCase::_setInternalShortCutVariablesModel()
{
  _nvar = _getNVar();
  _nfeq = _getNFeq();
  _neq  = _getNeq(); // reset as it depends on nech and Model
}
/**
 * Assign the values to local variables used as shortcuts
 * @return 1 if the number of active sample is zero
 */
int KrigingSystemSimpleCase::_setInternalShortCutVariablesNeigh(VectorInt& nbgh)
{
  _nech = nbgh.size();
  _neq  = _getNeq();
  return (_nech <= 0);
}
void KrigingSystemSimpleCase::_setInternalShortCutVariablesGeneral()
{
  _ndim = getNDim();
  _setInternalShortCutVariablesModel();
}
void KrigingSystemSimpleCase::_rhsDump(KrigingAlgebraSimpleCase& algebra)
{
  mestitle(0, "RHS of Kriging matrix");
  if (_nech > 0) message("Number of active samples    = %d\n", _nech);
  message("Total number of equations   = %d\n", _neq);
  message("Number of right-hand sides  = %d\n", 1);
  _dumpOptions();
  algebra.dumpRHS();
}

void KrigingSystemSimpleCase::_wgtDump(KrigingAlgebraSimpleCase& algebra)
{
  /* Header */
  mestitle(0, "(Co-) Kriging weights");
  algebra.dumpWGT();

  /* Auxiliary results for Drift */
  mestitle(0, "Drift or Mean Information");
  algebra.dumpAux();
}

/****************************************************************************/
/*!
 **  Calculate the final estimation and storage
 **
 ** \param[in] status   Kriging error status
 **
 *****************************************************************************/
void KrigingSystemSimpleCase::_estimateCalcul(int status, KrigingAlgebraSimpleCase& algebra, int iechout) const
{

  if (_flagEst)
    _estimateEstim(status, algebra, iechout);

  /* Variance of the estimation error */

  if (_flagStd)
    _estimateStdv(status, iechout, algebra);

  /* Variance of the estimator */

  if (_flagVarZ != 0)
    _estimateVarZ(status, iechout, algebra);

  /* Kriging weights (stored in _dbin) */

  if (_flagWeights != 0)
  {
    for (int ivarCL = 0; ivarCL < 1; ivarCL++)
    {
      for (int jech = 0; jech < _nech; jech++)
      {
        if (status != 0) continue;
        double wgt = algebra.getLambda()->getValue(jech, ivarCL);
        int iech   = algebra.getSampleRanks()->at(0)[jech];
        if (_flagSet)
          _dbin->setArray(iech, _iptrWeights + ivarCL, wgt);
        else
          _dbin->updArray(iech, _iptrWeights, EOperator::ADD, wgt);
      }
    }
  }
}

void KrigingSystemSimpleCase::_neighCalcul(int status, const VectorDouble& tab, int iechout)
{
  int ntab = (int)tab.size();
  for (int i = 0; i < ntab; i++)
  {

    /* Store the parameter */

    double value = (status == 0) ? tab[i] : TEST;
    _dbout->setArray(iechout, _iptrNeigh + i, value);
  }

  if (OptDbg::query(EDbg::NBGH) && status == 0)
  {
    mestitle(0, "Neighborhood Parameters");

    message("Number of selected samples          = %d\n", (int)tab[0]);
    message("Maximum neighborhood distance       = %lf\n", tab[1]);
    message("Minimum neighborhood distance       = %lf\n", tab[2]);
    message("Number of non-empty sectors         = %d\n", (int)tab[3]);
    message("Number of consecutive empty sectors = %d\n", (int)tab[4]);
  }
}

/****************************************************************************/
/*!
 **  Establish the calculation of estimation
 **
 ** \param[in]  status  Kriging error code
 **
 *****************************************************************************/
void KrigingSystemSimpleCase::_estimateEstim(int status, KrigingAlgebraSimpleCase& algebra, int iechout) const
{
  VectorDouble& local = algebra.getEstimation();
  if (local.size() <= 0) return;
  if (status) local.fill(TEST);
  if (FFFF(local[0]))
  {
    messerr("KrigingSystemSimpleCase::_estimateEstim: local[0] is FFFF");
    return;
  }
  for (int ivarCL = 0; ivarCL < 1; ivarCL++)
    _dbout->setArray(iechout, _iptrEst + ivarCL, local[ivarCL]);
}

/****************************************************************************/
/*!
 **  Establish the calculation of standard deviation
 **
 ** \param[in]  status  Kriging error code
 **
 *****************************************************************************/
void KrigingSystemSimpleCase::_estimateStdv(int status, int iechout, KrigingAlgebraSimpleCase& algebra) const
{
  VectorDouble local = algebra.getStdv();
  if (local.size() <= 0) return;
  if (status) local.fill(TEST);
  for (int ivarCL = 0; ivarCL < 1; ivarCL++)
    _dbout->setArray(iechout, _iptrStd + ivarCL, local[ivarCL]);
}

/****************************************************************************/
/*!
 **  Establish the variance of the estimator
 **
 ** \param[in]  status  Kriging error code
 **
 *****************************************************************************/
void KrigingSystemSimpleCase::_estimateVarZ(int status, int iechout, KrigingAlgebraSimpleCase& algebra) const
{
  DECLARE_UNUSED(status)
  for (int ivarCL = 0; ivarCL < 1; ivarCL++)
    _dbout->setArray(iechout, _iptrVarZ + ivarCL, algebra.getVarianceZstar(ivarCL));
}


void KrigingSystemSimpleCase::updateLHS(KrigingAlgebraSimpleCase& algebra, ModelGeneric& model)
{

  auto* Sigma       = algebra.getSigma();
  auto* X           = algebra.getX();
  auto* Z           = algebra.getZ();
  auto* sampleRanks = algebra.getSampleRanks();
  auto* nbgh        = algebra.getNbgh();
  _dbin->getSampleRanksInPlace(sampleRanks, VectorInt(), *nbgh);
  _dbin->getValuesByRanksInPlace(Z, *sampleRanks, _means, !model.hasDrift());
  _setInternalShortCutVariablesNeigh(*nbgh);
  if (model.evalCovMatSymInPlaceFromIdx(*Sigma, _dbin, *sampleRanks, nullptr, false)) return;
  if (model.evalDriftMatByRanks(*X, _dbin, *sampleRanks, ECalcMember::LHS)) return;
  algebra.updateSampleRanks();
  // if (algebra.setData(Z, sampleRanks, _meansTarget)) return;
}
/**
 * Performs the last operations before launching the loop on Estimations
 * @return
 */
bool KrigingSystemSimpleCase::isReady()
{
  if (!_isCorrect()) return false;

  // Attach the Input and Output Db
  _neigh->attach(_dbin, _dbout);

  // Define the means of each variable
  _means = _model->getMeans();
  _algebra.setMeans(_means);
  // Possible adjust the means in case of presence of 'matLC'
  _meansTarget = _means;

  _model->optimizationPreProcessForData(_dbin);
  _model->manage(_dbin, _dbout);
  if (_neighUnique)
  {
    auto* nbgh = _algebra.getNbgh();
    _neigh->select(0, *nbgh);
    updateLHS(_algebra, *_model);
    _algebra.prepare();
  }

  // Perform some pre-calculation when variance of estimator is requested
  if (_flagStd)
  {
    auto* Sigma00 = _algebra.getSigma00();
    if (_model->evalCovMat0InPlace(*Sigma00, _dbout, 0, _krigopt)) return false;
  }

  _isReady = true;


  return _isReady;
}

/**
 * This method closes the use of a KrigingSystemSimpleCase sequence
 */
void KrigingSystemSimpleCase::conclusion()
{
  _model->optimizationPostProcess();
}

int KrigingSystemSimpleCase::estimate(int iechout,
                                      SpacePoint& pin,
                                      SpacePoint& pout,
                                      VectorDouble& tabwork,
                                      KrigingAlgebraSimpleCase& algebra,
                                      ModelGeneric& model,
                                      ANeigh* neigh)
{

  if (!_dbout->isActive(iechout)) return 0;
  if (!_isReady)
  {
    messerr("You must call 'isReady' before launching 'estimate'");
    return 1;
  }

  // In case of Image Neighborhood, the neighboring samples have already
  // been selected in isReady(). No need to compute them again.

  //  // Store the Rank of the Target sample

  int status = 0;

  /* Establish the Kriging R.H.S. */
  MatrixRectangular* Sigma0 = algebra.getSigma0();
  MatrixRectangular* X0     = algebra.getX0();
  if (iechout + 1 == OptDbg::getReference())
  {
    OptDbg::setCurrentIndex(iechout + 1);
    if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
    {
      mestitle(1, "Target location");
      db_sample_print(_dbout, iechout, 1, 0, 0, 0);
    }
  }

  if (_neighUnique)
  {
    if (iechout == 0)
    {
      neigh->displayDebug(algebra.getSampleRanks()->at(0));
      _setInternalShortCutVariablesNeigh(*algebra.getNbgh());
      status = 0;
    }
  }
  else
  {
    auto* nbgh = algebra.getNbgh();
    neigh->select(iechout, *nbgh);
    if (nbgh->size() <= 0)
      return 0;
    if (!_neigh->isUnchanged() || _neigh->getFlagContinuous() || OptDbg::force())
    {
      updateLHS(algebra,model);
      status = _setInternalShortCutVariablesNeigh(*nbgh);
    }
  }

  model.evalCovVecRHSInPlace(Sigma0->getViewOnColumnModify(0),
                             _dbout,
                             algebra.getSampleRanks()->at(0),
                             iechout,
                             _krigopt,
                             pin,
                             pout,
                             tabwork);

  model.evalDriftMatByTarget(*X0, _dbout, iechout, _krigopt);
  algebra.updateRHS();

  // Printout for debugging case

   //if (!_neigh->isUnchanged() || _neigh->getFlagContinuous() || OptDbg::force())
  if (OptDbg::force())
  {
    // LHS is not printed systematically... only when it has been modified
    if (iechout + 1 == OptDbg::getReference())
      if (OptDbg::query(EDbg::KRIGING)) algebra.dumpLHS();
  }

  if (iechout + 1 == OptDbg::getReference())
    if (OptDbg::query(EDbg::KRIGING))
    {
      _rhsDump(algebra);
      _wgtDump(algebra);
    }

  /* Perform the final estimation */
  //  If status is not zero, cancel the current Neighborhood search status
  // if (status) _neigh->setIsChanged();

  // Store the results in the output Db

  _estimateCalcul(status, algebra, iechout);

  // Final printout
  if (iechout + 1 == OptDbg::getReference())
    if (OptDbg::query(EDbg::RESULTS))
    {
      _dumpKrigingResults(status, iechout, &algebra);
    }
  return 0;
}
/*****************************************************************....*********/
/*!
 **  Print the results
 **
 ** \param[in] status   Kriging error status
 **
 *****************************************************************************/
void KrigingSystemSimpleCase::_dumpKrigingResults(int status, int iechout, KrigingAlgebraSimpleCase* algebra) const
{

  mestitle(0, "(Co-) Kriging results");
  message("Target Sample = %d\n", iechout + 1);

  /* Loop on the results */

  for (int ivar = 0; ivar < 1; ivar++)
  {

    // Printout of the Estimation

    message("Variable Z%-2d\n", ivar + 1);
    if (_iptrEst >= 0)
    {
      double value = (status == 0) ? _dbout->getArray(iechout, _iptrEst + ivar) : TEST;
      tab_printg(" - Estimate  = ", value);
      message("\n");
    }
    if (_iptrStd >= 0)
    {
      double value = (status == 0) ? _dbout->getArray(iechout, _iptrStd + ivar) : TEST;
      tab_printg(" - Std. Dev. = ", value);
      message("\n");
      tab_printg(" - Variance  = ", FFFF(value) ? TEST : value * value);

      value = algebra->getSigma00()->getValue(ivar, ivar);
      message("\n");
      tab_printg(" - Cov(h=0)  = ", value);
      message("\n");
    }
    if (_iptrVarZ >= 0)
    {
      double value = (status == 0) ? _dbout->getArray(iechout, _iptrVarZ + ivar) : TEST;
      tab_printg(" - Var(Z*)   = ", value);
      message("\n");
    }
  }
}

/**
 * Set the calculation options
 * @param iptrEst  UID for storing the estimation(s)
 * @param iptrStd  UID for storing the Standard deviations(s)
 * @param iptrVarZ UID for storing the Variance(s) of estimator
 * @param forceNoDual Force that the algebra is not using the Dual option
 * @return Error returned code
 * @remark If a term must not be calculated, its UID must be negative
 */
int KrigingSystemSimpleCase::updKrigOptEstim(int iptrEst,
                                             int iptrStd,
                                             int iptrVarZ,
                                             bool forceNoDual)
{
  _iptrEst  = iptrEst;
  _iptrStd  = iptrStd;
  _iptrVarZ = iptrVarZ;

  _flagEst  = _iptrEst >= 0;
  _flagStd  = (_iptrStd >= 0);
  _flagVarZ = (_iptrVarZ >= 0);

  _flagDataChanged = true;

  // Set the Dual option automatically if:
  // - no Variance estimation(s) is requested
  // - the Unique Neighborhood is used
  // - the Dual option is not forced OFF
  if (!forceNoDual && !_flagStd && !_flagVarZ)
  {
    if (_neigh != nullptr && _neigh->getType() == ENeigh::UNIQUE)
      _algebra.setDual(true);
  }

  return 0;
}

int KrigingSystemSimpleCase::setKrigOptDataWeights(int iptrWeights, bool flagSet)
{
  _isReady = false;
  int nvar = _getNVar();
  if (iptrWeights >= 0 && nvar > 1)
  {
    messerr("The storage of the weights is only coded for Monovariate case");
    return 1;
  }
  _iptrWeights = iptrWeights;
  _flagWeights = true;
  _flagSet     = flagSet;
  return 0;
}

int KrigingSystemSimpleCase::setKrigOptCalcul(const EKrigOpt& calcul)
{
  _isReady = false;
  _krigopt.setOptionCalcul(calcul);
  return 0;
}

int KrigingSystemSimpleCase::setKrigOptFlagGlobal(bool flag_global)
{
  _isReady = false;
  if (!flag_global) return 0;
  messerr("Global is not handled within KrigingSystemSimpleCase anymore");
  return 1;
}

/**
 * Ask for the specific calculation of Z * A-1 * Z
 * @param flag_lterm Flag for asking this specific calculation
 * @return
 * @remark The calculated value can be retrieved using _getLTerm() method
 */

int KrigingSystemSimpleCase::setKrigOptFlagLTerm(bool flag_lterm)
{
  _isReady   = false;
  _flagLTerm = flag_lterm;
  return 0;
}

bool KrigingSystemSimpleCase::_isCorrect()
{
  /****************************/
  /* Checking Space Dimension */
  /****************************/

  int ndim = 0;
  if (_dbin != nullptr)
  {
    if (ndim > 0 && ndim != _dbin->getNDim())
    {
      messerr("Incompatible Space Dimension of '_dbin'");
      return false;
    }
    ndim = _dbin->getNDim();
  }
  if (_dbout != nullptr)
  {
    if (ndim > 0 && ndim != _dbout->getNDim())
    {
      messerr("Incompatible Space Dimension of '_dbout'");
      return false;
    }
    ndim = _dbout->getNDim();
  }
  if (_model != nullptr)
  {
    if (ndim > 0 && ndim != (int)_model->getNDim())
    {
      messerr("Incompatible Space Dimension of '_ model'");
      return false;
    }
    ndim = _model->getNDim();
  }
  if (_neigh != nullptr)
  {
    if (ndim > 0 && ndim != (int)_neigh->getNDim())
    {
      messerr("Incompatible Space Dimension of '_neigh'");
      return false;
    }
    ndim = (int)_neigh->getNDim();
  }

  /****************************/
  /* Checking Variable Number */
  /****************************/

  int nvar = 0;
  if (_dbin != nullptr)
  {
    if (nvar > 0 && nvar != _dbin->getNLoc(ELoc::Z))
    {
      messerr("Incompatible Variable Number of '_dbin'");
      return false;
    }
    nvar = _dbin->getNLoc(ELoc::Z);
  }
  if (_model != nullptr)
  {
    if (nvar > 0 && nvar != _model->getNVar())
    {
      messerr("Incompatible Variable Number of '_ model'");
      return false;
    }
  }

  /***************************************************/
  /* Checking the Number of Covariances in the Model */
  /***************************************************/

  if (_model != nullptr)
  {
    if (_model->getCov() == nullptr)
    {
      messerr("The Model should contain some Covariances defined before Kriging");
      return false;
    }
  }

  /**************************************/
  /* Checking the Validity of the Model */
  /**************************************/

  if (_model != nullptr && !_model->isValid()) return false;

  /******************************************/
  /* Checking the Number of External Drifts */
  /******************************************/

  int nfex = 0;
  if (_model != nullptr)
  {
    if (nfex > 0 && nfex != _model->getNExtDrift())
    {
      messerr("Incompatible Number of External Drifts of '_model'");
      return false;
    }
    nfex = _model->getNExtDrift();
  }
  if (nfex > 0)
  {
    // When External Drifts are defined in the Model,
    // the other files must be consistent

    if (_dbout != nullptr)
    {
      if (nfex != _dbout->getNLoc(ELoc::F))
      {
        messerr("Incompatible Number of External Drifts:");
        messerr("- In 'Model' = %d", nfex);
        messerr("- In '_dbout' = %d", _dbout->getNLoc(ELoc::F));
        return false;
      }
    }
    if (_dbin != nullptr)
    {
      if (_dbin->getNLoc(ELoc::F) == 0)
      {
        if (migrateByLocator(_dbout, _dbin, ELoc::F)) return false;
        // Store the UID of the newly created variables to be deleted at the end of the process
        _dbinUidToBeDeleted = _dbin->getUIDsByLocator(ELoc::F);
      }
      if (nfex != _dbin->getNLoc(ELoc::F))
      {
        messerr("Incompatible Number of External Drifts:");
        messerr("- In 'Model' = %d", nfex);
        messerr("- In 'dbin' = %d", _dbin->getNLoc(ELoc::F));
        return false;
      }
    }
  }

  /*********************************/
  /* Calculate the field extension */
  /*********************************/

  if (_model != nullptr)
  {
    VectorDouble db_mini;
    VectorDouble db_maxi;
    db_mini.resize(ndim, TEST);
    db_maxi.resize(ndim, TEST);

    /* Input Db structure */

    if (_dbin != nullptr)
      _dbin->getExtensionInPlace(db_mini, db_maxi, true);

    /* Output Db structure */

    if (_dbout != nullptr)
      _dbout->getExtensionInPlace(db_mini, db_maxi, true);

    // Merge extensions

    _model->setField(VH::extensionDiagonal(db_mini, db_maxi));
  }

  /*****************************/
  /* Checking the Neighborhood */
  /*****************************/

  if (_neigh != nullptr)
  {
    if (_neigh->getType() == ENeigh::IMAGE)
    {
      messerr("The Image neighborhood may not be used in KrigingSystemSimpleCase anymore");
      messerr("Use 'krimage' instead");
      return false;
    }

    if (_neigh->getType() == ENeigh::UNIQUE && _neigh->getFlagXvalid() && nvar > 1)
    {
      messerr("The algorithm for Cross-Validation in Unique Neighborhood");
      messerr("is restricted to a single variable");
      return false;
    }
  }

  /******************************/
  /* Preparing Non-stationarity */
  /******************************/

  if (_model != nullptr)
  {
    if (!_preparNoStat()) return false;
  }

  return true;
}

bool KrigingSystemSimpleCase::_preparNoStat()
{
  if (!_flagNoStat) return true;
  ACov* cova = _model->_getCovModify();
  if (cova == nullptr)
  {
    messerr("Your Model should contain an ACov item");
    return false;
  }
  cova->manage(_dbin, _dbout);
  return true;
}

/**
 * Returns the coordinates of the neighboring samples
 * @return Array organized by Coordinate (minor) then by Sample (major)
 */
VectorVectorDouble KrigingSystemSimpleCase::getSampleCoordinates(KrigingAlgebraSimpleCase& algebra, int iechout) const
{
  VectorVectorDouble xyz(_ndim);
  for (int idim = 0; idim < _ndim; idim++)
  {
    xyz[idim].resize(_nech);
    for (int iech = 0; iech < _nech; iech++)
    {
      int jech = algebra.getSampleRanks()->at(0)[iech];
      if (jech >= 0)
        xyz[idim][iech] = _dbin->getCoordinate(jech, idim);
      else
        xyz[idim][iech] = _dbout->getCoordinate(iechout, idim);
    }
  }
  return xyz;
}

MatrixRectangular KrigingSystemSimpleCase::getWeights(KrigingAlgebraSimpleCase& algebra) const
{
  const MatrixRectangular* lambda = algebra.getLambda();
  if (lambda == nullptr) return MatrixRectangular();
  return *lambda;
}
MatrixRectangular KrigingSystemSimpleCase::getMu(KrigingAlgebraSimpleCase& algebra) const
{
  const MatrixRectangular* mu = algebra.getMu();
  if (mu == nullptr) return MatrixRectangular();
  return *mu;
}
