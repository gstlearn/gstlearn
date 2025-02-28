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
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "geoslib_define.h"
#include "geoslib_old_f.h"

#include "Estimation/KrigingSystem.hpp"

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
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Polynomials/Hermite.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Space/SpaceRN.hpp"
#include "Estimation/KrigingAlgebra.hpp"
#include "Estimation/KrigOpt.hpp"

#include <math.h>

KrigingSystem::KrigingSystem(Db* dbin,
                             Db* dbout,
                             const ModelGeneric* model,
                             ANeigh* neigh)
  : _dbin(dbin)
  , _dbout(dbout)
  , _model(nullptr)
  , _neigh(neigh)
  , _anam(nullptr)
  , _isReady(false)
  , _algebra()
  , _krigopt()
  , _sampleRanks()
  , _Sigma00()
  , _Sigma()
  , _X()
  , _Sigma0()
  , _X0()
  , _Z()
  , _means()
  , _meansTarget()
  , _iptrEst(-1)
  , _iptrStd(-1)
  , _iptrVarZ(-1)
  , _flagEst(false)
  , _flagStd(false)
  , _flagVarZ(false)
  , _flagDataChanged(false)
  , _calcul(EKrigOpt::POINT)
  , _iptrWeights(-1)
  , _flagWeights(false)
  , _flagSet(true)
  , _flagSimu(false)
  , _nbsimu(0)
  , _rankPGS(-1)
  , _ndiscs()
  , _xvalidEstim(true)
  , _xvalidStdev(true)
  , _xvalidVarZ(false)
  , _rankColCok()
  , _valuesColCok()
  , _flagBayes(false)
  , _priorMean()
  , _priorCov()
  , _postMean()
  , _postCov()
  , _postSimu()
  , _varCorrec()
  , _flagDGM(false)
  , _flagFactorKriging(false)
  , _nclasses(0)
  , _factorClass(0)
  , _matLC(nullptr)
  , _flagLTerm(false)
  , _flagAnam(false)
  , _flagNeighOnly(false)
  , _iptrNeigh(-1)
  , _iechOut(-1)
  , _ndim(0)
  , _nvar(0)
  , _nvarCL(0)
  , _nech(0)
  , _nfeq(0)
  , _neq(0)
  , _nred(0)
  , _nbgh()
  , _dbinUidToBeDeleted()
  , _dboutUidToBeDeleted()
  , _space(SpaceRN::create(2))
  , _p0()
  , _p1()
  , _p2()
  , _p0_memo()
  , _flagNoMatLC(true)
  , _flagVerr(false)
  , _flagNoStat(false)
  , _cova(nullptr)
{
  // _model is a copy of input model to allow modification (still used???)
  if (model != nullptr) _model = (ModelGeneric*) model->clone();

  // Store the pointer casting the input ModelGeneric* into Model*
  // in order to avoid too many dynamic casts in the code
  // _cova = model->getCovAnisoListModify();
  ACov* cov = _model->_getCovModify();
  _cova = dynamic_cast<CovAnisoList*>(cov);
  
  if (model != nullptr)
    _flagNoStat = _cova->isNoStat();

  // Reset the neighborhood
  if (neigh != nullptr)
    neigh->reset();

  // Define local constants
  _flagNoMatLC = _matLC == nullptr;
  _flagVerr    = _dbin->hasLocVariable(ELoc::V);

  _resetMemoryGeneral();
}

KrigingSystem::~KrigingSystem()
{
  // Turn OFF this option for future task
  OptDbg::setCurrentIndex(-1);

  // Clean elements from _dbin

  if (_dbin != nullptr)
  {
    if (!_dbinUidToBeDeleted.empty())
    {
      (void) _dbin->deleteColumnsByUID(_dbinUidToBeDeleted);
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

int KrigingSystem::_getNVar() const
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
  if (!_flagFactorKriging)
  {
    if (_dbin != nullptr)
    {
      if (nvar > 0 && nvar != _dbin->getNLoc(ELoc::Z))
      {
        messerr("Inconsistent number of Variables - Value is returned as 0");
        return 0;
      }
      nvar = _dbin->getNLoc(ELoc::Z);
    }
  }
  return nvar;
}

int KrigingSystem::_getNVarCL() const
{
  if (_flagNoMatLC) return _getNVar();
  return (int)_matLC->getNRows();
}

int KrigingSystem::_getNbfl() const
{
  if (_model == nullptr) return 0;
  return _model->getNDrift();
}

int KrigingSystem::_getNFeq() const
{
  if (_model == nullptr) return 0;
  return _model->getNDriftEquation();
}

int KrigingSystem::_getNeq() const
{
  int neq = _nvar * _nech + _nfeq;
  return neq;
}

void KrigingSystem::_resetMemoryGeneral()
{
  _setInternalShortCutVariablesGeneral();

  _space = SpaceRN::create(_ndim);
  _p0 = SpacePoint(_space);
  _p1 = SpacePoint(_space);
  _p2 = SpacePoint(_space);
  _p0_memo = SpacePoint(_space);
}

/*****************************************************************************/
/*!
 **  Checks if the number of samples is compatible with the number of
 **  drift equations
 **
 ** \return  Error: 1 if an error is found; 0 otherwise
 **
 *****************************************************************************/
bool KrigingSystem::_isAuthorized() const
{
  int ncov = getCovSize();
  int ndrift = getDriftSize();
  return ncov > 0 && ncov >= ndrift;
}

/****************************************************************************/
/*!
 **  Returns the additional variance for continuous moving neighborhood
 **
 ** \return  Variance multiplier for Continuous Option
 **
 ** \param[in]  rank1      Rank of the sample in the first Db
 ** \param[in]  rank2      Rank of the sample in the second Db
 ** \param[in]  eps        Distance tolerance
 **
 ** \remarks In the case of a neighborhood which is not MOVING (or undefined),
 ** \remarks this function systematically returns 0.
 **
 ** \remarks In this function (which is not used currently... should be in ACov)
 ** \remarks the value of MHS(ivar, ivar, iech, iech)  should be multiplied
 ** \remarks by _continuousMultiplier(iech, _iechOut)
 ** \remarks where _iechOut is the rank of the target sample in _dbout
 **
 *****************************************************************************/
double KrigingSystem::_continuousMultiplier(int rank1,int rank2, double eps)
{
  if (_neigh == nullptr) return (0.);
  if (_neigh->getType() != ENeigh::MOVING) return (0.);
  const NeighMoving* neighM = dynamic_cast<const NeighMoving*>(_neigh);
  VectorDouble dd(_ndim);

  /* Calculate the distance increment */

  for (int idim = 0; idim < _ndim; idim++)
    dd[idim] = _dbin->getCoordinate(rank1, idim) - _dbout->getCoordinate(rank2, idim);

  double dist = neighM->getBiPtDist()->getNormalizedDistance(dd);

  double var = 0.;
  if (dist > neighM->getDistCont())
  {
    if (ABS(1. - dist) < eps) dist = 1. - eps;
    var = (dist - neighM->getDistCont()) / (1. - dist);
    var = var * var;
  }
  return (var);
}

void KrigingSystem::_dumpOptions() const
{
  /* Kriging option */

  switch (_calcul.toEnum())
  {
    case EKrigOpt::E_POINT: message("Punctual Estimation\n"); break;

    case EKrigOpt::E_BLOCK:
      message("Block Estimation : Discretization = ");
      for (int idim = 0; idim < _ndim; idim++)
      {
        if (idim != 0) message(" x ");
        message("%d", _ndiscs[idim]);
      }
      message("\n");
      break;

    case EKrigOpt::E_DRIFT: message("Drift Estimation\n"); break;

    case EKrigOpt::E_DGM: message("Discrete Gaussian Model\n"); break;
  }
  message("\n");
}

void KrigingSystem::_rhsDump()
{
  mestitle(0, "RHS of Kriging matrix");
  if (_nech > 0) message("Number of active samples    = %d\n", _nech);
  message("Total number of equations   = %d\n", _neq);
  message("Number of right-hand sides  = %d\n", _nvarCL);
  _dumpOptions();
  _algebra.dumpRHS();
}

void KrigingSystem::_wgtDump()
{
  /* Header */
  mestitle(0, "(Co-) Kriging weights");
  _algebra.dumpWGT();

  /* Auxiliary results for Drift */
  mestitle(0, "Drift or Mean Information");
  _algebra.dumpAux();
}

/****************************************************************************/
/*!
 **  Calculate the final conditional simulation
 **
 ** \param[in]  status    Kriging error status
 **
 *****************************************************************************/
void KrigingSystem::_simulateCalcul(int status)
{
  int ecr = 0;
  for (int isimu = ecr = 0; isimu < _nbsimu; isimu++)
    for (int ivar = 0; ivar < _nvar; ivar++, ecr++)
    {
      double simu = 0.;
      if (status == 0)
      {
        if (_flagBayes)
          simu = _model->evalDriftVarCoef(_dbout, _iechOut, ivar, _postSimu.getColumn(isimu));

        int lec = 0;
        for (int jvar = 0; jvar < _nvar; jvar++)
          for (int iech = 0; iech < _nech; iech++)
          {
            int jech = _nbgh[iech];

            // Get the simulated difference at data point (Simu - Data)
            double diff =
              _dbin->getSimvar(ELoc::SIMU, jech, isimu, jvar, _rankPGS, _nbsimu, _nvar);
            if (FFFF(diff)) continue;

            // Get the kriging weight
            double wgt = _algebra.getLambda()->getValue(lec, ivar);
            lec++;

            // Calculate the Kriging of simulated differences: -sum {wgt * diff}
            simu -= wgt * diff;
          }
      }
      else
      {
        // In case of failure with KS, set the conditional simulation to the mean
        if (_nfeq > 0) simu = TEST;
      }

      /* Add the conditioning kriging to the NC simulation at target */
      _dbout->updSimvar(ELoc::SIMU, _iechOut, isimu, ivar, _rankPGS, _nbsimu, _nvar,
                        EOperator::ADD, simu);
    }
}

/****************************************************************************/
/*!
 **  Calculate the final estimation and storage
 **
 ** \param[in] status   Kriging error status
 **
 *****************************************************************************/
void KrigingSystem::_estimateCalcul(int status)
{
  if (_flagEst)
    _estimateEstim(status);

  /* Variance of the estimation error */

  if (_flagStd)
    _estimateStdv(status);

  /* Variance of the estimator */

  if (_flagVarZ != 0)
    _estimateVarZ(status);

  // Modification specific to Cross-validation options

  if (_neigh->getFlagXvalid())
  {
    for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
    {
      double valdat = _dbin->getZVariable(_iechOut, ivarCL);
      double estim  = (_flagEst) ? _dbout->getArray(_iechOut, _iptrEst + ivarCL) : TEST;
      double stdv   = (_flagStd) ? _dbout->getArray(_iechOut, _iptrStd + ivarCL) : TEST;

      // Modification of Estimation

      if (_flagEst)
      {
        if (_xvalidEstim)
        {
          double estloc = (FFFF(valdat)) ? TEST : estim - valdat;
          _dbout->setArray(_iechOut, _iptrEst + ivarCL, estloc);
        }
      }

      // Modification of Standard Deviation

      if (_flagStd)
      {
        if (_xvalidStdev)
        {
          stdv = (FFFF(estim) || FFFF(valdat) || stdv <= 0.) ? TEST : (estim - valdat) / stdv;
          _dbout->setArray(_iechOut, _iptrStd + ivarCL, stdv);
        }
      }
    }
  }

  /* Kriging weights (stored in _dbin) */

  if (_flagWeights != 0)
  {
    for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
    {
      for (int jech = 0; jech < _nech; jech++)
      {
        if (status != 0) continue;
        double wgt = _algebra.getLambda()->getValue(jech, ivarCL);
        int iech = _nbgh[jech];
        if (_flagSet)
          _dbin->setArray(iech, _iptrWeights + ivarCL, wgt);
        else
          _dbin->updArray(iech, _iptrWeights, EOperator::ADD, wgt);
      }
    }
  }
}

void KrigingSystem::_neighCalcul(int status, const VectorDouble& tab)
{
  int ntab = (int) tab.size();
  for (int i = 0; i < ntab; i++)
  {

    /* Store the parameter */

    double value = (status == 0) ? tab[i] : TEST;
    _dbout->setArray(_iechOut, _iptrNeigh + i, value);
  }

  if (OptDbg::query(EDbg::NBGH) && status == 0)
  {
    mestitle(0, "Neighborhood Parameters");

    message("Number of selected samples          = %d\n", (int) tab[0]);
    message("Maximum neighborhood distance       = %lf\n", tab[1]);
    message("Minimum neighborhood distance       = %lf\n", tab[2]);
    message("Number of non-empty sectors         = %d\n", (int) tab[3]);
    message("Number of consecutive empty sectors = %d\n", (int) tab[4]);
  }
}

/****************************************************************************/
/*!
 **  Establish the calculation of estimation
 **
 ** \param[in]  status  Kriging error code
 **
 *****************************************************************************/
void KrigingSystem::_estimateEstim(int status)
{
  VectorDouble local = _algebra.getEstimation();
  if (local.size() <= 0) return;
  if (status) local.fill(TEST);
  for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
    _dbout->setArray(_iechOut, _iptrEst + ivarCL, local[ivarCL]);
}

/****************************************************************************/
/*!
 **  Establish the calculation of standard deviation
 **
 ** \param[in]  status  Kriging error code
 **
 *****************************************************************************/
void KrigingSystem::_estimateStdv(int status)
{
  VectorDouble local = _algebra.getStdv();
  if (local.size() <= 0) return;
  if (status) local.fill(TEST);
  for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
    _dbout->setArray(_iechOut, _iptrStd + ivarCL, local[ivarCL]);
}

/****************************************************************************/
/*!
 **  Establish the variance of the estimator
 **
 ** \param[in]  status  Kriging error code
 **
 *****************************************************************************/
void KrigingSystem::_estimateVarZ(int status)
{
  VectorDouble local = _algebra.getVarianceZstar();
  if (local.size() <= 0) return;
  if (status) local.fill(TEST);
  for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
    _dbout->setArray(_iechOut, _iptrVarZ + ivarCL, local[ivarCL]);
}

int KrigingSystem::resetData()
{
  const CovCalcMode calcmode(ECalcMember::LHS);
  _sampleRanks = _dbin->getSampleRanks(VectorInt(), _nbgh);
  _Z           = _dbin->getValuesByRanks(_sampleRanks, _means, !_model->hasDrift());
  if (_model->evalCovMatSymByRanks(_Sigma, _dbin, _sampleRanks, -1, &calcmode, false)) return 1;
  if (_model->evalDriftMatByRanks(_X, _dbin, _sampleRanks, -1, ECalcMember::LHS)) return 1;

  if (! _isAuthorized()) return 1;
  _algebra.resetNewData();
  if (_algebra.setData(&_Z, &_sampleRanks, &_meansTarget)) return 1;
  if (_algebra.setLHS(&_Sigma, &_X)) return 1;

  return 0;
}

/**
 * Performs the last operations before launching the loop on Estimations
 * @return
 */
bool KrigingSystem::isReady()
{
  if (!_isCorrect()) return false;

  // Define the means of each variable
  _means       = _model->getMeans();
  // Possible adjust the means in case of presence of 'matLC'
  _meansTarget = _means;
  if (_matLC != nullptr) _meansTarget = _matLC->prodMatVec(_means);

  if ((_neigh != nullptr && _neigh->getType() == ENeigh::UNIQUE) || _flagBayes)
  {
    _sampleRanks = _dbin->getSampleRanks();
    _Z           = _dbin->getValuesByRanks(_sampleRanks, 
                                           _means, !_model->hasDrift());
    if (_algebra.setData(&_Z, &_sampleRanks, &_meansTarget)) return false;

    if (_flagBayes)
    {
      const CovCalcMode calcmode(ECalcMember::LHS);
      if (_model->evalCovMatSymByRanks(_Sigma, _dbin, _sampleRanks, -1, &calcmode, false)) return false;
      if (_model->evalDriftMatByRanks(_X, _dbin, _sampleRanks, -1, ECalcMember::LHS)) return false;
      if (_algebra.setLHS(&_Sigma, &_X)) return false;
    }
  }

  // Perform some pre-calculation when variance of estimator is requested
  if (_flagStd)
  {
    _iechOut = 0;
    if (_model->evalCov0MatByTargetInPlace(_Sigma00, _dbout, _iechOut, _krigopt)) return false;
    if (_algebra.setVariance(&_Sigma00)) return false;
  }

  // Attach the Input and Output Db
  _neigh->attach(_dbin, _dbout);

  // In Bayesian case, calculate the Posterior information
  if (_flagBayes)
    _bayesPreCalculations();

  _isReady = true;
  return _isReady;
}

/**
 * This method closes the use of a KrigingSystem sequence
 */
void KrigingSystem::conclusion()
{
  if (_cova != nullptr)
    _cova->optimizationPostProcess();
}

/**
 * Perform the Kriging of target
 *
 * @param iech_out Rank of the target
 * @return
 */
int KrigingSystem::estimate(int iech_out)
{
  if (! _isReady)
  {
    messerr("You must call 'isReady' before launching 'estimate'");
    return 1;
  }

  // In case of Image Neighborhood, the neighboring samples have already
  // been selected in isReady(). No need to compute them again.
  bool skipCalculAll = false;
  if (_neigh->getType() == ENeigh::IMAGE) skipCalculAll = true;

  // Particular cases:
  // - Cross-validation in Unique Neighborhood
  bool caseXvalidUnique = (_neigh->getType() == ENeigh::UNIQUE && _neigh->getFlagXvalid());

  // Store the Rank of the Target sample
  _iechOut = iech_out;

  int status = 0;
  if (skipCalculAll) goto label_store;

  if (! _dbout->isActive(_iechOut)) return 0;
  OptDbg::setCurrentIndex(_iechOut + 1);
  if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
  {
    if (_flagFactorKriging)
      message("\nProcessing Factor %d / %d\n", _cova->getActiveFactor(), _nclasses);

    mestitle(1, "Target location");
    if (_rankColCok.empty())
      db_sample_print(_dbout, _iechOut, 1, 0, 0, 0);
    else
      db_sample_print(_dbout, _iechOut, 1, 1, 0, 0);
  }

  // Elaborate the Neighborhood
  // For XValid in Unique Neighborhood, turn the Xvalid option OFF during neighborhood search
  if (caseXvalidUnique) _neigh->setFlagXvalid(false);
  _neigh->select(_iechOut, _nbgh);
  status = _setInternalShortCutVariablesNeigh();
  if (_flagNeighOnly) goto label_store;
  if (status) goto label_store;

  /* Establish the Kriging L.H.S. */

  if (!_neigh->isUnchanged() || _neigh->getFlagContinuous() || OptDbg::force())
  {
    status = resetData();
    if (status) goto label_store;
  }

  // Establish the pre-calculation involving the data information

  if (caseXvalidUnique) _neigh->setFlagXvalid(true);

  /* Establish the Kriging R.H.S. */

  if (caseXvalidUnique)
  {
    // For XValid in Unique Neighborhood:
    // - no need to define the RHS information (it will be extracted from LHS°)
    // - only define the indices of the cross_validated columns
    VectorInt xvalidEqs = _xvalidUniqueIndices();
    if (xvalidEqs.size() <= 0)
    {
      // The sample to be cross-validated is not valid, skip
      status = 1;
      goto label_store;
    }
    VectorInt xvalidVars = VH::sequence(_getNVar());
    if (_algebra.setXvalidUnique(&xvalidEqs, &xvalidVars)) return 1;
  }
  else
  {
    if (_model->evalCovMatByTarget(_Sigma0, _dbin, _dbout, _sampleRanks, iech_out, _krigopt, false)) return 1;
    if (_model->evalDriftMatByTarget(_X0, _dbout, iech_out, _krigopt)) return 1;
    if (_algebra.setRHS(&_Sigma0, &_X0)) return 1;
  };

  // Special patch for Colocated CoKriging
  if (!_rankColCok.empty())
  {
    if (_neigh->getType() == ENeigh::MOVING)
    {
      _updateForColCokMoving();
    }
    else
    {
      _valuesColCok = _dbout->getLocVariables(ELoc::Z, _iechOut);
      if (_X.empty()) VH::subtractInPlace(_valuesColCok, _means);
      if (_algebra.setColCokUnique(&_valuesColCok, &_rankColCok)) return 1;
    }
  }

  // Printout for debugging case

  if (!_neigh->isUnchanged() || _neigh->getFlagContinuous() || OptDbg::force())
  {
    // LHS is not printed systematically... only when it has been modified
    if (OptDbg::query(EDbg::KRIGING)) _algebra.dumpLHS();
  }

  if (OptDbg::query(EDbg::KRIGING))
  {
    _rhsDump();
    _wgtDump();
  }

  /* Perform the final estimation */

  label_store:
  // If status is not zero, cancel the current Neighborhood search status
  if (status) _neigh->setIsChanged();

  // Store the results in the output Db

  if (_flagNeighOnly)
  {
    VectorDouble tab = _neigh->summary(_iechOut);
    _neighCalcul(status, tab);
  }
  else 
  {
    // Unique Neighborhood case

    if (_flagSimu)
      _simulateCalcul(status);
    else 
      _estimateCalcul(status);
  }

  // Gaussian transform (optional)
  if (_flagAnam) _transformGaussianToRaw();

  // Final printout
  if (OptDbg::query(EDbg::RESULTS))
  {
    if (_flagSimu)
      _dumpSimulationResults(status);
    else
      _dumpKrigingResults(status);
  }
  return 0;
}

int KrigingSystem::_updateForColCokMoving()
{
  int nvar = (int)_sampleRanks.size();
  int nbfl = _X.getNCols();
  int nrhs = _Sigma0.getNCols();
  int ndim = _dbin->getNDim();

  // If the target coincides with a data point, do not do anything
  // (otherwise the new CoKriging system will be regular)
  VectorDouble coor = _dbout->getSampleCoordinates(_iechOut);
  for (int jech = 0, nech = _nbgh.size(); jech < nech; jech++)
  {
    int iech          = _nbgh[jech];
    bool flagCoincide = true;
    for (int idim = 0; idim < ndim && flagCoincide; idim++)
    {
      if (ABS(_dbin->getCoordinate(iech, idim)) > EPSILON3) flagCoincide = false;
    }
    if (flagCoincide) return 0;
  }

  // Prepare the vector of values from the Target File for Colocated option
  int nAdd = 0;
  VectorDouble newValues(nvar, TEST);
  for (int jvar = 0; jvar < nvar; jvar++)
  {
    int ivar = _rankColCok[jvar];
    if (ivar < 0 || ivar >= _dbout->getNLoc(ELoc::Z)) continue;
    double value = _dbout->getZVariable(_iechOut, ivar);
    if (FFFF(value)) continue;
    newValues[ivar] = (nbfl > 0) ? value : value - _means[ivar];
    nAdd++;
  }
  if (nAdd <= 0) return 0;

  int oldSize = (int)_Z.size();
  int newSize = oldSize + nAdd;

  // Create the indexing vector (>0 for actual samples, <0 for additional sample)
  // Indices are 1-based values (to allow negative and positive distinction) 
  VectorInt adds(newSize);
  int ecr = 0;
  for (int ivar = 0, lec = 0; ivar < nvar; ivar++)
  {
    for (int i = 0, n = (int)_sampleRanks[ivar].size(); i < n; i++, lec++)
      adds[ecr++] = 1 + lec;
    if (! FFFF(newValues[ivar])) adds[ecr++] = -1 - ivar;
  }

  // Update _sampleRanks
  VectorVectorInt newVVI(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    newVVI[ivar] = _sampleRanks[ivar];
    if (! FFFF(newValues[ivar])) newVVI[ivar].push_back(-1);
  }
  _sampleRanks = newVVI;

  // Update Z vector
  VectorDouble newZ = VectorDouble(newSize);
  for (int i = 0; i < newSize; i++)
    newZ[i] = (adds[i] > 0) ? _Z[adds[i] - 1] : newValues[-adds[i] - 1];
  _Z = newZ;

  // Update _Sigma (symmetric square matrix)
  MatrixSquareSymmetric newS = MatrixSquareSymmetric(newSize);
  for (int i = 0; i < newSize; i++)
    for (int j = 0; j <= i; j++)
    {
      double value = TEST;
      if (adds[i] > 0)
      {
        if (adds[j] > 0)
          value = _Sigma.getValue(adds[i] - 1, adds[j] - 1);
        else
          value = _Sigma0.getValue(adds[i] - 1, -adds[j] - 1);
      }
      else
      {
        if (adds[j] > 0)
          value = _Sigma0.getValue(adds[j] - 1, -adds[i] - 1);
        else
          value = _Sigma00.getValue(-adds[i] - 1, -adds[j] - 1);
      }
      newS.setValue(i, j, value);
    }
  _Sigma = newS;

  // Update X
  MatrixRectangular newX = MatrixRectangular(newSize, nbfl);
  for (int i = 0; i < newSize; i++)
    for (int j = 0; j < nbfl; j++)
    {
      double value;
      if (adds[i] > 0)
        value = _X.getValue(adds[i] - 1, j);
      else
        value = _X0.getValue(-adds[i] - 1, j);
      newX.setValue(i, j, value);
    }
  _X = newX;

  // Update Sigma0
  MatrixRectangular newS0 = MatrixRectangular(newSize, nrhs);
  for (int i = 0; i < newSize; i++)
    for (int j = 0; j < nrhs; j++)
    {
      double value;
      if (adds[i] > 0)
        value = _Sigma0.getValue(adds[i] - 1, j);
      else
        value = _Sigma00.getValue(-adds[i] - 1, j);
      newS0.setValue(i, j, value);
    }
  _Sigma0 = newS0;

  // Store the result of Colocated updating in Moving Neighborhood
  if (!_isAuthorized()) return 1;
  _algebra.resetNewData();
  if (_algebra.setData(&_Z, &_sampleRanks, &_meansTarget)) return 1;
  if (_algebra.setLHS(&_Sigma, &_X)) return 1;
  if (_algebra.setRHS(&_Sigma0, &_X0)) return 1;
  return 0;
}

/**
 * @brief Identify the list of equations involving the target sample
 *        within the vector of vector of data indices.
 *        THis is used to mask off equations in the cross_validation in Unique Neighborhood
 * 
 * @return VectorInt Vector of indices to be masked in the Co-Kriging system due to XValid
 */
VectorInt KrigingSystem::_xvalidUniqueIndices() const
{
  VectorInt ranks;
  int lec = 0;
  for (int ivar = 0, nvar = (int)_sampleRanks.size(); ivar < nvar; ivar++)
  {
    for (int i = 0, n = (int) _sampleRanks[ivar].size(); i < n; i++, lec++)
      if (_sampleRanks[ivar][i] == _iechOut) ranks.push_back(lec);
  }
  return ranks;
}

/****************************************************************************/
/*!
 **  Print the results
 **
 ** \param[in] status   Kriging error status
 **
 *****************************************************************************/
void KrigingSystem::_dumpKrigingResults(int status)
{
  if (_neigh->getFlagXvalid())
    mestitle(0, "Cross-validation results");
  else
    mestitle(0, "(Co-) Kriging results");
  message("Target Sample = %d\n", _iechOut + 1);

  /* Loop on the results */

  for (int ivar = 0; ivar < _nvarCL; ivar++)
  {
    if (_neigh->getFlagXvalid())
    {

      // Printout of Cross-Validation test

      message("Variable Z%-2d\n", ivar + 1);
      double estval = TEST;
      double esterr = TEST;
      double sterr  = TEST;
      double sigma  = TEST;

      if (_iptrEst >= 0)
      {
        double trueval = (status == 0) ? _dbin->getZVariable(_iechOut, ivar) : TEST;
        double estim   = (status == 0) ? _dbout->getArray(_iechOut, _iptrEst + ivar) : TEST;

        if (status == 0)
        {
          if (_xvalidEstim)
          {
            estval = estim + trueval;
            esterr = estim;
          }
          else
          {
            estval = estim;
            esterr = estim - trueval;
          }

          tab_printg(" - True value        = ", trueval);
          message("\n");
          tab_printg(" - Estimated value   = ", estval);
          message("\n");
          tab_printg(" - Estimation Error  = ", esterr);
          message("\n");
        }
      }
      if (_iptrStd >= 0)
      {
        double stdev = (status == 0) ? _dbout->getArray(_iechOut, _iptrStd + ivar) : TEST;

        if (status == 0)
        {
          if (_xvalidStdev)
          {
            sterr = stdev;
            sigma = esterr / stdev;
          }
          else
          {
            sigma = stdev;
            sterr = esterr / stdev;
          }

          tab_printg(" - Std. deviation    = ", sigma);
          message("\n");
          tab_printg(" - Normalized Error  = ", sterr);
         message("\n");
        }
      }
    }
    else
    {

      // Printout of the Estimation

      message("Variable Z%-2d\n", ivar + 1);
      if (_iptrEst >= 0)
      {
        double value = (status == 0) ? _dbout->getArray(_iechOut, _iptrEst + ivar) : TEST;
        tab_printg(" - Estimate  = ", value);
        message("\n");
      }
      if (_iptrStd >= 0)
      {
        double value = (status == 0) ? _dbout->getArray(_iechOut, _iptrStd + ivar) : TEST;
        tab_printg(" - Std. Dev. = ", value);
        message("\n");
        tab_printg(" - Variance  = ", FFFF(value) ? TEST : value * value);

        value = _Sigma00.getValue(ivar, ivar);
        message("\n");
        tab_printg(" - Cov(h=0)  = ", value);
        message("\n");
      }
      if (_iptrVarZ >= 0)
      {
        double value = (status == 0) ? _dbout->getArray(_iechOut, _iptrVarZ + ivar) : TEST;
        tab_printg(" - Var(Z*)   = ", value);
        message("\n");
      }
    }
  }
}

void KrigingSystem::_dumpSimulationResults(int status)
{
  mestitle(0, "Simulation results");

  /* Loop on the results */

  int ecr = 0;
  for (int isimu = 0; isimu < _nbsimu; isimu++)
    for (int ivar = 0; ivar < _nvar; ivar++, ecr++)
    {
      message("Simulation #%d of Z%-2d : ", isimu + 1, ivar + 1);
      double value = (status == 0) ? _dbout->getArray(_iechOut, _iptrEst + ecr) : TEST;
      tab_printg(" = ", value);
      message("\n");
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
int KrigingSystem::updKrigOptEstim(int iptrEst,
                                   int iptrStd,
                                   int iptrVarZ,
                                   bool forceNoDual)
{
  _iptrEst  = iptrEst;
  _iptrStd  = iptrStd;
  _iptrVarZ = iptrVarZ;

  _flagEst  = _iptrEst >= 0 || (_neigh->getFlagXvalid() && _iptrStd >= 0);
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

int KrigingSystem::updKrigOptNeighOnly(int iptrNeigh)
{
  _isReady = false;
  if (iptrNeigh < 0)
  {
    messerr("UID for storing Neighborhood variable must be defined");
    return 1;
  }
  _iptrNeigh = iptrNeigh;
  _flagNeighOnly = true;
  return 0;
}

int KrigingSystem::setKrigOptDataWeights(int iptrWeights, bool flagSet)
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
  _flagSet = flagSet;
  return 0;
}

int KrigingSystem::setKrigOptCalcul(const EKrigOpt& calcul,
                                    const VectorInt& ndiscs,
                                    bool flag_per_cell)
{
  _isReady = false;
  _calcul  = calcul;
  bool flagPerCell = false;
  DbGrid* dbgrid   = dynamic_cast<DbGrid*>(_dbout);

  if (_calcul == EKrigOpt::BLOCK)
  {
    if (dbgrid == nullptr)
    {
      messerr("Block Estimation is only possible for Grid '_dbout'");
      return 1;
    }

    // Block support is defined per sample
    if (flag_per_cell)
    {
      flagPerCell = true;
    }
    if (_neigh->getType() == ENeigh::CELL)
    {
      flagPerCell = true;
    }

    // Check that discretization is defined
    if (ndiscs.empty())
    {
      messerr("In case of BLOCK kriging, you must define the discretization coefficients");
      messerr("i.e. a vector (dimension equal Space Dimension) filled with positive numbers");
      return 1;
    }

    // Discretization is stored

    _ndiscs = ndiscs;
  }
  else
  {
    _ndiscs.clear();
  }

  // New style operation
  _krigopt.setKrigingOption(calcul, dbgrid, ndiscs, flagPerCell);
  return 0;
}

/**
 * Set the flag for performing Cross-Validation
 * @param flag_xvalid True if the Cross-Validation option is switched ON
 * @param flag_kfold  True if the KFold option is switch ON
 * @param optionXValidEstim True for Z*-Z; False for Z*
 * @param optionXValidStdev True for (Z*-Z)/S; False for S
 * @param optionXValidVarZ  True for Var(Z*)
 * @return
 *
 * @remark The KFold option requires a Code to be assigned to each Datum
 * @remark In the Neighborhood search for a Target sample, some other data points
 * @remark are discarded:
 * @remark - either the Target sample itself (Leave-One-Point-Out) if KFold is False
 * @remark - all samples with same code as Target if KFold is True
 */
int KrigingSystem::setKrigOptXValid(bool flag_xvalid,
                                    bool flag_kfold,
                                    bool optionXValidEstim,
                                    bool optionXValidStdev,
                                    bool optionXValidVarZ)
{
  _isReady = false;
  if (! flag_xvalid)
  {
    _neigh->setFlagXvalid(false);
    _neigh->setFlagKFold(false);
  }
  else
  {
    _neigh->setFlagXvalid(flag_xvalid);
    if (flag_kfold)
    {
      if (_neigh->getType() == ENeigh::UNIQUE)
      {
        messerr("K-FOLD is not available in Unique Neighborhood");
        return 1;
      }
      if (! _dbin->hasLocVariable(ELoc::C))
        messerr("The K-FOLD option is ignored as no Code is defined");
    }
    _neigh->setFlagKFold(flag_kfold);
  }
  _xvalidEstim = optionXValidEstim;
  _xvalidStdev = optionXValidStdev;
  _xvalidVarZ  = optionXValidVarZ;
  return 0;
}

/****************************************************************************/
/*!
 **  Check the consistency of the Colocation specification
 **
 ** \return  Error return code
 **
 ** \param[in]  rank_colcok   Array of ranks of colocated variables
 **
 ** \remarks The array 'rank_colcok' (if present) must be dimensioned
 ** \remarks to the number of variables in Dbin.
 ** \remarks Each element gives the rank of the colocated variable within Dbout
 ** \remarks or -1 if not colocated
 ** \remarks If the array 'rank_colcok' is absent, colocation option is OFF.
 **
 *****************************************************************************/
int KrigingSystem::setKrigOptColCok(const VectorInt& rank_colcok)
{
  if (rank_colcok.empty()) return 0;

  _isReady = false;
  int nvar    = _getNVar();
  _rankColCok = rank_colcok;
  _valuesColCok.resize(nvar);

  /* Loop on the ranks of the colocated variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int jvar = _rankColCok[ivar];
    if (jvar < 0) continue;
    if (jvar > _dbout->getNLoc(ELoc::Z))
    {
      messerr("Error in the Colocation array:");
      messerr("Input variable (#%d): rank of the colocated variable is %d",
              ivar + 1, jvar);
      messerr("But the Output file only contains %d attributes(s)",
              _dbout->getNColumn());
      return (1);
    }
  }
  return 0;
}

int KrigingSystem::setKrigOptBayes(bool flag_bayes,
                                   const VectorDouble& prior_mean,
                                   const MatrixSquareSymmetric& prior_cov)
{
  _isReady = false;
  int nfeq = _getNFeq();
  if (flag_bayes)
  {
    VectorDouble local_mean = prior_mean;
    MatrixSquareSymmetric local_cov  = prior_cov;

    if (local_mean.empty())
      local_mean.resize(nfeq, 0.);
    if (local_cov.empty())
    {
      local_cov.resetFromValue(nfeq, nfeq, 0.);
      for (int i = 0; i < nfeq; i++)
        local_cov.setValue(i,i, 1.);
    }
    if ((int) local_mean.size() != nfeq)
    {
      messerr("Size of argument 'prior_mean'(%d)",(int) local_mean.size());
      messerr("should be equal to the Number of Drift Equations(%d)",nfeq);
      return 1;
    }
    if ((int) local_cov.size() != nfeq * nfeq)
    {
      messerr("Size of argument 'prior_cov'(%d)",(int) local_cov.size());
      messerr("should be equal to the Number of Drift Equations (squared) (%d)",
              nfeq * nfeq);
      return 1;
    }
    if (_neigh->getType() != ENeigh::UNIQUE)
    {
      messerr("The Bayesian Estimation of the Drift Coefficients");
      messerr("is only available in Unique Neighborhood");
      return 1;
    }

    // Set the parameters
    _priorMean = local_mean;
    _priorCov  = local_cov;
    _varCorrec.resize(_nvarCL, _nvarCL);

    // Pass the Bayesian information to '_algebra'
    if (_algebra.setBayes(&_priorMean, &_priorCov)) return 1;
  }
  _flagBayes = flag_bayes;
  return 0;
}

/**
 * Define the output as Linear Combinations of the Input Variables
 * @param matLC Vector of Vectors of weights (see remarks)
 * @return
 * @remarks The first dimension of 'matLC' is the number of Output variables
 * @remarks The second dimension is the number of input Variables.
 */
int KrigingSystem::setKrigOptMatLC(const MatrixRectangular* matLC)
{
  if (matLC == nullptr) return 0;
  _isReady = false;
  int n1 = (int) matLC->getNRows();
  int n2 = (int) matLC->getNCols();

  if (n1 > _getNVar())
  {
    messerr("First dimension of 'matLC' (%d)",(int) n1);
    messerr("should be smaller than the number of variables in the model (%d)",
            _getNVar());
    return 1;
  }
  if (n2 != _getNVar())
  {
    messerr("Second dimension of 'matLC' (%d)",(int) n2);
    messerr("should be equal to the number of variables in the model (%d)",
            _getNVar());
    return 1;
  }
  _matLC = matLC;
  _flagNoMatLC = false;
  _resetMemoryGeneral();

  
  // New style assignment
  _krigopt.setMatLC(matLC, _getNVar());
  return 0;
}

int KrigingSystem::setKrigOptFlagSimu(bool flagSimu, int nbsimu, int rankPGS)
{
  _isReady = false;
 _flagSimu = flagSimu;
 _nbsimu = nbsimu;
 _rankPGS = rankPGS;
 _neigh->setFlagSimu(flagSimu);
 return 0;
}


int KrigingSystem::setKrigOptDGM(bool flag_dgm, double eps)
{
  _isReady = false;
  if (! flag_dgm)
  {
    _flagDGM = flag_dgm;
    return 0;
  }

  const Model* model = dynamic_cast<const Model*>(_model);
  if (model->getCovMinIRFOrder() != -1)
  {
    messerr("The option DGM is limited to Stationary Covariances");
    return 1;
  }
  if (_model->getNVar() != 1)
  {
    messerr("The DGM option is limited to the Monovariate case");
    return 1;
  }
  Model* model_old = _castInOldModel();
  if (model_old == nullptr) return 1;
  if (ABS(model_old->getTotalSill(0,0) - 1.) > eps)
  {
    messerr("The DGM option requires a Model with Total Sill equal to 1.");
    return 1;
  }
  _flagDGM = flag_dgm;

  // New style assignment
  _krigopt.setKrigingDGM(flag_dgm);
  return 0;
}

int KrigingSystem::setKrigOptFlagGlobal(bool flag_global)
{
  _isReady = false;
  if (!flag_global) return 0;
  messerr("Global is not handled within KrigingSystem anymore");
  return 1;
}

/**
 * Ask for the specific calculation of Z * A-1 * Z
 * @param flag_lterm Flag for asking this specific calculation
 * @return
 * @remark The calculated value can be retrieved using _getLTerm() method
 */

int KrigingSystem::setKrigOptFlagLTerm(bool flag_lterm)
{
  _isReady = false;
  _flagLTerm = flag_lterm;
  return 0;
}

Model* KrigingSystem::_castInOldModel()
{
  Model* model_old = dynamic_cast<Model*>(_model);
  if (model_old == nullptr)
  {
    messerr("This method is only implemented for Model(old_style)");
  }
  return model_old;
}

/**
 * Perform Gaussian Anamoprhosis kriging
 * @param anam Pointer to the AAnam structure
 * @return
 */
int KrigingSystem::setKrigOptAnamophosis(AAnam* anam)
{

  _isReady = false;
  int nvar = _getNVar();
  if (nvar != 1)
  {
    messerr("This procedure is limited to the monovariate case");
    return 1;
  }

  // Check that sill of the (monovariate) Model is smaller or equal to 1
  Model* model_old = _castInOldModel();
  if (model_old == nullptr) return 1;
  double total = model_old->getTotalSill(0, 0);
  if (total > 1.)
  {
    messerr("This procedure requires the Sill of the Model (%lf)",total);
    messerr("to be smaller than 1.");
    return 1;
  }
  _flagAnam = true;
  _anam = anam;
  return 0;
}

int KrigingSystem::setKrigOptFactorKriging(bool flag_factor_kriging)
{
  _isReady = false;
  if (! flag_factor_kriging)
  {
    _flagFactorKriging = false;
  }
  else
  {
    CovLMCAnamorphosis* covAnam = dynamic_cast<CovLMCAnamorphosis*>(_cova);
    if (covAnam == nullptr)
    {
      messerr("Your Model should contain a CovLMCAnamorphosis covariance item");
      return 1;
    }
    if (!covAnam->hasAnam())
    {
      messerr("You may not use this option as there is no Anamorphosis defined");
      return 1;
    }
    _flagFactorKriging = true;
  }
  return 0;
}

int KrigingSystem::updKrigOptIclass(int index_class, int nclasses)
{
  if (! _flagFactorKriging)
  {
    messerr("Setting the Class Index only makes sense if 'flagFactorKriging' is ON");
    messerr("Use 'setKrigOptFactorKriging()' beforehand");
    return 1;
  }
  CovLMCAnamorphosis* covAnam = dynamic_cast<CovLMCAnamorphosis*>(_cova);
  if (covAnam == nullptr)
  {
    messerr("Your Model should contain a CovLMCAnamorphosis covariance item");
    return 1;
  }
  covAnam->setActiveFactor(index_class);
  _nclasses    = nclasses;
  _factorClass = index_class;

  // Update C00 if the variance calculation is required
  if (_flagStd)
  {
    if (_model->evalCov0MatByTargetInPlace(_Sigma00, _dbout, 0, _krigopt)) return 1;
    if (_algebra.setVariance(&_Sigma00)) return 1;
  }

  // Cancel any already existing Neighborhood
  _neigh->setIsChanged();

  return 0;
}

bool KrigingSystem::_isCorrect()
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
  if (_dbin != nullptr && ! _flagSimu)
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

  if (_model != nullptr && ! _model->isValid()) return false;

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
    db_mini.resize(ndim,TEST);
    db_maxi.resize(ndim,TEST);

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
      messerr("The Image neighborhood may not be used in KrigingSystem anymore");
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

  /**************************/
  /* Checking cross-options */
  /**************************/

  if (_flagDGM && (_calcul == EKrigOpt::BLOCK || _calcul == EKrigOpt::DRIFT))
  {
    messerr("The DGM option is incompatible with 'Block' calculation option");
    return false;
  }
  return true;
}

bool KrigingSystem::_preparNoStat()
{
  if (!_flagNoStat) return true;
  CovAnisoList* cova = dynamic_cast<CovAnisoList*>(_cova);
  if (cova == nullptr)
  {
    messerr("Your Model should contain a CovAnisoList covariance item");
    return false;
  }
  cova->manage(_dbin, _dbout);
  return true;
}

/**
 * Returns the coordinates of the neighboring samples
 * @return Array organized by Coordinate (minor) then by Sample (major)
 */
VectorVectorDouble KrigingSystem::getSampleCoordinates() const
{
  VectorVectorDouble xyz(_ndim);
  for (int idim = 0; idim < _ndim; idim++)
  {
    xyz[idim].resize(_nech);
    for (int iech = 0; iech < _nech; iech++)
    {
      int jech = _nbgh[iech];
      if (jech >= 0)
        xyz[idim][iech] = _dbin->getCoordinate(jech, idim);
      else
        xyz[idim][iech] = _dbout->getCoordinate(_iechOut, idim);
    }
  }
  return xyz;
}

/****************************************************************************/
/*!
 **  Perform the Bayesian estimation of the Drift Coefficients
 **  (only in Unique Neighborhood)
 **
 *****************************************************************************/
int KrigingSystem::_bayesPreCalculations()
{
  if (_dbin == nullptr) return 1;

  /* Retreive posterior statements */

  _postCov  = *(_algebra.getPostCov());
  _postMean = _algebra.getPostMean();

  if (OptDbg::query(EDbg::BAYES))
  {
    mestitle(0, "Bayesian Drift coefficients");
    _algebra.dumpAux();
  }

  // Particular case of Simulation: Simulate several outcomes for posterior means

  if (_flagSimu) _bayesPreSimulate();

  _neigh->reset();
  return 0;
}

/****************************************************************************/
/*!
 **  Simulate the drift coefficients from the posterior distributions
 **
 *****************************************************************************/
void KrigingSystem::_bayesPreSimulate()
{
  if (_nfeq <= 0) return;
  int memo = law_get_random_seed();
  CholeskyDense postCovChol;

  // Dimension '_postSimu' to store simulated posterior mean
  _postSimu.resize(_nfeq, _nbsimu);

  /* Core allocation */

  MatrixRectangular rndmat(_nfeq, 1);
  MatrixRectangular simu(_nfeq, 1);

  /* Cholesky decomposition */

  if (postCovChol.setMatrix(&_postCov))
  {
    messerr("Error in the Cholesky Decomposition of the covariance matrix");
    messerr("The Drift coefficients have been set to their posterior mean");
    for (int isimu = 0; isimu < _nbsimu; isimu++)
      for (int il = 0; il < _nfeq; il++)
        _postSimu.setValue(il, isimu, _postMean[il]);
  }
  else
  {
    VectorDouble trimat = postCovChol.getLowerTriangle();
    for (int isimu = 0; isimu < _nbsimu; isimu++)
    {

      /* Draw a vector of gaussian independent values */

      for (int il = 0; il < _nfeq; il++)
        rndmat.setValue(il, 0, law_gaussian());

      /* Product of the Lower triangular matrix by the random vector */

      postCovChol.matProductInPlace(1, rndmat, simu);

      /* Add the mean */

      for (int il = 0; il < _nfeq; il++)
        _postSimu.setValue(il, isimu, simu.getValue(il, 0) + _postMean[il]);
    }
  }

  /* If DEBUG option is switched ON, the values are printed out */

  if (OptDbg::query(EDbg::BAYES))
  {
    mestitle(1, "Simulation of Drift Coefficients (for Bayesian Simulation)");
    message("Rank     Drift Coefficients\n");
    for (int isimu = 0; isimu < _nbsimu; isimu++)
    {
      message(" %3d ", isimu + 1);
      for (int il = 0; il < _nfeq; il++)
        message(" %lf", _postSimu.getValue(il, isimu));
      message("\n");
    }
  }
  law_set_random_seed(memo);
}

/****************************************************************************/
/*!
 **  Transform the Kriging results from gaussian to raw
 **
 ** \remark  This procedure is designed for the monovariate case
 ** \remark  It assumes that the kriging estimate and variance are already
 ** \remark  calculated
 **
 *****************************************************************************/
void KrigingSystem::_transformGaussianToRaw()
{
  if (_anam == nullptr) return;
  const AnamHermite *anam_hermite = dynamic_cast<const AnamHermite*>(_anam);

  /* Get the estimation */

  double est = _dbout->getArray(_iechOut, _iptrEst);

  /* Get the variance of the kriging error */

  double std = _dbout->getArray(_iechOut, _iptrStd);

  /* Calculate the conditional expectation */

  double condexp = hermiteCondExpElement(est, std, anam_hermite->getPsiHns());
  _dbout->setArray(_iechOut, _iptrEst, condexp);

  /* Calculate the conditional variance */

  double condvar = hermiteCondStdElement(est, std, anam_hermite->getPsiHns());
  _dbout->setArray(_iechOut, _iptrStd, condvar);
}

void KrigingSystem::_setInternalShortCutVariablesModel()
{
  _nvar = _getNVar();
  _nfeq = _getNFeq();
  _neq  =  _getNeq(); // reset as it depends on nech and Model
}
/**
 * Assign the values to local variables used as shortcuts
 * @return 1 if the number of active sample is zero
 */
int KrigingSystem::_setInternalShortCutVariablesNeigh()
{
  _nech = getNech();
  _neq  = _getNeq();
  return (_nech <= 0);
}
void KrigingSystem::_setInternalShortCutVariablesGeneral()
{
  _ndim   = getNDim();
  _nvarCL = _getNVarCL();
  _setInternalShortCutVariablesModel();
}
MatrixRectangular KrigingSystem::getWeights() const
{
  const MatrixRectangular* lambda = _algebra.getLambda();
  if (lambda == nullptr) return MatrixRectangular();
  return *lambda;
}
MatrixRectangular KrigingSystem::getMu() const
{
  const MatrixRectangular* mu = _algebra.getMu();
  if (mu == nullptr) return MatrixRectangular();
  return *mu;
}
