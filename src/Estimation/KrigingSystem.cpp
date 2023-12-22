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
#include "geoslib_old_f.h"

#include "Estimation/KrigingSystem.hpp"

#include "Enum/EKrigOpt.hpp"
#include "Enum/ECalcMember.hpp"
#include "Enum/ELoc.hpp"
#include "Enum/ENeigh.hpp"

#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/PtrGeos.hpp"
#include "Model/Model.hpp"
#include "Model/ANoStat.hpp"
#include "Neigh/ANeigh.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighImage.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Covariances/ACovAnisoList.hpp"
#include "Polynomials/Hermite.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Calculators/CalcMigrate.hpp"
#include "Space/SpaceRN.hpp"

#include <math.h>

#define FF(ib,il)         (ff [(il) * shift + (ib)])
#define FF0(ib,iv)        (ff0 [(ib) + _nfeq * (iv)])
#define SIGMA(ib,jb)      (sigma[(jb)    * shift + (ib)])
#define SMEAN(i,isimu)    (_postSimu[isimu][i])
#define IND(iech, ivar)   ((iech) + (ivar) * _nech)
#define IJVAR(ivar,jvar)  ((ivar) + (jvar) * _nvar)

KrigingSystem::KrigingSystem(Db* dbin,
                             Db* dbout,
                             const Model* model,
                             ANeigh* neigh)
    : _dbin(dbin),
      _dbout(dbout),
      _modelInit(nullptr),
      _neigh(neigh),
      _anam(nullptr),
      _isReady(false),
      _model(nullptr),
      _optimEnabled(false),
      _iptrEst(-1),
      _iptrStd(-1),
      _iptrVarZ(-1),
      _flagEst(false),
      _flagStd(false),
      _flagVarZ(false),
      _flagGlobal(false),
      _flagDataChanged(false),
      _calcul(EKrigOpt::POINT),
      _iptrWeights(-1),
      _flagWeights(false),
      _flagSet(true),
      _flagSimu(false),
      _nbsimu(0),
      _rankPGS(-1),
      _flagCode(false),
      _flagPerCell(false),
      _ndiscNumber(0),
      _ndiscs(),
      _disc1(),
      _disc2(),
      _xvalidEstim(true),
      _xvalidStdev(true),
      _xvalidVarZ(false),
      _rankColCok(),
      _flagBayes(false),
      _seedForBayes(123456),
      _priorMean(),
      _priorCov(),
      _postMean(),
      _postCov(),
      _postSimu(),
      _varCorrec(),
      _modelSimple(nullptr),
      _flagDGM(false),
      _flagFactorKriging(false),
      _nclasses(0),
      _matCL(),
      _flagLTerm(false),
      _lterm(0.),
      _flagAnam(false),
      _seedForImage(14351),
      _dbaux(),
      _flagKeypairWeights(false),
      _flagNeighOnly(false),
      _iptrNeigh(-1),
      _iechOut(-1),
      _ndim(0),
      _nvar(0),
      _nvarCL(0),
      _nech(0),
      _nbfl(0),
      _nfeq(0),
      _nfex(0),
      _neq(0),
      _nred(0),
      _flagIsotopic(true),
      _nbgh(),
      _flag(),
      _covtab(),
      _drftab(),
      _lhs(),
      _lhsinv(),
      _rhs(),
      _wgt(),
      _zam(),
      _zext(),
      _var0(),
      _dbinUidToBeDeleted(),
      _dboutUidToBeDeleted(),
      _space(2), // empty constructor does not exist. Anyhow it will be overwritten next.
      _p0(),
      _p1(),
      _p2(),
      _p0_memo(),
      _flagNoMatCL(true),
      _flagVerr(false)
{
  // _modelInit is a copy of the input model (const) to allow modifying it
  if (model != nullptr)
    _modelInit = model->clone();

  // Set the current Model to _modelInit
  _model = _modelInit;

  // Reset the neighborhood
  if (neigh != nullptr)
  {
    neigh->reset();
  }

  // Define local constants
  _flagNoMatCL = _isMatCLempty();
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
      (void) _dbout->deleteColumnsByUID(_dboutUidToBeDeleted);
    }
  }

  // Clean the auxiliary file for Image

  if (_dbaux != nullptr)
  {
    delete _dbaux;
    _dbaux = nullptr;
  }

  // Clean auxiliary Model (if available)

  if (_modelSimple != nullptr)
  {
    delete _modelSimple;
    _modelSimple = nullptr;
  }

  // Clean elements from _modelInit

  if (_modelInit != nullptr)
  {
    if (_modelInit->isNoStat())
    {
      const ANoStat *nostat = _modelInit->getNoStat();

      // Detach the Input Db
      if (_dbin != nullptr) nostat->detachFromDb(_dbin, 1);

      // Detach the output Db
      if (_dbout != nullptr) nostat->detachFromDb(_dbout, 2);
    }

    delete _modelInit;
    _modelInit = nullptr;
  }
}

int KrigingSystem::_getNVar() const
{
  int nvar = 0;
  if (_model != nullptr)
  {
    if (nvar > 0 && nvar != _model->getVariableNumber())
    {
      messerr("Inconsistent number of Variables - Value is returned as 0");
      return 0;
    }
    nvar = _model->getVariableNumber();
  }

  // In the case of factor kriging, the number of Z-variables in the Data file
  // does not give the number of variables. Check should be avoided
  if (!_flagFactorKriging)
  {
    if (_dbin != nullptr)
    {
      if (nvar > 0 && nvar != _dbin->getLocNumber(ELoc::Z))
      {
        messerr("Inconsistent number of Variables - Value is returned as 0");
        return 0;
      }
      nvar = _dbin->getLocNumber(ELoc::Z);
    }
  }
  return nvar;
}

int KrigingSystem::_getNVarCL() const
{
  if (_flagNoMatCL)
    return _getNVar();
  else
    return (int) _matCL.size();
}

int KrigingSystem::_getNbfl() const
{
  if (_model == nullptr) return 0;
  return _model->getDriftNumber();
}

int KrigingSystem::_getNFeq() const
{
  if (_model == nullptr) return 0;
  return _model->getDriftEquationNumber();
}

/**
 * Returns the number of samples in the neighborhood
 * @return Number of samples
 * @remark If collocated option is switch ON, samples are alreadycount (rank < 0)
 */
int KrigingSystem::getNech() const
{
  if (_dbin == nullptr) return 0;
  return (int) _nbgh.size();
}

int KrigingSystem::getNDim() const
{
  if (_dbin == nullptr) return 0;
  return _dbin->getNDim();
}

int KrigingSystem::_getNFex() const
{
  if (_model == nullptr) return 0;
  return _model->getExternalDriftNumber();
}

int KrigingSystem::getNeq() const
{
  int neq = _nvar * _nech + _nfeq;
  return neq;
}

int KrigingSystem::_getNDisc() const
{
  return _ndiscNumber;
}

void KrigingSystem::_resetMemoryPerNeigh()
{
  _flag.resize(_neq);
  _lhs.resize (_neq, _neq);
  _rhs.resize (_neq, _nvarCL);
  _wgt.resize (_neq, _nvarCL);
  _zam.resize (_neq, 1);
  _zext.resize(_neq, 1); // in fact, it should be allocated to _nred <= _neq (unknown at this stage)
}

void KrigingSystem::_resetMemoryGeneral()
{
  _setInternalShortCutVariablesGeneral();

  _covtab.reset(_nvar, _nvar);
  _drftab.resize(_nbfl);
  _var0.reset(_nvarCL, _nvarCL);

  _space = SpaceRN(_ndim);
  _p0 = SpacePoint(&_space);
  _p1 = SpacePoint(&_space);
  _p2 = SpacePoint(&_space);
  _p0_memo = SpacePoint(&_space);
}

/****************************************************************************/
/*!
 **  Returns the coordinate of the data (at rank if rank >= 0)
 **  or of the target (if rank < 0)
 **
 ** \param[in]  loc_rank   Rank of the sample
 ** \param[in]  idim   Rank of the coordinate
 **
 *****************************************************************************/
double KrigingSystem::_getIdim(int loc_rank, int idim) const
{
  if (loc_rank >= 0)
  {
    return _dbin->getCoordinate(loc_rank, idim);
  }
  else
  {
    return _dbout->getCoordinate(_iechOut, idim);
  }
}

/****************************************************************************/
/*!
 **  Returns the value of the external drift "rank" (if rank >= 0)
 **  or of the target (if rank < 0)
 **
 ** \param[in]  rank   Rank of the sample
 ** \param[in]  ibfl   Rank of the external drift
 **
 *****************************************************************************/
double KrigingSystem::_getFext(int rank, int ibfl) const
{
  if (rank >= 0)
  {
    return _dbin->getLocVariable(ELoc::F,rank, ibfl);
  }
  else
  {
    return _dbout->getLocVariable(ELoc::F,_iechOut, ibfl);
  }
}

/****************************************************************************/
/*!
 **  Returns the value of the variable (at rank if rank >= 0)
 **  or of the target (if rank < 0)
 **
 ** \param[in]  rank   Rank of the sample
 ** \param[in]  ivar   Rank of the variable
 **
 ** \remarks   In case of simulation, the variable of the first simulation
 ** \remarks   is systematically returned. This has no influence on the rest
 ** \remarks   of the calculations
 **
 *****************************************************************************/
double KrigingSystem::_getIvar(int rank, int ivar) const
{
  if (rank >= 0)
  {

    // Variable in the Input file

    if (! _flagSimu)

      // Particular case of simulations

      return _dbin->getLocVariable(ELoc::Z,rank, ivar);

    else

      // Case of the traditional kriging based on Z-variables

      return _dbin->getSimvar(ELoc::SIMU, rank, 0, ivar, 0, 1, 0);
  }
  else
  {

    // Variable in the Output file: colocated case

    int jvar = (_rankColCok.empty()) ? -1 : _rankColCok[ivar];
    if (IFFFF(jvar))
      return TEST;
    else
      return _dbout->getArray(_iechOut, jvar);
  }
}

/****************************************************************************/
/*!
 **  Returns the value of the measurement error (at rank if rank >= 0)
 **  or of the target (if rank < 0)
 **
 ** \param[in]  rank     Rank of the sample
 ** \param[in]  ivar     Rank of the variable
 **
 *****************************************************************************/
double KrigingSystem::_getVerr(int rank, int ivar) const
{
  if (rank >= 0)
  {
    return _dbin->getLocVariable(ELoc::V,rank, ivar);
  }
  else
  {
    return _dbout->getLocVariable(ELoc::V,_iechOut, ivar);
  }
}
double KrigingSystem::_getMean(int ivarCL) const
{
  double value = 0.;
  if (_flagNoMatCL)
  {
    value = _model->getMean(ivarCL);
  }
  else
  {
    for (int ivar = 0; ivar < _nvar; ivar++)
      value += _matCL[ivarCL][ivar] * _model->getMean(ivar);
  }
  return value;
}

double KrigingSystem::_getDriftCoef(int ivar, int il, int ib) const
{
  return _model->getDriftCoef(ivar, il, ib);
}

void KrigingSystem::_setFlag(int iech, int ivar, int value)
{
  _flag[iech + ivar * _nech]= value;
}

int KrigingSystem::_getFlag(int iech, int ivar)
{
  return _flag[iech + ivar * _nech];
}

/****************************************************************************/
/*!
 **  Define the array flag to convert from isotropic to heterotopic case
 **  Stores the reduced number of equations in member '_nred'
 **
 *****************************************************************************/
void KrigingSystem::_flagDefine()
{
  for (int i = 0; i < _neq; i++) _flag[i] = 1;

  /* Check on the coordinates */

  for (int iech = 0; iech < _nech; iech++)
  {
    int nbgh_iech = _nbgh[iech];
    bool valid = true;
    for (int idim = 0; idim < _ndim; idim++)
      if (FFFF(_getIdim(nbgh_iech, idim))) valid = false;
    if (! valid)
      for (int ivar = 0; ivar < _nvar; ivar++)
        _setFlag(iech,ivar,0);
  }

  /* Check on the data values */

  for (int iech = 0; iech < _nech; iech++)
  {
    int nbgh_iech = _nbgh[iech];
    for (int ivar = 0; ivar < _nvar; ivar++)
      if (FFFF(_getIvar(nbgh_iech, ivar)))
        _setFlag(iech,ivar,0);
  }

  /* Check on the external drifts */

  if (_nfex > 0)
  {
    for (int iech = 0; iech < _nech; iech++)
    {
      int nbgh_iech = _nbgh[iech];
      for (int ibfl = 0; ibfl < _nfex; ibfl++)
        if (FFFF(_getFext(nbgh_iech, ibfl)))
          for (int ivar = 0; ivar < _nvar; ivar++)
            _setFlag(iech, ivar, 0);
    }
  }

  /* Check on the drift */

  for (int ib = 0; ib < _nfeq; ib++)
  {
    int valid = 0;
    for (int il = 0; il < _nbfl; il++)
      for (int ivar = 0; ivar < _nvar; ivar++)
      {
        if (_getDriftCoef(ivar, il, ib) == 0.) continue;
        for (int iech = 0; iech < _nech; iech++)
          if (!FFFF(_getIvar(_nbgh[iech], ivar))) valid++;
      }
    if (valid <= 0) _setFlag(_nech+ib, _nvar-1, 0);
  }

  /* Calculate the new number of equations */

  int count = 0;
  for (int i = 0; i < _neq; i++)
  {
    if (_flag[i] != 0) count++;
  }
  _nred = count;
  _flagIsotopic = (_nred == _neq);
  return;
}

/*****************************************************************************/
/*!
 **  Checks if the number of samples is compatible with the number of
 **  drift equations
 **
 ** \return  Error: 1 if an error is found; 0 otherwise
 **
 *****************************************************************************/
bool KrigingSystem::_isAuthorized()
{
  /* Preliminary check */

  if (_nech * _nvar < _nfeq) return false;

  /* Check that enough information is present */

  int n_cov = 0;
  for (int i = 0; i < _nvar * _nech; i++)
    n_cov += _flag[i];

  int n_drf = 0;
  for (int i = 0; i < _nfeq; i++)
    n_drf += _flag[i + _nvar * _nech];

  if (n_cov <= 0 || n_cov < n_drf) return false;

  return true;
}

/**
 * Modify the covariance before calling the covariance evaluation
 * This makes sense only when non-stationarity is defined
 * @param icas1     1 for dbin and 2 for Dbout
 * @param iech1     Rank of the first sample (or -1 for target)
 * @param icas2     1 for dbin and 2 for Dbout
 * @param iech2     Rank of the second sample (or -1 for target)
 */
void KrigingSystem::_covUpdate(int icas1, int iech1, int icas2, int iech2)
{
  const ANoStat *nostat = _model->getNoStat();
  nostat->updateModel(_model, icas1, iech1, icas2, iech2);
}

/**
 * Module for calculating the covariance internally (for 0 distance)
 * @param mode      CovCalcMode structure
 */
void KrigingSystem::_covtab0Calcul(const CovCalcMode* mode)
{
  _model->eval0MatInPlace(_covtab, mode);
}


/**
 * Module for calculating the covariance internally
 * It is called for LHS (iech1>=0 && iech2>=0), RHS (iech1>=0 && iech2=-1) and VAR (iech1=-1 && iech2=-1)
 * @param iech1     Rank of the first sample (or -1 for target)
 * @param iech2     Rank of the second sample (or -1 for target)
 * @param mode      CovCalcMode structure
 */
void KrigingSystem::_covtabCalcul(int iech1,
                                  int iech2,
                                  const CovCalcMode* mode)
{
  if (_optimEnabled)
  {
    _model->evalMatOptimInPlace(iech1, iech2, _covtab, mode);
  }
  else
  {
    if (iech1 >= 0)
      _dbin->getSampleCoordinatesAsSPInPlace(iech1, _p1);
    else
      _p1 = _p0;
    SpacePoint p2;
    if (iech2 >= 0)
      _dbin->getSampleCoordinatesAsSPInPlace(iech2, _p2);
    else
      _p2 = _p0;
    _model->evalMatInPlace(_p1, _p2, _covtab, mode);
  }
}

void KrigingSystem::_covCvvCalcul(const CovCalcMode* mode)
{
  VectorVectorDouble d1 = _getDISC1s();
  VectorVectorDouble d2 = _getDISC2s();

  for (int ivar = 0; ivar < _nvar; ivar++)
    for (int jvar = 0; jvar < _nvar; jvar++)
      _covtab(ivar,jvar) += _model->evalAverageIncrToIncr(d1, d2, ivar, jvar, mode);
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

int KrigingSystem::_drftabCalcul(const ECalcMember &member, int iech)
{
  const Db* db;
  int jech;
  if (iech >= 0)
  {
    db = _dbin;
    jech = iech;
  }
  else
  {
    db = _dbout;
    jech = _iechOut;
  }
  VectorDouble drft = _model->evalDriftVec(db, jech, member);
  for (int il = 0; il < (int) drft.size(); il++)
  {
    if (FFFF(drft[il])) return 1;
    _drftab[il] = drft[il];
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Establish the kriging L.H.S.
 **
 *****************************************************************************/
void KrigingSystem::_lhsCalcul()
{
  _lhs.fill(0.);

  /* Establish the covariance part */

  for (int iech = 0; iech < _nech; iech++)
  {
    for (int jech = 0; jech < _nech; jech++)
    {
      _covtab.fill(0.);
      if (_model->isNoStat()) _covUpdate(1, _nbgh[iech], 1, _nbgh[jech]);
      if (iech == jech && _model->isStationary())
        _covtab0Calcul(&_calcModeLHS);
      else
        _covtabCalcul(_nbgh[iech], _nbgh[jech], &_calcModeLHS);

      for (int ivar = 0; ivar < _nvar; ivar++)
        for (int jvar = 0; jvar < _nvar; jvar++)
        {
          _setLHS(iech,ivar,jech,jvar,_getCOVTAB(ivar, jvar));

          /* Correction due to measurement errors */

          if (_flagVerr)
          {
            double verr = 0.;
            if (_flagCode)
            {
              int code1 = (int) _dbin->getLocVariable(ELoc::C,_nbgh[iech],0);
              int code2 = (int) _dbin->getLocVariable(ELoc::C,_nbgh[jech],0);
              if (code1 != 0 && code2 != 0 && code1 == code2)
                verr = _dbin->getLocVariable(ELoc::V,_nbgh[iech], 0);
            }
            else
            {
              if (iech == jech && ivar == jvar)
              {
                verr = _dbin->getLocVariable(ELoc::V,_nbgh[iech], ivar);

                if (_neigh->getFlagContinuous())
                {
                  // In the case of continuous Kriging, we must update the LHS
                  // by considering the distance between data and target

                  double cref = _getLHS(iech, ivar, jech, jvar);
                  verr = cref * _continuousMultiplier(_nbgh[iech], _iechOut);
                }
              }
            }
            if (!FFFF(verr) && verr > 0) _addLHS(iech,ivar,jech,jvar,verr);
          }
        }
    }
  }

  /* Establish the drift part */

  if (_nfeq <= 0 || _nbfl <= 0) return;
  for (int iech = 0; iech < _nech; iech++)
  {
    (void) _drftabCalcul(ECalcMember::LHS, _nbgh[iech]);
    for (int ivar = 0; ivar < _nvar; ivar++)
      for (int ib = 0; ib < _nfeq; ib++)
      {
        double value = 0.;
        for (int il = 0; il < _nbfl; il++)
          value += _drftab[il] * _getDriftCoef(ivar, il, ib);
        _setLHS(iech,ivar,ib,_nvar,value,true);
        _setLHS(ib,_nvar,iech,ivar,value,true);
      }
  }
  return;
}

void KrigingSystem::_lhsIsoToHetero()
{
  if (_flagIsotopic) return;
  int ecri = 0;
  int ecrj = 0;

  ecri = 0;
  for (int i = 0; i < _neq; i++)
  {
    if (_flag[i] == 0) continue;
    ecrj = 0;
    for (int j = 0; j < _neq; j++)
    {
      if (_flag[j] == 0) continue;
      _lhs.setValue(ecri, ecrj, _lhs.getValue(i,j));
      ecri++;
    }
    ecrj++;
  }
  return;
}

VectorInt KrigingSystem::_getRelativePosition()
{
  VectorInt rel(_neq);
  int j = 0;
  for (int i = 0; i < _neq; i++)
  {
    if (! _flag.empty() && _flag[i])
      rel[j++] = i + 1;
    else
      rel[j++] = i + 1;
  }
  return rel;
}

void KrigingSystem::_lhsDump(int nbypas)
{
  VectorInt rel = _getRelativePosition();
  int npass = (_nred - 1) / nbypas + 1;

  /* General Header */

  mestitle(0, "LHS of Kriging matrix (compressed)");
  if (_nech > 0) message("Number of active samples    = %d\n", _nech);
  message("Total number of equations   = %d\n", _neq);
  message("Reduced number of equations = %d\n", _nred);

  /* Loop on the passes */

  for (int ipass = 0; ipass < npass; ipass++)
  {
    int ideb = ipass * nbypas;
    int ifin = MIN(_nred, ideb + nbypas);
    message("\n");

    /* Header line */

    tab_prints(NULL, "Rank");
    tab_prints(NULL, "    ");
    for (int j = ideb; j < ifin; j++)
      tab_printi(NULL, j + 1);
    message("\n");

    /* Flag line */

    if (! _flag.empty())
    {
      tab_prints(NULL, "    ");
      tab_prints(NULL, "Flag");
      for (int j = ideb; j < ifin; j++)
        tab_printi(NULL, rel[j]);
      message("\n");
    }

    /* Matrix lines */

    for (int i = 0; i < _nred; i++)
    {
      tab_printi(NULL, i + 1);
      tab_printi(NULL, rel[i]);
      for (int j = ideb; j < ifin; j++)
        tab_printg(NULL, _lhs.getValue(i,j));
      message("\n");
    }
  }
  return;
}

int KrigingSystem::_lhsInvert()
{
  // Duplicate the whole direct matrix
  _lhsinv = _lhs;

  /* Invert the L.H.S. matrix */

  int rank = _lhsinv.invert();
  if (rank > 0)
  {
    messerr("When estimating Target Site #%d",_iechOut+1);
    messerr("The Kriging Matrix (%d,%d) is singular - Rank = %d",
            _nred, _nred, rank);
    messerr("One of the usual reason is the presence of duplicates");
    return 1;
  }
  return 0;
}

void KrigingSystem::_rhsStore(int iech)
{
  if (_flagNoMatCL)
  {
    for (int ivar = 0; ivar < _nvar; ivar++)
      for (int jvar = 0; jvar < _nvar; jvar++)
        _setRHS(iech, ivar, jvar, _getCOVTAB(ivar, jvar));
  }
  else
  {
    for (int jvarCL = 0; jvarCL < _nvarCL; jvarCL++)
      for (int ivar = 0; ivar < _nvar; ivar++)
      {
        double value = 0.;
        for (int jvar = 0; jvar < _nvar; jvar++)
          value += _matCL[jvarCL][jvar] * _getCOVTAB(ivar, jvar);
        _setRHS(iech, ivar, jvarCL, value);
      }
  }
}

/****************************************************************************/
/*!
 **  Establish the kriging R.H.S (for Point Estimation)
 **
 *****************************************************************************/
void KrigingSystem::_rhsCalculPoint()
{
  if (_optimEnabled)
    _model->getCovAnisoList()->optimizationSetTarget(_p0);

  for (int iech = 0; iech < _nech; iech++)
  {
    if (_model->isNoStat()) _covUpdate(1, _nbgh[iech], 2, _iechOut);

    _covtab.fill(0.);
    _covtabCalcul(_nbgh[iech], -1, &_calcModeRHS);
    _rhsStore(iech);
  }
}

/****************************************************************************/
/*!
 **  Establish the kriging R.H.S (block)
 **
 *****************************************************************************/
void KrigingSystem::_rhsCalculBlock()
{
  // The Block calculation needs:
  // - to memorize the location of the target (center) before its randomization
  // - define a new matrix to cumulate '_covtab'  calculations
  _p0_memo = _p0;
  MatrixSquareGeneral covcum(_covtab);

  for (int iech = 0; iech < _nech; iech++)
  {

    if (_model->isNoStat()) _covUpdate(1, _nbgh[iech], 2, _iechOut);

    covcum.fill(0.);
    if (_flagPerCell) _blockDiscretize();
    int nscale = _getNDisc();

    for (int i = 0; i < nscale; i++)
    {
      // calculate the Local covariance between data and randomized target
      _p0 = _p0_memo;
      _p0.move(_getDISC1Vec(i));
      if (_optimEnabled)
        _model->getCovAnisoList()->optimizationSetTarget(_p0);
      _covtabCalcul(_nbgh[iech], -1, &_calcModeRHS);

      // Cumulate the Local covariance to '_covtab'
      covcum.addMatrix(_covtab);
    }

    // Normalization
    _covtab.copyElements(covcum);
    if (nscale > 1)
      _covtab.prodScalar(1. / (double) nscale);

    _rhsStore(iech);
  }
}

/****************************************************************************/
/*!
 **  Establish the kriging R.H.S (Drift case)
 **
 *****************************************************************************/
void KrigingSystem::_rhsCalculDrift()
{
  if (_optimEnabled)
    _model->getCovAnisoList()->optimizationSetTarget(_p0);

  for (int iech = 0; iech < _nech; iech++)
  {
    _covtab.fill(0.);
    _rhsStore(iech);
  }
}

/****************************************************************************/
/*!
 **  Establish the kriging R.H.S (DGM case)
 **
 *****************************************************************************/
void KrigingSystem::_rhsCalculDGM()
{
  if (_optimEnabled)
    _model->getCovAnisoList()->optimizationSetTarget(_p0);

  for (int iech = 0; iech < _nech; iech++)
  {
    if (_model->isNoStat()) _covUpdate(1, _nbgh[iech], 2, _iechOut);

    _covtab.fill(0.);
    _covtabCalcul(_nbgh[iech], -1, &_calcModeRHS);
    _rhsStore(iech);
  }
}

/****************************************************************************/
/*!
 **  Establish the kriging R.H.S
 **
 ** \remarks When 'matCL' is provided, 'nvar' stands for the first dimension of
 ** \remarks the matrix 'matCL' (its second dimension is equal to model->getNVar()).
 ** \remarks Otherwise nvar designates model->getNVar()
 **
 *****************************************************************************/
int KrigingSystem::_rhsCalcul()
{
  _dbout->getSampleCoordinatesAsSPInPlace(_iechOut, _p0);

  /* Establish the covariance part */

  switch (_calcul.toEnum())
  {
    case EKrigOpt::E_POINT:
    {
      _rhsCalculPoint();
      break;
    }

    case EKrigOpt::E_BLOCK:
    {
      _rhsCalculBlock();
      break;
    }

    case EKrigOpt::E_DRIFT:
    {
      _rhsCalculDrift();
      break;
    }

    case EKrigOpt::E_DGM:
    {
      _rhsCalculDGM();
      break;
    }
  }

  /* Establish the drift part */

  if (_nfeq <= 0) return 0;

  if (_drftabCalcul(ECalcMember::RHS, -1)) return 1;
  if (_flagNoMatCL)
  {
    for (int ivar = 0; ivar < _nvar; ivar++)
      for (int ib = 0; ib < _nfeq; ib++)
      {
        double value = 0.;
        for (int il = 0; il < _nbfl; il++)
          value += _drftab[il] * _getDriftCoef(ivar, il, ib);
        _setRHS(ib,_nvar,ivar,value,true);
      }
  }
  else
  {
    for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
    {
      int ib = 0;
      for (int jvar = 0; jvar < _nvar; jvar++)
        for (int jl = 0; jl < _nbfl; jl++, ib++)
        {
          double value = 0.;
          for (int il = 0; il < _nbfl; il++)
            value += _drftab[il] * _getDriftCoef(jvar, il, ib);
          value *= _matCL[ivarCL][jvar];
          _setRHS(ib,_nvar,ivarCL,value,true);
        }
    }
  }
  return 0;
}

void KrigingSystem::_rhsIsoToHetero()
{
  if (_flagIsotopic) return;
  int ecr_rhs = 0;
  for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
  {
    ecr_rhs = 0;
    for (int i = 0; i < _neq; i++)
    {
      if (_flag[i] == 0) continue;
      _rhs.setValue(ecr_rhs++, ivarCL, _rhs.getValue(i, ivarCL));
    }
  }
  return;
}

void KrigingSystem::_rhsDump()
{
  VectorInt rel = _getRelativePosition();

  /* General Header */

  mestitle(0, "RHS of Kriging matrix (compressed)");
  if (_nech > 0) message("Number of active samples    = %d\n", _nech);
  message("Total number of equations   = %d\n", _neq);
  message("Reduced number of equations = %d\n", _nred);
  message("Number of right-hand sides  = %d\n", _nvarCL);

  /* Kriging option */

  switch (_calcul.toEnum())
  {
    case EKrigOpt::E_POINT:
      message("Punctual Estimation\n");
      break;

    case EKrigOpt::E_BLOCK:
      message("Block Estimation : Discretization = ");
      for (int idim = 0; idim < _ndim; idim++)
      {
        if (idim != 0) message(" x ");
        message("%d", _ndiscs[idim]);
      }
      message("\n");
      break;

    case EKrigOpt::E_DRIFT:
      message("Drift Estimation\n");
      break;

    case EKrigOpt::E_DGM:
      message("Discrete Gaussian Model\n");
      break;
  }
  message("\n");

  /* Header line */

  tab_prints(NULL, "Rank");
  if (! _flag.empty()) tab_prints(NULL, "Flag");
  for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
    tab_printi(NULL, ivarCL + 1);
  message("\n");

  /* Matrix lines */

  for (int i = 0; i < _nred; i++)
  {
    tab_printi(NULL, i + 1);
    if (! _flag.empty()) tab_printi(NULL, rel[i]);
    for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
      tab_printg(NULL, _rhs.getValue(i,ivarCL));
    message("\n");
  }
  return;
}

void KrigingSystem::_wgtCalcul()
{
  _wgt.prodMatrix(_lhsinv, _rhs);
}

void KrigingSystem::_wgtDump(int status)
{
  int ndisc  = _getNDisc();
  VectorDouble sum(_nvarCL);
  char string[20];

  /* Header */

  mestitle(0, "(Co-) Kriging weights");
  const DbGrid* dbgrid = dynamic_cast<const DbGrid*>(_dbout);

  /* First line */

  tab_prints(NULL, "Rank");
  for (int idim = 0; idim < _ndim; idim++)
  {
    String strloc = getLocatorName(ELoc::X, idim);
    tab_prints(NULL, strloc.c_str());
  }
  if (_dbin->hasLocVariable(ELoc::C)) tab_prints(NULL, "Code");
  if (_dbin->getLocNumber(ELoc::V) > 0)
    tab_prints(NULL, "Err.");
  if (ndisc > 0)
    for (int idim = 0; idim < _ndim; idim++)
    {
      (void) gslSPrintf(string, "Size%d", idim + 1);
      tab_prints(NULL, string);
    }
  tab_prints(NULL, "Data");
  for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
  {
    (void) gslSPrintf(string, "Z%d*", ivarCL + 1);
    tab_prints(NULL, string);
  }
  message("\n");

  /* Display the information and the weights */

  int lec = 0;
  int cumflag = 0;
  for (int jvarCL = 0; jvarCL < _nvarCL; jvarCL++)
  {
    if (_nvarCL > 1) message("Using variable Z%-2d\n", jvarCL + 1);

    /* Loop on the samples */

    for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
      sum[ivarCL] = 0.;
    for (int iech = 0; iech < _nech; iech++, lec++)
    {
      int flag_value = (! _flag.empty()) ? _flag[lec] : 1;
      tab_printi(NULL, iech + 1);
      for (int idim = 0; idim < _ndim; idim++)
        tab_printg(NULL, _getIdim(_nbgh[iech], idim));
      if (_dbin->hasLocVariable(ELoc::C))
        tab_printg(NULL, _dbin->getLocVariable(ELoc::C,_nbgh[iech],0));
      if (_dbin->getLocNumber(ELoc::V) > 0)
        tab_printg(NULL, _getVerr(_nbgh[iech], (_flagCode) ? 0 : jvarCL));
      if (ndisc > 0)
      {
        for (int idim = 0; idim < _ndim; idim++)
          if (! _flagPerCell)
            tab_printg(NULL, dbgrid->getDX(idim));
          else
            tab_printg(NULL, dbgrid->getLocVariable(ELoc::BLEX,_nbgh[iech], idim));
      }
      if (_rankPGS < 0)
        tab_printg(NULL, _getIvar(_nbgh[iech], jvarCL));
      else
        tab_prints(NULL,  "    ");

      for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
      {
        double value = (! _wgt.isEmpty() && status == 0 && flag_value) ?
            _wgt.getValue(cumflag, ivarCL) : TEST;
        if (!FFFF(value)) sum[ivarCL] += value;
        tab_printg(NULL, value);
      }
      if (flag_value) cumflag++;
      message("\n");
    }

    int number = 1 + _ndim + 1;
    if (_dbin->getLocNumber(ELoc::V) > 0) number++;
    if (ndisc > 0) number += _ndim;
    tab_prints(NULL, "Sum of weights", number, EJustify::LEFT);
    for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
    {
      double value = (status == 0) ? sum[ivarCL] : TEST;
      tab_printg(NULL, value);
    }
    message("\n");
  }
  if (_nfeq <= 0) return;

  /* Header */

  mestitle(0, "Drift coefficients");

  /* First line */

  tab_prints(NULL, "Rank");
  tab_prints(NULL, "Lagrange");
  tab_prints(NULL, "Coeff");
  message("\n");

  /* Loop on the drift coefficients */

  cumflag = _nred - _nfeq;
  for (int ib = 0; ib < _nfeq; ib++)
  {
    int iwgt = ib + cumflag;
    tab_printi(NULL, ib + 1);
    tab_printg(NULL, (status == 0) ? _wgt.getValue(iwgt,0) : TEST);
    if (_flagSimu)
      tab_printg(NULL, 0.);
    else
      tab_printg(NULL, (status == 0) ? _zam.getValue(iwgt,0) : TEST);
    message("\n");
  }
  return;
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
      if (_nfeq <= 0) simu = _getMean(ivar);

      if (status == 0)
      {
        if (_flagBayes)
          simu = _model->evalDriftCoef(_dbout, _iechOut, ivar, _postSimu[isimu].data());

        int lec = ivar * _nred;
        for (int jvar = 0; jvar < _nvar; jvar++)
          for (int iech = 0; iech < _nech; iech++)
          {
            int jech = _nbgh[iech];

            double mean = 0.;
            if (_nfeq <= 0) mean = _getMean(jvar);
            if (_flagBayes)
              mean = _model->evalDriftCoef(_dbin, jech, jvar,_postSimu[isimu].data());
            double data = _dbin->getSimvar(ELoc::SIMU, jech, isimu, ivar, _rankPGS, _nbsimu, _nvar);
            simu -= _wgt.getValue(lec++,0) * (data + mean);
          }

        if (OptDbg::query(EDbg::KRIGING))
        {
          double value = _dbout->getArray(_iechOut, _iptrEst + ecr);
          message("Non-conditional simulation #%d = %lf\n", isimu + 1, value);
          message("Kriged difference = %lf\n", -simu);
          message("Conditional simulation #%d = %lf\n", isimu + 1,
                  value + simu);
        }
      }
      else
      {
        // In case of failure with KS, set the conditional simulation to the mean
        if (_nfeq > 0) simu = TEST;
      }

      /* Add the conditioning kriging to the NC simulation at target */
      _dbout->updSimvar(ELoc::SIMU, _iechOut, isimu, ivar, _rankPGS, _nbsimu, _nvar,
                       0, simu);
    }

  return;
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
      double valdat = _dbin->getLocVariable(ELoc::Z,_iechOut, ivarCL);
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
      for (int i = 0; i < _nech; i++)
      {
        if (status != 0) continue;
        double wgt = _wgt.getValue(i,ivarCL);
        int iech = _nbgh[i];
        if (_flagSet)
          _dbin->setArray(iech, _iptrWeights + ivarCL, wgt);
        else
          _dbin->updArray(iech, _iptrWeights, 0, wgt);
      }
    }
  }
  return;
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
  return;
}

void KrigingSystem::_estimateCalculImage(int status)
{
  if (! _flagEst) return;

  const DbGrid* dbgrid = dynamic_cast<const DbGrid*>(_dbout);

  VectorInt indn0(_ndim);
  VectorInt indg0(_ndim);
  VectorInt indnl(_ndim);
  VectorInt indgl(_ndim);

  _dbaux->rankToIndice(_dbaux->getSampleNumber()/2, indn0);
  dbgrid->rankToIndice(_iechOut, indg0);
  VH::subtractInPlace(indg0, indn0);

  for (int ivar = 0; ivar < _nvar; ivar++)
  {
    double estim = 0.;
    if (_nfeq <= 0) estim = _getMean(ivar);

    if (status == 0)
    {
      /* Loop on the neighboring points */

      int ecr = ivar * _neq;
      for (int jvar = 0; jvar < _nvar; jvar++)
      {
        for (int iech = 0; iech < _nech; iech++)
        {
          _dbaux->rankToIndice(_nbgh[iech], indnl);
          for (int idim = 0; idim < _ndim; idim++)
            indgl[idim] = dbgrid->getMirrorIndex(idim, indg0[idim] + indnl[idim]);
          int jech = dbgrid->indiceToRank(indgl);
          double data = dbgrid->getLocVariable(ELoc::Z,jech, jvar);
          if (FFFF(data) || FFFF(estim))
          {
            estim = TEST;
          }
          else
          {
            if (_nfeq <= 0) data -= _getMean(jvar);
            estim += data * _wgt.getValue(ecr++,0);
          }
        }
      }
      _dbout->setArray(_iechOut, _iptrEst + ivar, estim);
    }
    else
    {
      // In case of failure with KS, set the result to mean
      if (_nfeq > 0) estim = TEST;
      _dbout->setArray(_iechOut, _iptrEst + ivar, estim);
    }
  }
  return;
}

void KrigingSystem::_estimateCalculXvalidUnique(int /*status*/)
{
  int iech  = _iechOut;
  int iiech = _getFlagAddress(iech, 0);

  // Do not process as this sample is either masked or its variable undefined
  if (iiech < 0) return;

  double valdat = _dbin->getLocVariable(ELoc::Z,iech, 0);
  if (! FFFF(valdat))
  {
    double variance = 1. / _getLHSINV(iiech, 0, iiech, 0);
    double stdv = sqrt(variance);

    /* Perform the estimation */

    double valest = 0.;
    if (_nfeq <= 0) valest = _getMean(0);
    for (int jech = 0; jech < _dbin->getSampleNumber(); jech++)
    {
      int jjech = _getFlagAddress(jech, 0);
      if (jjech < 0) continue;
      if (iiech != jjech)
        valest -= _getLHSINV(iiech,0,jjech,0) * variance * (_dbin->getLocVariable(ELoc::Z, jech, 0) - _getMean(0));
      jjech++;
    }

    if (_flagEst)
    {
      double valloc = valest;
      if (_xvalidEstim) valloc -= valdat;
      _dbin->setArray(iech, _iptrEst, valloc);
    }
    if (_flagStd)
    {
      if (_xvalidStdev) stdv = (valest - valdat) / stdv;
      _dbin->setArray(iech, _iptrStd, stdv);
    }
    if (_flagVarZ)
    {
      _dbin->setArray(iech, _iptrVarZ, TEST);
    }
  }
}

/****************************************************************************/
/*!
 **  Establish the constant term for the variance calculation
 **
 *****************************************************************************/
void KrigingSystem::_variance0()
{
  _dbout->getSampleCoordinatesAsSPInPlace(_iechOut, _p0);
  if (_optimEnabled)
    _model->getCovAnisoList()->optimizationSetTarget(_p0);

  if (_model->isNoStat()) _covUpdate(2, _iechOut, 2, _iechOut);

  _covtab.fill(0.);
  switch (_calcul.toEnum())
  {
    case EKrigOpt::E_POINT:
      _covtab0Calcul(&_calcModeVAR);
      break;

    case EKrigOpt::E_BLOCK:
      _covCvvCalcul(&_calcModeVAR);
      break;

    case EKrigOpt::E_DRIFT:
      break;

    case EKrigOpt::E_DGM:
      _covtab0Calcul(&_calcModeVAR);
      break;
  }

  /* Storage */

  if (_flagNoMatCL)
  {
    for (int ivar = 0; ivar < _nvar; ivar++)
      for (int jvar = 0; jvar < _nvar; jvar++)
        _setVAR0(ivar, jvar, _getCOVTAB(ivar, jvar));
  }
  else
  {
    for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
      for (int jvarCL = 0; jvarCL < _nvarCL; jvarCL++)
      {
        double value = 0.;
        for (int ivar = 0; ivar < _nvar; ivar++)
          for (int jvar = 0; jvar < _nvar; jvar++)
            value += _matCL[ivarCL][ivar] * _getCOVTAB(ivar, jvar) * _matCL[jvarCL][jvar];
        _setVAR0(ivarCL, jvarCL, value);
      }
  }
  return;
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
  MatrixRectangular estims(_nvarCL, 1);

  // Calculate the solution
  if (status == 0)
    estims.prodTMatrix(_rhs, _zam);

  // Loop for writing the estimation

  for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
  {
    double estim0 = 0.;
    if (_nfeq <= 0)
      estim0 = _getMean(ivarCL);
    if (_flagBayes)
      estim0 = _model->evalDriftCoef(_dbout, _iechOut, ivarCL, _postMean.data());

    if (status == 0)
      _dbout->setArray(_iechOut, _iptrEst + ivarCL, estims.getValue(ivarCL,0) + estim0);
    else
      _dbout->setArray(_iechOut, _iptrEst + ivarCL, TEST);
  }
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
  MatrixRectangular vars(_nvarCL,1);

  // Calculate the solution
  if (status == 0)
    vars.prodTMatrix(_rhs, _wgt);

  // Loop for writing the estimation

  for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
  {
    if (status == 0)
    {
      // The Variance at origin must be updated
      // - in Non-stationary case
      // - for Block Estimation, when the block size if defined per cell
      if (_model->isNoStat() || _flagPerCell) _variance0();

      double var = _getVAR0(ivarCL, ivarCL);
      if (_flagBayes) var += _varCorrec[ivarCL * _nvarCL + ivarCL];

      var -= vars.getValue(ivarCL,0);

      double stdv = 0.;
      if (var > 0)
        stdv = sqrt(var);

      _dbout->setArray(_iechOut, _iptrStd + ivarCL, stdv);
    }
    else
      _dbout->setArray(_iechOut, _iptrStd + ivarCL, TEST);
  }
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
  int cumflag = _nred - _nfeq;

  for (int ivarCL = 0; ivarCL < _nvarCL; ivarCL++)
  {
    if (status == 0)
    {
      double varZ = VH::innerProduct(_rhs.getColumn(ivarCL),
                                     _wgt.getColumn(ivarCL), cumflag);
      _dbout->setArray(_iechOut, _iptrVarZ + ivarCL, varZ);
    }
    else
      _dbout->setArray(_iechOut, _iptrVarZ + ivarCL, TEST);
  }
}

/****************************************************************************/
/*!
 **  Define the array flag[] and the kriging L.H.S.
 **
 ** \return The status for the preparation
 **
 *****************************************************************************/
int KrigingSystem::_prepar()
{
  // Resize the internal working arrays

  _resetMemoryPerNeigh();

  /* Define the array flag */

  _flagDefine();

  /* Check if the number of points is compatible with the model */

  if (! _isAuthorized()) return 1;

  /* Establish the Kriging L.H.S. */

  _lhsCalcul();
  _lhsIsoToHetero();
  if (OptDbg::query(EDbg::KRIGING)) _lhsDump();

  /* Invert the L.H.S. matrix */

  if (_lhsInvert()) return 1;

  return 0;
}

/****************************************************************************/
/*!
 **  Extract the valid data
 **  Operate the product by the inverse covariance matrix
 **
 **  Warning; the array 'zext' is dimensioned to '_neq' instead of 'nred'
 **
 *****************************************************************************/
void KrigingSystem::_dualCalcul()
{
  _zext.fill(0.);

  /* Extract the data */

  int ecr = 0;
  for (int ivar = 0; ivar < _nvar; ivar++)
  {
    for (int iech = 0; iech < _nech; iech++)
    {
      if (! _getFLAG(iech, ivar)) continue;
      double mean = 0.;
      if (_nfeq <= 0)
        mean = _getMean(ivar);
      if (_flagBayes)
        mean = _model->evalDriftCoef(_dbout, _iechOut, ivar, _postMean.data());
      _zext.setValue(ecr, 0, _getIvar(_nbgh[iech], ivar) - mean);
      ecr++;
    }
  }

  /* Operate the product : Z * A-1 */

  _zam.prodMatrix(_lhsinv, _zext);

  /* Operate the product : Z * A-1 * Z */

  if (_flagLTerm)
  {
    MatrixSquareGeneral lterMat(1);
    lterMat.prodTMatrix(_zam, _zext);
    _lterm = lterMat.getValue(0,0);
  }

  // Turn back the flag to OFF in order to avoid provoking
  // the _dualCalcul() calculations again

  _flagDataChanged = false;

  return;
}

/**
 * Performs the last operations before launching the loop on Estimations
 * @return
 */
bool KrigingSystem::isReady()
{
  if (! _isCorrect()) return false;

  // Define the calculation modes
  _calcModeLHS = CovCalcMode(ECalcMember::LHS);
  _calcModeLHS.setActiveCovList(_model->getAllActiveCovList(), true);
  _calcModeRHS = CovCalcMode(ECalcMember::RHS);
  _calcModeRHS.setActiveCovList(_model->getActiveCovList(), _model->isAllActiveCovList());
  _calcModeVAR = CovCalcMode(ECalcMember::VAR);
  _calcModeVAR.setActiveCovList(_model->getActiveCovList(), _model->isAllActiveCovList());

  // Perform some pre-calculation when variance of estimator is requested
  if (_flagStd)
  {
    _iechOut = 0;
    _variance0();
  }

  if (_neigh->getType() == ENeigh::IMAGE)
  {
    // Perform some preliminary work in the case of Image
    const NeighImage* neighI = dynamic_cast<const NeighImage*>(_neigh);
    if (_prepareForImage(neighI)) return false;

    // Prepare the projection of data on different covariances of the Model
    if (_optimEnabled)
    {
      if (_model != nullptr && _dbaux != nullptr)
        _model->getCovAnisoList()->optimizationPreProcess(_dbaux->getSamplesAsSP());
    }

    // Setup the Kriging pre-processings
    if (_prepareForImageKriging(_dbaux, neighI)) return false;
  }
  else
  {
    if (_optimEnabled)
    {
      // Prepare the projection of data on different covariances of the Model
      if (_model != nullptr && _dbin != nullptr)
        _model->getCovAnisoList()->optimizationPreProcess(_dbin->getSamplesAsSP());

      if (_flagBayes && _modelSimple != nullptr && _dbin != nullptr)
        _modelSimple->getCovAnisoList()->optimizationPreProcess(_dbin->getSamplesAsSP());
    }
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
  if (_model != nullptr)
    _model->getCovAnisoList()->optimizationPostProcess();
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

  // In case of Cross-validation in Unique Neighborhood, do not establish the RHS
  bool caseXvalidUnique = false;
  if (_neigh->getType() == ENeigh::UNIQUE &&
      _neigh->getFlagXvalid()) caseXvalidUnique = true;

  // Store the Rank of the Target sample
  _iechOut = iech_out;

  int status = 0;
  if (skipCalculAll) goto label_store;

  if (! _dbout->isActive(_iechOut)) return 0;
  OptDbg::setCurrentIndex(_iechOut + 1);
  if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
  {
    if (_flagFactorKriging)
      message("\nProcessing Factor %d / %d\n",_model->getActiveFactor(), _nclasses);

    mestitle(1, "Target location");
    db_sample_print(_dbout, _iechOut, 1, 0, 0);
  }

  // Elaborate the Neighborhood

  if (caseXvalidUnique) _neigh->setFlagXvalid(false);
  _nbgh = _neigh->select(_iechOut);
  status = _setInternalShortCutVariablesNeigh();
  if (_flagNeighOnly) goto label_store;
  if (status) goto label_store;

  /* Establish the Kriging L.H.S. */

  if (!_neigh->isUnchanged() || _neigh->getFlagContinuous() || OptDbg::force())
  {
    if (_flagBayes) _setLocalModel(_modelSimple);
    status = _prepar();
    if (status) goto label_store;
    if (_flagBayes) _setLocalModel(_modelInit);
  }

  // Establish the pre-calculation involving the data information

  if (!_neigh->isUnchanged() || _neigh->getFlagContinuous() || _flagDataChanged || OptDbg::force())
  {
    _dualCalcul();
  }

  if (caseXvalidUnique) _neigh->setFlagXvalid(true);

  /* Establish the Kriging R.H.S. */

  if (caseXvalidUnique) goto label_store;

  if (_flagBayes) _setLocalModel(_modelSimple);
  _rhsCalcul();
  if (_flagBayes) _setLocalModel(_modelInit);

  if (status != 0) goto label_store;
  _rhsIsoToHetero();
  if (OptDbg::query(EDbg::KRIGING)) _rhsDump();

  /* Derive the kriging weights */

  if (_flagStd || _flagVarZ || _flagSimu || _flagWeights || _flagKeypairWeights)
    _wgtCalcul();
  if (OptDbg::query(EDbg::KRIGING)) _wgtDump(status);

  // Optional Save of the Kriging weights
  if (_flagKeypairWeights) _saveWeights(status);

  /* Perform the final estimation */

  label_store:
  // If status is not zero, cancel the current Neighborhood search status
  if (status) _neigh->setIsChanged();

  // Correct the Variance in Bayesian case

  if (_flagBayes) _bayesCorrectVariance();

  // Store the results in the output Db

  if (_flagNeighOnly)
  {
    VectorDouble tab = _neigh->summary(_iechOut);
    _neighCalcul(status, tab);
  }
  else if (_neigh->getType() == ENeigh::IMAGE)
  {
    // Image Neighborhood case

    _estimateCalculImage(status);
  }
  else if (_neigh->getType() == ENeigh::UNIQUE)
  {
    // Unique Neighborhood case

    if (_neigh->getFlagXvalid())
      _estimateCalculXvalidUnique(status);

    else if (_flagSimu)
      _simulateCalcul(status);

    else if (! _flagGlobal)
      _estimateCalcul(status);
  }
  else
  {
    // Moving Neighborhood case

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
      _simulateDump(status);
    else
      _krigingDump(status);
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Print the results
 **
 ** \param[in] status   Kriging error status
 **
 *****************************************************************************/
void KrigingSystem::_krigingDump(int status)
{
  if (_neigh->getFlagXvalid())
    mestitle(0, "Cross-validation results");
  else
    mestitle(0, "(Co-) Kriging results");
  message("Target Sample = %d\n", _iechOut + 1);

  /* Loop on the results */

  for (int ivar = 0; ivar < _nvar; ivar++)
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
        double trueval = (status == 0) ? _dbin->getLocVariable(ELoc::Z,_iechOut, ivar) : TEST;
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
        value = (status == 0) ? _var0(ivar, ivar) : TEST;
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
  return;
}

void KrigingSystem::_simulateDump(int status)
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
  return;
}

/**
 * Set the calculation options
 * @param iptrEst  UID for storing the estimation(s)
 * @param iptrStd  UID for storing the Standard deviations(s)
 * @param iptrVarZ UID for storing the Variance(s) of estimator
 * @return
 * @remark If a term must not be calculated, its UID must be negative
 */
int KrigingSystem::updKrigOptEstim(int iptrEst, int iptrStd, int iptrVarZ)
{
  _iptrEst = iptrEst;
  _iptrStd = iptrStd;
  _iptrVarZ = iptrVarZ;

  _flagEst = _iptrEst >= 0 || (_neigh->getFlagXvalid() && _iptrStd >= 0);
  _flagStd = (_iptrStd >= 0);
  _flagVarZ = (_iptrVarZ >= 0);

  _flagDataChanged = true;

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

void KrigingSystem::_blockDiscretize()
{
  const DbGrid* dbgrid = dynamic_cast<const DbGrid*>(_dbout);
  _disc1 = dbgrid->getDiscretizedBlock(_ndiscs, _iechOut, _flagPerCell, false, 1234546);
  _disc2 = dbgrid->getDiscretizedBlock(_ndiscs, _iechOut, _flagPerCell, true, 1234546);
}

int KrigingSystem::setKrigOptCalcul(const EKrigOpt& calcul,
                                    const VectorInt& ndiscs,
                                    bool flag_per_cell)
{
  _isReady = false;
  _calcul = calcul;
  _flagPerCell = false;
  if (_calcul == EKrigOpt::BLOCK)
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(_dbout);
    if (dbgrid == nullptr)
    {
      messerr("Block Estimation is only possible for Grid '_dbout'");
      return 1;
    }

    // Block support is defined per sample
    if (flag_per_cell)
    {
      _flagPerCell = true;
    }
    if (_neigh->getType() == ENeigh::CELL)
    {
      _flagPerCell = true;
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
    _ndiscNumber = VH::product(_ndiscs);
    _disc1.resize(_ndiscNumber);
    _disc2.resize(_ndiscNumber);
    for (int i = 0; i < _ndiscNumber; i++)
    {
      _disc1[i].resize(getNDim());
      _disc2[i].resize(getNDim());
    }

    // For constant discretization, calculate discretization coordinates

    if (! _flagPerCell) _blockDiscretize();
  }
  else
  {
    _ndiscs.clear();
    _ndiscNumber = 0;
  }
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
 ** \remarks In input, the numbering in ; rank_colcok' starts from 1
 ** \remarks In output, the numbering starts from 0
 **
 *****************************************************************************/
int KrigingSystem::setKrigOptColCok(const VectorInt& rank_colcok)
{
  if (rank_colcok.empty()) return 0;

  _isReady = false;
  _rankColCok = rank_colcok;
  int nvar = _getNVar();
  _neigh->setRankColCok(rank_colcok);

  /* Loop on the ranks of the colocated variables */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    int jvar = rank_colcok[ivar];
    if (IFFFF(jvar)) continue;
    if (jvar > _dbout->getColumnNumber())
    {
      messerr("Error in the Colocation array:");
      messerr("Input variable (#%d): rank of the colocated variable is %d",
              ivar + 1, jvar);
      messerr("But the Output file only contains %d attributes(s)",
              _dbout->getColumnNumber());
      return (1);
    }
  }
  return 0;
}

int KrigingSystem::setKrigOptBayes(bool flag_bayes,
                                   const VectorDouble& prior_mean,
                                   const VectorDouble& prior_cov,
                                   int seed)
{
  _isReady = false;
  int nfeq = _getNFeq();
  if (flag_bayes)
  {
    VectorDouble local_mean = prior_mean;
    VectorDouble local_cov  = prior_cov;

    if (local_mean.empty())
      local_mean.resize(nfeq, 0.);
    if (local_cov.empty())
    {
      local_cov.resize(nfeq * nfeq);
      for (int i = 0; i < nfeq; i++)
        for (int j = 0; j < nfeq; j++)
          local_cov[i + j*nfeq] = (i == j) ? 1. : 0.;
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
    if (!is_matrix_symmetric(nfeq, local_cov.data(), 1))
    {
      messerr("Argument 'prior_cov' should be a symmetric matrix");
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
    _varCorrec.resize(nfeq * nfeq);

    // Duplicate the Model and suppress any Drift component

    _modelSimple = _modelInit->clone();
    _modelSimple->delAllDrifts();
  }
  _flagBayes = flag_bayes;
  _seedForBayes = seed;
  return 0;
}

int KrigingSystem::setKrigOptImage(int seed)
{
  _seedForImage = seed;
  return 0;
}

/**
 * Define the output as Linear Combinations of the Input Variables
 * @param matCL Vector of Vectors of weights (see remarks)
 * @return
 * @remarks The first dimension of 'matCL' is the number of Output variables
 * @remarks The second dimension is the number of input Variables.
 */
int KrigingSystem::setKrigOptMatCL(const VectorVectorDouble& matCL)
{
  if (_flagNoMatCL) return 0;
  _isReady = false;
  int n1 = (int) matCL.size();
  int n2 = (int) matCL[0].size();

  if (n1 > _getNVar())
  {
    messerr("First dimension of 'matCL' (%d)",(int) n1);
    messerr("should be smaller than the number of variables in the model (%d)",
            _getNVar());
    return 1;
  }
  if (n2 != _getNVar())
  {
    messerr("Second dimension of 'matCL' (%d)",(int) n2);
    messerr("should be equal to the number of variables in the model (%d)",
            _getNVar());
    return 1;
  }
  _matCL = matCL;
  return 0;
}

/**
 * Define the option for Kriging By Profile
 * @param flag_code True if the option is switched ON
 * @return
 *
 * @remark When the Option for Kriging By Profile is ON, the covariance term in LHS
 * @remark is incremented by V(iech) when both samples 'iech' and 'jech' have
 * @remark the same Code.
 */
int KrigingSystem::setKrigoptCode(bool flag_code)
{
  _isReady = false;
  if (flag_code)
  {
    if (! _dbin->hasLocVariable(ELoc::C) || _dbin->getLocNumber(ELoc::V) != 1)
    {
      messerr("This method requires variables CODE and V to be defined");
      return 1;
    }
  }
  _flagCode = flag_code;
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

/**
 * Switch the option for saving the Kriging Weights using Keypair mechanism
 * @param flag_save Value of the switch
 */
int KrigingSystem::setKrigOptSaveWeights(bool flag_save)
{
  _isReady = false;
  _flagKeypairWeights = flag_save;
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

  if (_model->getCovaMinIRFOrder() != -1)
  {
    messerr("The option DGM is limited to Stationary Covariances");
    return 1;
  }
  if (_model->getVariableNumber() != 1)
  {
    messerr("The DGM option is limited to the Monovariate case");
    return 1;
  }
  if (ABS(_model->getTotalSill(0,0) - 1.) > eps)
  {
    messerr("The DGM option requires a Model with Total Sill equal to 1.");
    return 1;
  }
  _flagDGM = flag_dgm;
  return 0;
}

int KrigingSystem::setKrigOptFlagGlobal(bool flag_global)
{
  _isReady = false;
  if (flag_global)
  {
    if (! _dbout->isGrid())
    {
      messerr("The Global Option needs '_dbout' to be a Grid");
      return 1;
    }
    if (_neigh->getType() != ENeigh::UNIQUE)
    {
      messerr("The Global Option requires Unique Neighborhood");
      return 1;
    }
  }
  _flagGlobal = flag_global;
  return 0;
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
  double total = _model->getTotalSill(0, 0);
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
    if (! _model->hasAnam())
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
  _model->setActiveFactor(index_class);
  _nclasses = nclasses;

  // Update C00 if the variance calculation is required
  if (_flagStd) _variance0();

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
    if (ndim > 0 && ndim != _model->getDimensionNumber())
    {
      messerr("Incompatible Space Dimension of '_ model'");
      return false;
    }
    ndim = _model->getDimensionNumber();
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
    if (nvar > 0 && nvar != _dbin->getLocNumber(ELoc::Z))
    {
      messerr("Incompatible Variable Number of '_dbin'");
      return false;
    }
    nvar = _dbin->getLocNumber(ELoc::Z);
  }
  if (_model != nullptr)
  {
    if (nvar > 0 && nvar != _model->getVariableNumber())
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
    if (_model->getCovaNumber() <= 0)
    {
      messerr("The Model should contain some Covariances defined before Kriging");
      return false;
    }
  }

  /**************************************/
  /* Checking the Validity of the Model */
  /**************************************/

  if (! _model->isValid()) return false;

  /******************************************/
  /* Checking the Number of External Drifts */
  /******************************************/

  int nfex = 0;
  if (_model != nullptr)
  {
    if (nfex > 0 && nfex != _model->getExternalDriftNumber())
    {
      messerr("Incompatible Number of External Drifts of '_model'");
      return false;
    }
    nfex = _model->getExternalDriftNumber();
  }
  if (nfex > 0)
  {
    // When External Drifts are defined in the Model,
    // the other files must be consistent

    if (_dbout != nullptr)
    {
      if (nfex != _dbout->getLocNumber(ELoc::F))
      {
        messerr("Incompatible Number of External Drifts:");
        messerr("- In 'Model' = %d", nfex);
        messerr("- In '_dbout' = %d", _dbout->getLocNumber(ELoc::F));
        return false;
      }
    }
    if (_dbin != nullptr)
    {
      if (_dbin->getLocNumber(ELoc::F) == 0)
      {
        if (migrateByLocator(_dbout, _dbin, ELoc::F)) return false;
        // Store the UID of the newly created variables to be deleted at the end of the process
        _dbinUidToBeDeleted = _dbin->getUIDsByLocator(ELoc::F);
      }
      if (nfex != _dbin->getLocNumber(ELoc::F))
      {
        messerr("Incompatible Number of External Drifts:");
        messerr("- In 'Model' = %d", nfex);
        messerr("- In 'dbin' = %d", _dbin->getLocNumber(ELoc::F));
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
      db_extension(_dbin, db_mini, db_maxi, true);

    /* Output Db structure */

    if (_dbout != nullptr)
      db_extension(_dbout, db_mini, db_maxi, true);

    // Merge extensions

    _model->setField(VH::extensionDiagonal(db_mini, db_maxi));
  }

  /*****************************/
  /* Checking the Neighborhood */
  /*****************************/

  if (_neigh != nullptr)
  {
    if (_neigh->getType() == ENeigh::IMAGE && (_dbout == nullptr || ! _dbout->isGrid()))
    {
      messerr("The Image neighborhood can only be used when the output Db is a grid");
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
    if (_model->isNoStat())
    {
      const ANoStat *nostat = _model->getNoStat();

      if (_dbin != nullptr)
      {
        // Attach the Input Db
        if (nostat->attachToDb(_dbin, 1)) return false;
      }

      if (_dbout != nullptr)
      {
        // Attach the Output Db
        if (nostat->attachToDb(_dbout, 2)) return false;
      }

      // Discard optimization in the non-stationary case
      _model->setOptimEnabled(false);
      _optimEnabled = false;
    }
    else
    {
      _optimEnabled = _model->isOptimEnabled();
    }
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

// Shortcuts for addressing functions

int KrigingSystem::_getFLAG(int iech, int ivar) const
{
  int ind = IND(iech, ivar);
  return _flag[ind];
}
double KrigingSystem::_getCOVTAB(int ivar,int jvar) const
{
  return _covtab.getValue(ivar, jvar);
}
double KrigingSystem::_getRHS(int iech, int ivar, int jvCL) const
{
  int ind  = IND(iech, ivar);
  return _rhs.getValue(ind, jvCL);
}
void KrigingSystem::_setRHS(int iech, int ivar, int jvCL, double value, bool isForDrift)
{
  int ind  = IND(iech,ivar);
  _rhs.setValue(ind, jvCL, value);
}
double KrigingSystem::_getRHSC(int i, int jvCL) const
{
  return _rhs.getValue(i, jvCL);
}
double KrigingSystem::_getWGTC(int i,int jvCL) const
{
  return _wgt.getValue(i, jvCL);
}

/**
 * Setting the LHS element
 * @param iech Rank of the first sample
 * @param ivar Rank of the first variable
 * @param jech Rank of the second sample
 * @param jvar Rank of the second variable
 * @param value Assigned value
 * @param isForDrift True is used to set the drift element
 *
 * @remark: When used for setting a Drift element (say in 'i'), then:
 * @remark: 'iech' is set for 'ib' (which must be within [0,nfeq[) and 'ivar' is set to 'nvar'
 */
void KrigingSystem::_setLHS(int iech, int ivar, int jech, int jvar, double value, bool isForDrift)
{
  int indi = IND(iech, ivar);
  int indj = IND(jech, jvar);
  _lhs.setValue(indi, indj, value);
}
void KrigingSystem::_addLHS(int iech, int ivar, int jech, int jvar, double value)
{
  int indi = IND(iech, ivar);
  int indj = IND(jech, jvar);
  _lhs.setValue(indi, indj, _lhs.getValue(indi, indj) + value);
}
double KrigingSystem::_getLHS(int iech, int ivar, int jech, int jvar) const
{
  int indi = IND(iech, ivar);
  int indj = IND(jech, jvar);
  return _lhs.getValue(indi, indj);
}
double KrigingSystem::_getLHSC(int i, int j) const
{
  return _lhs.getValue(i,j);
}
double KrigingSystem::_getLHSINV(int iech, int ivar, int jech, int jvar) const
{
  int indi = IND(iech, ivar);
  int indj = IND(jech, jvar);
  return _lhsinv.getValue(indi, indj);
}
double KrigingSystem::_getDISC1(int idisc, int idim) const
{
  return _disc1[idisc][idim];
}
VectorDouble KrigingSystem::_getDISC1Vec(int idisc) const
{
  VectorDouble vec(_ndim);
  for (int idim = 0; idim < _ndim; idim++)
    vec[idim] = _disc1[idisc][idim];
  return vec;
}
double KrigingSystem::_getDISC2(int idisc,int idim) const
{
  return _disc2[idisc][idim];
}
VectorDouble KrigingSystem::_getDISC2Vec(int idisc) const
{
  VectorDouble vec(_ndim);
  for (int idim = 0; idim < _ndim; idim++)
    vec[idim] = _disc2[idisc][idim];
  return vec;
}
VectorVectorDouble KrigingSystem::_getDISC1s() const
{
  VectorVectorDouble vecvec(_ndiscNumber);
  for (int idisc = 0; idisc < _ndiscNumber; idisc++)
    vecvec[idisc] = _getDISC1Vec(idisc);
  return vecvec;
}
VectorVectorDouble KrigingSystem::_getDISC2s() const
{
  VectorVectorDouble vecvec(_ndiscNumber);
  for (int idisc = 0; idisc < _ndiscNumber; idisc++)
    vecvec[idisc] = _getDISC2Vec(idisc);
  return vecvec;
}
double KrigingSystem::_getVAR0(int ivCL, int jvCL) const
{
  return _var0(ivCL, jvCL);
}
void KrigingSystem::_setVAR0(int ivCL, int jvCL, double value)
{
  _var0(ivCL, jvCL) = value;
}

void KrigingSystem::_checkAddress(const String& title,
                                  const String& theme,
                                  int ival,
                                  int nval) const
{
  if (ival < 0)
    messageAbort("Error in %s: %s (%d) may not be negative",
                 title.c_str(),theme.c_str(),ival);
  if (ival >= nval)
    messageAbort("Error in %s: %s (%d) may not be larger or equal than %d",
                 title.c_str(),theme.c_str(),ival,nval);
}

bool KrigingSystem::_prepareForImage(const NeighImage* neighI)
{
  if (! _dbout->isGrid()) return 1;
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(_dbout);
  double seuil = 1. / neighI->getSkip();

  /* Core allocation */

  VectorInt nx(_ndim);
  int nech = 1;
  for (int i=0; i<_ndim; i++)
  {
    nx[i] = 2 * neighI->getImageRadius(i) + 1;
    nech *= nx[i];
  }

  law_set_random_seed(_seedForImage);
  VectorBool sel(nech);
  for (int iech = 0; iech < nech; iech++)
     sel[iech] = (law_uniform(0., 1.) < seuil) ? true : false;

  VectorDouble tab(nech * _nvar);
  int iecr = 0;
    for (int ivar = 0; ivar < _nvar; ivar++)
      for (int iech = 0; iech < nech; iech++)
        tab[iecr++] = (sel[iech]) ? 0. : TEST;

  /* Create the grid */

  _dbaux = DbGrid::create(nx, dbgrid->getDXs(), dbgrid->getX0s(), dbgrid->getAngles());
  _dbaux->addColumns(tab, "Test", ELoc::Z);

  /* Shift the origin */

  VectorDouble coor(_ndim);
  _dbaux->rankToCoordinatesInPlace(nech/2, coor);
  for (int i=0; i<_ndim; i++) _dbaux->setX0(i, _dbaux->getX0(i) - coor[i]);
  if (db_grid_define_coordinates(_dbaux)) return 1;

  return 0;
}

bool KrigingSystem::_prepareForImageKriging(Db* dbaux, const NeighImage* neighI)
{
  // Save pointers to previous Data Base (must be restored at the end)
  Db* dbin_loc  = _dbin;
  Db* dbout_loc = _dbout;

  _dbin  = dbaux;
  _dbout = dbaux;
  int error = 1;

  /* Prepare the neighborhood (mimicking the Unique neighborhood) */

  SpaceRN space(_ndim);
  NeighUnique neighU = NeighUnique::create(false, &space);
  neighU.attach(dbaux, dbaux);

  _iechOut = dbaux->getSampleNumber() / 2;
  _nbgh = neighU.select(_iechOut);
  bool status = _setInternalShortCutVariablesNeigh();

  /* Establish the L.H.S. */

  status = _prepar();
  if (status) goto label_end;
  _dualCalcul();

  /* Establish the R.H.S. */

  _rhsCalcul();
  if (status != 0) goto label_end;
  _rhsIsoToHetero();
  if (OptDbg::query(EDbg::KRIGING)) _rhsDump();

  /* Derive the Kriging weights (always necessary) */

  _wgtCalcul();
  if (OptDbg::query(EDbg::KRIGING)) _wgtDump(status);

  error = 0;

  label_end:
  _dbin  = dbin_loc;
  _dbout = dbout_loc;
  return error;
}

void KrigingSystem::_saveWeights(int status)
{
  if (status != 0) return;
  VectorDouble values(5);
  values[0] = _iechOut;

  /* Loop on the output variables */

  int lec = 0;
  int cumflag = 0;
  for (int jvar = 0; jvar < _nvar; jvar++)
  {
    values[1] = jvar;

    /* Loop on the input samples */

    for (int iech = 0; iech < _nech; iech++, lec++)
    {
      int flag_value = (! _flag.empty()) ? _flag[lec] : 1;
      if (flag_value)
      {
        values[2] = _nbgh[iech];

        /* Loop on the input variables */

        for (int ivar = 0; ivar < _nvar; ivar++)
        {
          int iwgt = _nred * ivar + cumflag;
          double wgtloc = (! _wgt.isEmpty() && flag_value) ? _wgt.getValue(iwgt,0) : TEST;
          if (!FFFF(wgtloc))
          {
            values[3] = ivar;
            values[4] = wgtloc;
            app_keypair("KrigingWeights", 1, 1, 5, values.data());
          }
        }
        if (flag_value) cumflag++;
      }
    }
  }
  return;
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
      xyz[idim][iech] = _getIdim(_nbgh[iech], idim);
  }
  return xyz;
}

VectorDouble KrigingSystem::getSampleData() const
{
  VectorDouble zval(_nred, 0.);

  /* Extract the data */

  int ecr = 0;
  for (int ivar = 0; ivar < _nvar; ivar++)
  {
    for (int iech = 0; iech < _nech; iech++)
    {
      if (!_getFLAG(iech, ivar)) continue;
      zval[ecr++] = _getIvar(_nbgh[iech], ivar);
    }
  }
  return zval;
}

VectorDouble KrigingSystem::getRHSC(int ivar) const
{
  return _rhs.getColumn(ivar);
}

VectorDouble KrigingSystem::getZamC() const
{
  return _zam.getColumn(0);
}

/****************************************************************************/
/*!
 **  Perform the Bayesian estimation of the Drift Coefficients
 **  (only in Unique Neighborhood)
 **
 *****************************************************************************/
int KrigingSystem::_bayesPreCalculations()
{
  _iechOut = _dbin->getSampleNumber() / 2;

  // Elaborate the (Unique) Neighborhood
  _nbgh = _neigh->select(_iechOut);
  if (_setInternalShortCutVariablesNeigh()) return 1;

  /* Prepare the Kriging matrix (without correction) */

  _flagBayes = false;
  if (_prepar()) return 1;
  _flagBayes = true;
  int shift = _nred - _nfeq;

  /* Complementary core allocation */

  VectorDouble ff(shift * _nfeq);
  VectorDouble smu(_nfeq);
  VectorDouble sigma(shift * shift);
  VectorDouble vars(shift);

  // Create the array of variables

  int ind = 0;
  for (int iech = 0; iech < _dbin->getSampleNumber(); iech++)
  {
    if (! _dbin->isActive(iech)) continue;
    for (int ivar = 0; ivar < _nvar; ivar++)
    {
      double value = _dbin->getLocVariable(ELoc::Z,_nbgh[iech], ivar);
      if (FFFF(value)) continue;
      vars[ind++] = value;
    }
  }

  /* Copy DCOV into S and DMEAN into RMEAN */

  _postCov  = _priorCov;
  _postMean = _priorMean;

  /* Establish the drift array FF */

  for (int il = 0; il < _nfeq; il++)
    for (int ib = 0; ib < shift; ib++)
      FF(ib,il) = _getLHSC(ib, shift + il);

  /* Calculate S-1 */

  if (matrix_invert(_postCov.data(), _nfeq, -1)) return 1;

  /* Calculate: SMU = S-1 * MEAN */

  matrix_product_safe(_nfeq, _nfeq, 1, _postCov.data(), _priorMean.data(), smu.data());

  /* Covariance matrix SIGMA */

  for (int ib = 0; ib < shift; ib++)
    for (int jb = 0; jb < shift; jb++)
      SIGMA(ib,jb) = _getLHSC(ib, jb);

  /* Calculate SIGMA-1 */

  if (matrix_invert(sigma.data(), shift, -1)) return 1;

  /* Inverse of posterior covariance matrix: SC-1 = FFt * SIGMA-1 * FF + S-1 */

  int ecr = 0;
  for (int il = 0; il < _nfeq; il++)
    for (int jl = 0; jl < _nfeq; jl++)
    {
      double value = 0.;
      for (int ib=0; ib<shift; ib++)
        for (int jb=0; jb<shift; jb++)
          value += FF(ib,il) * SIGMA(ib,jb) * FF(jb,jl);
      _postCov[ecr++] += value;
    }

  /* Calculating: SMU = FFt * SIGMA-1 * Z + S-1 * MU */

  for (int il = 0; il < _nfeq; il++)
  {
    double value = 0.;
    for (int ib=0; ib<shift; ib++)
      for (int jb=0; jb<shift; jb++)
        value += FF(ib,il) * SIGMA(ib,jb) * vars[jb];
    smu[il] += value;
  }

  /* Posterior mean: _postMean = SC * SMU */

  if (matrix_invert(_postCov.data(), _nfeq, -1)) return 1;
  matrix_product_safe(_nfeq, _nfeq, 1, _postCov.data(), smu.data(), _postMean.data());

  if (OptDbg::query(EDbg::BAYES))
  {
    mestitle(0, "Bayesian Drift coefficients");
    print_matrix("Prior Mean", 0, 1, _nfeq, 1, NULL, _priorMean.data());
    print_matrix("Prior Variance-Covariance", 0, 1, _nfeq, _nfeq, NULL,
                 _priorCov.data());
    print_matrix("Posterior Mean", 0, 1, _nfeq, 1, NULL, _postMean.data());
    print_matrix("Posterior Variance-Covariance", 0, 1, _nfeq, _nfeq, NULL,
                 _postCov.data());
    message("\n");
  }

  // Particular case of Simulation: Simulate several outcomes for posterior means

  if (_flagSimu) _bayesPreSimulate();

  return 0;
}

/****************************************************************************/
/*!
 **  Correct the array of Variance in Bayesian case
 **
 *****************************************************************************/
void KrigingSystem::_bayesCorrectVariance()
{
  /* Establish the Drift matrix */

  if (_nbfl <= 0 || _nfeq <= 0) return;

  VectorDouble ff0(_nvar * _nfeq);
  for (int ivar = 0; ivar < _nvar; ivar++)
    for (int ib = 0; ib < _nfeq; ib++)
    {
      double value = 0.;
      for (int il = 0; il < _nbfl; il++)
        value += _drftab[il] * _getDriftCoef(ivar, il, ib);
      FF0(ib,ivar) = value;
    }

  /* Correct the arrays */

  int ecr = 0;
  for (int ivar = 0; ivar < _nvar; ivar++)
    for (int jvar = 0; jvar < _nvar; jvar++)
    {
      double value = 0.;
      int lec = 0;
      for (int il = 0; il < _nfeq; il++)
        for (int jl = 0; jl < _nfeq; jl++)
        {
          value += FF0(il,ivar) * _postCov[lec] * FF0(jl, jvar);
          lec++;
        }
      _varCorrec[ecr] = value;
      ecr++;
    }
  return;
}

/****************************************************************************/
/*!
 **  Simulate the drift coefficients from the posterior distributions
 **
 *****************************************************************************/
void KrigingSystem::_bayesPreSimulate()
{
  if (_nfeq <= 0) return;
  int nftri = _nfeq * (_nfeq + 1) / 2;
  int memo = law_get_random_seed();

  // Dimension '_smean' to store simulated posterior mean
  _postSimu.resize(_nbsimu);
  for (int isimu = 0; isimu < _nbsimu; isimu++) _postSimu[isimu].resize(_nfeq);

  /* Core allocation */

  VectorDouble trimat(nftri);
  VectorDouble rndmat(_nfeq);
  VectorDouble simu(_nfeq);
  // The array _postCov is duplicated as the copy is destroyed by matrix_cholesky_decompose
  VectorDouble rcov = _postCov;

  /* Cholesky decomposition */

  int rank = matrix_cholesky_decompose(rcov.data(), trimat.data(), _nfeq);

  if (rank > 0)
  {
    messerr("Error in the Cholesky Decomposition of the covariance matrix");
    messerr("Rank of the Matrix = %d", rank);
    messerr("The Drift coefficients have been set to their posterior mean");
    for (int isimu = 0; isimu < _nbsimu; isimu++)
      for (int il = 0; il < _nfeq; il++)
        SMEAN(il,isimu) = _postMean[il];
  }
  else
  {
    for (int isimu = 0; isimu < _nbsimu; isimu++)
    {

      /* Draw a vector of gaussian independent values */

      for (int il = 0; il < _nfeq; il++) rndmat[il] = law_gaussian();

      /* Product of the Lower triangular matrix by the random vector */

      matrix_cholesky_product(1, _nfeq, 1, trimat.data(), rndmat.data(), simu.data());

      /* Add the mean */

      for (int il = 0; il < _nfeq; il++)
        SMEAN(il,isimu) = simu[il] + _postMean[il];
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
        message(" %lf", SMEAN(il, isimu));
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

int KrigingSystem::_getFlagAddress(int iech0, int ivar0)
{
  int rank = 0;
  for (int ivar = 0; ivar < _nvar; ivar++)
    for (int iech = 0; iech < (int) _dbin->getSampleNumber(); iech++)
    {
      bool found = (ivar == ivar0 && iech == iech0);
      if (!_dbin->isActive(iech) || !_dbin->isIsotopic(iech))
      {
        if (found) return -1;
        continue;
      }
      if (found) return rank;
      rank++;
    }
  return rank;
}

bool KrigingSystem::_isMatCLempty() const
{
  if (_matCL.empty()) return true;
  int size = (int) _matCL.size();
  if (size == 1 && _matCL.at(0).empty()) return true;
  return false;
}

/**
 * This (internal) method is used to modify the Model locally
 * Note: It also modifies the shortcut variables consequently
 * @param model Pointer to the new model
 */
void KrigingSystem::_setLocalModel(Model* model)
{
  _model = model;
  _setInternalShortCutVariablesModel();
}

void KrigingSystem::_setInternalShortCutVariablesModel()
{
  _nvar = _getNVar();
  _nbfl = _getNbfl();
  _nfeq = _getNFeq();
  _nfex = _getNFex();
  _neq  = getNeq(); // reset as it depends on nech and Model
}
/**
 * Assign the values to local variables used as shortcuts
 * @return 1 if the number of active sample is zero
 */
int KrigingSystem::_setInternalShortCutVariablesNeigh()
{
  _nech = getNech();
  _neq  = getNeq();
  return (_nech <= 0);
}
void KrigingSystem::_setInternalShortCutVariablesGeneral()
{
  _ndim = getNDim();
  _nvarCL = _getNVarCL();
  _setInternalShortCutVariablesModel();
}
