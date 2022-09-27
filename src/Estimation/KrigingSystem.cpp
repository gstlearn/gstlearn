/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_old_f.h"

#include "Enum/EKrigOpt.hpp"
#include "Enum/ECalcMember.hpp"
#include "Enum/ELoc.hpp"
#include "Enum/ENeigh.hpp"

#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/PtrGeos.hpp"
#include "Model/Model.hpp"
#include "Model/ANoStat.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighImage.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighWork.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/Law.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Covariances/ACovAnisoList.hpp"
#include "Polynomials/Hermite.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Calculators/CalcMigrate.hpp"

#include <math.h>

#define FF(ib,il)         (ff [(il) * shift + (ib)])
#define FF0(ib,iv)        (ff0 [(ib) + nfeq * (iv)])
#define SIGMA(ib,jb)      (sigma[(jb)    * shift + (ib)])
#define SMEAN(i,isimu)    (_postSimu[isimu][i])

KrigingSystem::KrigingSystem(Db* dbin,
                             Db* dbout,
                             const Model* model,
                             ANeighParam* neighParam)
    : _dbin(dbin),
      _dbout(dbout),
      _modelInit(nullptr),
      _neighParam(neighParam),
      _anam(nullptr),
      _isReady(false),
      _model(nullptr),
      _iptrEst(-1),
      _iptrStd(-1),
      _iptrVarZ(-1),
      _flagEst(false),
      _flagStd(false),
      _flagVarZ(false),
      _flagGlobal(false),
      _flagDataChanged(false),
      _calcul(EKrigOpt::PONCTUAL),
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
      _priorMean(),
      _priorCov(),
      _postMean(),
      _postCov(),
      _postSimu(),
      _varCorrec(),
      _modelSimple(nullptr),
      _flagDGM(false),
      _rDGM(1.),
      _flagFactorKriging(false),
      _nclasses(0),
      _matCL(),
      _flagLTerm(false),
      _lterm(0.),
      _flagAnam(false),
      _dbaux(),
      _flagKeypairWeights(false),
      _flagNeighOnly(false),
      _iptrNeigh(-1),
      _iechOut(-1),
      _nred(0),
      _flagCheckAddress(false),
      _nbghWork(_dbin, _neighParam),
      _nbgh(),
      _flag(),
      _covtab(),
      _drftab(),
      _lhs(),
      _lhsinv(),
      _rhs(),
      _wgt(),
      _zam(),
      _var0(),
      _dbinUidToBeDeleted(),
      _dboutUidToBeDeleted()
{
  // _modelInit is a copy of the input model (const) to allow modifying it
  if (model != nullptr) _modelInit = model->clone();

  // Set the current Model to _modelInit
  _model = _modelInit;

  _resetMemoryGeneral();
}

KrigingSystem::~KrigingSystem()
{
  // Turn OFF this option for future task

  OptDbg::setIndex(0);

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

  // Reset elements in _neighParam

  _neighParam->setFlagXvalid(false);
  _neighParam->setFlagKFold(false);
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
  if (_dbin != nullptr)
  {
    if (nvar > 0 && nvar != _dbin->getVariableNumber())
    {
      messerr("Inconsistent number of Variables - Value is returned as 0");
      return 0;
    }
    nvar = _dbin->getVariableNumber();
  }
  return nvar;
}

int KrigingSystem::_getNVarCL() const
{
  if (_isMatCLempty())
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
  int nech = getNech();
  int nvar = _getNVar();
  int nfeq = _getNFeq();
  int neq = nvar * nech + nfeq;
  return neq;
}

int KrigingSystem::_getNDisc() const
{
  return _ndiscNumber;
}

void KrigingSystem::_resetMemoryPerNeigh()
{
  int nvarCL = _getNVarCL();
  int neq    = getNeq();
  _flag.resize(neq);
  _lhs.resize (neq * neq);
  _rhs.resize (neq * nvarCL);
  _wgt.resize (neq * nvarCL);
  _zam.resize (neq);
}

void KrigingSystem::_resetMemoryGeneral()
{
  int nvar   = _getNVar();
  int nvarCL = _getNVarCL();
  int nbfl   = _getNbfl();
  _covtab.resize(nvar * nvar);
  _drftab.resize(nbfl);
  _var0.resize(nvarCL * nvarCL);
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

/**
 * Update the vector of distances between two samples of two Dbs
 * @param loc_rank1 First Sample Rank in input Db (>=0) or Target (<0)
 * @param loc_rank2 Second Sample Rank in input Db (>=0) or Target (<0)
 * @param dd Vector of distances
 */
void KrigingSystem::_getDistance(int loc_rank1, int loc_rank2, VectorDouble& dd) const
{
  const Db* db1;
  const Db* db2;
  int rank1, rank2;

  if (loc_rank1 >= 0)
  {
    db1 = _dbin;
    rank1 = loc_rank1;
  }
  else
  {
    db1 = _dbout;
    rank1 = _iechOut;
  }

  if (loc_rank2 >= 0)
  {
    db2 = _dbin;
    rank2 = loc_rank2;
  }
  else
  {
    db2 = _dbout;
    rank2 = _iechOut;
  }

  db1->getDistanceVec(rank1, rank2, dd, db2);
  return;
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
    return _dbin->getExternalDrift(rank, ibfl);
  }
  else
  {
    return _dbout->getExternalDrift(_iechOut, ibfl);
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

      return _dbin->getVariable(rank, ivar);

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
    return _dbin->getVarianceError(rank, ivar);
  }
  else
  {
    return _dbout->getVarianceError(_iechOut, ivar);
  }
}
double KrigingSystem::_getMean(int ivarCL) const
{
  double value = 0.;
  if (_isMatCLempty())
  {
    value = _model->getMean(ivarCL);
  }
  else
  {
    int nvar = _getNVar();
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      value += _matCL[ivarCL][ivar] * _model->getMean(ivar);
    }
  }
  return value;
}

double KrigingSystem::_getCoefDrift(int ivar, int il, int ib) const
{
  return _model->getCoefDrift(ivar, il, ib);
}

void KrigingSystem::_setFlag(int iech, int ivar, int value)
{
  int nech = getNech();
  _flag[iech + ivar * nech]= value;
}

int KrigingSystem::_getFlag(int iech, int ivar)
{
  int nech = getNech();
  return _flag[iech + ivar * nech];
}

/****************************************************************************/
/*!
 **  Define the array flag to convert from isotropic to heterotopic case
 **  Stores the reduced number of equations in member '_nred'
 **
 *****************************************************************************/
void KrigingSystem::_flagDefine()
{
  int nech = getNech();
  int nvar = _getNVar();
  int nfeq = _getNFeq();
  int nbfl = _getNbfl();
  int next = _getNFex();
  int neq  = getNeq();
  int ndim = getNDim();
  for (int i = 0; i < neq; i++) _flag[i] = 1;

  /* Check on the coordinates */

  for (int iech = 0; iech < nech; iech++)
  {
    bool valid = true;
    for (int idim = 0; idim < ndim; idim++)
      if (FFFF(_getIdim(_nbgh[iech], idim))) valid = false;
    if (! valid)
      for (int ivar = 0; ivar < nvar; ivar++)
        _setFlag(iech,ivar,0);
  }

  /* Check on the data values */

  for (int iech = 0; iech < nech; iech++)
    for (int ivar = 0; ivar < nvar; ivar++)
      if (FFFF(_getIvar(_nbgh[iech], ivar)))
        _setFlag(iech,ivar,0);

  /* Check on the external drifts */

  if (next > 0)
  {
    for (int iech = 0; iech < nech; iech++)
      for (int ibfl = 0; ibfl < next; ibfl++)
        if (FFFF(_getFext(_nbgh[iech], ibfl)))
          for (int ivar = 0; ivar < nvar; ivar++)
            _setFlag(iech, ivar, 0);
  }

  /* Check on the drift */

  for (int ib = 0; ib < nfeq; ib++)
  {
    int valid = 0;
    for (int il = 0; il < nbfl; il++)
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        if (_getCoefDrift(ivar, il, ib) == 0.) continue;
        for (int iech = 0; iech < nech; iech++)
          if (!FFFF(_getIvar(_nbgh[iech], ivar))) valid++;
      }
    _setFlag(nech+ib, nvar-1, (valid > 0));
  }

  /* Calculate the new number of equations */

  int count = 0;
  for (int i = 0; i < neq; i++)
  {
    if (_flag[i] != 0) count++;
  }
  _nred = count;
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
  int nech = getNech();
  int nvar = _getNVar();
  int nfeq = _getNFeq();

  /* Preliminary check */

  if (nech * nvar < nfeq) return false;

  /* Check that enough information is present */

  int n_cov = 0;
  for (int i = 0; i < nvar * nech; i++)
    n_cov += _flag[i];

  int n_drf = 0;
  for (int i = 0; i < nfeq; i++)
    n_drf += _flag[i + nvar * nech];

  if (n_cov <= 0 || n_cov < n_drf) return false;

  return true;
}

void KrigingSystem::_covtabInit()
{
  int nvar = _getNVar();
  for (int i = 0; i < nvar * nvar; i++) _covtab[i] = 0.;
}

/**
 * Module for calculating the covariance internally
 * It is called for LHS (iech1>=0 && iech2>=0), RHS (iech1>=0 && iech2=-1) and VAR (iech1=-1 && iech2=-1)
 * @param member    Type of usage (LHS, RHS or VAR)
 * @param mode      CovCalcMode structure
 * @param iech1     Rank of the first sample (or -1)
 * @param iech2     Rank of the second sample (or -1)
 * @param d1
 */
void KrigingSystem::_covtabCalcul(const ECalcMember &member,
                                  const CovCalcMode& mode,
                                  int iech1,
                                  int iech2,
                                  const VectorDouble& d1)
{
  // Load the non-stationary parameters if needed

  if (_model->isNoStat())
  {
    const ANoStat *nostat = _model->getNoStat();
    int jech1 = 0;
    int jech2 = 0;
    int icas1 = 0;
    int icas2 = 0;

    switch (member.getValue())
    {
      case ECalcMember::E_LHS:
        icas1 = (iech1 >= 0) ? 1 : 2;
        jech1 = (iech1 >= 0) ? iech1 : _iechOut;
        icas2 = (iech2 >= 0) ? 1 : 2;
        jech2 = (iech2 >= 0) ? iech2 : _iechOut;
        break;

      case ECalcMember::E_RHS:
        icas1 = (iech1 >= 0) ? 1 : 2;
        jech1 = (iech1 >= 0) ? iech1 : _iechOut;
        icas2 = 2;
        jech2 = _iechOut;
        break;

      case ECalcMember::E_VAR:
        icas1 = 2;
        jech1 = _iechOut;
        icas2 = 2;
        jech2 = _iechOut;
        break;
    }
    nostat->updateModel(_model, icas1, jech1, icas2, jech2);
  }

  // Evaluate the Model

  MatrixSquareGeneral mat;
  bool flagSameData = (iech1 == iech2 && iech1 > 0);
  if (flagSameData)
    mat = _model->eval0Nvar(mode);
  else
    mat = _model->evalNvarIpas(1., d1, VectorDouble(), mode);

  // Modify the Model (DGM case). Only provided for the monovariate case

  if (_flagDGM) _covtabModifyDGM(member, iech1, iech2, d1, mat);

  // Expand the Model to all terms of the LHS of the Kriging System

  int nvar = _getNVar();
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      double value = mat.getValue(ivar, jvar);
      _addCOVTAB(ivar,jvar,value);
    }
  return;
}

/**
 * This function modifies the Covariance returned array between 'iech1' and 'iech2'
 * due to the presence of 'flagDGM' flag.
 * This is programmed in the monovariate case only
 * @param member Type of usage (LHS, RHS or VAR)
 * @param iech1  Rank of the first sample (or <0 for target)
 * @param iech2  Rank of the second sample (or <0 for target)
 * @param d1     Vector of increment vector between samples
 * @param mat    Covariance matrix table
 */
void KrigingSystem::_covtabModifyDGM(const ECalcMember &member,
                                     int iech1,
                                     int iech2,
                                     const VectorDouble& d1,
                                     MatrixSquareGeneral& mat)
{
  double covn = 0.;
  double dist = ut_vector_norm(d1);
  if (member == ECalcMember::LHS)
  {
    if (iech1 >= 0 && iech1 == iech2)
      covn = 1.;
    else if (dist <= 0.)
      covn = _rDGM * _rDGM; // Samples at zero distance (due to centering) but different
    else
      covn = _rDGM * _rDGM * mat.getValue(0);
  }
  else if (member == ECalcMember::RHS)
  {
    covn = _rDGM * mat.getValue(0);
  }
  else
  {
    covn = 1.;
  }
  mat.setValue(0, covn);
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
  if (_neighParam == nullptr) return (0.);
  if (_neighParam->getType() != ENeigh::MOVING) return (0.);
  const NeighMoving* neighM = dynamic_cast<const NeighMoving*>(_neighParam);
  int ndim = getNDim();
  VectorDouble dd(ndim);

  /* Calculate the distance increment */

  for (int idim = 0; idim < ndim; idim++)
    dd[idim] = _dbin->getCoordinate(rank1, idim) - _dbout->getCoordinate(rank2, idim);

  /* Anisotropic neighborhood */

  if (neighM->getFlagAniso())
  {

    /* Rotated anisotropy ellipsoid */

    if (neighM->getFlagRotation())
      matrix_product_safe(1, ndim, ndim, dd.data(), neighM->getAnisoRotMats().data(),
                          dd.data());
    for (int idim = 0; idim < ndim; idim++)
      dd[idim] /= neighM->getAnisoCoeff(idim);
  }

  /* Calculate the distance */

  double dist;
  matrix_product(1, ndim, 1, dd.data(), dd.data(), &dist);
  dist = sqrt(dist) / neighM->getRadius();
  double var = 0.;
  if (dist > neighM->getDistCont())
  {
    if (ABS(1. - dist) < eps) dist = 1. - eps;
    var = (dist - neighM->getDistCont()) / (1. - dist);
    var = var * var;
  }
  return (var);
}

void KrigingSystem::_drftabCalcul(const ECalcMember &member, int iech)
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
    _drftab[il] = drft[il];
  return;
}

/****************************************************************************/
/*!
 **  Establish the kriging L.H.S.
 **
 *****************************************************************************/
void KrigingSystem::_lhsCalcul()
{
  int nech   = getNech();
  int nvar   = _getNVar();
  int nfeq   = _getNFeq();
  int nbfl   = _getNbfl();
  int neq    = getNeq();
  int ndim   = getNDim();
  for (int i = 0; i < neq * neq; i++) _lhs[i] = 0.;
  VectorDouble d1(ndim);

  CovCalcMode mode(ECalcMember::LHS);

  /* Establish the covariance part */

  for (int iech = 0; iech < nech; iech++)
    for (int jech = 0; jech < nech; jech++)
    {
      _covtabInit();
      _getDistance(_nbgh[jech], _nbgh[iech], d1);
      _covtabCalcul(ECalcMember::LHS, mode, _nbgh[iech], _nbgh[jech], d1);

      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar < nvar; jvar++)
        {
          _setLHS(iech,ivar,jech,jvar,_getCOVTAB(ivar, jvar));

          /* Correction due to measurement errors */

          double verr = 0.;
          if (_flagCode)
          {
            int code1 = (int) _dbin->getCode(_nbgh[iech]);
            int code2 = (int) _dbin->getCode(_nbgh[jech]);
            if (code1 != 0 && code2 != 0 && code1 == code2)
              verr = _dbin->getVarianceError(_nbgh[iech], 0);
          }
          else
          {
            if (iech == jech)
            {
              verr = _dbin->getVarianceError(_nbgh[iech], ivar);

              if (_neighParam->getFlagContinuous())
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

  /* Establish the drift part */

  if (nfeq <= 0 || nbfl <= 0) return;
  for (int iech = 0; iech < nech; iech++)
  {
    _drftabCalcul(ECalcMember::LHS, _nbgh[iech]);
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int ib = 0; ib < nfeq; ib++)
      {
        double value = 0.;
        for (int il = 0; il < nbfl; il++)
          value += _drftab[il] * _getCoefDrift(ivar, il, ib);
        _setLHS(iech,ivar,ib,nvar,value,true);
        _setLHS(ib,nvar,iech,ivar,value,true);
      }
  }
  return;
}

void KrigingSystem::_lhsIsoToHetero()

{
  int neq = getNeq();
  int lec_lhs = 0;
  int ecr_lhs = 0;
  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
    {
      if (_flag[i] != 0 && _flag[j] != 0) _lhs[ecr_lhs++] = _lhs[lec_lhs];
      lec_lhs++;
    }
  return;
}

VectorInt KrigingSystem::_getRelativePosition()
{
  int neq = getNeq();
  VectorInt rel(neq);
  int j = 0;
  for (int i = 0; i < neq; i++)
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
  int nech = getNech();
  int neq = getNeq();
  VectorInt rel = _getRelativePosition();
  int npass = (_nred - 1) / nbypas + 1;

  /* General Header */

  mestitle(0, "LHS of Kriging matrix (compressed)");
  if (nech > 0) message("Number of active samples    = %d\n", nech);
  message("Total number of equations   = %d\n", neq);
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
        tab_printg(NULL, _lhs[(i) + _nred * (j)]);
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

  int rank = matrix_invert(_lhsinv.data(), _nred, -1);
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
  double nscale = 1;
  int nech   = getNech();
  int nvar   = _getNVar();
  int nvarCL = _getNVarCL();
  int nbfl   = _getNbfl();
  int nfeq   = _getNFeq();
  int ndim   = getNDim();
  int ndisc  = _getNDisc();
  VectorDouble d1(ndim);

  CovCalcMode mode(ECalcMember::RHS);

  /* Establish the covariance part */

  for (int iech = 0; iech < nech; iech++)
  {
    _covtabInit();
    switch (_calcul.toEnum())
    {
      case EKrigOpt::E_PONCTUAL:
        nscale = 1;
        _getDistance(-1, _nbgh[iech], d1);
        _covtabCalcul(ECalcMember::RHS, mode, _nbgh[iech], -1, d1);
        break;

      case EKrigOpt::E_BLOCK:
        if (_flagPerCell) _blockDiscretize();
        nscale = ndisc;
        for (int i = 0; i < nscale; i++)
        {
          _getDistance(-1, _nbgh[iech], d1);
          for (int idim = 0; idim < ndim; idim++) d1[idim] += _getDISC1(i, idim);
          _covtabCalcul(ECalcMember::RHS, mode, _nbgh[iech], -1, d1);
        }
        break;

      case EKrigOpt::E_DRIFT:
        nscale = 1;
        _covtabInit();
        break;
    }

    /* Normalization */

    if (nscale != 1)
    {
      double ratio = 1. / (double) nscale;
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar < nvar; jvar++)
          _prodCOVTAB(ivar,jvar,ratio);
    }

    /* Storage */

    if (_isMatCLempty())
    {
      for (int ivar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar < nvar; jvar++)
          _setRHS(iech,ivar,jvar,_getCOVTAB(ivar,jvar));
    }
    else
    {
      for (int jvarCL = 0; jvarCL < nvarCL; jvarCL++)
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          double value = 0.;
          for (int jvar = 0; jvar < nvar; jvar++)
            value += _matCL[jvarCL][jvar] * _getCOVTAB(ivar, jvar);
          _setRHS(iech,ivar,jvarCL,value);
        }
    }
  }

  /* Establish the drift part */

  if (nfeq <= 0) return 0;

  _drftabCalcul(ECalcMember::RHS, -1);
  for (int il = 0; il < nbfl; il++)
    if (FFFF(_drftab[il])) return 1;

  if (_isMatCLempty())
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int ib = 0; ib < nfeq; ib++)
      {
        double value = 0.;
        for (int il = 0; il < nbfl; il++)
          value += _drftab[il] * _getCoefDrift(ivar, il, ib);
        _setRHS(ib,nvar,ivar,value,true);
      }
  }
  else
  {
    for (int ivarCL = 0; ivarCL < nvarCL; ivarCL++)
    {
      int ib = 0;
      for (int jvar = 0; jvar < nvar; jvar++)
        for (int jl = 0; jl < nbfl; jl++, ib++)
        {
          double value = 0.;
          for (int il = 0; il < nbfl; il++)
            value += _drftab[il] * _getCoefDrift(jvar, il, ib);
          value *= _matCL[ivarCL][jvar];
          _setRHS(ib,nvar,ivarCL,value,true);
        }
    }
  }
  return 0;
}

void KrigingSystem::_rhsIsoToHetero()
{
  int nvarCL = _getNVarCL();
  int neq    = getNeq();

  int lec_rhs = 0;
  int ecr_rhs = 0;
  for (int ivarCL = 0; ivarCL < nvarCL; ivarCL++)
    for (int i = 0; i < neq; i++, lec_rhs++)
    {
      if (_flag[i] == 0) continue;
      _rhs[ecr_rhs++] = _rhs[lec_rhs];
    }
  return;
}

void KrigingSystem::_rhsDump()
{
  int nech   = getNech();
  int neq    = getNeq();
  int nvarCL = _getNVarCL();
  int ndim   = getNDim();
  VectorInt rel = _getRelativePosition();

  /* General Header */

  mestitle(0, "RHS of Kriging matrix (compressed)");
  if (nech > 0) message("Number of active samples    = %d\n", nech);
  message("Total number of equations   = %d\n", neq);
  message("Reduced number of equations = %d\n", _nred);
  message("Number of right-hand sides  = %d\n", nvarCL);

  /* Kriging option */

  switch (_calcul.toEnum())
  {
    case EKrigOpt::E_PONCTUAL:
      message("Punctual Estimation\n");
      break;

    case EKrigOpt::E_BLOCK:
      message("Block Estimation : Discretization = ");
      for (int idim = 0; idim < ndim; idim++)
      {
        if (idim != 0) message(" x ");
        message("%d", _ndiscs[idim]);
      }
      message("\n");
      break;

    case EKrigOpt::E_DRIFT:
      message("Drift Estimation\n");
      break;
  }
  message("\n");

  /* Header line */

  tab_prints(NULL, "Rank");
  if (! _flag.empty()) tab_prints(NULL, "Flag");
  for (int ivarCL = 0; ivarCL < nvarCL; ivarCL++)
    tab_printi(NULL, ivarCL + 1);
  message("\n");

  /* Matrix lines */

  for (int i = 0; i < _nred; i++)
  {
    tab_printi(NULL, i + 1);
    if (! _flag.empty()) tab_printi(NULL, rel[i]);
    for (int ivarCL = 0; ivarCL < nvarCL; ivarCL++)
      tab_printg(NULL, _getRHSC(i,ivarCL));
    message("\n");
  }
  return;
}

void KrigingSystem::_wgtCalcul()
{
  int nvarCL = _getNVarCL();
  matrix_product(_nred, _nred, nvarCL, _lhsinv.data(), _rhs.data(), _wgt.data());
}

void KrigingSystem::_wgtDump(int status)
{
  int nech   = getNech();
  int ndim   = getNDim();
  int nvarCL = _getNVarCL();
  int nfeq   = _getNFeq();
  int ndisc  = _getNDisc();
  VectorDouble sum(nvarCL);
  char string[20];

  /* Header */

  mestitle(0, "(Co-) Kriging weights");
  const DbGrid* dbgrid = dynamic_cast<const DbGrid*>(_dbout);

  /* First line */

  tab_prints(NULL, "Rank");
  for (int idim = 0; idim < ndim; idim++)
  {
    String strloc = getLocatorName(ELoc::X, idim);
    tab_prints(NULL, strloc.c_str());
  }
  if (_dbin->hasCode()) tab_prints(NULL, "Code");
  if (_dbin->getVarianceErrorNumber() > 0)
    tab_prints(NULL, "Err.");
  if (ndisc > 0)
    for (int idim = 0; idim < ndim; idim++)
    {
      (void) gslSPrintf(string, "Size%d", idim + 1);
      tab_prints(NULL, string);
    }
  tab_prints(NULL, "Data");
  for (int ivarCL = 0; ivarCL < nvarCL; ivarCL++)
  {
    (void) gslSPrintf(string, "Z%d*", ivarCL + 1);
    tab_prints(NULL, string);
  }
  message("\n");

  /* Display the information and the weights */

  int lec = 0;
  int cumflag = 0;
  for (int jvarCL = 0; jvarCL < nvarCL; jvarCL++)
  {
    if (nvarCL > 1) message("Using variable Z%-2d\n", jvarCL + 1);

    /* Loop on the samples */

    for (int ivarCL = 0; ivarCL < nvarCL; ivarCL++)
      sum[ivarCL] = 0.;
    for (int iech = 0; iech < nech; iech++, lec++)
    {
      int flag_value = (! _flag.empty()) ? _flag[lec] : 1;
      tab_printi(NULL, iech + 1);
      for (int idim = 0; idim < ndim; idim++)
        tab_printg(NULL, _getIdim(_nbgh[iech], idim));
      if (_dbin->hasCode())
        tab_printg(NULL, _dbin->getCode(_nbgh[iech]));
      if (_dbin->getVarianceErrorNumber() > 0)
        tab_printg(NULL, _getVerr(_nbgh[iech], (_flagCode) ? 0 : jvarCL));
      if (ndisc > 0)
      {
        for (int idim = 0; idim < ndim; idim++)
          if (! _flagPerCell)
            tab_printg(NULL, dbgrid->getDX(idim));
          else
            tab_printg(NULL, dbgrid->getBlockExtension(_nbgh[iech], idim));
      }
      if (_rankPGS < 0)
        tab_printg(NULL, _getIvar(_nbgh[iech], jvarCL));
      else
        tab_prints(NULL,  "    ");

      for (int ivarCL = 0; ivarCL < nvarCL; ivarCL++)
      {
        double value = (! _wgt.empty() && status == 0 && flag_value) ? _getWGTC(cumflag,ivarCL) : TEST;
        if (!FFFF(value)) sum[ivarCL] += value;
        tab_printg(NULL, value);
      }
      if (flag_value) cumflag++;
      message("\n");
    }

    int number = 1 + ndim + 1;
    if (_dbin->getVarianceErrorNumber() > 0) number++;
    if (ndisc > 0) number += ndim;
    tab_prints(NULL, "Sum of weights", number, EJustify::LEFT);
    for (int ivarCL = 0; ivarCL < nvarCL; ivarCL++)
    {
      double value = (status == 0) ? sum[ivarCL] : TEST;
      tab_printg(NULL, value);
    }
    message("\n");
  }
  if (nfeq <= 0) return;

  /* Header */

  mestitle(0, "Drift coefficients");

  /* First line */

  tab_prints(NULL, "Rank");
  tab_prints(NULL, "Lagrange");
  tab_prints(NULL, "Coeff");
  message("\n");

  /* Loop on the drift coefficients */

  cumflag = _nred - nfeq;
  for (int ib = 0; ib < nfeq; ib++)
  {
    int iwgt = ib + cumflag;
    tab_printi(NULL, ib + 1);
    tab_printg(NULL, (status == 0) ? _wgt[iwgt] : TEST);
    if (_flagSimu)
      tab_printg(NULL, 0.);
    else
      tab_printg(NULL, (status == 0) ? _zam[iwgt] : TEST);
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
  int nech = getNech();
  int nvar = _getNVar();
  int nfeq = _getNFeq();

  /* Simulation */

  int ecr = 0;
  for (int isimu = ecr = 0; isimu < _nbsimu; isimu++)
    for (int ivar = 0; ivar < nvar; ivar++, ecr++)
    {
      double simu = 0.;
      if (nfeq <= 0) simu = _getMean(ivar);

      if (status == 0)
      {
        if (_flagBayes)
          simu = _model->_evalDriftCoef(_dbout, _iechOut, ivar, _postSimu[isimu].data());

        int lec = ivar * _nred;
        for (int jvar = 0; jvar < nvar; jvar++)
          for (int iech = 0; iech < nech; iech++)
          {
            int jech = _nbgh[iech];

            double mean = 0.;
            if (nfeq <= 0) mean = _getMean(jvar);
            if (_flagBayes)
              mean = _model->_evalDriftCoef(_dbin, jech, jvar,_postSimu[isimu].data());
            double data = _dbin->getSimvar(ELoc::SIMU, jech, isimu, ivar, _rankPGS, _nbsimu, nvar);
            simu -= _wgt[lec++] * (data + mean);
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
        if (nfeq > 0) simu = TEST;
      }

      /* Add the conditioning kriging to the NC simulation at target */
      _dbout->updSimvar(ELoc::SIMU, _iechOut, isimu, ivar, _rankPGS, _nbsimu, nvar,
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
  int nech   = getNech();
  int nfeq   = _getNFeq();
  int nvarCL = _getNVarCL();

  /* Estimation */

  if (_flagEst)
  {
    for (int ivarCL = 0; ivarCL < nvarCL; ivarCL++)
    {
      double estim = 0.;
      if (nfeq <= 0) estim = _getMean(ivarCL);

      if (status == 0 && (_nred > 0 || nfeq <= 0 || _flagBayes))
      {
        if (_flagBayes)
          estim = _model->_evalDriftCoef(_dbout, _iechOut, ivarCL,
                                         _postMean.data());
        for (int i = 0; i < _nred; i++)
          estim += _getRHSC(i,ivarCL) * _zam[i];

        if (_neighParam->getFlagXvalid() && _xvalidEstim)
        {
          double valdat = _dbin->getVariable(_iechOut, ivarCL);
          estim = (FFFF(valdat)) ? TEST : estim - valdat;
        }
      }
      else
      {
        // In case of failure with KS, set the result to mean
        if (nfeq > 0) estim = TEST;
      }
      _dbout->setArray(_iechOut, _iptrEst + ivarCL, estim);
    }
  }

  /* Variance of the estimation error */

  if (_flagStd)
  {
    for (int ivarCL = 0; ivarCL < nvarCL; ivarCL++)
    {
      double stdv = TEST;
      if (status == 0 && (_nred > 0 || nfeq <= 0 || _flagBayes))
      {
        stdv = _variance(ivarCL, ivarCL);
        if (stdv < 0) stdv = 0.;
        stdv = sqrt(stdv);

        if (_neighParam->getFlagXvalid() && _xvalidStdev)
        {
          double estim = _dbout->getArray(_iechOut, _iptrEst + ivarCL);
          stdv = (FFFF(estim) || stdv <= 0.) ? TEST : estim / stdv;
        }
      }
      _dbout->setArray(_iechOut, _iptrStd + ivarCL, stdv);
    }
  }

  /* Variance of the estimator */

  if (_flagVarZ != 0)
  {
    for (int ivarCL = 0; ivarCL < nvarCL; ivarCL++)
    {
      double varZ = TEST;
      if (status == 0 && (_nred > 0 || nfeq <= 0))
        varZ = _estimateVarZ(ivarCL, ivarCL);
      _dbout->setArray(_iechOut, _iptrVarZ + ivarCL, varZ);
    }
  }

  /* Kriging weights (stored in _dbin) */

  if (_flagWeights != 0)
  {
    for (int ivarCL = 0; ivarCL < nvarCL; ivarCL++)
    {
      int lec = ivarCL * _nred;
      for (int i = 0; i < nech; i++)
      {
        if (status != 0) continue;
        double wgt = _wgt[lec + i];
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

  int ndim   = getNDim();
  int nech   = getNech();
  int nfeq   = _getNFeq();
  int neq    = getNeq();
  int nvar   = _getNVar();
  const DbGrid* dbgrid = dynamic_cast<const DbGrid*>(_dbout);

  VectorInt indn0(ndim);
  VectorInt indg0(ndim);
  VectorInt indnl(ndim);
  VectorInt indgl(ndim);

  _dbaux->rankToIndice(_dbaux->getSampleNumber()/2, indn0);
  dbgrid->rankToIndice(_iechOut, indg0);

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    double estim = 0.;
    if (nfeq <= 0) estim = _getMean(ivar);

    if (status == 0)
    {
      /* Loop on the neighboring points */

      int ecr = ivar * neq;
      for (int jvar = 0; jvar < nvar; jvar++)
      {
        for (int iech = 0; iech < nech; iech++)
        {
          if (FFFF(dbgrid->getVariable(iech, 0))) continue;
          _dbaux->rankToIndice(_nbgh[iech], indnl);
          for (int idim = 0; idim < ndim; idim++)
          {
            indgl[idim] = indg0[idim] - indn0[idim] + indnl[idim];
            indgl[idim] = dbgrid->getMirrorIndex(idim, indgl[idim]);
          }
          int jech = dbgrid->indiceToRank(indgl);
          double data = dbgrid->getVariable(jech, jvar);
          if (FFFF(data) || FFFF(estim))
          {
            estim = TEST;
          }
          else
          {
            if (nfeq <= 0) data -= _getMean(jvar);
            estim += data * _wgt[ecr];
            ecr++;
          }
        }
      }
      _dbout->setArray(_iechOut, _iptrEst + ivar, estim);
    }
    else
    {
      // In case of failure with KS, set the result to mean
      if (nfeq > 0) estim = TEST;
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

  double trueval = _dbin->getVariable(iech, 0);
  if (! FFFF(trueval))
  {
    double variance = 1. / _getLHSINV(iiech, 0, iiech, 0);
    double stdv = sqrt(variance);

    /* Perform the estimation */

    double value = (_xvalidEstim) ? -trueval : 0.;
    for (int jech = 0; jech < _dbin->getSampleNumber(); jech++)
    {
      int jjech = _getFlagAddress(jech, 0);
      if (jjech < 0) continue;
      if (iiech != jjech)
        value -= _getLHSINV(iiech,0,jjech,0) * variance * _dbin->getVariable(jech, 0);
      jjech++;
    }

    if (_flagEst) _dbin->setArray(iech, _iptrEst, value);
    if (_flagStd)
    {
      if (_xvalidStdev) stdv = value / stdv;
      _dbin->setArray(iech, _iptrStd, stdv);
    }
    if (_flagVarZ)
    {
      // TODO Do something
      _dbin->setArray(iech, _iptrVarZ, TEST);
    }
  }
}

/****************************************************************************/
/*!
 **  Establish the variance of the estimator
 **
 ** \param[in]  ivarCL    Rank of the target variable
 ** \param[in]  jvarCL  Rank of the auxiliary variable
 **
 *****************************************************************************/
double KrigingSystem::_estimateVarZ(int ivarCL, int jvarCL)
{
  int nfeq = _getNFeq();
  int cumflag = _nred - nfeq;

  double var = 0.;
  for (int i = 0; i < _nred; i++)
  {
    double signe = (i < cumflag) ? 1. : -1.;
    var += signe * _getRHSC(i,jvarCL) * _getWGTC(i,ivarCL);
  }
  return var;
}

/****************************************************************************/
/*!
 **  Establish the constant term for the variance calculation
 **
 *****************************************************************************/
void KrigingSystem::_variance0()
{
  int nscale = 1;
  int nvar   = _getNVar();
  int nvarCL = _getNVarCL();
  int ndim   = getNDim();
  int ndisc  = _getNDisc();
  VectorDouble d1(ndim, 0.);

  CovCalcMode mode(ECalcMember::VAR);

  _covtabInit();
  switch (_calcul.toEnum())
  {
    case EKrigOpt::E_PONCTUAL:
      nscale = 1;
      _covtabCalcul(ECalcMember::VAR, mode, -1, -1, d1);
      break;

    case EKrigOpt::E_BLOCK:
      nscale = ndisc;
      for (int i = 0; i < nscale; i++)
        for (int j = 0; j < nscale; j++)
        {
          for (int idim = 0; idim < ndim; idim++)
            d1[idim] = _getDISC1(i, idim) - _getDISC2(j,idim);
          _covtabCalcul(ECalcMember::VAR, mode, -1, -1, d1);
       }
      nscale = nscale * nscale;
      break;

    case EKrigOpt::E_DRIFT:
      nscale = 1;
      break;
  }

  /* Normalization */

  double ratio = 1. / (double) nscale;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      _prodCOVTAB(ivar,jvar,ratio);

  /* Storage */

  if (_isMatCLempty())
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        _setVAR0(ivar,jvar,_getCOVTAB(ivar,jvar));
  }
  else
  {
    for (int ivarCL = 0; ivarCL < nvarCL; ivarCL++)
      for (int jvarCL = 0; jvarCL < nvarCL; jvarCL++)
      {
        double value = 0.;
        for (int ivar = 0; ivar < nvar; ivar++)
          for (int jvar = 0; jvar < nvar; jvar++)
            value += _matCL[ivarCL][ivar] * _getCOVTAB(ivar, jvar) * _matCL[jvarCL][jvar];
        _setVAR0(ivarCL,jvarCL,value);
      }
  }
  return;
}

/****************************************************************************/
/*!
 **  Establish the calculation of variance or standard deviation
 **
 ** \param[in]  ivarCL  Rank of the target variable
 ** \param[in]  jvarCL  Rank of the auxiliary variable
 **
 *****************************************************************************/
double KrigingSystem::_variance(int ivarCL, int jvarCL)
{
  int nvarCL = _getNVarCL();

  // The Variance at origin must be updated
  // - in Non-stationary case
  // - for Block Estimation, when the block size if defined per cell
  if (_model->isNoStat() || _flagPerCell) _variance0();

  double var = _getVAR0(ivarCL, jvarCL);
  if (_flagBayes) var += _varCorrec[jvarCL * nvarCL + ivarCL];
  for (int i = 0; i < _nred; i++)
    var -= _getRHSC(i,jvarCL) * _getWGTC(i,ivarCL);

  return var;
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
 *****************************************************************************/
void KrigingSystem::_dual()
{
  int nech = getNech();
  int nvar = _getNVar();
  int nfeq = _getNFeq();

  /* Set the whole array to 0 */

  VectorDouble zext(_nred);
  for (int i = 0; i < _nred; i++) zext[i] = 0.;

  /* Extract the data */

  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    for (int iech = 0; iech < nech; iech++)
    {
      if (! _getFLAG(iech, ivar)) continue;
      double mean = 0.;
      if (nfeq <= 0)
        mean = _getMean(ivar);
      if (_flagBayes)
        mean = _model->_evalDriftCoef(_dbout, _iechOut, ivar,
                                      _postMean.data());
      zext[ecr++] = _getIvar(_nbgh[iech], ivar) - mean;
    }
  }

  /* Operate the product : Z * A-1 */

  matrix_product(_nred, _nred, 1, _lhsinv.data(), zext.data(), _zam.data());

  /* Operate the product : Z * A-1 * Z */

  if (_flagLTerm)
    matrix_product(1, _nred, 1, _zam.data(), zext.data(), &_lterm);

  // Turn back the flag to OFF in order to avoid provoking
  // the _dual() calculations again

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

  // Perform some pre-calculation when variance of estimator is requested
  if (_flagStd)
    _variance0();

  // Perform some preliminary work in the case of Image
  const NeighImage* neighI = dynamic_cast<const NeighImage*>(_neighParam);
  if (neighI != nullptr)
  {
    if (_prepareForImage(neighI)) return false;
  }

  // In Bayesian case, calculate the Posterior information
  if (_flagBayes)
  {
    _bayesPreCalculations();
  }

  _isReady = true;
  return _isReady;
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
  // been selected (in isReady(). No need to compute them again.
  bool skipCalculAll = false;
  if (_neighParam->getType().toEnum() == ENeigh::E_IMAGE) skipCalculAll = true;

  // In case of Cross-validation in Unique Neighborhood, do not establish the RHS
  bool caseXvalidUnique = false;
  if (_neighParam->getType().toEnum() == ENeigh::E_UNIQUE &&
      _neighParam->getFlagXvalid()) caseXvalidUnique = true;

  // Store the Rank of the Target sample
  _iechOut = iech_out;

  int status = 0;
  if (skipCalculAll) goto label_store;

  if (! _dbout->isActive(_iechOut)) return 0;
  OptDbg::setIndex(_iechOut + 1);
  if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
  {
    if (_flagFactorKriging)
      message("\nProcessing Factor %d / %d\n",_model->getAnamIClass(), _nclasses);

    mestitle(1, "Target location");
    db_sample_print(_dbout, _iechOut, 1, 0, 0);
  }

  // Elaborate the Neighborhood

  if (caseXvalidUnique) _neighParam->setFlagXvalid(false);
  _nbgh = _nbghWork.select(_dbout,_iechOut,_rankColCok);
  status = (getNech() <= 0);
  if (_flagNeighOnly) goto label_store;
  if (status) goto label_store;

  /* Establish the Kriging L.H.S. */

  if (!_nbghWork.isUnchanged() || _neighParam->getFlagContinuous() ||
      OptDbg::force())
  {
    if (_flagBayes) _model = _modelSimple;
    status = _prepar();
    if (_flagBayes) _model = _modelInit;
  }

  // Establish the pre-calculation involving the data information

  if (!_nbghWork.isUnchanged() || _neighParam->getFlagContinuous() ||
      _flagDataChanged || OptDbg::force())
  {
    _dual();
  }

  if (caseXvalidUnique) _neighParam->setFlagXvalid(true);

  /* Establish the Kriging R.H.S. */

  if (caseXvalidUnique) goto label_store;

  if (_flagBayes) _model = _modelSimple;
  _rhsCalcul();
  if (_flagBayes) _model = _modelInit;

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

  // Correct the Variance in Bayesian case

  if (_flagBayes) _bayesCorrectVariance();

  // Store the results in the output Db

  if (_flagNeighOnly)
  {
    VectorDouble tab = _nbghWork.summary(_dbout, _iechOut, _rankColCok);
    _neighCalcul(status, tab);

  }
  else if (_neighParam->getType().toEnum() == ENeigh::E_IMAGE)
  {
    // Image Neighborhood case

    _estimateCalculImage(status);
  }
  else if (_neighParam->getType().toEnum() == ENeigh::E_UNIQUE)
  {
    // Unique Neighborhood case

    if (_neighParam->getFlagXvalid())
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
  if (_neighParam->getFlagXvalid())
    mestitle(0, "Cross-validation results");
  else
    mestitle(0, "(Co-) Kriging results");
  message("Target Sample = %d\n", _iechOut + 1);

  /* Loop on the results */

  int nvar = _getNVar();
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    if (_neighParam->getFlagXvalid())
    {

      // Printout of Cross-Validation test

      message("Variable Z%-2d\n", ivar + 1);
      double estval = TEST;
      double esterr = TEST;
      double sterr  = TEST;
      double sigma  = TEST;

      if (_iptrEst >= 0)
      {
        double trueval = (status == 0) ? _dbin->getVariable(_iechOut, ivar) : TEST;
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
        value = (status == 0) ? _var0[ivar + nvar * ivar] : TEST;
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
  int nvar = _getNVar();

  mestitle(0, "Simulation results");

  /* Loop on the results */

  int ecr = 0;
  for (int isimu = 0; isimu < _nbsimu; isimu++)
    for (int ivar = 0; ivar < nvar; ivar++, ecr++)
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

  _flagEst = _iptrEst >= 0 || (_neighParam->getFlagXvalid() && _iptrStd >= 0);
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
  int ndim = getNDim();
  int ntot = _getNDisc();
  int memo = law_get_random_seed();
  law_set_random_seed(1234546);
  const DbGrid* dbgrid = dynamic_cast<const DbGrid*>(_dbout);

  /* Loop on the discretization points */

  for (int i = 0; i < ntot; i++)
  {
    int jech = i;
    int nval = ntot;
    for (int idim = ndim - 1; idim >= 0; idim--)
    {
      double taille =
          (! _flagPerCell) ? dbgrid->getDX(idim) : _dbout->getBlockExtension(_iechOut, idim);
      int nd = _ndiscs[idim];
      nval /= nd;
      int j = jech / nval;
      jech -= j * nval;
      double local = taille * ((j + 0.5) / nd - 0.5);
      _setDISC1(i,idim, local);
      _setDISC2(i,idim, local + taille * law_uniform(-0.5, 0.5) / (double) nd);
    }
  }
  law_set_random_seed(memo);
}

int KrigingSystem::setKrigOptCalcul(const EKrigOpt& calcul,
                                    const VectorInt& ndiscs,
                                    bool flag_per_cell)
{
  _isReady = false;
  int ndim = getNDim();
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

    if (flag_per_cell)
    {
      _flagPerCell = true;
    }
    NeighMoving* neighM = dynamic_cast<NeighMoving*>(_neighParam);
    if (neighM != nullptr)
    {
      _flagPerCell = neighM->getForceWithinBlock();
    }

    // Discretization is stored

    _ndiscs = ndiscs;
    _ndiscNumber = ut_vector_prod(_ndiscs);
    _disc1.resize(_ndiscNumber * ndim);
    _disc2.resize(_ndiscNumber * ndim);

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
    _neighParam->setFlagXvalid(false);
    _neighParam->setFlagKFold(false);
  }
  else
  {
    _neighParam->setFlagXvalid(flag_xvalid);
    if (flag_kfold)
    {
      if (_neighParam->getType() == ENeigh::UNIQUE)
      {
        messerr("K-FOLD is not available in Unique Neighborhood");
        return 1;
      }
      if (! _dbin->hasCode())
        messerr("The K-FOLD option is ignored as no Code is defined");
    }
    _neighParam->setFlagKFold(flag_kfold);
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
                                   const VectorDouble& prior_cov)
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
    const NeighUnique* neighU = dynamic_cast<const NeighUnique*>(_neighParam);
    if (neighU == nullptr)
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
  if (_isMatCLempty()) return 0;
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
    if (! _dbin->hasCode() || _dbin->getVarianceErrorNumber() != 1)
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
 _nbghWork.setFlagSimu(flagSimu);
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

int KrigingSystem::setKrigOptDGM(bool flag_dgm, double rcoeff, double eps)
{
  _isReady = false;
  if (! flag_dgm)
  {
    _flagDGM = flag_dgm;
    _rDGM = 1.;
    return 0;
  }

  if (_model->getMinOrder() != -1)
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
  if (rcoeff <= 0. || rcoeff >= 1.)
  {
    messerr("For DGM option, the change of support coefficient must be in ]0.,1[ (%lf)",rcoeff);
    return 1;
  }
  _flagDGM = flag_dgm;
  _rDGM = rcoeff;
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
    if (_neighParam->getType() != ENeigh::UNIQUE)
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
  _model->setAnamIClass(index_class);
  _nclasses = nclasses;

  // Update C00 if the variance calculation is required
  if (_flagStd) _variance0();

  return 0;
}

/**
 * This function switches ON the systematic test of addresses before running
 * @param flagCheckAddress True if addresses must be systematically checked
 * @remark When turned ON, this option slows the process.
 * @remark It should only be used for Debugging purposes.
 */
int KrigingSystem::updKrigOptCheckAddress(bool flagCheckAddress)
{
  _flagCheckAddress = flagCheckAddress;
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
  if (_neighParam != nullptr)
  {
    if (ndim > 0 && ndim != _neighParam->getNDim())
    {
      messerr("Incompatible Space Dimension of '_neighParam'");
      return false;
    }
    ndim = _neighParam->getNDim();
  }

  /****************************/
  /* Checking Variable Number */
  /****************************/

  int nvar = 0;
  if (_dbin != nullptr && ! _flagSimu)
  {
    if (nvar > 0 && nvar != _dbin->getVariableNumber())
    {
      messerr("Incompatible Variable Number of '_dbin'");
      return false;
    }
    nvar = _dbin->getVariableNumber();
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
      if (nfex != _dbout->getExternalDriftNumber())
      {
        messerr("Incompatible Number of External Drifts:");
        messerr("- In 'Model' = %d", nfex);
        messerr("- In '_dbout' = %d", _dbout->getExternalDriftNumber());
        return false;
      }
    }
    if (_dbin != nullptr)
    {
      if (_dbin->getExternalDriftNumber() == 0)
      {
        if (migrateByLocator(_dbout, _dbin, ELoc::F)) return false;
        // Store the UID of the newly created variables to be deleted at the end of the process
        _dbinUidToBeDeleted = _dbin->getUIDsByLocator(ELoc::F);
      }
      if (nfex != _dbin->getExternalDriftNumber())
      {
        messerr("Incompatible Number of External Drifts:");
        messerr("- In 'Model' = %d", nfex);
        messerr("- In 'dbin' = %d", _dbin->getExternalDriftNumber());
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

    _model->setField(ut_vector_extension_diagonal(db_mini, db_maxi));
  }

  /*****************************/
  /* Checking the Neighborhood */
  /*****************************/

  if (_neighParam != nullptr)
  {
    if (_neighParam->getType() == ENeigh::IMAGE && (_dbout == nullptr || ! _dbout->isGrid()))
    {
      messerr("The Image neighborhood can only be used when the output Db is a grid");
      return false;
    }

    if (_neighParam->getType() == ENeigh::UNIQUE &&
        _neighParam->getFlagXvalid() && nvar > 1)
    {
     messerr("The algorithm for Cross-Validation in Unique Neighborhood");
     messerr("is restricted to a single variable");
     return false;
    }
  }

  /******************************/
  /* Preparing Non-stationarity */
  /******************************/

  if (_model != nullptr && _model->isNoStat())
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
  }

  /**************************/
  /* Checking cross-options */
  /**************************/

  if (_flagDGM && _calcul != EKrigOpt::PONCTUAL)
  {
    messerr("The DGM option is incompatible with 'Block' calculation option");
    return false;
  }
  return true;
}

// Shortcuts for addressing functions

int KrigingSystem::_IND(int iech, int ivar, int nech) const
{
  return   ((iech) + (ivar) * nech);
}
int KrigingSystem::_getFLAG(int iech, int ivar) const
{
  int nech = getNech();
  int ind  = _IND(iech, ivar, nech);
  if (_flagCheckAddress)
  {
    _checkAddress("_getFLAG","iech",iech,getNech());
    _checkAddress("_getFLAG","ivar",ivar,_getNVar());
    _checkAddress("_getFLAG","address",ind,(int) _flag.size());
  }
  return _flag[ind];
}
double KrigingSystem::_getCOVTAB(int ivar,int jvar) const
{
  int nvar = _getNVar();
  int iad = (jvar) * nvar + (ivar);
  if (_flagCheckAddress)
  {
    _checkAddress("_getCOVTAB","ivar",ivar,_getNVar());
    _checkAddress("_getCOVTAB","jvar",jvar,_getNVar());
    _checkAddress("_getCOVTAB","address",iad,(int) _covtab.size());
  }
  return _covtab[iad];
}
void KrigingSystem::_setCOVTAB(int ivar,int jvar,double value)
{
  int nvar = _getNVar();
  int iad = (jvar) * nvar + (ivar);
  if (_flagCheckAddress)
  {
    _checkAddress("_setCOVTAB","ivar",ivar,_getNVar());
    _checkAddress("_setCOVTAB","jvar",jvar,_getNVar());
    _checkAddress("_setCOVTAB","address",iad,(int) _covtab.size());
  }
  _covtab[iad] = value;
}
void KrigingSystem::_addCOVTAB(int ivar,int jvar,double value)
{
  int nvar = _getNVar();
  int iad = (jvar) * nvar + (ivar);
  if (_flagCheckAddress)
   {
     _checkAddress("_addCOVTAB","ivar",ivar,_getNVar());
     _checkAddress("_addCOVTAB","jvar",jvar,_getNVar());
     _checkAddress("_addCOVTAB","address",iad,(int) _covtab.size());
   }
  _covtab[iad] += value;
}
void KrigingSystem::_prodCOVTAB(int ivar,int jvar,double value)
{
  int nvar = _getNVar();
  int iad = (jvar) * nvar + (ivar);
  if (_flagCheckAddress)
   {
     _checkAddress("_prodCOVTAB","ivar",ivar,_getNVar());
     _checkAddress("_prodCOVTAB","jvar",jvar,_getNVar());
     _checkAddress("_prodCOVTAB","address",iad,(int) _covtab.size());
   }
  _covtab[iad] *= value;
}
double KrigingSystem::_getRHS(int iech, int ivar, int jvCL) const
{
  int nech = getNech();
  int neq  = getNeq();
  int ind  = _IND(iech, ivar, nech);
  int iad  = ind + neq * (jvCL);
  if (_flagCheckAddress)
  {
    _checkAddress("_getRHS","iech",iech, getNech());
    _checkAddress("_getRHS","ivar",ivar,_getNVar());
    _checkAddress("_getRHS","jvCL",jvCL,_getNVarCL());
    _checkAddress("_getRHS","address",iad,(int) _rhs.size());
  }
  return _rhs[iad];
}
void KrigingSystem::_setRHS(int iech, int ivar, int jvCL, double value, bool isForDrift)
{
  int nech = getNech();
  int neq = getNeq();
  int ind = _IND(iech,ivar,nech);
  int iad = ind + neq * (jvCL);
  if (_flagCheckAddress)
  {
    if (isForDrift)
    {
      if (ivar == _getNVar())
      {
        _checkAddress("_setRHS","ib",iech,_getNFeq());
      }
      else
      {
        _checkAddress("_setRHS", "iech", iech, getNech());
        _checkAddress("_setRHS", "ivar", ivar, _getNVar());
      }
    }
    else
    {
      _checkAddress("_setRHS", "iech", iech, getNech());
      _checkAddress("_setRHS", "ivar", ivar, _getNVar());
    }
    _checkAddress("_setRHS", "jvCL", jvCL, _getNVarCL());
    _checkAddress("_setRHS","address",iad, (int) _rhs.size());
  }
  _rhs[iad] = value;
}
double KrigingSystem::_getRHSC(int i, int jvCL) const
{
  int iad = (i) + _nred * (jvCL);
  if (_flagCheckAddress)
  {
    _checkAddress("_getRHSC","i",i,_nred);
    _checkAddress("_getRHSC","jvCL",jvCL,_getNVarCL());
    _checkAddress("_getRHSC","address",iad, (int) _rhs.size());
  }
  return _rhs[iad];
}
double KrigingSystem::_getWGTC(int i,int jvCL) const
{
  int iad = (i) + _nred * (jvCL);
  if (_flagCheckAddress)
  {
    _checkAddress("_getWGTC","i",i,_nred);
    _checkAddress("_getWGTC","jvCL",jvCL,_getNVarCL());
    _checkAddress("_getWGTC","address",iad, (int) _wgt.size());
  }
  return _wgt[iad];
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
  int nech = getNech();
  int neq  = getNeq();
  int indi = _IND(iech, ivar, nech);
  int indj = _IND(jech, jvar, nech);
  int iad  = indi + neq * indj;

  if (_flagCheckAddress)
  {
    if (isForDrift)
    {
      if (ivar == _getNVar())
      {
        _checkAddress("_setLHS","ib",iech,_getNFeq());
      }
      else
      {
        _checkAddress("_setLHS","iech",iech,getNech());
        _checkAddress("_setLHS","ivar",ivar,_getNVar());
      }
      if (jvar == _getNVar())
      {
        _checkAddress("_setLHS","jb",jech,_getNFeq());
      }
      else
      {
        _checkAddress("_setLHS","jech",jech,getNech());
        _checkAddress("_setLHS","jvar",jvar,_getNVar());
      }
    }
    else
    {
      _checkAddress("_setLHS","iech",iech,getNech());
      _checkAddress("_setLHS","ivar",ivar,_getNVar());
      _checkAddress("_setLHS","jech",jech,getNech());
      _checkAddress("_setLHS","jvar",jvar,_getNVar());
    }
    _checkAddress("_setLHS","address",iad, (int) _lhs.size());
  }
  _lhs [iad] = value;
}
void KrigingSystem::_addLHS(int iech, int ivar, int jech, int jvar, double value)
{
  int nech = getNech();
  int neq  = getNeq();
  int indi = _IND(iech, ivar, nech);
  int indj = _IND(jech, jvar, nech);
  int iad  = indi + neq * indj;
  if (_flagCheckAddress)
  {
    _checkAddress("_addLHS","iech",iech,getNech());
    _checkAddress("_addLHS","ivar",ivar,_getNVar());
    _checkAddress("_addLHS","jech",jech,getNech());
    _checkAddress("_addLHS","jvar",jvar,_getNVar());
    _checkAddress("_addLHS","address",iad, (int) _lhs.size());
  }
  _lhs[iad] += value;
}
void KrigingSystem::_prodLHS(int iech, int ivar, int jech, int jvar, double value)
{
  int nech = getNech();
  int neq  = getNeq();
  int indi = _IND(iech, ivar, nech);
  int indj = _IND(jech, jvar, nech);
  int iad  = indi + neq * indj;

  if (_flagCheckAddress)
  {
    _checkAddress("_prodLHS","iech",iech,getNech());
    _checkAddress("_prodLHS","ivar",ivar,_getNVar());
    _checkAddress("_prodLHS","jech",jech,getNech());
    _checkAddress("_prodLHS","jvar",jvar,_getNVar());
    _checkAddress("_prodLHS","address",iad, (int) _lhs.size());
  }
  _lhs[iad] *= value;
}
double KrigingSystem::_getLHS(int iech, int ivar, int jech, int jvar) const
{
  int nech = getNech();
  int neq  = getNeq();
  int indi = _IND(iech, ivar, nech);
  int indj = _IND(jech, jvar, nech);
  int iad  = indi + neq * indj;

  if (_flagCheckAddress)
  {
    _checkAddress("_getLHS","iech",iech,getNech());
    _checkAddress("_getLHS","ivar",ivar,_getNVar());
    _checkAddress("_getLHS","jech",jech,getNech());
    _checkAddress("_getLHS","jvar",jvar,_getNVar());
    _checkAddress("_getLHS","address",iad, (int) _lhs.size());
  }
  return _lhs[iad];
}

double KrigingSystem::_getLHSC(int i, int j) const
{
  int iad = (i) + _nred * (j);
  if (_flagCheckAddress)
  {
    _checkAddress("_getLHSC","i",i,_nred);
    _checkAddress("_getLHSC","j",j,_nred);
    _checkAddress("_getLHSC","address",iad, (int) _lhs.size());
  }
  return _lhs[iad];
}
double KrigingSystem::_getLHSINV(int iech, int ivar, int jech, int jvar) const
{
  int nech = getNech();
  int neq  = getNeq();
  int indi = _IND(iech, ivar, nech);
  int indj = _IND(jech, jvar, nech);
  int iad  = indi + neq * indj;

  if (_flagCheckAddress)
  {
    _checkAddress("_getLHSINV","iech",iech,getNech());
    _checkAddress("_getLHSINV","ivar",ivar,_getNVar());
    _checkAddress("_getLHSINV","jech",jech,getNech());
    _checkAddress("_getLHSINV","jvar",jvar,_getNVar());
    _checkAddress("_getLHSINV","address",iad, (int) _lhsinv.size());
  }
  return _lhsinv[iad];
}
double KrigingSystem::_getDISC1(int idisc, int idim) const
{
  int iad = (idim) * _ndiscNumber + (idisc);
  if (_flagCheckAddress)
  {
    _checkAddress("_getDISC1","idisc",idisc,_ndiscNumber);
    _checkAddress("_getDISC1","idim",idim,getNDim());
    _checkAddress("_getDISC1","address",iad, (int) _disc1.size());
  }
  return _disc1[iad];
}
void KrigingSystem::_setDISC1(int idisc, int idim, double value)
{
  int iad = (idim) * _ndiscNumber + (idisc);
  if (_flagCheckAddress)
  {
    _checkAddress("_setDISC1","idisc",idisc,_ndiscNumber);
    _checkAddress("_setDISC1","idim",idim,getNDim());
    _checkAddress("_setDISC1","address",iad, (int) _disc1.size());
  }
  _disc1[iad] = value;
}
double KrigingSystem::_getDISC2(int idisc,int idim) const
{
  int iad = (idim) * _ndiscNumber + (idisc);
  if (_flagCheckAddress)
  {
    _checkAddress("_getDISC2","idisc",idisc,_ndiscNumber);
    _checkAddress("_getDISC2","idim",idim,getNDim());
    _checkAddress("_getDISC2","address",iad, (int) _disc2.size());
  }
  return _disc2[iad];
}
void KrigingSystem::_setDISC2(int idisc,int idim, double value)
{
  int iad = (idim) * _ndiscNumber + (idisc);
  if (_flagCheckAddress)
  {
    _checkAddress("_setDISC2","idisc",idisc,_ndiscNumber);
    _checkAddress("_setDISC2","idim",idim,getNDim());
    _checkAddress("_setDISC2","address",iad, (int) _disc2.size());
  }
  _disc2[iad] = value;
}
double KrigingSystem::_getVAR0(int ivCL, int jvCL) const
{
  int nvarCL = _getNVarCL();
  int iad    = (jvCL) + nvarCL * (ivCL);

  if (_flagCheckAddress)
  {
    _checkAddress("_getVAR0","ivCL",ivCL,_getNVarCL());
    _checkAddress("_getVAR0","jvCL",jvCL,_getNVarCL());
    _checkAddress("_getVAR0","address",iad, (int) _var0.size());
  }
  return _var0[iad];
}
void KrigingSystem::_setVAR0(int ivCL, int jvCL, double value)
{
  int nvarCL = _getNVarCL();
  int iad    = (jvCL) + nvarCL * (ivCL);

  if (_flagCheckAddress)
  {
    _checkAddress("_setVAR0","ivCL",ivCL,_getNVarCL());
    _checkAddress("_setVAR0","jvCL",jvCL,_getNVarCL());
    _checkAddress("_setVAR0","address",iad, (int) _var0.size());
  }
  _var0[iad] = value;
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
  int ndim = getNDim();
  int nvar = _getNVar();
  double seuil = 1. / neighI->getSkip();

  /* Core allocation */

  VectorInt nx(ndim);
  int nech = 1;
  for (int i=0; i<ndim; i++)
  {
    nx[i] = 2 * neighI->getImageRadius(i) + 1;
    nech *= nx[i];
  }

  law_set_random_seed(12345); // TODO to be suppressed when in KrigingSystem
  VectorBool sel(nech);
  for (int iech = 0; iech < nech; iech++)
     sel[iech] = (law_uniform(0., 1.) < seuil) ? true : false;

  VectorDouble tab(nech * nvar);
  int iecr = 0;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int iech = 0; iech < nech; iech++)
        tab[iecr++] = (sel[iech]) ? 0. : TEST;

  /* Create the grid */

  _dbaux = DbGrid::create(nx, dbgrid->getDXs(), dbgrid->getX0s(), dbgrid->getAngles());
  _dbaux->addColumns(tab, "Test", ELoc::Z);

  /* Shift the origin */

  VectorDouble coor(ndim);
  _dbaux->rankToCoordinateInPlace(nech/2, coor);
  for (int i=0; i<ndim; i++) _dbaux->setX0(i, _dbaux->getX0(i) - coor[i]);
  if (db_grid_define_coordinates(_dbaux)) return 1;

  // Setup the Kriging pre-processings

  if (_prepareForImageKriging(_dbaux)) return 1;

  return 0;
}

bool KrigingSystem::_prepareForImageKriging(Db* dbaux)
{
  // Save pointers to previous Data Base (must be restored at the end)
  Db* dbin_loc  = _dbin;
  Db* dbout_loc = _dbout;

  _dbin  = dbaux;
  _dbout = dbaux;
  int error = 1;
  int ndim  = getNDim();

  /* Prepare the neighborhood (mimicking the Unique neighborhood) */

  NeighUnique* neighU = NeighUnique::create(ndim, false);
  _nbghWork.initialize(dbaux, neighU);

  _iechOut = dbaux->getSampleNumber() / 2;
  _nbgh = _nbghWork.select(_dbout,_iechOut,VectorInt());
  int nech = getNech();
  bool status = (nech <= 0);

  /* Establish the L.H.S. */

  status = _prepar();
  if (status) goto label_end;
  _dual();

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
  _nbghWork.initialize(_dbin, _neighParam);
  return error;
}

void KrigingSystem::_saveWeights(int status)
{
  if (status != 0) return;

  int nech = getNech();
  int nvar = _getNVar();
  VectorDouble values(5);
  values[0] = _iechOut;

  /* Loop on the output variables */

  int lec = 0;
  int cumflag = 0;
  for (int jvar = 0; jvar < nvar; jvar++)
  {
    values[1] = jvar;

    /* Loop on the input samples */

    for (int iech = 0; iech < nech; iech++, lec++)
    {
      int flag_value = (! _flag.empty()) ? _flag[lec] : 1;
      if (flag_value)
      {
        values[2] = _nbgh[iech];

        /* Loop on the input variables */

        for (int ivar = 0; ivar < nvar; ivar++)
        {
          int iwgt = _nred * ivar + cumflag;
          double wgtloc = (! _wgt.empty() && flag_value) ? _wgt[iwgt] : TEST;
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
VectorDouble KrigingSystem::getSampleCoordinates() const
{
  int ndim = getNDim();
  int nech = getNech();

  VectorDouble xyz(ndim * nech);
  int ecr = 0;
  for (int idim = 0; idim < ndim; idim++)
    for (int iech = 0; iech < nech; iech++)
      xyz[ecr++] = _getIdim(_nbgh[iech], idim);
  return xyz;
}

VectorDouble KrigingSystem::getSampleData() const
{
  VectorDouble zext(_nred);
  int nvar = _getNVar();
  int nech = getNech();
  for (int i = 0; i < _nred; i++) zext[i] = 0.;

  /* Extract the data */

  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    for (int iech = 0; iech < nech; iech++)
    {
      if (!_getFLAG(iech, ivar)) continue;
      zext[ecr++] = _getIvar(_nbgh[iech], ivar);
    }
  }
  return zext;
}

VectorDouble KrigingSystem::getRHSC(int ivar) const
{
  VectorDouble rhs(_nred);

  for (int i = 0; i < _nred; i++)
    rhs[i] = _getRHSC(i, ivar);
  return rhs;
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
  int nfeq = _getNFeq();
  int nvar = _getNVar();

  // Elaborate the (Unique) Neighborhood
  _nbgh = _nbghWork.select(_dbout,_iechOut,_rankColCok);
  if (getNech() <= 0) return 1;

  /* Prepare the Kriging matrix (without correction) */

  _flagBayes = false;
  if (_prepar()) return 1;
  _flagBayes = true;
  int shift = _nred - nfeq;

  /* Complementary core allocation */

  VectorDouble ff(shift * nfeq);
  VectorDouble smu(nfeq);
  VectorDouble sigma(shift * shift);
  VectorDouble vars(shift);

  // Create the array of variables

  int ind = 0;
  for (int iech = 0; iech < _dbin->getSampleNumber(); iech++)
  {
    if (! _dbin->isActive(iech)) continue;
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      double value = _dbin->getVariable(_nbgh[iech], ivar);
      if (FFFF(value)) continue;
      vars[ind++] = value;
    }
  }

  /* Copy DCOV into S and DMEAN into RMEAN */

  _postCov  = _priorCov;
  _postMean = _priorMean;

  /* Establish the drift array FF */

  for (int il = 0; il < nfeq; il++)
    for (int ib = 0; ib < shift; ib++)
      FF(ib,il) = _getLHSC(ib, shift + il);

  /* Calculate S-1 */

  if (matrix_invert(_postCov.data(), nfeq, -1)) return 1;

  /* Calculate: SMU = S-1 * MEAN */

  matrix_product(nfeq, nfeq, 1, _postCov.data(), _priorMean.data(), smu.data());

  /* Covariance matrix SIGMA */

  for (int ib = 0; ib < shift; ib++)
    for (int jb = 0; jb < shift; jb++)
      SIGMA(ib,jb) = _getLHSC(ib, jb);

  /* Calculate SIGMA-1 */

  if (matrix_invert(sigma.data(), shift, -1)) return 1;

  /* Inverse of posterior covariance matrix: SC-1 = FFt * SIGMA-1 * FF + S-1 */

  int ecr = 0;
  for (int il = 0; il < nfeq; il++)
    for (int jl = 0; jl < nfeq; jl++)
    {
      double value = 0.;
      for (int ib=0; ib<shift; ib++)
        for (int jb=0; jb<shift; jb++)
          value += FF(ib,il) * SIGMA(ib,jb) * FF(jb,jl);
      _postCov[ecr++] += value;
    }

  /* Calculating: SMU = FFt * SIGMA-1 * Z + S-1 * MU */

  for (int il = 0; il < nfeq; il++)
  {
    double value = 0.;
    for (int ib=0; ib<shift; ib++)
      for (int jb=0; jb<shift; jb++)
        value += FF(ib,il) * SIGMA(ib,jb) * vars[jb];
    smu[il] += value;
  }

  /* Posterior mean: _postMean = SC * SMU */

  if (matrix_invert(_postCov.data(), nfeq, -1)) return 1;
  matrix_product(nfeq, nfeq, 1, _postCov.data(), smu.data(), _postMean.data());

  if (OptDbg::query(EDbg::BAYES))
  {
    mestitle(0, "Bayesian Drift coefficients");
    print_matrix("Prior Mean", 0, 1, nfeq, 1, NULL, _priorMean.data());
    print_matrix("Prior Variance-Covariance", 0, 1, nfeq, nfeq, NULL,
                 _priorCov.data());
    print_matrix("Posterior Mean", 0, 1, nfeq, 1, NULL, _postMean.data());
    print_matrix("Posterior Variance-Covariance", 0, 1, nfeq, nfeq, NULL,
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
  int nvar = _getNVar();
  int nfeq = _getNFeq();
  int nbfl = _getNbfl();

  /* Establish the Drift matrix */

  if (nbfl <= 0 || nfeq <= 0) return;

  VectorDouble ff0(nvar * nfeq);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int ib = 0; ib < nfeq; ib++)
    {
      double value = 0.;
      for (int il = 0; il < nbfl; il++)
        value += _drftab[il] * _getCoefDrift(ivar, il, ib);
      FF0(ib,ivar) = value;
    }

  /* Correct the arrays */

  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      double value = 0.;
      int lec = 0;
      for (int il = 0; il < nfeq; il++)
        for (int jl = 0; jl < nfeq; jl++)
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
  int nfeq = _getNFeq();
  if (nfeq <= 0) return;
  int nftri = nfeq * (nfeq + 1) / 2;
  int memo = law_get_random_seed();

  // Dimension '_smean' to store simulated posterior mean
  _postSimu.resize(_nbsimu);
  for (int isimu = 0; isimu < _nbsimu; isimu++) _postSimu[isimu].resize(nfeq);

  /* Core allocation */

  VectorDouble trimat(nftri);
  VectorDouble rndmat(nfeq);
  VectorDouble simu(nfeq);
  // The array _postCov is duplicated as the copy is destroyed by matrix_cholesky_decompose
  VectorDouble rcov = _postCov;

  /* Cholesky decomposition */

  int rank = matrix_cholesky_decompose(rcov.data(), trimat.data(), nfeq);

  if (rank > 0)
  {
    messerr("Error in the Cholesky Decomposition of the covariance matrix");
    messerr("Rank of the Matrix = %d", rank);
    messerr("The Drift coefficients have been set to their posterior mean");
    for (int isimu = 0; isimu < _nbsimu; isimu++)
      for (int il = 0; il < nfeq; il++)
        SMEAN(il,isimu) = _postMean[il];
  }
  else
  {
    for (int isimu = 0; isimu < _nbsimu; isimu++)
    {

      /* Draw a vector of gaussian independent values */

      for (int il = 0; il < nfeq; il++) rndmat[il] = law_gaussian();

      /* Product of the Lower triangular matrix by the random vector */

      matrix_cholesky_product(1, nfeq, 1, trimat.data(), rndmat.data(), simu.data());

      /* Add the mean */

      for (int il = 0; il < nfeq; il++)
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
      for (int il = 0; il < nfeq; il++)
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

  double condexp = hermiteCondExpElement(est, std, anam_hermite->getPsiHn());
  _dbout->setArray(_iechOut, _iptrEst, condexp);

  /* Calculate the conditional variance */

  double condvar = hermiteCondStdElement(est, std, anam_hermite->getPsiHn());
  _dbout->setArray(_iechOut, _iptrStd, condvar);
}

int KrigingSystem::_getFlagAddress(int iech0, int ivar0)
{
  int nvar = _getNVar();

  int rank = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
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
  if (size == 1 && _matCL[0].empty()) return true;
  return false;
}
