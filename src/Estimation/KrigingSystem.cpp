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
#include "Estimation/KrigingSystem.hpp"

#include "geoslib_old_f.h"

#include "Enum/EKrigOpt.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/ELoc.hpp"
#include "Db/PtrGeos.hpp"
#include "Model/Model.hpp"
#include "Model/ANoStat.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighWork.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/Law.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Covariances/ECalcMember.hpp"
#include "Covariances/ACovAnisoList.hpp"

#include <math.h>

KrigingSystem::KrigingSystem(Db* dbin,
                             Db* dbout,
                             Model* model,
                             const ANeighParam* neighParam)
    : _dbin(dbin),
      _dbout(dbout),
      _model(model),
      _neighParam(neighParam),
      _isReady(false),
      _iptrEst(-1),
      _iptrStd(-1),
      _iptrVarZ(-1),
      _flagEst(false),
      _flagStd(false),
      _flagVarZ(false),
      _flagSimu(false),
      _calcul(EKrigOpt::PONCTUAL),
      _flagCode(false),
      _discreteMode(0),
      _ndisc(0),
      _ndiscs(),
      _disc1(),
      _disc2(),
      _xvalidEstim(false),
      _xvalidStdev(false),
      _rankColCok(),
      _flagBayes(false),
      _rmean(),
      _flagDGM(false),
      _supportCoeff(1.),
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
      _var0()
{
  _resetMemoryGeneral();
}

KrigingSystem::~KrigingSystem()
{
  OptDbg::setIndex(0); // Turn OFF this option for future task
}

int KrigingSystem::_getNVar() const
{
  if (_model == nullptr) return 0;
  return _model->getVariableNumber();
}

int KrigingSystem::_getNVarCL() const
{
  if (_matCL.empty())
    return _getNVar();
  else
    return _matCL.size();
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

int KrigingSystem::_getNech() const
{
  if (_dbin == nullptr) return 0;
  return _nbgh.size();
}

int KrigingSystem::_getNDim() const
{
  if (_dbin == nullptr) return 0;
  return _dbin->getNDim();
}

int KrigingSystem::_getNFex() const
{
  if (_model == nullptr) return 0;
  return _model->getExternalDriftNumber();
}

int KrigingSystem::_getNeq() const
{
  int nech = _getNech();
  int nvar = _getNVar();
  int nfeq = _getNFeq();
  int neq = nvar * nech + nfeq;
  if (! _rankColCok.empty()) neq += nvar;
  return neq;
}

int KrigingSystem::_getNDisc() const
{
  return _ndisc;
}

void KrigingSystem::_resetMemoryPerNeigh()
{
  int nvarCL = _getNVarCL();
  int neq    = _getNeq();
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
 **  or of the target (at iech_out if rank < 0)
 **
 ** \param[in]  loc_rank   Rank of the sample
 ** \param[in]  idim   Rank of the coordinate
 ** \param[in]  iech_out Rank of the target (if rank < 0)
 **
 *****************************************************************************/
double KrigingSystem::_getIdim(int loc_rank, int idim, int iech_out) const
{
  if (loc_rank >= 0)
  {
    return _dbin->getCoordinate(loc_rank, idim);
  }
  else
  {
    return _dbout->getCoordinate(iech_out, idim);
  }
}

/****************************************************************************/
/*!
 **  Returns the value of the external drift "rank" (if rank >= 0)
 **  or of the target (at iech_out if rank < 0)
 **
 ** \param[in]  rank   Rank of the sample
 ** \param[in]  ibfl   Rank of the external drift
 ** \param[in]  iech_out Rank of the target (used when rank < 0)
 **
 *****************************************************************************/
double KrigingSystem::_getFext(int rank, int ibfl, int iech_out) const
{
  if (rank >= 0)
  {
    return _dbin->getExternalDrift(rank, ibfl);
  }
  else
  {
    return _dbout->getExternalDrift(iech_out, ibfl);
  }
}

/****************************************************************************/
/*!
 **  Returns the value of the variable (at rank if rank >= 0)
 **  or of the target (at iech_out if rank < 0)
 **
 ** \param[in]  rank   Rank of the sample
 ** \param[in]  ivar   Rank of the variable
 ** \param[in]  iech_out Rank of the target (if rank < 0)
 **
 ** \remarks   In case of simulation, the variable of the first simulation
 ** \remarks   is systematically returned. This has no influence on the rest
 ** \remarks   of the calculations
 **
 *****************************************************************************/
double KrigingSystem::_getIvar(int rank,
                               int ivar,
                               int iech_out) const
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
    if (jvar < 0)
      return TEST;
    else
      return _dbout->getArray(iech_out, jvar);
  }
}

/****************************************************************************/
/*!
 **  Returns the value of the measurement error (at rank if rank >= 0)
 **  or of the target (at iech_out if rank < 0)
 **
 ** \param[in]  rank     Rank of the sample
 ** \param[in]  ivar     Rank of the variable
 ** \param[in]  iech_out Rank of the target
 **
 *****************************************************************************/
double KrigingSystem::_getVerr(int rank, int ivar, int iech_out) const
{
  if (rank >= 0)
  {
    return _dbin->getVarianceError(rank, ivar);
  }
  else
  {
    return _dbout->getVarianceError(iech_out, ivar);
  }
}
double KrigingSystem::_getMean(int ivarCL)
{
  double value = 0.;
  if (_matCL.empty())
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

void KrigingSystem::_setFlag(int iech, int ivar, int value)
{
  int nech = _getNech();
  _flag[iech + ivar * nech]= value;
}

int KrigingSystem::_getFlag(int iech, int ivar)
{
  int nech = _getNech();
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
  int nech = _getNech();
  int nvar = _getNVar();
  int nfeq = _getNFeq();
  int nbfl = _getNbfl();
  int next = _getNFex();
  int neq  = _getNeq();
  int ndim = _getNDim();
  for (int i = 0; i < neq; i++) _flag[i] = 1;

  /* Check on the coordinates */

  for (int iech = 0; iech < nech; iech++)
  {
    bool valid = true;
    for (int idim = 0; idim < _dbin->getNDim(); idim++)
      if (FFFF(_getIdim(_nbgh[iech], idim))) valid = false;
    if (! valid)
      for (int ivar = 0; ivar < nvar; ivar++)
        _setFlag(iech,ivar,0);
  }

  /* Check on the data values */

  for (int iech = 0; iech < nech; iech++)
    for (int ivar = 0; ivar < nvar; ivar++)
      if (FFFF(_getIvar(_nbgh[iech], ivar))) _setFlag(iech,ivar,0);

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
        if (_model->getCoefDrift(ivar, il, ib) == 0.) continue;
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
 ** \return  Error returned code: 1 if an error is found; 0 otherwise
 **
 *****************************************************************************/
bool KrigingSystem::_isAuthorized()
{
  int nech = _getNech();
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

void KrigingSystem::_covtabInit(bool flag_init)
{
  if (! flag_init) return;
  int nvar = _getNVar();
  for (int i = 0; i < nvar * nvar; i++) _covtab[i] = 0.;
}

void KrigingSystem::_covtabCalcul(bool flag_init,
                                  const CovCalcMode& mode,
                                  int iech1,
                                  int iech2,
                                  VectorDouble d1)
{
  // Load the non-stationary parameters if needed

  if (_model->isNoStat())
  {
    const ANoStat *nostat = _model->getNoStat();

    int jech1 = iech1;
    int icas1 = 1;
    if (iech1 < 0)
    {
      icas1 = 1;
      jech1 = _iechOut;
    }

    int jech2 = iech2;
    int icas2 = 1;
    if (iech2 < 0)
    {
      icas2 = 2;
      jech2 = _iechOut;
    }
    nostat->updateModel(_model, icas1, jech1, icas2, jech2);
  }

  // Evaluate the Model

  MatrixSquareGeneral mat = _model->getCovAnisoList()->evalNvarIpas(1., d1, VectorDouble(), mode);

  int nvar = _getNVar();
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      double value = mat.getValue(ivar, jvar);
      if (flag_init)
        _setCOVTAB(ivar,jvar,value);
      else
        _addCOVTAB(ivar,jvar,value);
    }
  return;
}

/****************************************************************************/
/*!
 **  Returns the additional variance for continuous moving neighborhood
 **
 ** \return  Variance mutipplier for Continuous Option
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
  int ndim = _getNDim();
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

void KrigingSystem::_drftabCalcul(const ECalcMember &member, const Db* db, int iech)
{
  VectorDouble drft = _model->evalDriftVec(db, iech, member);
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
  int nech   = _getNech();
  int nvar   = _getNVar();
  int nfeq   = _getNFeq();
  int nbfl   = _getNbfl();
  int neq    = _getNeq();
  int ndim   = _getNDim();
  for (int i = 0; i < neq * neq; i++) _lhs[i] = 0.;
  VectorDouble d1(ndim);

  CovCalcMode mode(ECalcMember::LHS);
  mode.update(0, 0, ECalcMember::LHS, -1, 0, 1);

  /* Establish the covariance part */

  for (int iech = 0; iech < nech; iech++)
    for (int jech = 0; jech < nech; jech++)
    {
      _covtabInit(true);
      for (int idim = 0; idim < ndim; idim++)
        d1[idim] = (_getIdim(_nbgh[jech], idim) - _getIdim(_nbgh[iech], idim));
      _covtabCalcul(false, mode, _nbgh[iech], _nbgh[jech], d1);

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

          // Correction in the DGM case

          if (_flagDGM)
          {
            if (iech != jech)
            _prodLHS(iech,ivar,jech,jvar,_supportCoeff * _supportCoeff);
          }
        }
    }

  /* Establish the drift part */

  if (nfeq <= 0 || nbfl <= 0) return;
  for (int iech = 0; iech < nech; iech++)
  {
    if (_nbgh[iech] >= 0)
      _drftabCalcul(ECalcMember::LHS, _dbin, _nbgh[iech]);
    else
      _drftabCalcul(ECalcMember::LHS, _dbout, _iechOut);

    for (int ivar = 0; ivar < nvar; ivar++)
      for (int ib = 0; ib < nfeq; ib++)
      {
        double value = 0.;
        for (int il = 0; il < nbfl; il++)
          value += _drftab[il] * _model->getCoefDrift(ivar, il, ib);
        _setLHS(iech,ivar,ib,nvar,value,true);
        _setLHS(ib,nvar,iech,ivar,value,true);
      }
  }
  return;
}

void KrigingSystem::_lhsIsoToHetero()

{
  int neq = _getNeq();
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
  int neq = _getNeq();
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
  int nech = _getNech();
  int neq = _getNeq();
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

  if (matrix_invert(_lhsinv.data(), _nred, -1))
  {
    messerr("The Kriging Matrix (%d,%d) is singular", _nred, _nred);
    messerr("One of the usual reason is the presence of duplicates");
    return 1;
  }
  return 0;
}

/****************************************************************************/
/*!
 **  Establish the kriging R.H.S
 **
 ** \param[in]  rankRandom Rank of the discretization for DGM case
 **
 ** \remarks When 'matCL' is provided, 'nvar' stands for the first dimension of
 ** \remarks the matrix 'matCL' (its second dimension is equal to model->getNVar()).
 ** \remarks Otherwise nvar designates model->getNVar()
 **
 *****************************************************************************/
int KrigingSystem::_rhsCalcul(int rankRandom)
{
  double nscale = 1;
  int nech   = _getNech();
  int nvar   = _getNVar();
  int nvarCL = _getNVarCL();
  int nbfl   = _getNbfl();
  int nfeq   = _getNFeq();
  int ndim   = _getNDim();
  int neq    = _getNeq();
  int ndisc  = _getNDisc();
  VectorDouble d1(ndim);

  CovCalcMode mode(ECalcMember::RHS);
  mode.update(0, 0, ECalcMember::RHS, -1, 0, 1);

  /* Establish the covariance part */

  for (int iech = 0; iech < nech; iech++)
  {
    _covtabInit(true);
    switch (_calcul.toEnum())
    {
      case EKrigOpt::E_PONCTUAL:
        nscale = 1;
        for (int idim = 0; idim < ndim; idim++)
        {
          d1[idim] = (_dbout->getCoordinate(_iechOut, idim) - _getIdim(_nbgh[iech], idim));
          // In the Point-Block Model, Randomize the target
          if (rankRandom >= 0 && ! _disc1.empty()) d1[idim] += _getDISC1(rankRandom, idim);
        }
        _covtabCalcul(false, mode, _nbgh[iech], -1, d1);
        break;

      case EKrigOpt::E_BLOCK:
        if (_discreteMode == 2) _blockDiscretize();
        nscale = ndisc;
        for (int i = 0; i < nscale; i++)
        {
          for (int idim = 0; idim < ndim; idim++)
            d1[idim] = (_dbout->getCoordinate(_iechOut, idim) - _getIdim(_nbgh[iech], idim)
                        + _getDISC1(i, idim));
          _covtabCalcul(false, mode, _nbgh[iech], -1, d1);
        }
        break;

      case EKrigOpt::E_DRIFT:
        nscale = 1;
        _covtabInit(true);
        break;
    }

    /* Normalization */

    double ratio = 1. / (double) nscale;
    if (_flagDGM) ratio *= _supportCoeff;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        _prodCOVTAB(ivar,jvar,ratio);

    /* Storage */

    if (_matCL.empty())
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

  _drftabCalcul(ECalcMember::RHS, _dbout, _iechOut);
  for (int il = 0; il < nbfl; il++)
    if (FFFF(_drftab[il])) return 1;

  if (_matCL.empty())
  {
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int ib = 0; ib < nfeq; ib++)
      {
        double value = 0.;
        for (int il = 0; il < nbfl; il++)
          value += _drftab[il] * _model->getCoefDrift(ivar, il, ib);
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
            value += _drftab[il] * _model->getCoefDrift(jvar, il, ib);
          value *= _matCL[ivarCL][jvar];
          _setRHS(ib,nvar,ivarCL,value,true);
        }
    }
  }
  return 0;
}

void KrigingSystem::_rhsIsoToHetero()
{
  int nvar = _getNVar();
  int neq  = _getNeq();

  int lec_rhs = 0;
  int ecr_rhs = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int i = 0; i < neq; i++, lec_rhs++)
    {
      if (_flag[i] == 0) continue;
      _rhs[ecr_rhs++] = _rhs[lec_rhs];
    }
  return;
}

void KrigingSystem::_rhsDump()
{
  int nech   = _getNech();
  int neq    = _getNeq();
  int nvarCL = _getNVarCL();
  int ndim   = _getNDim();
  int ndisc  = _getNDisc();
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
        message("%d", ndisc);
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
  int nech   = _getNech();
  int ndim   = _getNDim();
  int nvarCL = _getNVarCL();
  int nfeq   = _getNFeq();
  int ndisc  = _getNDisc();
  VectorDouble sum(nvarCL);
  char string[20];

  /* Header */

  mestitle(0, "(Co-) Kriging weights");

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
          tab_printg(NULL, _dbin->getBlockExtension(_nbgh[iech], idim));
      }
        tab_printg(NULL, _getIvar(_nbgh[iech], jvarCL));

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
  tab_prints(NULL, "Coeff");
  message("\n");

  /* Loop on the drift coefficients */

  cumflag = _nred - nfeq;
  for (int ib = 0; ib < nfeq; ib++)
  {
    int iwgt = ib + cumflag;
    tab_printi(NULL, ib + 1);
    double value = (status == 0) ? _zam[iwgt] : TEST;
    tab_printg(NULL, value);
    message("\n");
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
  int nfeq   = _getNFeq();
  int nvarCL = _getNVarCL();

  /* Estimation */

  if (_flagEst != 0)
  {
    for (int ivarCL = 0; ivarCL < nvarCL; ivarCL++)
    {
      double estim = 0.;
      if (nfeq <= 0) estim = _getMean(ivarCL);

      if (status == 0 && (_nred > 0 || nfeq <= 0 || _flagBayes))
      {
        if (_flagBayes)
          estim = _model->_evalDriftCoef(_dbout, _iechOut, ivarCL, _rmean.data());
        for (int i = 0; i < _nred; i++)
          estim += _getRHSC(i,ivarCL) * _zam[i];

        if (_neighParam->getFlagXvalid() && _flagEst > 0)
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

  if (_flagStd != 0)
  {
    for (int ivarCL = 0; ivarCL < nvarCL; ivarCL++)
    {
      double stdv = TEST;
      if (status == 0 && (_nred > 0 || nfeq <= 0 || _flagBayes))
      {
        stdv = _variance(ivarCL, ivarCL);
        if (stdv < 0) stdv = 0.;
        stdv = sqrt(stdv);

        if (_neighParam->getFlagXvalid() && _flagStd > 0)
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
  return;
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
  int ndim   = _getNDim();
  int ndisc  = _getNDisc();
  VectorDouble d1(ndim, 0.);

  CovCalcMode mode(ECalcMember::VAR);

  _covtabInit(true);
  switch (_calcul.toEnum())
  {
    case EKrigOpt::E_PONCTUAL:
      nscale = 1;
      _covtabCalcul(false, mode, -1, -1, d1);
      break;

    case EKrigOpt::E_BLOCK:
      nscale = ndisc;
      for (int i = 0; i < nscale; i++)
        for (int j = 0; j < nscale; j++)
        {
          for (int idim = 0; idim < ndim; idim++)
            d1[idim] = _getDISC1(i, idim) - _getDISC2(j,idim);
          _covtabCalcul(false, mode, -1, -1, d1);
       }
      nscale = nscale * nscale;
      break;

    case EKrigOpt::E_DRIFT:
      nscale = 1;
      break;
  }

  /* Normalization */

  double ratio = 1. / (double) nscale;
  if (_flagDGM) ratio *= _supportCoeff;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      _prodCOVTAB(ivar,jvar,ratio);

  /* Storage */

  if (_matCL.empty())
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
 ** \param[in]  varb    Matrix of variances
 **
 *****************************************************************************/
double KrigingSystem::_variance(int ivarCL,
                                int jvarCL,
                                const double* varb)
{
  int nvarCL = _getNVarCL();

  // In non-stationary case, the Variance at origin must be updated
  if (_model->isNoStat()) _variance0();

  double var = _getVAR0(ivarCL, jvarCL);
  if (_flagBayes) var += varb[jvarCL * nvarCL + ivarCL];
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
 ** \param[in]  flagLterm True if the LTerm must be calculated
 **
 ** \param[out] lterm  Product Z*C-1*Z
 **                    (only produced if FLAG_LTERM is true)
 **
 *****************************************************************************/
void KrigingSystem::_dual(bool flagLterm,
                          double *lterm)
{
  int nech = _getNech();
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
        mean = _model->getMean(ivar);
      if (_flagBayes)
        mean = _model->_evalDriftCoef(_dbout, _iechOut, ivar, _rmean.data());
      zext[ecr++] = _getIvar(_nbgh[iech], ivar) - mean;
    }
  }

  /* Operate the product : Z * A-1 */

  matrix_product(_nred, _nred, 1, _lhsinv.data(), zext.data(), _zam.data());

  /* Operate the product : Z * A-1 * Z */

  if (flagLterm)
    matrix_product(1, _nred, 1, _zam.data(), zext.data(), lterm);

  return;
}

/**
 * Performs the last operations before launching the loop on Estimations
 * @return
 */
bool KrigingSystem::isReady()
{
  if (! _isCorrect()) return false;
  if (_flagStd)
  {
    _variance0();
  }
  _isReady = true;
  return _isReady;
}

/**
 * Perform the Kriging of target iech_out
 *
 * @param iech_out Rank of the target
 * @return
 */
int KrigingSystem::estimate(int iech_out)
{
  double ldum;
  if (! _isReady)
  {
    messerr("You must call 'isReady' before launching 'estimate'");
    return 1;
  }

  // Store the Rank of the Target sample
  _iechOut = iech_out;

  int status = 0;
  OptDbg::setIndex(iech_out + 1);
  if (! _dbout->isActive(_iechOut)) return status;

  if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH) || OptDbg::query(EDbg::RESULTS))
  {
    mestitle(1, "Target location");
    db_sample_print(_dbout, iech_out, 1, 0, 0);
  }

  // Check the Neighborhood

  _nbgh = _nbghWork.select(_dbout,iech_out,_rankColCok);
  int nech = _getNech();
  status = (nech <= 0);

  /* Establish the Kriging L.H.S. */

  if (! _nbghWork.isUnchanged() || _neighParam->getFlagContinuous()
      || OptDbg::force())
  {
    status = _prepar();
    if (status) goto label_store;
    _dual(false, &ldum);
  }

  /* Establish the Kriging R.H.S. */

  _rhsCalcul();
  if (status) goto label_store;
  _rhsIsoToHetero();

  if (OptDbg::query(EDbg::KRIGING)) _rhsDump();

  /* Derive the kriging weights */

  if (_flagStd || _flagVarZ) _wgtCalcul();
  if (OptDbg::query(EDbg::KRIGING)) _wgtDump(status);

  /* Perform the estimation */

  label_store:
  _estimateCalcul(status);
  if (OptDbg::query(EDbg::RESULTS)) _krigingDump(status);
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
  int nvar = _getNVar();

  if (_neighParam->getFlagXvalid())
    mestitle(0, "Cross-validation results");
  else
    mestitle(0, "(Co-) Kriging results");
  message("Target Sample = %d\n", _iechOut + 1);

  /* Loop on the results */

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    if (_neighParam->getFlagXvalid())
    {
      message("Variable Z%-2d\n", ivar + 1);
      if (_iptrEst >= 0)
      {
        double trueval = (status == 0) ? _dbin->getVariable(_iechOut, ivar) : TEST;
        double estim = (status == 0) ? _dbout->getArray(_iechOut, _iptrEst + ivar) : TEST;

        double estval;
        double esterr;
        if (_xvalidEstim)
        {
          estval = (status == 0) ? estim + trueval : TEST;
          esterr = (status == 0) ? estim : TEST;
        }
        else
        {
          estval = (status == 0) ? estim : TEST;
          esterr = (status == 0) ? estim - trueval : TEST;
        }

        tab_printg(" - True value        = ", trueval);
        message("\n");
        tab_printg(" - Estimated value   = ", estval);
        message("\n");
        tab_printg(" - Estimation Error  = ", esterr);
        message("\n");

        if (_iptrStd >= 0)
        {
          double stdev = (status == 0) ? _dbout->getArray(_iechOut, _iptrStd + ivar) : TEST;
          double sterr;
          double sigma;

          if (_xvalidStdev)
          {
            sterr = stdev;
            sigma = (status == 0) ? esterr / stdev : TEST;
          }
          else
          {
            sigma = stdev;
            sterr = (status == 0) ? esterr / stdev : TEST;
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

int KrigingSystem::setKrigOptEstim(int iptrEst, int iptrStd, int iptrVarZ)
{
  _iptrEst = iptrEst;
  _iptrStd = iptrStd;
  _iptrVarZ = iptrVarZ;

  _flagEst = _iptrEst >= 0 || (_neighParam->getFlagXvalid() && _iptrStd >= 0);
  _flagStd = (_iptrStd >= 0);
  _flagVarZ = (_iptrVarZ >= 0);

  return 0;
}

void KrigingSystem::_blockDiscretize()
{
  int ndim = _getNDim();
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
          (_discreteMode == 1) ? dbgrid->getDX(idim) : _dbout->getBlockExtension(_iechOut, idim);
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
                                    const VectorInt& ndiscs)
{
  _calcul = calcul;
  if (_calcul == EKrigOpt::BLOCK)
  {
    _ndiscs = ndiscs;

    // Discretization calculated from Block extension variable (1) or Grid mesh (0)

    _discreteMode = 0;
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(_dbout);
    if (_dbout->getBlockExtensionNumber() > 0)
      _discreteMode = 2;
    else
    {
      if (dbgrid != nullptr) _discreteMode = 1;
    }
    if (_discreteMode <= 0)
    {
      messerr("Block estimation is impossible as Block is not defined");
      return 1;
    }

    // Core dimensioning

    int ndisc = _getNDisc();
    _disc1.resize(ndisc);
    _disc2.resize(ndisc);

    // For constant discretization, calculate the discretization coordinate immediately

    if (_discreteMode == 1) _blockDiscretize();
   }
  return 0;
}

int KrigingSystem::setKrigOptXValid(bool optionXValidEstim,
                                    bool optionXValidStdev)
{
  if (! _neighParam->getFlagXvalid()) return 0;
  _xvalidEstim = optionXValidEstim;
  _xvalidStdev = optionXValidStdev;
  return 0;
}
int KrigingSystem::setKrigOptColCok(const VectorInt& rank_colcok)
{
  if (rank_colcok.empty()) return 0;

  _rankColCok = rank_colcok;
  int ivar, jvar;
  int nvar = _getNVar();

  /* Loop on the ranks of the colocated variables */

  for (int ivar = 0; ivar < _getNVar(); ivar++)
  {
    int jvar = rank_colcok[ivar];
    if (IFFFF(jvar)) jvar = 0;
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

int KrigingSystem::setKrigOptBayes(bool flag_bayes)
{
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
  if (matCL.empty()) return 0;
  int n1 = matCL.size();
  int n2 = matCL[0].size();

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

void KrigingSystem::setKrigOptFlagSimu(bool flagSimu)
{
 _flagSimu = flagSimu;
 _nbghWork.setFlagSimu(flagSimu);
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

  /******************************************/
  /* Checking the Number of External Drifts */
  /******************************************/

  int nfex = 0;
  if (_dbout != nullptr)
  {
    if (nfex > 0 && nfex != _dbout->getExternalDriftNumber())
    {
      messerr("Incompatible Number of External Drifts of '_dbout'");
      return false;
    }
    nfex = _dbout->getExternalDriftNumber();
  }
  if (_dbin != nullptr)
  {
    if (nfex > 0)
    {
      if (_dbin->getExternalDriftNumber() != 0)
      {
        if (nfex != _dbin->getExternalDriftNumber())
        {
          messerr("Incompatible Number of External Drifts of '_dbin'");
          return false;
        }
        nfex = _dbin->getExternalDriftNumber();
      }
      else
      {
        if (_dbout == nullptr || ! _dbout->isGrid())
        {
          messerr("External Drift is not defined on '_dbin'");
          messerr("It cannot be interpolated from '_dbout' as it is not a Grid");
          return false;
        }
      }
    }
  }
  if (_model != nullptr)
  {
    if (nfex > 0 && nfex != _model->getExternalDriftNumber())
    {
      messerr("Incompatible NUmber of External Drifts of '_model'");
      return false;
    }
    nfex = _model->getExternalDriftNumber();
  }

  /*********************************/
  /* Calculate the field extension */
  /*********************************/

  VectorDouble dbin_mini;
  VectorDouble dbin_maxi;
  VectorDouble dbout_mini;
  VectorDouble dbout_maxi;

  if (_model != nullptr)
  {

    /* Input Db structure */

    if (_dbin != nullptr)
    {
      dbin_mini.resize(ndim,TEST);
      dbin_maxi.resize(ndim,TEST);
      if (db_extension(_dbin, dbin_mini.data(), dbin_maxi.data(), nullptr)) return false;
    }

    /* Output Db structure */

    if (_dbout != nullptr)
    {
      dbout_mini.resize(ndim,TEST);
      dbout_maxi.resize(ndim,TEST);
      if (db_extension(_dbout, dbout_mini.data(), dbout_maxi.data(), nullptr)) return false;
    }

    // Merge extensions

    if (_dbin != nullptr && _dbout != nullptr)
      _model->setField(
          ut_merge_extension(ndim, dbin_mini.data(), dbin_maxi.data(),
                             dbout_mini.data(), dbout_maxi.data()));
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
  return true;
}

// Shortcuts for addressing functions

int KrigingSystem::_IND(int iech, int ivar, int nech)
{
  return   ((iech) + (ivar) * nech);
}
int KrigingSystem::_getFLAG(int iech, int ivar)
{
  if (_flagCheckAddress)
  {
    _checkAddress("_getFLAG","iech",iech,_getNech());
    _checkAddress("_getFLAG","ivar",ivar,_getNVar());
  }
  int nech = _getNech();
  int ind  = _IND(iech, ivar, nech);
  return _flag[ind];
}
double KrigingSystem::_getCOVTAB(int ivar,int jvar)
{
  if (_flagCheckAddress)
  {
    _checkAddress("_getCOVTAB","ivar",ivar,_getNVar());
    _checkAddress("_getCOVTAB","jvar",jvar,_getNVar());
  }
  int nvar = _getNVar();
  return _covtab[(jvar) * nvar + (ivar)];
}
void KrigingSystem::_setCOVTAB(int ivar,int jvar,double value)
{
  if (_flagCheckAddress)
  {
    _checkAddress("_setCOVTAB","ivar",ivar,_getNVar());
    _checkAddress("_setCOVTAB","jvar",jvar,_getNVar());
  }
  int nvar = _getNVar();
  _covtab[(jvar) * nvar + (ivar)] = value;
}
void KrigingSystem::_addCOVTAB(int ivar,int jvar,double value)
{
  if (_flagCheckAddress)
   {
     _checkAddress("_addCOVTAB","ivar",ivar,_getNVar());
     _checkAddress("_addCOVTAB","jvar",jvar,_getNVar());
   }
  int nvar = _getNVar();
  _covtab[(jvar) * nvar + (ivar)] += value;
}
void KrigingSystem::_prodCOVTAB(int ivar,int jvar,double value)
{
  if (_flagCheckAddress)
   {
     _checkAddress("_prodCOVTAB","ivar",ivar,_getNVar());
     _checkAddress("_prodCOVTAB","jvar",jvar,_getNVar());
   }
  int nvar = _getNVar();
  _covtab[(jvar) * nvar + (ivar)] *= value;
}
double KrigingSystem::_getRHS(int iech, int ivar, int jvCL)
{
  if (_flagCheckAddress)
  {
    _checkAddress("_getRHS","iech",iech,_getNech());
    _checkAddress("_getRHS","ivar",ivar,_getNVar());
    _checkAddress("_getRHS","jvCL",jvCL,_getNVarCL());
  }
  int nech = _getNech();
  int neq = _getNeq();
  int ind  = _IND(iech, ivar, nech);
  return _rhs [ind + neq * (jvCL)];
}
void KrigingSystem::_setRHS(int iech, int ivar, int jvCL, double value, bool isForDrift)
{
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
        _checkAddress("_setRHS", "iech", iech, _getNech());
        _checkAddress("_setRHS", "ivar", ivar, _getNVar());
      }
    }
    else
    {
      _checkAddress("_setRHS", "iech", iech, _getNech());
      _checkAddress("_setRHS", "ivar", ivar, _getNVar());
    }
    _checkAddress("_setRHS", "jvCL", jvCL, _getNVarCL());
  }
  int nech = _getNech();
  int neq = _getNeq();
  int ind = _IND(iech,ivar,nech);
  _rhs [ind + neq * (jvCL)] = value;
}
double KrigingSystem::_getRHSC(int i, int jvCL)
{
  if (_flagCheckAddress)
  {
    _checkAddress("_getRHSC","i",i,_nred);
    _checkAddress("_getRHSC","jvCL",jvCL,_getNVarCL());
  }
  return _rhs [(i) + _nred * (jvCL)];
}
double KrigingSystem::_getWGTC(int i,int jvCL)
{
  if (_flagCheckAddress)
  {
    _checkAddress("_getWGTC","i",i,_nred);
    _checkAddress("_getWGTC","jvCL",jvCL,_getNVarCL());
  }
  return _wgt [(i) + _nred * (jvCL)];
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
        _checkAddress("_setLHS","iech",iech,_getNech());
        _checkAddress("_setLHS","ivar",ivar,_getNVar());
      }
      if (jvar == _getNVar())
      {
        _checkAddress("_setLHS","jb",jech,_getNFeq());
      }
      else
      {
        _checkAddress("_setLHS","jech",jech,_getNech());
        _checkAddress("_setLHS","jvar",jvar,_getNVar());
      }
    }
    else
    {
      _checkAddress("_setLHS","iech",iech,_getNech());
      _checkAddress("_setLHS","ivar",ivar,_getNVar());
      _checkAddress("_setLHS","jech",jech,_getNech());
      _checkAddress("_setLHS","jvar",jvar,_getNVar());
    }
  }
  int nech = _getNech();
  int neq = _getNeq();
  int indi = _IND(iech, ivar, nech);
  int indj = _IND(jech, jvar, nech);
  _lhs [indi + neq * indj] = value;
}
void KrigingSystem::_addLHS(int iech, int ivar, int jech, int jvar, double value)
{
  if (_flagCheckAddress)
  {
    _checkAddress("_addLHS","iech",iech,_getNech());
    _checkAddress("_addLHS","ivar",ivar,_getNVar());
    _checkAddress("_addLHS","jech",jech,_getNech());
    _checkAddress("_addLHS","jvar",jvar,_getNVar());
  }
  int nech = _getNech();
  int neq = _getNeq();
  int indi = _IND(iech, ivar, nech);
  int indj = _IND(jech, jvar, nech);
  _lhs [indi + neq * indj] += value;
}
void KrigingSystem::_prodLHS(int iech, int ivar, int jech, int jvar, double value)
{
  if (_flagCheckAddress)
  {
    _checkAddress("_prodLHS","iech",iech,_getNech());
    _checkAddress("_prodLHS","ivar",ivar,_getNVar());
    _checkAddress("_prodLHS","jech",jech,_getNech());
    _checkAddress("_prodLHS","jvar",jvar,_getNVar());
  }
  int nech = _getNech();
  int neq = _getNeq();
  int indi = _IND(iech, ivar, nech);
  int indj = _IND(jech, jvar, nech);
  _lhs [indi + neq * indj] *= value;
}
double KrigingSystem::_getLHS(int iech, int ivar, int jech, int jvar)
{
  if (_flagCheckAddress)
  {
    _checkAddress("_getLHS","iech",iech,_getNech());
    _checkAddress("_getLHS","ivar",ivar,_getNVar());
    _checkAddress("_getLHS","jech",jech,_getNech());
    _checkAddress("_getLHS","jvar",jvar,_getNVar());
  }
  int nech = _getNech();
  int neq = _getNeq();
  int indi = _IND(iech, ivar, nech);
  int indj = _IND(jech, jvar, nech);
  return _lhs [indi + neq * indj];
}
double KrigingSystem::_getDISC1(int idisc, int idim)
{
  if (_flagCheckAddress)
  {
    _checkAddress("_getDISC1","idisc",idisc,_ndisc);
    _checkAddress("_getDISC1","idim",idim,_getNDim());
  }
  return _disc1[(idim) * _ndisc + (idisc)];
}
void KrigingSystem::_setDISC1(int idisc, int idim, double value)
{
  if (_flagCheckAddress)
  {
    _checkAddress("_setDISC1","idisc",idisc,_ndisc);
    _checkAddress("_setDISC1","idim",idim,_getNDim());
  }
  _disc1[(idim) * _ndisc + (idisc)] = value;
}
double KrigingSystem::_getDISC2(int idisc,int idim)
{
  if (_flagCheckAddress)
  {
    _checkAddress("_getDISC2","idisc",idisc,_ndisc);
    _checkAddress("_getDISC2","idim",idim,_getNDim());
  }
  return _disc2[(idim) * _ndisc + (idisc)];
}
void KrigingSystem::_setDISC2(int idisc,int idim, double value)
{
  if (_flagCheckAddress)
  {
    _checkAddress("_setDISC2","idisc",idisc,_ndisc);
    _checkAddress("_setDISC2","idim",idim,_getNDim());
  }
  _disc2[(idim) * _ndisc + (idisc)] = value;
}
double KrigingSystem::_getVAR0(int ivCL, int jvCL)
{
  if (_flagCheckAddress)
  {
    _checkAddress("_getVAR0","ivCL",ivCL,_getNVarCL());
    _checkAddress("_getVAR0","jvCL",jvCL,_getNVarCL());
  }
  int nvarCL = _getNVarCL();
  return _var0[(jvCL) + nvarCL * (ivCL)];
}
void KrigingSystem::_setVAR0(int ivCL, int jvCL, double value)
{
  if (_flagCheckAddress)
  {
    _checkAddress("_setVAR0","ivCL",ivCL,_getNVarCL());
    _checkAddress("_setVAR0","jvCL",jvCL,_getNVarCL());
  }
  int nvarCL = _getNVarCL();
  _var0[(jvCL) + nvarCL * (ivCL)] = value;
}
void KrigingSystem::_checkAddress(const String& title,const String& theme,int ival,int nval)
{
  if (ival < 0)
    messageAbort("Error in %s: %s (%d) may not be negative",
                 title.c_str(),theme.c_str(),ival);
  if (ival >= nval)
    messageAbort("Error in %s: %s (%d) may not be larger or equal than %d",
                 title.c_str(),theme.c_str(),ival,nval);
}
