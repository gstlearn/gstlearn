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
#include "Db/ELoc.hpp"
#include "Db/PtrGeos.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighWork.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/OptDbg.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Covariances/ECalcMember.hpp"

#include <math.h>

#define IND(iech,ivar)    ((iech) + (ivar) * nech)
#define FLAG(iech,ivar)   (_flag[IND(iech,ivar)])
#define COVTAB(ivar,jvar) (_covtab[(jvar) * nvar + (ivar)])
#define RHS(i,iv,jv)      (_rhs [IND(i,iv) + neq * (jv)])
#define LHS(i,iv,j,jv)    (_lhs [IND(i,iv) + neq * IND(j,jv)])
#define DISC(i,idim)      (_disc[(idim) * ndisc + (i)])

KrigingSystem::KrigingSystem(const Db* dbin,
                             Db* dbout,
                             const Model* model,
                             const ANeighParam* neighParam,
                             bool flagSimu)
    : _dbin(dbin),
      _dbout(dbout),
      _model(model),
      _neighParam(neighParam),
      _iptrEst(-1),
      _iptrStd(-1),
      _iptrVarZ(-1),
      _flagEst(false),
      _flagStd(false),
      _flagVarZ(false),
      _calcul(EKrigOpt::PONCTUAL),
      _flagCode(false),
      _disc(),
      _xvalidEstim(false),
      _xvalidStdev(false),
      _rankColCok(),
      _flagBayes(false),
      _rmean(),
      _flagDGM(false),
      _supportCoeff(1.),
      _flagWeights(false),
      _iechOut(-1),
      _nred(0),
      _nbghWork(_dbin, _neighParam, flagSimu),
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
  int nvar   = _getNVar();
  int ndrift = _getNDrift();
  _covtab.resize(nvar * nvar);
  _drftab.resize(ndrift);
  _var0.resize(nvar * nvar);
}

//KrigingSystem::KrigingSystem(const KrigingSystem &m)
//    : _dbin(m._dbin),
//      _dbout(m._dbout),
//      _model(m._model),
//      _neighParam(m._neighParam),
//      _iptrEst(m._iptrEst),
//      _iptrStd(m._iptrStd),
//      _iptrVarZ(m._iptrVarZ),
//      _flagEst(m._flagEst),
//      _flagStd(m._flagStd),
//      _flagVarZ(m._flagVarZ),
//      _calcul(m._calcul),
//      _flagCode(m._flagCode),
//      _disc(m._disc),
//      _xvalidEstim(m._xvalidEstim),
//      _xvalidStdev(m._xvalidStdev),
//      _rankColCok(m._rankColCok),
//      _flagBayes(m._flagBayes),
//      _rmean(m._rmean),
//      _flagDGM(m._flagDGM),
//      _supportCoeff(m._supportCoeff),
//      _flagWeights(m._flagWeights),
//      _iechOut(m._iechOut),
//      _nred(m._nred),
//      _nbghw(m._nbghw),
//      _nbgh(m._nbgh),
//      _flag(m._flag),
//      _covtab(m._covtab),
//      _drftab(m._drftab),
//      _lhs(m._lhs),
//      _lhsinv(m._lhsinv),
//      _rhs(m._rhs),
//      _wgt(m._wgt),
//      _zam(m._zam),
//      _var0(m._var0)
//{
//
//}
//
//KrigingSystem& KrigingSystem::operator=(const KrigingSystem &m)
//{
//  if (this != &m)
//  {
//    _dbin = m._dbin;
//    _dbout = m._dbout;
//    _model = m._model;
//    _neighParam = m._neighParam;
//    _iptrEst = m._iptrEst;
//    _iptrStd = m._iptrStd;
//    _iptrVarZ = m._iptrVarZ;
//    _flagEst = m._flagEst;
//    _flagStd = m._flagStd;
//    _flagVarZ = m._flagVarZ;
//    _calcul = m._calcul;
//    _flagCode = m._flagCode;
//    _disc = m._disc;
//    _xvalidEstim = m._xvalidEstim;
//    _xvalidStdev = m._xvalidStdev;
//    _rankColCok = m._rankColCok;
//    _flagBayes = m._flagBayes;
//    _rmean = m._rmean;
//    _flagDGM = m._flagDGM;
//    _supportCoeff = m._supportCoeff;
//    _flagWeights = m._flagWeights;
//    _iechOut = m._iechOut;
//    _nred = m._nred;
//    _nbghw = m._nbghw;
//    _nbgh = m._nbgh;
//    _flag = m._flag;
//    _covtab = m._covtab;
//    _drftab = m._drftab;
//    _lhs = m._lhs;
//    _lhsinv = m._lhsinv;
//    _rhs = m._rhs;
//    _wgt = m._wgt;
//    _zam = m._zam;
//    _var0 = m._var0;
//  }
//  return *this;
//}

KrigingSystem::~KrigingSystem()
{

}

int KrigingSystem::_getNVar() const
{
  if (_model == nullptr) return 0;
  return _model->getVariableNumber();
}

int KrigingSystem::_getNDrift() const
{
  if (_model == nullptr) return 0;
  return _model->getDriftNumber();
}

int KrigingSystem::_getNBfl() const
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
  int nbfl = _getNBfl();
  int neq = nvar * nech + nbfl;
  return neq;
}

int KrigingSystem::_getNDisc() const
{
  int ndim = _getNDim();
  int ndisc = (_disc.size() / ndim);
  return ndisc;
}

void KrigingSystem::_reset()
{
  int nvar   = _getNVar();
  int neq    = _getNeq();
  _lhs.resize(neq * neq);
  _rhs.resize(neq * nvar);
  _wgt.resize(neq * nvar);
  _zam.resize(neq);
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
 ** \param[in]  flagSimu: Case of Simulation
 **
 ** \remarks   In case of simulation, the variable of the first simulation
 ** \remarks   is systematically returned. This has no influence on the rest
 ** \remarks   of the calculations
 **
 *****************************************************************************/
double KrigingSystem::_getIvar(int rank,
                               int ivar,
                               int iech_out,
                               bool flagSimu) const
{
  if (rank >= 0)
  {

    // Variable in the Input file

    if (! flagSimu)

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
  int nbfl = _getNBfl();
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
      if (FFFF(_getIdim(_nbgh[iech], ivar))) _setFlag(iech,ivar,0);

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


  for (int ib = 0; ib < nbfl; ib++)
  {
    int valid = 0;
    for (int il = 0; il < _model->getDriftNumber(); il++)
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        if (_model->getCoefDrift(ivar, il, ib) == 0.) continue;
        for (int iech = 0; iech < nech; iech++)
          if (!FFFF(_getIvar(_nbgh[iech], ivar))) valid++;
      }
    _setFlag(nech+ib,nvar-1,(valid > 0));
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
  int nbfl = _getNBfl();

  /* Preliminary check */

  if (nech * nvar < nbfl) return false;

  /* Check that enough information is present */

  int n_cov = 0;
  for (int i = 0; i < nvar * nech; i++)
    n_cov += _flag[i];

  int n_drf = 0;
  for (int i = 0; i < nbfl; i++)
    n_drf += _flag[i + nvar * nech];

  if (n_cov <= 0 || n_cov < n_drf) return false;

  return true;
}

void KrigingSystem::_covtabInit(bool flag_init)
{
  int nvar = _getNVar();
  if (! flag_init) return;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        COVTAB(ivar,jvar)= 0.;
}

void KrigingSystem::_covtabCalcul(bool flag_init,
                                  const CovCalcMode& mode,
                                  int rank1,
                                  int rank2,
                                  VectorDouble d1)
{
  // Load the non-stationary parameters if needed

//  if (_model->isNoStat()) model_nostat_update(covint, model);

  // Evaluate the Model

  MatrixSquareGeneral mat = _model->getCovAnisoList()->evalNvarIpas(1., d1, VectorDouble(), mode);

  int nvar = _getNVar();
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
    {
      double value = mat.getValue(ivar, jvar);
      if (flag_init)
        COVTAB(ivar,jvar) = value;
      else
        COVTAB(ivar,jvar) += value;
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
  int nbfl   = _getNBfl();
  int ndrift = _getNDrift();
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
          LHS(iech,ivar,jech,jvar) = COVTAB(ivar, jvar);

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

                double cref = LHS(iech, ivar, jech, jvar);
                verr = cref * _continuousMultiplier(_nbgh[iech], _iechOut);
              }
            }
          }
          if (!FFFF(verr) && verr > 0) LHS(iech,ivar,jech,jvar) += verr;

          // Correction in the DGM case

          if (_flagDGM)
          {
            if (iech != jech)
            LHS(iech,ivar,jech,jvar) *= (_supportCoeff * _supportCoeff);
          }
        }
    }

  /* Establish the drift part */

  if (nbfl <= 0 || ndrift <= 0) return;
  for (int iech = 0; iech < nech; iech++)
  {
    if (_nbgh[iech] >= 0)
      _drftabCalcul(ECalcMember::LHS, _dbin, _nbgh[iech]);
    else
      _drftabCalcul(ECalcMember::LHS, _dbout, _iechOut);

    for (int ivar = 0; ivar < nvar; ivar++)
      for (int ib = 0; ib < nbfl; ib++)
      {
        double value = 0.;
        for (int il = 0; il < ndrift; il++)
          value += _drftab[il] * _model->getCoefDrift(ivar, il, ib);
        LHS(iech,ivar,ib,nvar) = value;
        LHS(ib,nvar,iech,ivar) = value;
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

  // Shrink it
  int neq = _getNeq();
  if (_nred < neq)
    _lhsinv.resize(_nred * _nred);

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
  int ndrift = _getNDrift();
  int nfeq   = _getNFex();
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
          // The next option is plugged for the case of target randomization
          // for the case of Point-Block Model
          if (rankRandom >= 0 && ! _disc.empty()) d1[idim] += DISC(rankRandom, idim);
        }
        _covtabCalcul(false, mode, _nbgh[iech], -1, d1);
        break;

      case EKrigOpt::E_BLOCK:
        nscale = ndisc;
        for (int i = 0; i < nscale; i++)
        {
          for (int idim = 0; idim < ndim; idim++)
            d1[idim] = (_dbout->getCoordinate(_iechOut, idim) - _getIdim(_nbgh[iech], idim)
                        + DISC(i, idim));
          _covtabCalcul(false, mode, _nbgh[iech], -1, d1);
        }
        break;

      case EKrigOpt::E_DRIFT:
        nscale = 1;
        for (int ivar = 0; ivar < nvar; ivar++)
          for (int jvar = 0; jvar < nvar; jvar++)
            COVTAB(ivar,jvar) = 0;
        break;
    }

    /* Normalization */

    double ratio = 1. / (double) nscale;
    if (_flagDGM) ratio *= _supportCoeff;
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar < nvar; jvar++)
        COVTAB(ivar,jvar) *= ratio;

    for (int jvar = 0; jvar < nvar; jvar++)
      for (int ivar = 0; ivar < nvar; ivar++)
        RHS(iech,ivar,jvar) = COVTAB(ivar, jvar);
  }

  /* Establish the drift part */

  if (nfeq <= 0) return 0;

  _drftabCalcul(ECalcMember::RHS, _dbout, _iechOut);
  for (int il = 0; il < ndrift; il++)
    if (FFFF(_drftab[il])) return 1;

  for (int ivar = 0; ivar < nvar; ivar++)
    for (int ib = 0; ib < nfeq; ib++)
    {
      double value = 0.;
      for (int il = 0; il < ndrift; il++)
        value += _drftab[il] * _model->getCoefDrift(ivar, il, ib);
      RHS(ib,nvar,ivar) = value;
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
  int nech  = _getNech();
  int neq   = _getNeq();
  int nvar  = _getNVar();
  int ndim  = _getNDim();
  int ndisc = _getNDisc();
  VectorInt rel = _getRelativePosition();

  /* General Header */

  mestitle(0, "RHS of Kriging matrix (compressed)");
  if (nech > 0) message("Number of active samples    = %d\n", nech);
  message("Total number of equations   = %d\n", neq);
  message("Reduced number of equations = %d\n", _nred);
  message("Number of right-hand sides  = %d\n", nvar);

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
  for (int ivar = 0; ivar < nvar; ivar++)
    tab_printi(NULL, ivar + 1);
  message("\n");

  /* Matrix lines */

  for (int i = 0; i < _nred; i++)
  {
    tab_printi(NULL, i + 1);
    if (! _flag.empty()) tab_printi(NULL, rel[i]);
    for (int ivar = 0; ivar < nvar; ivar++)
      tab_printg(NULL, _rhs[(i) + _nred * (ivar)]);
    message("\n");
  }
  return;
}

void KrigingSystem::_wgtDump(int status)
{
  if (! _flagWeights) return;

  int nech  = _getNech();
  int ndim  = _getNDim();
  int nvar  = _getNVar();
  int nfeq  = _getNDrift();
  int ndisc = _getNDisc();
  VectorDouble sum(nvar);
  char string[20];

  /* Calculate the weights */

  VectorDouble wgt(_rhs.size());
  matrix_product(_nred, _nred, nvar, _lhsinv.data(), _rhs.data(), wgt.data());

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
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    (void) gslSPrintf(string, "Z%d*", ivar + 1);
    tab_prints(NULL, string);
  }
  message("\n");

  /* Display the information and the weights */

  int lec = 0;
  int cumflag = 0;
  for (int jvar = 0; jvar < nvar; jvar++)
  {
    if (nvar > 1) message("Using variable Z%-2d\n", jvar + 1);

    /* Loop on the samples */

    for (int ivar = 0; ivar < nvar; ivar++)
      sum[ivar] = 0.;
    for (int iech = 0; iech < nech; iech++, lec++)
    {
      int flag_value = (! _flag.empty()) ? _flag[lec] : 1;
      tab_printi(NULL, iech + 1);
      for (int idim = 0; idim < ndim; idim++)
        tab_printg(NULL, _getIdim(_nbgh[iech], idim));
      if (_dbin->hasCode())
        tab_printg(NULL, _dbin->getCode(_nbgh[iech]));
      if (_dbin->getVarianceErrorNumber() > 0)
        tab_printg(NULL, _getVerr(_nbgh[iech], (_flagCode) ? 0 : jvar));
      if (ndisc > 0)
      {
        for (int idim = 0; idim < ndim; idim++)
          tab_printg(NULL, _dbin->getBlockExtension(_nbgh[iech], idim));
      }
        tab_printg(NULL, _getIvar(_nbgh[iech], jvar));

      for (int ivar = 0; ivar < nvar; ivar++)
      {
        int iwgt = _nred * ivar + cumflag;
        double value = (! wgt.empty() && status == 0 && flag_value) ? wgt[iwgt] : TEST;
        if (!FFFF(value)) sum[ivar] += value;
        tab_printg(NULL, value);
      }
      if (flag_value) cumflag++;
      message("\n");
    }

    int number = 1 + ndim + 1;
    if (_dbin->getVarianceErrorNumber() > 0) number++;
    if (ndisc > 0) number += ndim;
    tab_prints(NULL, "Sum of weights", number, EJustify::LEFT);
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      double value = (status == 0) ? sum[ivar] : TEST;
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
  int nfeq = _getNBfl();
  int nvar = _getNVar();

  /* Estimation */

  if (_flagEst != 0)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      double estim = 0.;
      if (nfeq <= 0) estim = _model->getMean(ivar);

      if (status == 0 && (_nred > 0 || nfeq <= 0 || _flagBayes))
      {
        if (_flagBayes)
          estim = _model->_evalDriftCoef(_dbout, _iechOut, ivar, _rmean.data());
        for (int i = 0; i < _nred; i++)
          estim += _rhs[i + _nred * ivar] * _zam[i];

        if (_neighParam->getFlagXvalid() && _flagEst > 0)
        {
          double valdat = _dbin->getVariable(_iechOut, ivar);
          estim = (FFFF(valdat)) ? TEST : estim - valdat;
        }
      }
      else
      {
        // In case of failure with KS, set the result to mean
        if (nfeq > 0) estim = TEST;
      }
      _dbout->setArray(_iechOut, _iptrEst + ivar, estim);
    }
  }

  /* Variance of the estimation error */

  if (_flagStd != 0)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      double stdv = TEST;
      if (status == 0 && (_nred > 0 || nfeq <= 0 || _flagBayes))
      {
        stdv = _variance(ivar, ivar);
        if (stdv < 0) stdv = 0.;
        stdv = sqrt(stdv);

        if (_neighParam->getFlagXvalid() && _flagStd > 0)
        {
          double estim = _dbout->getArray(_iechOut, _iptrEst + ivar);
          stdv = (FFFF(estim) || stdv <= 0.) ? TEST : estim / stdv;
        }
      }
      _dbout->setArray(_iechOut, _iptrStd + ivar, stdv);
    }
  }

  /* Variance of the estimator */

  if (_flagVarZ != 0)
  {
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      double varZ = TEST;
      if (status == 0 && (_nred > 0 || nfeq <= 0))
        varZ = _estimateVarZ(ivar, ivar);
      _dbout->setArray(_iechOut, _iptrVarZ + ivar, varZ);
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Establish the variance of the estimator
 **
 ** \param[in]  ivar    Rank of the target variable
 ** \param[in]  jvar    Rank of the auxiliary variable
 **
 *****************************************************************************/
double KrigingSystem::_estimateVarZ(int ivar, int jvar)
{
  int nfeq = _getNDrift();
  int cumflag = _nred - nfeq;

  double var = 0.;
  for (int i = 0; i < _nred; i++)
  {
    double signe = (i < cumflag) ? 1. : -1.;
    var += signe * _rhs[i + _nred * jvar] * _wgt[i + _nred * ivar];
  }
  return var;
}

/****************************************************************************/
/*!
 **  Establish the calculation of variance or standard deviation
 **
 ** \param[in]  ivar    Rank of the target variable
 ** \param[in]  jvar    Rank of the auxiliary variable
 ** \param[in]  varb    Matrix of variances
 **
 *****************************************************************************/
double KrigingSystem::_variance(int ivar,
                                int jvar,
                                const double* varb)
{
  int nvar = _getNVar();

  // In the stationary case, var0 has already been calculated
  // Otherwise if must be derived here

//  if (model->isNoStat() || model->hasExternalCov())
//  {
//    COVINT.setIcas1(2);
//    COVINT.setIcas2(2);
//    COVINT.setIech1(iech_out);
//    COVINT.setIech2(iech_out);
//    COVINT.setNdim(model->getDimensionNumber());
//    model_variance0_nostat(model, KOPTION, &COVINT, covtab, var0);
//  }

  double var = _var0[jvar * nvar + ivar];
  if (_flagBayes) var += varb[jvar * nvar + ivar];
  for (int i = 0; i < _nred; i++)
    var -= _rhs[i + _nred * jvar] * _wgt[i + _nred * ivar];

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

  _reset();

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
  int nfeq = _getNBfl();

  /* Set the whole array to 0 */

  VectorDouble zext(_nred);
  for (int i = 0; i < _nred; i++) zext[i] = 0.;

  /* Extract the data */

  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    for (int iech = 0; iech < nech; iech++)
    {
      if (! FLAG(iech, ivar)) continue;
      double mean = 0.;
      if (nfeq <= 0) mean = _model->getMean(ivar);
      if (_flagBayes)
        mean = _model->_evalDriftCoef(_dbout, _iechOut, ivar, _rmean.data());
      zext[ecr++] = _getIvar(_nbgh[iech], ivar) - mean;
    }
  }

  /* Operate the product : Z * A-1 */

  matrix_product(_nred, _nred, 1, _lhs.data(), zext.data(), _zam.data());

  /* Operate the product : Z * A-1 * Z */

  if (flagLterm)
    matrix_product(1, _nred, 1, _zam.data(), zext.data(), lterm);

  return;
}

/**
 * Perform the Kriging of target iech_out
 *
 * @param iech_out Rank of the target
 * @param verbose  Verbose option
 * @return
 */
int KrigingSystem::estimate(int iech_out, bool verbose)
{
  double ldum;

  // Store the Rank of the Target sample
  _iechOut = iech_out;

  int status = 0;
  OptDbg::setIndex(iech_out + 1);
  if (! _dbout->isActive(_iechOut)) return status;

  if (OptDbg::query(EDbg::KRIGING) || OptDbg::query(EDbg::NBGH)
      || OptDbg::query(EDbg::RESULTS))
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

  /* Establish the kriging R.H.S. */

  _rhsCalcul();
  if (status) goto label_store;
  _rhsIsoToHetero();

  if (OptDbg::query(EDbg::KRIGING)) _rhsDump();

  /* Derive the kriging weights */

  if (_flagWeights) _wgtDump(status);

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

int KrigingSystem::setKrigOptCalcul(const EKrigOpt& calcul,
                                    const VectorDouble& disc)
{
  _calcul = calcul;
  if (_calcul == EKrigOpt::BLOCK)
  {
    _disc = disc;
  }
  return 0;
}

int KrigingSystem::setKrigOptXValid(bool optionXValidEstim,
                                    bool optionXValidStdev)
{
  if (! _neighParam->getFlagXvalid())
  {
    messerr("This method should not be used here");
    messerr("You are not running a Cross-Validation Option");
    return 1;
  }
  return 0;
}
int KrigingSystem::setKrigOptColCok(const VectorInt& rank_colcok)
{
  _rankColCok = rank_colcok;
  return 0;
}
int KrigingSystem::setKrigOptBayes(bool flag_bayes)
{
  _flagBayes = flag_bayes;
}
