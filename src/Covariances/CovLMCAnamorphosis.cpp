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
#include "Covariances/ACov.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Enum/EAnam.hpp"
#include "Enum/ECalcMember.hpp"

#include "Space/ASpace.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"

#include <math.h>

CovLMCAnamorphosis::CovLMCAnamorphosis(const AAnam* anam,
                                       const VectorInt& strcnt,
                                       const CovContext& ctxt)
  : CovAnisoList(ctxt)
  , _activeFactor(0)
  , _anamStrCount()
  , _anam(anam)
{
  init(strcnt);
}

CovLMCAnamorphosis::CovLMCAnamorphosis(const CovAnisoList& lmc,
                                       const AAnam* anam,
                                       const VectorInt& strcnt)
  : CovAnisoList(lmc)
  , _activeFactor(0)
  , _anamStrCount()
  , _anam(anam)
{
  init(strcnt);
}

CovLMCAnamorphosis::CovLMCAnamorphosis(const CovLMCAnamorphosis& r)
  : CovAnisoList(r)
  , _activeFactor(r._activeFactor)
  , _anamStrCount(r._anamStrCount)
  , _anam(r._anam)
{
}

CovLMCAnamorphosis& CovLMCAnamorphosis::operator=(const CovLMCAnamorphosis& r)
{
  if (this != &r)
  {
    CovAnisoList::operator=(r);
    _activeFactor = r._activeFactor;
    _anamStrCount = r._anamStrCount;
    _anam         = r._anam;
  }
  return *this;
}

void CovLMCAnamorphosis::_loadAndAddEvalCovMatBiPointInPlace(
  MatrixSquareGeneral& mat,
  const SpacePoint& p1,
  const SpacePoint& p2,
  const CovCalcMode* mode) const
{
  // TODO: cannot replace by CovAnisoList???
  ACov::_loadAndAddEvalCovMatBiPointInPlace(mat, p1, p2, mode); 
}

void CovLMCAnamorphosis::_addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                                     const SpacePoint &pwork1,
                                                     const SpacePoint &pwork2,
                                                     const CovCalcMode *mode) const
{
  // TODO: cannot replace by CovAnisoList???
  ACov::_addEvalCovMatBiPointInPlace(mat, pwork1, pwork2, mode);
}

CovLMCAnamorphosis::~CovLMCAnamorphosis()
{
}

int CovLMCAnamorphosis::init(const VectorInt& anam_strcnt)
{
  for (auto &e: _covAnisos)
  {
    e->setOptimEnabled(false);
  }
  if (_anam == nullptr)
  {
    messerr("You must define 'anam' beforehand");
    return 1;
  }
  EAnam type = _anam->getType();
  if (type != EAnam::HERMITIAN && type != EAnam::DISCRETE_IR && type != EAnam::DISCRETE_DD)
  {
    messerr("Unknown Anamorphosis Definition for Model Transformation");
    messerr("It must be either 'HERMITIAN' or 'DISCRETE_IR' or 'DISCRETE_DD'");
    return 1;
  }
  if (type == EAnam::DISCRETE_IR)
  {
    int nfact = _anam->getNFactor();
    if ((int) anam_strcnt.size() != nfact)
    {
      messerr("Argument 'anam_strcnt' must be dimensioned to the number of factors (%d)",nfact);
      return 1;
    }
    int ncov = getNCov();
    for (int i=0; i<nfact; i++)
    {
      if (anam_strcnt[i] < 0 || anam_strcnt[i] >= ncov)
      {
        messerr("Argument 'anam_strcnt' must contain ranks of covariances of each factor");
        messerr("This value (%d) should lie within [1,ncov[ where ncov=%d",anam_strcnt[i],ncov);
        return 1;
      }
    }
    _anamStrCount = anam_strcnt;
  }

  return 0;
}

String CovLMCAnamorphosis::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << CovAnisoList::toString(strfmt);

  sstr << _anam->toString(strfmt);

  int iclass = getActiveFactor();
  if (iclass == -1)
    sstr << "Option switch to Raw Variable" << std::endl;
  else if (iclass > 0)
    sstr << "Active Factor: Rank" << iclass << std::endl;

  return sstr.str();
}

double CovLMCAnamorphosis::eval0(int ivar,
                                 int jvar,
                                 const CovCalcMode* mode) const
{
  if (_anam == nullptr) return TEST;

  const CovCalcMode* modeloc;
  if (mode == nullptr)
  {
    modeloc = new CovCalcMode();
  }
  else
  {
    modeloc = mode;
  }

  EAnam type = _anam->getType();
  double value = TEST;
  if (type == EAnam::HERMITIAN)
  {
    value = _evalHermite0(ivar, jvar, modeloc);
  }

  if (type == EAnam::DISCRETE_DD)
  {
    value = _evalDiscreteDD0(ivar, jvar, modeloc);
  }

  if (type == EAnam::DISCRETE_IR)
  {
    value = _evalDiscreteIR0(ivar, jvar, modeloc);
  }

  if (mode == nullptr) delete modeloc;

  return value;
}

double CovLMCAnamorphosis::eval(const SpacePoint& p1,
                                const SpacePoint& p2,
                                int ivar,
                                int jvar,
                                const CovCalcMode* mode) const
{
  if (_anam == nullptr) return TEST;

  const CovCalcMode* modeloc;
  if (mode == nullptr)
  {
    modeloc = new CovCalcMode();
  }
  else
  {
    modeloc = mode;
  }

  EAnam type = _anam->getType();
  double value = TEST;
  if (type == EAnam::HERMITIAN)
  {
    value = _evalHermite(ivar, jvar, p1, p2, modeloc);
  }

  if (type == EAnam::DISCRETE_DD)
  {
    value = _evalDiscreteDD(ivar, jvar, p1, p2, modeloc);
  }

  if (type == EAnam::DISCRETE_IR)
  {
    value = _evalDiscreteIR(ivar, jvar, p1, p2, modeloc);
  }

  if (mode == nullptr) delete modeloc;

  return value;
}

double CovLMCAnamorphosis::_evalHermite(int ivar,
                                        int jvar,
                                        const SpacePoint &p1,
                                        const SpacePoint &p2,
                                        const CovCalcMode *mode) const
{
  const AnamHermite *anamH = dynamic_cast<const AnamHermite*>(_anam);

  CovCalcMode modeloc(*mode);
  modeloc.setAsVario(false);

  double rho = 1.;
  if (getDistance(p1, p2) > 0.)
    rho = CovAnisoList::eval(p1, p2, ivar, jvar, &modeloc);
  double r = 1.;
  if (anamH->isChangeSupportDefined()) r = anamH->getRCoef();

  double cov = TEST;
  int iclass = getActiveFactor();

  if (iclass == 0)
  {

    // For the Gaussian variable

    cov = rho;
    if (mode->getAsVario()) cov = 1. - cov;
  }
  else if (iclass == -1)
  {

    // For the Raw variable

    cov = 0.;
    double rhon = 1.;
    double rn = 1.;
    double val = 0.;
    for (int jclass = 1; jclass < getAnamNClass(); jclass++)
    {
      rhon *= rho;
      rn   *= r;
      double psin = anamH->getPsiHn(jclass);
      val = rhon;
      if (mode->getAsVario()) val = 1. - val;
      switch (mode->getMember().getValue())
      {
        case ECalcMember::E_LHS:
          cov += psin * psin * val;
          break;

        case ECalcMember::E_RHS:
          cov += psin * psin * val / rn;
          break;

        case ECalcMember::E_VAR:
          cov += psin * psin * val;
          break;
      }
    }
  }
  else
  {

    // For the factor 'iclass'

    double rhon = pow(rho, (double) iclass);
    double rn = pow(r, (double) iclass);
    switch (mode->getMember().getValue())
    {
      case ECalcMember::E_LHS:
        cov =  rn * rn * rhon;
        break;

      case ECalcMember::E_RHS:
        cov = rn * rhon;
        break;

      case ECalcMember::E_VAR:
        cov = rhon;
        break;
    }
    if (mode->getAsVario()) cov = 1. - cov;
  }

  return cov;
}

double CovLMCAnamorphosis::_evalHermite0(int ivar,
                                         int jvar,
                                         const CovCalcMode* mode) const
{
  const AnamHermite *anamH = dynamic_cast<const AnamHermite*>(_anam);
  int iclass = getActiveFactor();

  double r = 1.;
  if (anamH->isChangeSupportDefined()) r = anamH->getRCoef();

  double cov = 0.;
  if (iclass == 0)
  {
    // For the Gaussian variable

    cov = CovAnisoList::eval0(ivar, jvar, mode);
  }
  else if (iclass == -1)
  {

    // For the Raw variable

    cov = 0.;
    double rn = 1.;
    for (int jclass = 1; jclass < getAnamNClass(); jclass++)
    {
      rn   *= r;
      double psin = anamH->getPsiHn(jclass);
      switch (mode->getMember().getValue())
      {
        case ECalcMember::E_LHS:
          cov += psin * psin / (rn * rn);
          break;

        case ECalcMember::E_RHS:
          cov += psin * psin  / rn;
          break;

        case ECalcMember::E_VAR:
          cov += psin * psin;
          break;
      }
    }
  }
  else
  {

    // For the factor 'iclass'

    cov = 1.;
  }

  return cov;
}

double CovLMCAnamorphosis::_evalDiscreteDD(int ivar,
                                           int jvar,
                                           const SpacePoint& p1,
                                           const SpacePoint& p2,
                                           const CovCalcMode* mode) const
{
  const AnamDiscreteDD *anamDD = dynamic_cast<const AnamDiscreteDD*>(_anam);
  int iclass = getActiveFactor();

  double gamma = 0.;
  if (getDistance(p1, p2) > 0.)
  {
    gamma = CovAnisoList::eval(p1, p1, ivar, jvar, mode) -
            CovAnisoList::eval(p1, p2, ivar, jvar, mode);
  }

  if (iclass == 0)
  {
    // Structure for the whole discretized variables

   double cov = 0.;
   for (int jclass = 1; jclass < getAnamNClass(); jclass++)
   {
     double li  = anamDD->getDDStatLambda(iclass);
     double csi = anamDD->getDDStatCnorm(iclass);
     double mui = anamDD->getDDStatMul(iclass);

     double coeff = 0.;
     switch (mode->getMember().getValue())
     {
       case ECalcMember::E_LHS:
         coeff = csi * csi;
         break;

       case ECalcMember::E_RHS:
         coeff = csi * csi / mui;
         break;

       case ECalcMember::E_VAR:
         coeff =  csi * csi;
         break;
     }
     cov += coeff * exp(-li * gamma);
   }
   return cov;
  }

  // Structure for the factor 'iclass'

  double li  = anamDD->getDDStatLambda(iclass);
  double mui = anamDD->getDDStatMul(iclass);

  double coeff = 0.;
  switch (mode->getMember().getValue())
  {
    case ECalcMember::E_LHS: return 1.; break;
    case ECalcMember::E_RHS: return mui; break;
    case ECalcMember::E_VAR: return 1.; break;
  }
  return coeff * exp(-li * gamma);
}

double CovLMCAnamorphosis::_evalDiscreteDD0(int /*ivar*/,
                                            int /*jvar*/,
                                            const CovCalcMode* mode) const
{
  if (mode == nullptr)
    messageAbort("In _evalHermite, mode MUST be defined");
  const AnamDiscreteDD *anamDD = dynamic_cast<const AnamDiscreteDD*>(_anam);
  int iclass = getActiveFactor();

  double cov = TEST;
  if (iclass == 0)
  {
    // Structure for the whole discretized variable

    cov = 0.;
    for (int jclass = 1; jclass < getAnamNClass(); jclass++)
    {
      double csi = anamDD->getDDStatCnorm(iclass);
      double mui = anamDD->getDDStatMul(iclass);

      double coeff = 0.;
      switch (mode->getMember().getValue())
      {
        case ECalcMember::E_LHS:
          coeff = csi * csi;
          break;
        case ECalcMember::E_RHS:
          coeff = csi * csi / mui;
          break;
        case ECalcMember::E_VAR:
          coeff = csi * csi;
          break;
      }
      cov += coeff;
    }
  }
  else
  {
    double mui = anamDD->getDDStatMul(iclass);

    double coeff = 0.;
    switch (mode->getMember().getValue())
    {
      case ECalcMember::E_LHS:
        coeff = 1.;
        break;

      case ECalcMember::E_RHS:
        coeff = mui;
        break;

      case ECalcMember::E_VAR:
        coeff = 1.;
        break;
    }
    cov = coeff;
  }
  return cov;
}

void CovLMCAnamorphosis::_transformCovCalcModeIR(CovCalcMode* mode, int iclass) const
{
  int from = 0;
  if (iclass > 0) from = _anamStrCount[iclass-1];
  mode->setActiveCovListFromInterval(from, _anamStrCount[iclass]);
}

double CovLMCAnamorphosis::_evalDiscreteIR(int ivar,
                                           int jvar,
                                           const SpacePoint& p1,
                                           const SpacePoint& p2,
                                           const CovCalcMode* mode) const
{
  if (mode == nullptr)
    messageAbort("In _evalHermite, mode MUST be defined");
  const AnamDiscreteIR *anamIR = dynamic_cast<const AnamDiscreteIR*>(_anam);
  int iclass = getActiveFactor();
  CovCalcMode modeloc(*mode);

  double r = 1.;
  bool flag_support = anamIR->isChangeSupportDefined();
  if (flag_support) r = anamIR->getRCoef();

  if (iclass == 0)
  {

    // Structure for the whole discretized variable

    double cov = 0.;
    double cov1 = 0.;
    double cov2 = 1.;
    for (int jclass = 1; jclass < getAnamNClass(); jclass++)
    {
      double bi = anamIR->getIRStatB(jclass);
      cov1 = cov2;
      _transformCovCalcModeIR(&modeloc, iclass);
      cov2 = pow(1. + CovAnisoList::eval(p1, p2, ivar, jvar, &modeloc) * anamIR->getIRStatR(jclass),r);
      cov += bi * bi * (cov2 - cov1);
    }
    return cov;
  }

  // Structure for the factor 'iclass´

  _transformCovCalcModeIR(&modeloc, iclass - 1);
  double cov1 = pow(1. + CovAnisoList::eval(p1, p2, ivar, jvar, &modeloc) *
                           anamIR->getIRStatR(iclass - 1),
                    r);
  _transformCovCalcModeIR(&modeloc, iclass);
  double cov2 = pow(1. + CovAnisoList::eval(p1, p2, ivar, jvar, &modeloc) *
                           anamIR->getIRStatR(iclass),
                    r);
  return (cov2 - cov1);
}

double CovLMCAnamorphosis::_evalDiscreteIR0(int /*ivar*/,
                                            int /*jvar*/,
                                            const CovCalcMode* mode) const
{
  if (mode == nullptr)
    messageAbort("In _evalHermite, mode MUST be defined");
  const AnamDiscreteIR *anamIR = dynamic_cast<const AnamDiscreteIR*>(_anam);
  int iclass = getActiveFactor();
  CovCalcMode modeloc(*mode);

  double r = 1.;
  if (anamIR->isChangeSupportDefined()) r = anamIR->getRCoef();

  if (iclass == 0)
  {

    // Structure for the whole discretized variable

    double cov = 0.;
    for (int jclass = 1; jclass < getAnamNClass(); jclass++)
    {
      double bi   = anamIR->getIRStatB(jclass);
      double cov2 = pow(anamIR->getIRStatR(jclass), r);
      cov += bi * bi * cov2;
    }
    return cov;
  }

  // Structure for the factor 'iclass´
  return pow(anamIR->getIRStatR(iclass - 1), r);
}

void CovLMCAnamorphosis::setActiveFactor(int anam_iclass)
{
  if (anam_iclass != 0 && anam_iclass > _anam->getNFactor())
  {
    messerr("The rank of the active factor (%d) is incorrect", anam_iclass);
    messerr("It should lie between 1 and the number of factors (%d)",
            _anam->getNFactor() - 1);
    messerr("or be set to 0 to estimate the whole discretized grade");
    messerr("The rank is set back to 0 (Gaussian Variable)");
    return;
  }
  _activeFactor = anam_iclass;
}

EAnam CovLMCAnamorphosis::getAnamType() const
{
  if (_anam == nullptr) return EAnam::UNKNOWN;
  return _anam->getType();
}

void CovLMCAnamorphosis::addCovAniso(const CovAniso* cov)
{
  // In this context, check that the Covariance is monovariate

  if (cov->getNVariables() != 1)
  {
    messerr("You can only add Monovariate Covariances in 'CovLMCAnamorphosis' "
            "object");
    messerr("Operation bypassed");
    return;
  }
  CovAnisoList::addCov(cov);
}
