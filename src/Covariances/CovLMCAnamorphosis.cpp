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
#include "Space/ASpace.hpp"
#include "Basic/AException.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Covariances/CovLMC.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/EConvType.hpp"
#include "Covariances/EConvDir.hpp"
#include "Covariances/ECalcMember.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Anamorphosis/AnamDiscrete.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Anamorphosis/EAnam.hpp"

#include <math.h>

CovLMCAnamorphosis::CovLMCAnamorphosis(const AAnam* anam,
                                       const VectorInt& strcnt,
                                       const ASpace* space)
    : CovLMC(space),
      _anamIClass(0),
      _anamStrCount(),
      _anam(anam)
{
  init(strcnt);
}

CovLMCAnamorphosis::CovLMCAnamorphosis(const CovLMC* lmc,
                                       const AAnam* anam,
                                       const VectorInt& strcnt)
    : CovLMC(*lmc),
      _anamIClass(0),
      _anamStrCount(),
      _anam(anam)
{
  init(strcnt);
}

CovLMCAnamorphosis::CovLMCAnamorphosis(const CovLMCAnamorphosis &r)
    : CovLMC(r),
      _anamIClass(r._anamIClass),
      _anamStrCount(r._anamStrCount),
      _anam(r._anam)
{
}

CovLMCAnamorphosis& CovLMCAnamorphosis::operator=(const CovLMCAnamorphosis &r)
{
  if (this != &r)
  {
    CovLMC::operator=(r);
    _anamIClass = r._anamIClass;
    _anamStrCount = r._anamStrCount;
    _anam = r._anam;
  }
  return *this;
}

CovLMCAnamorphosis::~CovLMCAnamorphosis()
{
}

int CovLMCAnamorphosis::init(const VectorInt& anam_strcnt)
{
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
    int ncov = getCovNumber();
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

  sstr << ACovAnisoList::toString(strfmt);

  sstr << _anam->toString(strfmt);

  return sstr.str();
}

double CovLMCAnamorphosis::eval0(int ivar,
                                 int jvar,
                                 const CovCalcMode& mode) const
{
  if (_anam == nullptr) return TEST;
  EAnam type = _anam->getType();
  if (type == EAnam::HERMITIAN)
  {
    return _evalHermite0(ivar, jvar, mode);
  }

  if (type == EAnam::DISCRETE_DD)
  {
    return _evalDiscreteDD0(ivar, jvar, mode);
  }

  if (type == EAnam::DISCRETE_IR)
  {
    return _evalDiscreteIR0(ivar, jvar, mode);
  }

  return TEST;
}

double CovLMCAnamorphosis::eval(int ivar,
                                int jvar,
                                const SpacePoint& p1,
                                const SpacePoint& p2,
                                const CovCalcMode& mode) const
{
  if (_anam == nullptr) return TEST;
  EAnam type = _anam->getType();
  if (type == EAnam::HERMITIAN)
  {
    return _evalHermite(ivar, jvar, p1, p2, mode);
  }

  if (type == EAnam::DISCRETE_DD)
  {
    return _evalDiscreteDD(ivar, jvar, p1, p2, mode);
  }

  if (type == EAnam::DISCRETE_IR)
  {
    return _evalDiscreteIR(ivar, jvar, p1, p2, mode);
  }

  return TEST;
}

double CovLMCAnamorphosis::_evalHermite(int ivar,
                                        int jvar,
                                        const SpacePoint& p1,
                                        const SpacePoint& p2,
                                        const CovCalcMode& mode) const
{
  const AnamHermite *anamH = dynamic_cast<const AnamHermite*>(_anam);

  double rho = 1.;
  if (getDistance(p1, p2) > 0.)
    rho = CovLMC::eval(ivar, jvar, p1, p2, mode);
  double r = 1.;
  if (anamH->isChangeSupportDefined()) r = anamH->getRCoef();

  double cov = TEST;
  int iclass = getAnamIClass();
  if (iclass == 0)
  {

    // For the whole discretized variable
    cov = 0.;
    double rhon = 1.;
    double rn = 1.;
    for (int jclass = 1; jclass < getAnamNClass(); jclass++)
    {
      rhon *= rho;
      rn *= r;
      double psin = anamH->getPsiHn(jclass);
      switch (mode.getMember().getValue())
      {
        case ECalcMember::E_LHS:
          cov += psin * psin * rn * rn * rhon;
          break;

        case ECalcMember::E_RHS:
          cov += psin * psin * rn * rhon;
          break;

        case ECalcMember::E_VAR:
          cov += psin * psin * rhon;
          break;
      }
    }
  }
  else
  {

    // For the given factor 'iclass'
    double rhon = pow(rho, (double) iclass);
    double rn = pow(r, (double) iclass);
    switch (mode.getMember().getValue())
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
  }
  return cov;
}

double CovLMCAnamorphosis::_evalHermite0(int /*ivar*/,
                                         int /*jvar*/,
                                         const CovCalcMode& mode) const
{
  const AnamHermite *anamH = dynamic_cast<const AnamHermite*>(_anam);
  int iclass = getAnamIClass();

  if (mode.getMember().getValue() != ECalcMember::E_LHS)
    messageAbort("CovLMCAnamorphosis eval0");

  double cov = 1.;
  if (iclass == 0)
  {
    // For the whole discretized variables
    cov = 0.;
    for (int jclass = 1; jclass < getAnamNClass(); jclass++)
    {
      double psin = anamH->getPsiHn(jclass);
      cov += psin * psin;
    }
  }
  return cov;
}

double CovLMCAnamorphosis::_evalDiscreteDD(int ivar,
                                           int jvar,
                                           const SpacePoint& p1,
                                           const SpacePoint& p2,
                                           const CovCalcMode& mode) const
{
  const AnamDiscreteDD *anamDD = dynamic_cast<const AnamDiscreteDD*>(_anam);
  int iclass = getAnamIClass();

  CovCalcMode modeloc(mode);
  modeloc.setAsVario(true);

  double gamma = 0.;
  double dist2 = getDistance(p1, p2);
  if (dist2 > 0.)
    gamma = CovLMC::eval(ivar, jvar, p1, p2, modeloc);

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
     switch (mode.getMember().getValue())
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
  else
  {
    // Structure for the factor 'iclass'

    double li  = anamDD->getDDStatLambda(iclass);
    double mui = anamDD->getDDStatMul(iclass);

    double coeff = 0.;
    switch (mode.getMember().getValue())
    {
      case ECalcMember::E_LHS:
        return 1.;
        break;
      case ECalcMember::E_RHS:
        return mui;
        break;
      case ECalcMember::E_VAR:
        return 1.;
        break;
    }
    return coeff * exp(-li * gamma);
  }
  return TEST;
}

double CovLMCAnamorphosis::_evalDiscreteDD0(int /*ivar*/,
                                            int /*jvar*/,
                                            const CovCalcMode& mode) const
{
  const AnamDiscreteDD *anamDD = dynamic_cast<const AnamDiscreteDD*>(_anam);
  int iclass = getAnamIClass();

  double cov = TEST;
  if (iclass == 0)
  {
    // Structure for the whole discretized variable

    double cov = 0.;
    for (int jclass = 1; jclass < getAnamNClass(); jclass++)
    {
      double csi = anamDD->getDDStatCnorm(iclass);
      double mui = anamDD->getDDStatMul(iclass);

      double coeff = 0.;
      switch (mode.getMember().getValue())
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
    switch (mode.getMember().getValue())
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

void CovLMCAnamorphosis::_transformCovCalcModeIR(CovCalcMode& mode, int iclass) const
{
  mode.setAllCovFiltered(getCovNumber(), true);
  int from = 0;
  if (iclass > 0) from = _anamStrCount[iclass-1];
  for (int i = from; i < _anamStrCount[iclass]; i++)
    mode.setCovFiltered(i, false);
}

double CovLMCAnamorphosis::_evalDiscreteIR(int ivar,
                                           int jvar,
                                           const SpacePoint& p1,
                                           const SpacePoint& p2,
                                           const CovCalcMode& mode) const
{
  const AnamDiscreteIR *anamIR = dynamic_cast<const AnamDiscreteIR*>(_anam);
  int iclass = getAnamIClass();
  CovCalcMode mode_loc(mode);

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
      _transformCovCalcModeIR(mode_loc, iclass);
      cov2 = pow(1. + eval(ivar, jvar, p1, p2, mode_loc) * anamIR->getIRStatR(jclass),r);
      cov += bi * bi * (cov2 - cov1);
    }
    return cov;
  }
  else
  {

    // Structure for the factor 'iclass´

    _transformCovCalcModeIR(mode_loc, iclass - 1);
    double cov1 = pow(1. + eval(ivar, jvar, p1, p2, mode_loc) *
                      anamIR->getIRStatR(iclass - 1), r);
    _transformCovCalcModeIR(mode_loc, iclass);
    double cov2 = pow(1. + eval(ivar, jvar, p1, p2, mode_loc) *
                      anamIR->getIRStatR(iclass), r);
    return (cov2 - cov1);
  }
  return TEST;
}

double CovLMCAnamorphosis::_evalDiscreteIR0(int /*ivar*/,
                                            int /*jvar*/,
                                            const CovCalcMode& mode) const
{
  const AnamDiscreteIR *anamIR = dynamic_cast<const AnamDiscreteIR*>(_anam);
  int iclass = getAnamIClass();
  CovCalcMode mode_loc(mode);

  double r = 1.;
  if (anamIR->isChangeSupportDefined()) r = anamIR->getRCoef();

  if (iclass == 0)
  {

    // Structure for the whole discretized variable

    double cov = 0.;
    for (int jclass = 1; jclass < getAnamNClass(); jclass++)
    {
      double bi = anamIR->getIRStatB(jclass);
      double cov2 = pow(anamIR->getIRStatR(jclass), r);
      cov += bi * bi * cov2;
    }
    return cov;
  }
  else
  {

    // Structure for the factor 'iclass´

    return pow(anamIR->getIRStatR(iclass - 1), r);
  }
  return TEST;
}

int CovLMCAnamorphosis::setAnamIClass(int anam_iclass)
{
  if (! (anam_iclass == 0 || anam_iclass <= _anam->getNFactor()))
  {
    messerr("The rank of the active factor (%d) is incorrect", anam_iclass);
    messerr("It should lie between 1 and the number of factors (%d)", _anam->getNFactor() - 1);
    messerr("or be set to 0 to estimate the whole discretized grade");
    return 1;
  }
  _anamIClass = anam_iclass;
  return 0;
}

const EAnam CovLMCAnamorphosis::getAnamType() const
{
  if (_anam == nullptr) return EAnam::UNKNOWN;
  return _anam->getType();
}

void CovLMCAnamorphosis::addCov(const CovAniso* cov)
{
  // In this context, check that the Covariance is monovariate

  if (cov->getNVariables() != 1)
  {
    messerr("You can only add Monovariate Covariances in 'CovLMCAnamorphosis' object");
    messerr("Operation bypassed");
    return;
  }
  CovLMC::addCov(cov);
}
