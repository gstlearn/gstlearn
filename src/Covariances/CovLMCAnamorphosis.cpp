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
#include "geoslib_f.h"

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

  /* Check if the distance is zero */

  int iclass = getAnamIClass();
  bool flag_support = anamH->isChangeSupportDefined();
  double rn = 1.;
  if (flag_support) rn = pow(anamH->getRCoef(), (double) iclass);

  double cov = 0.;
  double dist2 = getDistance(p1, p2);
  if (dist2 <= 0.)
  {
    cov = 1.;
  }
  else
  {
    cov = CovLMC::eval(ivar, jvar, p1, p2, mode);
    cov = pow(cov, (double) iclass);
  }
  if (flag_support) cov *= rn * rn;
  return cov;
}

double CovLMCAnamorphosis::_evalHermite0(int /*ivar*/,
                                         int /*jvar*/,
                                         const CovCalcMode& /*mode*/) const
{
  return 1.;
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

  double cov = 0.;
  double dist2 = getDistance(p1, p2);
  if (dist2 <= 0.)
  {
    cov = 1.;
  }
  else
  {
    double gamma = CovLMC::eval(ivar, jvar, p1, p2, mode);
    double li = anamDD->getDDStatLambda(iclass);
    cov = exp(- li * gamma);
  }
  return cov;
}

double CovLMCAnamorphosis::_evalDiscreteDD0(int /*ivar*/,
                                            int /*jvar*/,
                                            const CovCalcMode& /*mode*/) const
{
  return 1.;
}

double CovLMCAnamorphosis::_evalDiscreteIR(int ivar,
                                           int jvar,
                                           const SpacePoint& p1,
                                           const SpacePoint& p2,
                                           const CovCalcMode& mode) const
{
  int iclass = getAnamIClass();

  // Highlight the only covariances of the Model which are valid for current class

  CovCalcMode mode_loc(mode);
  mode_loc.setAllCovFiltered(getCovNumber(), true);
  int from = 0;
  if (iclass > 0) from = _anamStrCount[iclass-1];
  for (int i = from; i < _anamStrCount[iclass]; i++)
    mode_loc.setCovFiltered(i, false);

  double cov = 0.;
  double dist2 = getDistance(p1, p2);
  if (dist2 <= 0)
  {
    cov = 1.;
  }
  else
  {
    cov = eval(ivar, jvar, p1, p2, mode_loc);
  }
  return cov;
}

double CovLMCAnamorphosis::_evalDiscreteIR0(int /*ivar*/,
                                            int /*jvar*/,
                                            const CovCalcMode& /*mode*/) const
{
  return 1.;
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
