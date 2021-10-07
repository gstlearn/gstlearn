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
#include "Covariances/CovLMGradient.hpp"

#include "Space/ASpace.hpp"
#include "Basic/AException.hpp"
#include "Space/SpacePoint.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/ACovAnisoList.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "geoslib_f.h"

CovLMGradient::CovLMGradient(const ASpace* space)
: ACovAnisoList(space)
{
}

CovLMGradient::CovLMGradient(const CovLMGradient &r)
: ACovAnisoList(r)
{
}

CovLMGradient& CovLMGradient::operator=(const CovLMGradient &r)
{
  if (this != &r)
  {
    ACovAnisoList::operator=(r);
  }
  return *this;
}

CovLMGradient::~CovLMGradient()
{
  /// TODO : Delete pointers ?
}

IClonable* CovLMGradient::clone() const
{
  return new CovLMGradient(*this);
}

void CovLMGradient::evalZAndGradients(const SpacePoint& p1,
                                      const SpacePoint& p2,
                                      double* covVal,
                                      VectorDouble& covGp,
                                      VectorDouble& covGg,
                                      const CovCalcMode& mode,
                                      bool flagGrad) const
{
  _initGradients(covVal, covGp, covGg, flagGrad);

  for (unsigned int i = 0, n = getCovNumber(); i < n; i++)
  {
    evalZAndGradients(p1, p2, covVal, covGp, covGg, mode, flagGrad);
  }
}

void CovLMGradient::evalZAndGradients(const VectorDouble& vec,
                                      double* covVal,
                                      VectorDouble& covGp,
                                      VectorDouble& covGg,
                                      const CovCalcMode& mode,
                                      bool flagGrad) const
{
  _initGradients(covVal, covGp, covGg, flagGrad);

  /// TODO : Not true whatever the space
  SpacePoint p1(getOrigin());
  SpacePoint p2(getOrigin());
  p2.move(vec);

  for (unsigned int i = 0, n = getCovNumber(); i < n; i++)
  {
    evalZAndGradients(p1, p2, covVal, covGp, covGg, mode, flagGrad);
  }
}

void CovLMGradient::addCov(const CovAniso* cov)
{
  const ACovGradient* covgrad = dynamic_cast<const ACovGradient*>(cov);
  if (covgrad == nullptr)
  {
    messerr("This covariance cannot be added");
    return;
  }
  ACovAnisoList::addCov(cov);
}

void CovLMGradient::_initGradients(double* covVal,
                                   VectorDouble& covGp,
                                   VectorDouble& covGg,
                                   bool flagGrad) const
{
  (*covVal) = 0.;
  for (int i = 0; i < 3; i++)
    covGp[i] = 0.;
  if (flagGrad) for (int i = 0; i < 9; i++)
    covGg[i] = 0.;
}

