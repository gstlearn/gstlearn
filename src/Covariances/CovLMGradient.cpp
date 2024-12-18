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
#include "Space/ASpace.hpp"
#include "Space/SpacePoint.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "Covariances/ACovAnisoList.hpp"
#include "Covariances/CovGradientFunctional.hpp"

CovLMGradient::CovLMGradient(const ASpace* space)
: ACovAnisoList(space)
{
}

CovLMGradient::CovLMGradient(const CovLMGradient &r)
: ACovAnisoList(r)
{

}

CovLMGradient::CovLMGradient(const ACovAnisoList& r)
    : ACovAnisoList()
{

  for (int icov = r.getCovaNumber()-1; icov >= 0; icov--)
  {
    const CovAniso *cov = r.getCova(icov);
    if (!cov->hasCovDerivative())
    {
      messerr("The covariance %s is not compatible with Gradients",
              cov->getCovName().c_str());
    }
    else
    {
      CovGradientFunctional* newcov = new CovGradientFunctional(*cov);
      addCov(dynamic_cast<const CovAniso*>(newcov));
      delete newcov;
    }
  }
  for (auto &e: _covs)
  {
    e->setOptimEnabled(false);
  }
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


void CovLMGradient::evalZAndGradients(const SpacePoint& p1,
                                      const SpacePoint& p2,
                                      double& covVal,
                                      VectorDouble& covGp,
                                      VectorDouble& covGG,
                                      const CovCalcMode* mode,
                                      bool flagGrad) const
{
  _initGradients(covVal, covGp, covGG, flagGrad);

  for (unsigned int i = 0, n = getCovaNumber(); i < n; i++)
  {
    ACovGradient* covloc = dynamic_cast<ACovGradient *>(_covs[i]);
    if (covloc != nullptr)
      covloc->evalZAndGradients(p1, p2, covVal, covGp, covGG, mode, flagGrad);
  }
}

void CovLMGradient::evalZAndGradients(const VectorDouble& vec,
                                      double& covVal,
                                      VectorDouble& covGp,
                                      VectorDouble& covGG,
                                      const CovCalcMode* mode,
                                      bool flagGrad) const
{
  /// TODO : Not true whatever the space
  SpacePoint p1(getOrigin(),-1);
  SpacePoint p2(getOrigin(),-1);
  p2.move(vec);

  evalZAndGradients(p1, p2, covVal, covGp, covGG, mode, flagGrad);
}

void CovLMGradient::addCov(const CovAniso* cov)
{
  // TODO This should be checked for some cases of Gradient (probably non numerical)
  const ACovGradient* covgrad = dynamic_cast<const ACovGradient*>(cov);
  if (covgrad == nullptr)
  {
    messerr("This covariance cannot be added");
    return;
  }
  cov->setOptimEnabled(false);
  ACovAnisoList::addCov(cov);

}

/**
 * Initialize the resulting arrays (coded explicitely for dimension 3)
 * @param covVal  Point covariance
 * @param covGp   Vector for Point-Gradient covariance
 * @param covGG   Vector for Gradient-Gradient covariance
 * @param flagGrad True if Gradient must be calculated
 */
void CovLMGradient::_initGradients(double& covVal,
                                   VectorDouble& covGp,
                                   VectorDouble& covGG,
                                   bool flagGrad)
{
  covVal = 0.;
  for (int i = 0; i < 3; i++) covGp[i] = 0.;
  if (flagGrad)
    for (int i = 0; i < 9; i++) covGG[i] = 0.;
}
