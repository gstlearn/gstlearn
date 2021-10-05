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
#include "Covariances/CovStable.hpp"

#include <math.h>
#include "Covariances/CovContext.hpp"

CovStable::CovStable(const CovContext& ctxt)
: ACovFunc(ECov::STABLE, ctxt)
{
  setParam(1);
}

CovStable::CovStable(const CovStable &r)
: ACovFunc(r)
{
}

CovStable& CovStable::operator=(const CovStable &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovStable::~CovStable()
{
}

double CovStable::getScadef() const
{
  return pow(3., 1. / getParam());
}

double CovStable::_evaluateCov(double h) const
{
  double cov = 1.;
  if (h > 0) cov = exp(-pow(h, getParam()));

  return (cov);
}

