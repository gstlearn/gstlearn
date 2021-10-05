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
#include "Covariances/CovCosExp.hpp"

#include "math.h"
#include "Covariances/CovContext.hpp"

CovCosExp::CovCosExp(const CovContext& ctxt)
: ACovFunc(ECov::COSEXP, ctxt)
{
  setParam(1);
}

CovCosExp::CovCosExp(const CovCosExp &r)
: ACovFunc(r)
{
}

CovCosExp& CovCosExp::operator=(const CovCosExp &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovCosExp::~CovCosExp()
{
}

double CovCosExp::getScadef() const
{
  return (2.995732);
}

double CovCosExp::_evaluateCov(double h) const
{
  double cov = 1.;
  if (h > 100) return (0.);
  cov = exp(-h);
  double h2 = h / getParam();
  cov *= cos(2. * GV_PI * h2);
  return (cov);
}
