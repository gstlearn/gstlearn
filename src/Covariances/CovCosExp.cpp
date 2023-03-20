/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
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
