/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Covariances/CovCosinus.hpp"

#include <math.h>
#include "Covariances/CovContext.hpp"

CovCosinus::CovCosinus(const CovContext& ctxt)
: ACovFunc(ECov::COSINUS, ctxt)
{
}

CovCosinus::CovCosinus(const CovCosinus &r)
: ACovFunc(r)
{
}

CovCosinus& CovCosinus::operator=(const CovCosinus &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovCosinus::~CovCosinus()
{
}

double CovCosinus::_evaluateCov(double h) const
{
  double cov = cos(2. * GV_PI * h);
  return (cov);
}

