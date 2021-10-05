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

