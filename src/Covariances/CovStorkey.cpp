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
#include "Covariances/CovStorkey.hpp"

#include <math.h>
#include "Covariances/CovContext.hpp"

CovStorkey::CovStorkey(const CovContext& ctxt)
: ACovFunc(ECov::STORKEY, ctxt)
{
}

CovStorkey::CovStorkey(const CovStorkey &r)
: ACovFunc(r)
{
}

CovStorkey& CovStorkey::operator=(const CovStorkey &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovStorkey::~CovStorkey()
{
}

double CovStorkey::_evaluateCov(double h) const
{
  double cov = 0.;
  double pi2 = 2. * GV_PI;
  if (h < 1)
    cov = (2. * (1. - h) * (1. + cos(pi2 * h) / 2.) + 3 / pi2 * sin(pi2 * h)) / 3.;
  return (cov);
}

