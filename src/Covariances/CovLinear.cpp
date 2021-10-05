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
#include "Covariances/CovLinear.hpp"

#include "Covariances/CovContext.hpp"

CovLinear::CovLinear(const CovContext& ctxt)
: ACovFunc(ECov::LINEAR, ctxt)
{
}

CovLinear::CovLinear(const CovLinear &r)
: ACovFunc(r)
{
}

CovLinear& CovLinear::operator=(const CovLinear &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovLinear::~CovLinear()
{
}

double CovLinear::_evaluateCov(double h) const
{
  double cov, r;

  r = getContext().getField();

  int ndim = getContext().getNDim();
  if (ndim == 1)
    cov = r - h;
  else if (ndim == 2)
    cov = r * GV_PI / 2 - h;
  else
    cov = r * 2 - h;

  return (cov);
}
