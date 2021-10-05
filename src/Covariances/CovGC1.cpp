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
#include "Covariances/CovGC1.hpp"

#include "Covariances/CovContext.hpp"

CovGC1::CovGC1(const CovContext& ctxt)
: ACovFunc(ECov::ORDER1_GC, ctxt)
{
}

CovGC1::CovGC1(const CovGC1 &r)
: ACovFunc(r)
{
}

CovGC1& CovGC1::operator=(const CovGC1 &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovGC1::~CovGC1()
{
}

double CovGC1::_evaluateCov(double h) const
{
  double cov;
  double r = getContext().getField();
  int ndim = getContext().getNDim();

  if (ndim == 1)
    cov = r - h;
  else if (ndim == 2)
    cov = r * GV_PI / 2 - h;
  else
    cov = r * 2 - h;

  return (cov);
}
