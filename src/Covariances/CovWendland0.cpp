/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Covariances/CovWendland0.hpp"

#include "Covariances/CovContext.hpp"

CovWendland0::CovWendland0(const CovContext& ctxt)
: ACovFunc(ECov::WENDLAND0, ctxt)
{
}

CovWendland0::CovWendland0(const CovWendland0 &r)
: ACovFunc(r)
{
}

CovWendland0& CovWendland0::operator=(const CovWendland0 &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovWendland0::~CovWendland0()
{
}

double CovWendland0::_evaluateCov(double h) const
{
  double cov = 0.;
  double h2 = h * h;
  if (h < 1) cov = 1. - 2*h + h2;
  return (cov);
}

