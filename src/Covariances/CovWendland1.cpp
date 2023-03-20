/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Covariances/CovWendland1.hpp"

#include "Covariances/CovContext.hpp"

CovWendland1::CovWendland1(const CovContext& ctxt)
: ACovFunc(ECov::WENDLAND1, ctxt)
{
}

CovWendland1::CovWendland1(const CovWendland1 &r)
: ACovFunc(r)
{
}

CovWendland1& CovWendland1::operator=(const CovWendland1 &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovWendland1::~CovWendland1()
{
}

double CovWendland1::_evaluateCov(double h) const
{
  double cov = 0.;
  double h2 = h * h;
  if (h < 1) cov = 1. - h2 * (10. - h * (20. - h * (15. - h * 4)));
  return (cov);
}

