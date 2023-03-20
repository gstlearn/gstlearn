/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "Covariances/CovNugget.hpp"

#include "Covariances/CovContext.hpp"

CovNugget::CovNugget(const CovContext& ctxt)
: ACovFunc(ECov::NUGGET, ctxt)
{
}

CovNugget::CovNugget(const CovNugget &r)
: ACovFunc(r)
{
}

CovNugget& CovNugget::operator=(const CovNugget &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovNugget::~CovNugget()
{
}

double CovNugget::_evaluateCov(double h) const
{
  double cov = 0.;
  if (ABS(h) < 1.e-10) cov = 1.;
  return (cov);
}

String CovNugget::getFormula() const
{
  return "C(h)=\\delta(h)";
}
