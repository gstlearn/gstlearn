/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Covariances/CovCauchyGen.hpp"
#include "Covariances/CovContext.hpp"

#include <cmath>

namespace gstlrn
{
CovCauchyGen::CovCauchyGen(const CovContext& ctxt)
  : ACovFunc(ECov::CAUCHY, ctxt)
{
  setParam(1);
  setAlpha(2.0);
}

CovCauchyGen::CovCauchyGen(const CovCauchyGen& r)
  : ACovFunc(r)
  , _alpha(r._alpha)
{
}

CovCauchyGen& CovCauchyGen::operator=(const CovCauchyGen& r)
{
  if (this != &r)
  {
    ACovFunc::operator=(r);
    _alpha = r._alpha;
  }
  return *this;
}

CovCauchyGen::~CovCauchyGen()
{
}

bool CovCauchyGen::setAlpha(double alpha)
{
    if((alpha > 0)&&(alpha <= 2.0)){
        _alpha = alpha;
        return true;
    }
  _alpha = 0.0;
  return false;
}

double CovCauchyGen::getScadef() const
{
  return pow(pow(20., 1. / getParam()) - 1., 1. / getAlpha());
}

double CovCauchyGen::_evaluateCov(double h) const
{
  double cov = 1. / pow(1. + pow(abs(h), getAlpha()), getParam());
  return (cov);
}

String CovCauchyGen::getFormula() const
{
  return "C(h)=\\frac{1}{\\left( 1+ \\frac{h^\\alpha}{a_t^\\alpha} \\right)^\\beta";
}
} // namespace gstlrn