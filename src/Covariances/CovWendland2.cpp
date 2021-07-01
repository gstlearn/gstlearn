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
#include "Covariances/CovWendland2.hpp"

#include "Covariances/CovContext.hpp"

CovWendland2::CovWendland2(const CovContext& ctxt)
: ACovFunc(COV_WENDLAND2, ctxt)
{
}

CovWendland2::CovWendland2(const CovWendland2 &r)
: ACovFunc(r)
{
}

CovWendland2& CovWendland2::operator=(const CovWendland2 &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovWendland2::~CovWendland2()
{
}

double CovWendland2::_evaluateCov(double h) const
{
  double cov = 0.;
  double h2 = h * h;
  if (h < 1)
    cov = 1
        - h2 * ((28. / 3.)
            - h2 * (70.
                - h * ((448. / 3.) - h * (140. - h * (64 - h * (35. / 3.))))));
  return (cov);
}

