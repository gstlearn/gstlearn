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
#include "Covariances/CovP8.hpp"

#include "Covariances/CovContext.hpp"

CovP8::CovP8(const CovContext& ctxt)
: ACovFunc(ECov::P8, ctxt)
{
}

CovP8::CovP8(const CovP8 &r)
: ACovFunc(r)
{
}

CovP8& CovP8::operator=(const CovP8 &r)
{
  if (this != &r)
  {
    ACovFunc::operator =(r);
  }
  return *this;
}

CovP8::~CovP8()
{
}

double CovP8::_evaluateCov(double h) const
{
  double cov = 0.;
  double h2 = h * h;
  if (h < 1)
    cov = 3.
        - h2 * (28.
            - h2 * (210. - h * (448. - h * (420. - h * (192. - h * 35.)))));

  return (cov);
}

double CovP8::_evaluateCovDerivate(int degree, double h) const
{
  double h2, res;

  res = 0.;
  h2 = h * h;
  if (h > 1) return res;

  switch (degree)
  {
    case 1:
      res = -h
          * (56. - h2
              * (840. - h * (2240. - h * (2520. - h * (1344. - h * 280.)))));
      break;

    case 2:
      res = -56.
          + h2 * (2520. - h * (8960. - h * (12600. - h * (8064. - h * 1960.))));
      break;

    case 3:
      res = h * (5040. - h * (26880. - h * (50400. - h * (40320. - h * 11760.))));
      break;

    case 4:
      res = 5040. - h * (53760. - h * (151200. - h * (161280. - h * 58800.)));
      break;
  }
  return (res);
}
