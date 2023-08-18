/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Covariances/CovWendland2.hpp"

#include "Covariances/CovContext.hpp"

CovWendland2::CovWendland2(const CovContext& ctxt)
: ACovFunc(ECov::WENDLAND2, ctxt)
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

double CovWendland2::_evaluateCovDerivative(int degree, double h) const
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
