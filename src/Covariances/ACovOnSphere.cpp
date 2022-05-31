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
#include "Covariances/ACovOnSphere.hpp"

#include <math.h>

ACovOnSphere::ACovOnSphere()
{
}

ACovOnSphere::~ACovOnSphere()
{
}

double ACovOnSphere::evalCovOnSphere0(double scale, int degree) const
{
  double s0 = 0.;
  for (int i = 0; i < degree; i++)
  {
    s0 += evalSpectrumOnSphere(i-1, scale);
  }
  return s0;
}

double ACovOnSphere::evalCovOnSphere(double alpha, double scale, int degree) const
{
  double s = 0.;
  double u0 = 1.;
  double u2;
  double calpha = cos(alpha);
  double u1 = calpha;
  for (int i = 1; i < (degree + 2); i++)
  {
    u2 = 1. / (i + 1) * ((2 * i + 1) * calpha * u1 - i * u0);
    s += u0 * evalSpectrumOnSphere(i-1, scale);
    u0 = u1;
    u1 = u2;
  }
  return s;
}
