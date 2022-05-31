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
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/Vector.hpp"

/* Covariance basic function expressed on a Sphere
 * */

class GSTLEARN_EXPORT ACovOnSphere
{
public:
  ACovOnSphere();
  virtual ~ACovOnSphere();

  virtual double evalSpectrumOnSphere(int degree, double scale) const = 0;

  double evalCovOnSphere0(double scale, int degree = 50) const;
  double evalCovOnSphere(double alpha, double scale, int degree = 50) const;
};
