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
#include "Basic/AFunctional.hpp"

class GSTLEARN_EXPORT FunctionalSpirale : public AFunctional
{
public:
  FunctionalSpirale();
  FunctionalSpirale(double a, double b, double c, double d, double sx, double sy);
  FunctionalSpirale(const FunctionalSpirale &m);
  FunctionalSpirale& operator=(const FunctionalSpirale &m);
  virtual ~FunctionalSpirale();

  virtual double getFunctionValue(const VectorDouble& coor) const override;

  VectorVectorDouble getFunctionVectors(const VectorDouble& coor) const;

private:
  double _linearCombination(double x, double y, double a, double b) const;

private:
  double _a;
  double _b;
  double _c;
  double _d;
  double _sx;
  double _sy;
};
