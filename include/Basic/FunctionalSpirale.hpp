/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
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
