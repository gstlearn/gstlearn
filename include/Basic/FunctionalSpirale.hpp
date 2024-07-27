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
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/AFunctional.hpp"

class MatrixSquareGeneral;
class CovAniso;

class GSTLEARN_EXPORT FunctionalSpirale : public AFunctional
{
public:
  FunctionalSpirale();
  FunctionalSpirale(double a, double b, double c, double d, double sx, double sy);
  FunctionalSpirale(const FunctionalSpirale &m);
  FunctionalSpirale& operator=(const FunctionalSpirale &m);
  virtual ~FunctionalSpirale();

  virtual double getFunctionValue(const VectorDouble& coor) const override;

  MatrixSquareGeneral getFunctionMatrix(const VectorDouble& coor) const;
  VectorVectorDouble getFunctionVectors(const Db *db, const CovAniso* cova) const;

private:
  static double _linearCombination(double x, double y, double a, double b);

private:
  double _a;
  double _b;
  double _c;
  double _d;
  double _xcenter;
  double _ycenter;
};
