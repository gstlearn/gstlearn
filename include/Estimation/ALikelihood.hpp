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

#include "Basic/VectorNumT.hpp"
#include "Estimation/AModelOptim.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Basic/NamingConvention.hpp"

namespace gstlrn
{
class Db;
class ModelGeneric;

class GSTLEARN_EXPORT ALikelihood: public AModelOptim
{
public:
  ALikelihood(ModelGeneric* model,
              const Db* db,
              bool reml = false);
  ALikelihood(const ALikelihood& r);
  ALikelihood& operator=(const ALikelihood& r);
  virtual ~ALikelihood();

  double computeCost(bool verbose = false) override;
  double computeLogLikelihood(bool verbose = false);
  VectorDouble getBeta() const { return _beta; }

protected:
  void _initLikelihood(bool verbose = false);

private:
  virtual void _updateModel(bool verbose = false)
  {
    DECLARE_UNUSED(verbose);
  }
  virtual void _computeCm1X()           = 0;
  virtual void _computeCm1Yc()          = 0;
  virtual double _computeLogDet() const = 0;
  virtual void _init(bool verbose = false)
  {
    DECLARE_UNUSED(verbose);
  }

protected:
  const Db* _db;
  VectorDouble _Y;  // Vector of multivariate data
  VectorDouble _Yc; // Centered data
  MatrixDense _X;   // Matrix of drifts
  VectorDouble _beta;
  MatrixDense _Cm1X;
  VectorDouble _Cm1Yc;
  MatrixSymmetric _XtCm1X; // X^T * C^{-1} * X
  bool _reml;
  Id _nDrift;
};
} // namespace gstlrn