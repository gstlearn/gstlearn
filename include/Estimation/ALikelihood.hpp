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
#include "Matrix/MatrixDense.hpp"
#include "Estimation/AModelOptim.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Basic/NamingConvention.hpp"

class Db;
class ModelGeneric;

class GSTLEARN_EXPORT ALikelihood: public AModelOptim
{
public:
  ALikelihood(ModelGeneric* model,
              const Db* db);
  ALikelihood(const ALikelihood& r);
  ALikelihood& operator=(const ALikelihood& r);
  virtual ~ALikelihood();

  void init(bool verbose = false);
  double computeCost(bool verbose = false) override;
  double computeLogLikelihood(bool verbose = false);

private:
  virtual void _updateModel(bool verbose = false)
  {
    DECLARE_UNUSED(verbose);
  }
  virtual void _computeCm1X()           = 0;
  virtual void _computeCm1Z()           = 0;
  virtual double _computeLogDet() const = 0;
  virtual void _init(bool verbose = false)
  {
    DECLARE_UNUSED(verbose);
  }

protected:
  const Db* _db;
  VectorDouble _Y; // Vector of multivariate data
  MatrixDense _X;  // Matrix of drifts
  VectorDouble _beta;
  MatrixDense _Cm1X;
  VectorDouble _Cm1Y;
};
