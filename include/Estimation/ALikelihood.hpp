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

#include "Basic/NamingConvention.hpp"
#include "Basic/VectorNumT.hpp"
#include "Estimation/AModelOptim.hpp"
#include "LinearOp/ASimulableMatrix.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class Db;
  class ModelGeneric;

  class GSTLEARN_EXPORT ALikelihood: public AModelOptim, public ASimulableMatrix
  {
  public:
    ALikelihood(ModelGeneric* model, const Db* db, bool reml = false);

    double computeCost(bool flagPrint = false, bool verbose = false) override;
    double computeLogLikelihood(bool flagPrint = false, bool verbose = false);

    VectorDouble getBeta() const { return _beta; }

    void initLikelihood(bool verbose = false);
    void updateModel(bool verbose = false);
    double computeLogDet(Id nMC = 1) const override;
    Id getSize() const override;

  protected:
    void _initLikelihoodForOptim(bool verbose = false);
    Id _addSimulateToDest(const constvect whitenoise, vect outv) const override;
    Id _addToDest(constvect inv, vect outv) const override;
    bool _calculateBeta(bool verbose = false);

  private:
    virtual void _solveQ(constvect inv, vect outv) const = 0;

    virtual void _updateModel(bool verbose = false) { DECLARE_UNUSED(verbose); }

    virtual void _computeCm1X() = 0;
    virtual void _computeCm1Yc() = 0;
    virtual double _computeLogDet() const = 0;

    virtual void _init(bool verbose = false) { DECLARE_UNUSED(verbose); }

  protected:
    const Db* _db;
    VectorDouble _Z; // Vector of multivariate data (raw)
    VectorDouble _Y; // Vector of multivariate data (Gaussian)
    VectorDouble _Yc; // Centered data
    MatrixDense _X; // Matrix of drifts
    VectorDouble _beta;
    MatrixDense _Cm1X;
    VectorDouble _Cm1Yc;
    MatrixSymmetric _XtCm1X; // X^T * C^{-1} * X
    bool _reml;
    Id _nDrift;
  };
} // namespace gstlrn
