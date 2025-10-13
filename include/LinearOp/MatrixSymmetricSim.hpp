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

#include "Matrix/MatrixSymmetric.hpp"
#include "gstlearn_export.hpp"

#include "LinearOp/ASimulable.hpp"
#include "LinearOp/ACholesky.hpp"
#include "Matrix/AMatrix.hpp"

namespace gstlrn
{
class AMatrix;

/**
 * Square Symmetric matrices
 */
class GSTLEARN_EXPORT MatrixSymmetricSim : public ASimulable
{
public:
  MatrixSymmetricSim(const AMatrix& m, bool inverse = true);
  MatrixSymmetricSim(const MatrixSymmetricSim &m) = delete;
  MatrixSymmetricSim& operator=(const MatrixSymmetricSim &m) = delete;
  virtual ~MatrixSymmetricSim();

  Id  getSize() const override;
  bool isEmpty() const { return _factor == nullptr; }
  double computeLogDet(Id nMC = 1) const override;
  
  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

protected:
  Id _addSimulateToDest(const constvect whitenoise, vect outv) const override;
  Id _addToDest(const constvect inv, vect outv) const override;

private:
  bool _inverse;
  ACholesky* _factor;
  const AMatrix& _mat;
  MatrixSymmetric _matSymConverted; // Used if the input matrix is not a MatrixSymmetric, just to avoid premature destruction
};
}