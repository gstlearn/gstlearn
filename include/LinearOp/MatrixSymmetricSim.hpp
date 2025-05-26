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

#include "LinearOp/ASimulable.hpp"
#include "LinearOp/ACholesky.hpp"
#include "Matrix/AMatrix.hpp"

class AMatrix;
class Cholesky;

/**
 * Square Symmetric matrices
 */
class GSTLEARN_EXPORT MatrixSymmetricSim : public ASimulable
{
public:
  MatrixSymmetricSim(const AMatrix* m, bool inverse = true);
  MatrixSymmetricSim(const MatrixSymmetricSim &m) = delete;
  MatrixSymmetricSim& operator=(const MatrixSymmetricSim &m) = delete;
  virtual ~MatrixSymmetricSim();

  const AMatrix* getMatrix() const;
  int  getSize() const override;
  bool isEmpty() const { return _factor == nullptr; }
  double computeLogDet(int nMC = 1) const override;
  
  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

protected:
  virtual int _addSimulateToDest(const constvect whitenoise, vect outv) const override;
  virtual int _addToDest(const constvect inv, vect outv) const override;

private:
  bool _inverse;
  ACholesky* _factor;
};
