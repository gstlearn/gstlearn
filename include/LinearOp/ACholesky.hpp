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
#include "Matrix/AMatrix.hpp"


namespace gstlrn
{ 
class AMatrix;
class MatrixDense;


class GSTLEARN_EXPORT ACholesky: public ASimulable
{
public:
  ACholesky(const AMatrix& mat);
  ACholesky(const ACholesky& m);
  ACholesky& operator=(const ACholesky& m);
  virtual ~ACholesky() {}

  Id getSize() const override { return _size; }
  Id solve(const constvect vecin, vect vecout) const;
  Id InvLtX(const constvect whitenoise, vect vecout) const;
  Id LtX(const constvect whitenoise, vect vecout) const;
  Id LX(const constvect whitenoise, vect vecout) const;
  Id InvLX(const constvect whitenoise, vect vecout) const;
  Id solveMatrix(const MatrixDense& b, MatrixDense& x) const;
  bool isReady() const { return _ready; }
  VectorDouble invLtX(const VectorDouble& vecin) const;
  VectorDouble LtX(const VectorDouble& vecin) const;
  VectorDouble LX(const VectorDouble& vecin) const;
  VectorDouble invLX(const VectorDouble& vecin) const; 
  VectorDouble solveX(const VectorDouble& vecin) const;
  
  virtual double computeLogDeterminant() const                    = 0;
  virtual Id addSolveX(const constvect vecin, vect vecout) const = 0;
  virtual Id addInvLtX(const constvect vecin, vect vecout) const = 0;
  virtual Id addLtX(const constvect vecin, vect vecout) const    = 0;
  virtual Id addLX(const constvect vecin, vect vecout) const     = 0;
  virtual Id addInvLX(const constvect vecin, vect vecout) const  = 0;

protected:
  void _setReady() const { _ready = true; _empty = false; }

private:
  Id _addToDest(const constvect vecin, vect vecout) const override;
  Id _addSimulateToDest(const constvect whitenoise, vect vecout) const override;
  double computeLogDet(Id nMC = 1) const override; //just for ASimulable interface

protected:
  Id _size;
  mutable bool _ready;
  mutable bool _empty;
};
}