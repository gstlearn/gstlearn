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

class MatrixRectangular;

class GSTLEARN_EXPORT ACholesky: public ASimulable
{
public:
  ACholesky(const AMatrix* mat);
  ACholesky(const ACholesky& m);
  ACholesky& operator=(const ACholesky& m);
  virtual ~ACholesky() {}

  int getSize() const override { return _size; }
  const AMatrix* getMatrix() const { return _mat; }

  int solve(const constvect vecin, vect vecout) const;
  int InvLtX(const constvect whitenoise, vect vecout) const;
  int LtX(const constvect whitenoise, vect vecout) const;
  int LX(const constvect whitenoise, vect vecout) const;
  int InvLX(const constvect whitenoise, vect vecout) const;
  int solveMatrix(const MatrixRectangular& b, MatrixRectangular& x) const;
  bool isReady() const { return _ready; }


  VectorDouble invLtX(const VectorDouble& vecin) const;
  VectorDouble LtX(const VectorDouble& vecin) const;
  VectorDouble LX(const VectorDouble& vecin) const;
  VectorDouble invLX(const VectorDouble& vecin) const;
  VectorDouble solveX(const VectorDouble& vecin) const;
  
  virtual double computeLogDeterminant() const                    = 0;
  virtual int addSolveX(const constvect vecin, vect vecout) const = 0;
  virtual int addInvLtX(const constvect vecin, vect vecout) const = 0;
  virtual int addLtX(const constvect vecin, vect vecout) const    = 0;
  virtual int addLX(const constvect vecin, vect vecout) const     = 0;
  virtual int addInvLX(const constvect vecin, vect vecout) const  = 0;

protected:
  void _setReady() const { _ready = true; }

private:
  int _addToDest(const constvect vecin, vect vecout) const override;
  int _addSimulateToDest(const constvect whitenoise, vect vecout) const override;

protected:
  const AMatrix* _mat; // Pointer to original matrix (not to be deleted)
  int _size;
  mutable bool _ready;
};
