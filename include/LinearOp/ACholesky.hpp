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

#include "LinearOp/ASimulable.hpp"
#include "gstlearn_export.hpp"
#include "Matrix/AMatrix.hpp"

class GSTLEARN_EXPORT ACholesky: public ASimulable
{
public:
  ACholesky(const AMatrix* mat, bool inverse = true);
  ACholesky(const ACholesky& m)            = delete;
  ACholesky& operator=(const ACholesky& m) = delete;
  virtual ~ACholesky() {}

  int getSize() const override { return _size; }
  int solve(const constvect vecin, vect vecout) const;

  virtual double computeLogDeterminant() const = 0;

private:
  int _addToDest(const constvect vecin, vect vecout) const override;
  int _addSimulateToDest(const constvect whitenoise, vect vecout) const override;

  virtual int _addSolveX(const constvect vecin, vect vecout) const = 0;
  virtual int _addInvLtX(const constvect vecin, vect vecout) const  = 0;
  virtual int _addLX(const constvect vecin, vect vecout) const  = 0;

protected:
  const AMatrix* _mat; // Pointer to original matrix (not to be deleted)
  bool _inverse;
  int _size;
};
