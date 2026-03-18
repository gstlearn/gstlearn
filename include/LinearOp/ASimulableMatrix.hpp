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
#include "LinearOp/CholeskySparse.hpp"
#include "gstlearn_export.hpp"
#include "LinearOp/ASimulable.hpp"

namespace gstlrn{

class MatrixSparse;
class GSTLEARN_EXPORT ASimulableMatrix : public virtual ASimulable
{
public:
  ASimulableMatrix();
  virtual ~ASimulableMatrix();
  double computeLogDet(Id nMC = 1) const override; 
  virtual const MatrixSparse& getQMat() const = 0;

private:
  const CholeskySparse& getChol() const;
#ifndef SWIG

protected:
  Id _addToDest(const constvect whitenoise, vect outv) const override;
  Id _addSimulateToDest(const constvect whitenoise, vect outv) const override;
#endif

private:
  mutable CholeskySparse* _chol; // when needed, e.g in a const method.
};
}