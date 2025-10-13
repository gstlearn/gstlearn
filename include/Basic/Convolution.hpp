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

#include "Basic/VectorNumT.hpp"

namespace gstlrn
{
class DbGrid;
class MatrixDense;

class GSTLEARN_EXPORT Convolution
{
public:
  Convolution(DbGrid* dbgrid = nullptr);
  Convolution(const Convolution& m);
  Convolution& operator=(const Convolution& m);
  virtual ~Convolution();

  Id ConvolveSparse(Id iatt,
                     const VectorVectorInt& ranks,
                     const MatrixDense& wgt,
                     const VectorDouble& means = VectorDouble(),
                     Id optionVerbose         = 0);
  Id ConvolveFFT(Id iatt,
                  Id nvar,
                  const DbGrid* marpat,
                  const VectorDouble& means = VectorDouble());

private:
  bool _isDbGridDefined() const;

private:
  DbGrid* _dbgrid; // Pointer to external DbGrid: do not delete
};
}