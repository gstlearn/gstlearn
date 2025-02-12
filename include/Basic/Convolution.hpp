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

class DbGrid;
class MatrixRectangular;

class GSTLEARN_EXPORT Convolution
{
public:
  Convolution(DbGrid* dbgrid = nullptr);
  Convolution(const Convolution& m);
  Convolution& operator=(const Convolution& m);
  virtual ~Convolution();

  int ConvolveSparse(int iatt,
                     const VectorVectorInt& ranks,
                     const MatrixRectangular& wgt,
                     const VectorDouble& means = VectorDouble());
  int ConvolveFFT(int iatt,
                  int nvar,
                  const DbGrid* marpat,
                  const VectorDouble& means = VectorDouble());

private:
  bool _isDbGridDefined() const;

private:
  DbGrid* _dbgrid; // Pointer to external DbGrid: do not delete
};
