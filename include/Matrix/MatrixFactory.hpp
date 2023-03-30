/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

class AMatrix;
class AMatrixSquare;

class GSTLEARN_EXPORT MatrixFactory {

public:
  MatrixFactory();
  virtual ~MatrixFactory();

  /// TODO : Use smartpointer
  static AMatrix* matProduct(const AMatrix* x, const AMatrix* y);
  static AMatrixSquare* matNorm(const AMatrixSquare* x, const AMatrix* y);
  static AMatrix* createIdentity(int nrow, bool sparse);
  static AMatrixSquare* createMatrixSquare(const AMatrixSquare* x,int nrow);
};
