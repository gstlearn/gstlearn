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
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"

class AMatrix;
class AMatrixSquare;

class GSTLEARN_EXPORT MatrixFactory {

public:
  /// TODO : Use smartpointer
  static AMatrix* prodMatMat(const AMatrix *x,
                             const AMatrix *y,
                             bool transposeX = false,
                             bool transposeY = false);
  template <typename T>
  static T* prodMatMat(const AMatrix *x,
                       const AMatrix *y,
                       bool transposeX = false,
                       bool transposeY = false);
  static AMatrixSquare* createMatrixSquare(const AMatrixSquare* x,int nrow);
  static AMatrix* createReduce(const AMatrix *x,
                               const VectorInt &selRows = VectorInt(),
                               const VectorInt &selCols = VectorInt(),
                               bool flagKeepRows = true,
                               bool flagKeepCols = true);
  static AMatrix* createReduceOne(const AMatrix *x,
                                  int selRow = -1,
                                  int selCol = -1,
                                  bool flagKeepRow = true,
                                  bool flagKeepCol = true);
  static AMatrix* createGlue(const AMatrix *a1,
                             const AMatrix *a2,
                             bool flagShiftRow,
                             bool flagShiftCol);
};

/****************************************************************************/
/*!
 **  Performs the product of two matrices: X * Y
 **
 ** \return Pointer to the newly created AMatrix matrix
 **
 ** \param[in]  x          First AMatrix matrix
 ** \param[in]  y          Second AMatrix matrix
 ** \param[in]  transposeX True if First matrix is transposed
 ** \param[in]  transposeY True if Second matrix is transposed
 **
 ** \remarks: To be called as follows:
 **      MatrixSparse* mat = MatrixFactory::prodMatMat<MatrixSparse>(x, y);
 **
 *****************************************************************************/
template<typename T>
T* MatrixFactory::prodMatMat(const AMatrix *x,
                             const AMatrix *y,
                             bool transposeX,
                             bool transposeY)
{
  T* res = new T();

  if (x->getNCols() != y->getNRows())
  {
    messerr("Incompatible dimensions when making product of two matrices");
    return res;
  }

  res->AMatrix::reset(x->getNRows(), y->getNCols(), 0., x->isFlagEigen());
  res->prodMatMatInPlace(x, y, transposeX, transposeY);

  return res;
}

