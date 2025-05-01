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
#include "Matrix/MatrixSymmetric.hpp"

class AMatrix;
class MatrixSquare;

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
  static MatrixSquare* createMatrixSquare(const MatrixSquare* x,int nrow);
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
 ** TODO : Why 2 methods for MatrixFactory::prodMatMat ?
 *****************************************************************************/
template<typename T>
T* MatrixFactory::prodMatMat(const AMatrix *x,
                             const AMatrix *y,
                             bool transposeX,
                             bool transposeY)
{
  T* res = new T(); /// TODO : if MatrixSparse => x or y 'eigen flag' is ignored

  int nxrows = (! transposeX) ? x->getNRows() : x->getNCols();
  int nxcols = (! transposeX) ? x->getNCols() : x->getNRows();
  int nyrows = (! transposeY) ? y->getNRows() : y->getNCols();
  int nycols = (! transposeY) ? y->getNCols() : y->getNRows();

  if (nxcols != nyrows)
  {
    messerr("Incompatible dimensions when making product of two matrices");
    return res;
  }

  res->reset(nxrows, nycols);
  res->prodMatMatInPlace(x, y, transposeX, transposeY);

  return res;
}

