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
#include "Matrix/AMatrix.hpp"
#include "Matrix/AMatrixSquare.hpp"
#include "Matrix/MatrixFactory.hpp"

#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixSparse.hpp"

#include "Basic/VectorHelper.hpp"
#include "Basic/AException.hpp"

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
 *****************************************************************************/
AMatrix* MatrixFactory::prodMatMat(const AMatrix *x,
                                   const AMatrix *y,
                                   bool transposeX,
                                   bool transposeY)
{
  if (x->getNCols() != y->getNRows())
  {
    my_throw("Incompatible dimensions when making product of two matrices");
  }

  const MatrixSparse* mxsparse = dynamic_cast<const MatrixSparse*>(x);
  const MatrixSparse* mysparse = dynamic_cast<const MatrixSparse*>(y);

  AMatrix* res = nullptr;
  if (mxsparse != nullptr && mysparse != nullptr)
  {
    // Case of a resulting Sparse matrix

    res = new MatrixSparse();
  }
  else
  {

    // Case of a resulting Dense matrix

    const MatrixSquareSymmetric* mxsym = dynamic_cast<const MatrixSquareSymmetric*>(x);
    const MatrixSquareSymmetric* mysym = dynamic_cast<const MatrixSquareSymmetric*>(y);

    if (x->getNRows() == y->getNCols())
    {

      // Case of a resulting Square matrix

      if (mxsym != nullptr || mysym != nullptr)
      {

        // Cas of a resulting Square Symmetric matrix

        res = new MatrixSquareSymmetric();
      }
      else
      {

        // Case of a resulting Square general matrix

        res = new MatrixSquareGeneral();
      }
    }
    else
    {

      // Case of a resulting Rectangular matrix

      res = new MatrixRectangular();
    }
  }

  res->reset(x->getNRows(), y->getNCols(), 0., x->isFlagEigen());
  res->prodMatMatInPlace(x, y, transposeX, transposeY);

  return res;
}

/****************************************************************************/
/*!
 **  Create a Matrix similar to the input one with a given row number
 **
 ** \return Pointer to the newly created AMatrix matrix
 **
 ** \param[in]  x          First AMatrix matrix
 ** \param[in]  nrow       Number of rows
 **
 *****************************************************************************/
AMatrixSquare* MatrixFactory::createMatrixSquare(const AMatrixSquare *x,
                                                 int nrow)
{
  /// TODO : use typeinfo
  const MatrixSquareGeneral*     mxsg  = dynamic_cast<const MatrixSquareGeneral*>(x);
  const MatrixSquareSymmetric*   mxsym = dynamic_cast<const MatrixSquareSymmetric*>(x);

  AMatrixSquare* res = nullptr;
  if (mxsg != nullptr)
  {
    res = new MatrixSquareGeneral(nrow);
  }
  else if (mxsym != nullptr)
  {
    res = new MatrixSquareSymmetric(nrow);
  }
  return res;
}

AMatrix* MatrixFactory::createReduce(const AMatrix *x,
                                     const VectorInt &validRows,
                                     const VectorInt &validCols)
{
  // Order and shrink the input vectors
  VectorInt localValidRows = VH::filter(validRows, 0, x->getNRows());
  VectorInt localValidCols = VH::filter(validCols, 0, x->getNCols());
  int newNRows = (int) localValidRows.size();
  int newNCols = (int) localValidCols.size();
  if (newNRows <= 0)
  {
    messerr("The new Matrix has no Row left");
    return nullptr;
  }
  if (newNCols <= 0)
  {
    messerr("The new Matrix has no Column left");
    return nullptr;
  }
  bool flagSame = (localValidRows == localValidCols);

  /// TODO : use typeinfo
  AMatrix* res = nullptr;
  const MatrixRectangular*        mxrg  = dynamic_cast<const MatrixRectangular*>(x);
  const MatrixSquareGeneral*      mxsg  = dynamic_cast<const MatrixSquareGeneral*>(x);
  const MatrixSquareSymmetric*    mxsym = dynamic_cast<const MatrixSquareSymmetric*>(x);

  if (mxrg != nullptr)
  {
    res = new MatrixRectangular(newNRows, newNCols);
  }
  else if (mxsg != nullptr)
  {
    if (flagSame)
      res = new MatrixSquareGeneral(newNRows);
    else
      res = new MatrixRectangular(newNRows, newNCols);
  }
  else if (mxsym != nullptr)
  {
    if (flagSame)
      res = new MatrixSquareSymmetric(newNRows);
    else
      res = new MatrixRectangular(newNRows, newNCols);
  }
  res->copyReduce(x, localValidRows, localValidCols);

  return res;
}

