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

MatrixFactory::MatrixFactory()
{
}

MatrixFactory::~MatrixFactory()
{
}

/****************************************************************************/
/*!
 **  Performs the product of two matrices: X * Y
 **
 ** \return Pointer to the newly created AMatrix matrix
 **
 ** \param[in]  x          First AMatrix matrix
 ** \param[in]  y          Second AMatrix matrix
 **
 *****************************************************************************/
AMatrix* MatrixFactory::matProduct(const AMatrix* x, const AMatrix* y)
{
  if (x->getNCols() != y->getNRows())
  {
    my_throw("Incompatible dimensions when making product of two matrices");
  }

  /// TODO : use typeinfo
  const MatrixSquareSymmetric*   mxsym = dynamic_cast<const MatrixSquareSymmetric*>(x);
  const MatrixSquareSymmetric*   mysym = dynamic_cast<const MatrixSquareSymmetric*>(y);

  AMatrix* res = nullptr;

  if (x->getNRows() == y->getNCols())
  {
    // Case of a resulting Square matrix

    if (mxsym != nullptr || mysym != nullptr)
    {
      res = new MatrixSquareSymmetric();
    }
    else
    {
      res = new MatrixSquareGeneral();
    }
  }
  else
  {
    res = new MatrixRectangular();
  }

  res->reset(x->getNRows(), y->getNCols());
  res->prodMatrix(*x, *y);

  return res;
}

/****************************************************************************/
/*!
 **  Performs the norm product of matrix V by matrix X: t(Y) * X * Y
 **
 ** \return Pointer to the newly created AMatrix matrix
 **
 ** \param[in]  x          AMatrixSquare matrix
 ** \param[in]  y          Second AMatrix matrix
 **
 *****************************************************************************/
AMatrixSquare* MatrixFactory::matNorm(const AMatrixSquare *x, const AMatrix *y)
{
  if (x->getNCols() != y->getNRows())
  {
    my_throw("Incompatible dimensions when making norm product of two matrices");
  }

  const MatrixSquareSymmetric*     mxsym = dynamic_cast<const MatrixSquareSymmetric*>(x);
  const MatrixSquareSymmetric*     mysym = dynamic_cast<const MatrixSquareSymmetric*>(y);

  AMatrixSquare* res = nullptr;

  if (mxsym != nullptr && mysym != nullptr)
  {
    res = new MatrixSquareSymmetric();
  }
  else
  {
    res = new MatrixSquareGeneral();
  }

  res->reset(y->getNCols(), y->getNCols());
  res->normMatrix(*x, *y);

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

