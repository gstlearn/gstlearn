/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Matrix/AMatrix.hpp"
#include "Matrix/AMatrixSquare.hpp"
#include "Matrix/MatrixFactory.hpp"

#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareDiagonal.hpp"
#include "Matrix/MatrixSquareDiagonalCst.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"

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
  const MatrixSquareDiagonal*    mxsd  = dynamic_cast<const MatrixSquareDiagonal*>(x);
  const MatrixSquareDiagonalCst* mxsdc = dynamic_cast<const MatrixSquareDiagonalCst*>(x);
  const MatrixSquareSymmetric*   mxsym = dynamic_cast<const MatrixSquareSymmetric*>(x);
  const MatrixSquareDiagonal*    mysd  = dynamic_cast<const MatrixSquareDiagonal*>(y);
  const MatrixSquareDiagonalCst* mysdc = dynamic_cast<const MatrixSquareDiagonalCst*>(y);
  const MatrixSquareSymmetric*   mysym = dynamic_cast<const MatrixSquareSymmetric*>(y);

  AMatrix* res = nullptr;

  if (x->getNRows() == y->getNCols())
  {
    // Case of a resulting Square matrix

    if (mxsdc != nullptr && mysdc != nullptr)
    {
      res = new MatrixSquareDiagonalCst();
    }
    else if ((mxsdc != nullptr || mxsd != nullptr)
          && (mysdc != nullptr || mysd != nullptr))
    {
      res = new MatrixSquareDiagonal();
    }
    else if ((mxsdc != nullptr || mxsd != nullptr || mxsym != nullptr)
          && (mysdc != nullptr || mysd != nullptr || mysym != nullptr))
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

  const MatrixSquareDiagonal*    mxsd  = dynamic_cast<const MatrixSquareDiagonal*>(x);
  const MatrixSquareDiagonalCst* mxsdc = dynamic_cast<const MatrixSquareDiagonalCst*>(x);
  const MatrixSquareSymmetric*     mxsym = dynamic_cast<const MatrixSquareSymmetric*>(x);
  const MatrixSquareDiagonal*    mysd  = dynamic_cast<const MatrixSquareDiagonal*>(y);
  const MatrixSquareDiagonalCst* mysdc = dynamic_cast<const MatrixSquareDiagonalCst*>(y);
  const MatrixSquareSymmetric*     mysym = dynamic_cast<const MatrixSquareSymmetric*>(y);

  AMatrixSquare* res = nullptr;

  if (mxsdc != nullptr && mysdc != nullptr)
  {
    res = new MatrixSquareDiagonalCst();
  }
  else if ((mxsdc != nullptr || mxsd != nullptr)
      && (mysdc != nullptr || mysd != nullptr))
  {
    res = new MatrixSquareDiagonal();
  }
  else if ((mxsdc != nullptr || mxsd != nullptr || mxsym != nullptr)
      && (mysdc != nullptr || mysd != nullptr || mysym != nullptr))
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

AMatrix* MatrixFactory::createIdentity(int nrow, bool sparse)
{
  return new MatrixSquareDiagonalCst(nrow, sparse);
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
  const MatrixSquareDiagonal*    mxsd  = dynamic_cast<const MatrixSquareDiagonal*>(x);
  const MatrixSquareDiagonalCst* mxsdc = dynamic_cast<const MatrixSquareDiagonalCst*>(x);
  const MatrixSquareSymmetric*   mxsym = dynamic_cast<const MatrixSquareSymmetric*>(x);

  AMatrixSquare* res = nullptr;
  if (mxsg != nullptr)
  {
    res = new MatrixSquareGeneral(nrow);
  }
  else if (mxsd != nullptr)
  {
    res = new MatrixSquareDiagonal(nrow);
  }
  else if (mxsdc != nullptr)
  {
    res = new MatrixSquareDiagonalCst(nrow);
  }
  else if (mxsym != nullptr)
  {
    res = new MatrixSquareSymmetric(nrow);
  }
  return res;
}

