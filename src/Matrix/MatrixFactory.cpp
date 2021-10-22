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
#include "Matrix/MatrixSDiag.hpp"
#include "Matrix/MatrixSDiagCst.hpp"
#include "Matrix/MatrixSGeneral.hpp"
#include "Matrix/MatrixSSym.hpp"

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
AMatrix* MatrixFactory::matProduct(const AMatrix* x,
                                     const AMatrix* y)
{
  if (x->getNCols() != y->getNRows())
  {
    my_throw("Incompatible dimensions when making product of two matrices");
  }

  /// TODO : use typeinfo
  const MatrixSDiag*    mxsd  = dynamic_cast<const MatrixSDiag*>(x);
  const MatrixSDiagCst* mxsdc = dynamic_cast<const MatrixSDiagCst*>(x);
  const MatrixSSym*     mxsym = dynamic_cast<const MatrixSSym*>(x);
  const MatrixSDiag*    mysd  = dynamic_cast<const MatrixSDiag*>(y);
  const MatrixSDiagCst* mysdc = dynamic_cast<const MatrixSDiagCst*>(y);
  const MatrixSSym*     mysym = dynamic_cast<const MatrixSSym*>(y);

  AMatrix* res = nullptr;

  if (x->getNRows() == y->getNCols())
  {
    // Case of a resulting Square matrix

    if (mxsdc != nullptr && mysdc != nullptr)
    {
      res = new MatrixSDiagCst();
    }
    else if ((mxsdc != nullptr || mxsd != nullptr)
          && (mysdc != nullptr || mysd != nullptr))
    {
      res = new MatrixSDiag();
    }
    else if ((mxsdc != nullptr || mxsd != nullptr || mxsym != nullptr)
          && (mysdc != nullptr || mysd != nullptr || mysym != nullptr))
    {
      res = new MatrixSSym();
    }
    else
    {
      res = new MatrixSGeneral();
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
AMatrixSquare* MatrixFactory::matNorm(const AMatrixSquare* x,
                                        const AMatrix* y)
{
  if (x->getNCols() != y->getNRows())
  {
    my_throw("Incompatible dimensions when making norm product of two matrices");
  }

  const MatrixSDiag*    mxsd  = dynamic_cast<const MatrixSDiag*>(x);
  const MatrixSDiagCst* mxsdc = dynamic_cast<const MatrixSDiagCst*>(x);
  const MatrixSSym*     mxsym = dynamic_cast<const MatrixSSym*>(x);
  const MatrixSDiag*    mysd  = dynamic_cast<const MatrixSDiag*>(y);
  const MatrixSDiagCst* mysdc = dynamic_cast<const MatrixSDiagCst*>(y);
  const MatrixSSym*     mysym = dynamic_cast<const MatrixSSym*>(y);

  AMatrixSquare* res = nullptr;

  if (mxsdc != nullptr && mysdc != nullptr)
  {
    res = new MatrixSDiagCst();
  }
  else if ((mxsdc != nullptr || mxsd != nullptr)
      && (mysdc != nullptr || mysd != nullptr))
  {
    res = new MatrixSDiag();
  }
  else if ((mxsdc != nullptr || mxsd != nullptr || mxsym != nullptr)
      && (mysdc != nullptr || mysd != nullptr || mysym != nullptr))
  {
    res = new MatrixSSym();
  }
  else
  {
    res = new MatrixSGeneral();
  }

  res->reset(y->getNCols(), y->getNCols());
  res->normMatrix(*x, *y);

  return res;

}

AMatrix* MatrixFactory::createIdentity(int nrow, bool sparse)
{
  return new MatrixSDiagCst(nrow, sparse);
}
