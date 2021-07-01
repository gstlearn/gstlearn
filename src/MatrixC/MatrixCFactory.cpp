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
#include "MatrixC/AMatrixC.hpp"
#include "MatrixC/AMatrixCSquare.hpp"
#include "MatrixC/MatrixCFactory.hpp"

#include "MatrixC/MatrixCRectangular.hpp"
#include "MatrixC/MatrixCSDiag.hpp"
#include "MatrixC/MatrixCSDiagCst.hpp"
#include "MatrixC/MatrixCSGeneral.hpp"
#include "MatrixC/MatrixCSSym.hpp"

#include "Basic/AException.hpp"

MatrixCFactory::MatrixCFactory()
{
}

MatrixCFactory::~MatrixCFactory()
{
}

/****************************************************************************/
/*!
 **  Performs the product of two matrices: X * Y
 **
 ** \return Pointer to the newly created AMatrixC matrix
 **
 ** \param[in]  x          First AMatrixC matrix
 ** \param[in]  y          Second AMatrixC matrix
 **
 *****************************************************************************/
AMatrixC* MatrixCFactory::matProduct(const AMatrixC* x,
                                     const AMatrixC* y)
{
  if (x->getNCols() != y->getNRows())
  {
    my_throw("Incompatible dimensions when making product of two matrices");
  }

  /// TODO : use typeinfo
  const MatrixCSDiag*    mxsd  = dynamic_cast<const MatrixCSDiag*>(x);
  const MatrixCSDiagCst* mxsdc = dynamic_cast<const MatrixCSDiagCst*>(x);
  const MatrixCSSym*     mxsym = dynamic_cast<const MatrixCSSym*>(x);
  const MatrixCSDiag*    mysd  = dynamic_cast<const MatrixCSDiag*>(y);
  const MatrixCSDiagCst* mysdc = dynamic_cast<const MatrixCSDiagCst*>(y);
  const MatrixCSSym*     mysym = dynamic_cast<const MatrixCSSym*>(y);

  AMatrixC* res = nullptr;

  if (x->getNRows() == y->getNCols())
  {
    // Case of a resulting Square matrix

    if (mxsdc != nullptr && mysdc != nullptr)
    {
      res = new MatrixCSDiagCst();
    }
    else if ((mxsdc != nullptr || mxsd != nullptr)
          && (mysdc != nullptr || mysd != nullptr))
    {
      res = new MatrixCSDiag();
    }
    else if ((mxsdc != nullptr || mxsd != nullptr || mxsym != nullptr)
          && (mysdc != nullptr || mysd != nullptr || mysym != nullptr))
    {
      res = new MatrixCSSym();
    }
    else
    {
      res = new MatrixCSGeneral();
    }
  }
  else
  {
    res = new MatrixCRectangular();
  }

  res->reset(x->getNRows(), y->getNCols());
  res->prodMatrix(*x, *y);

  return res;
}

/****************************************************************************/
/*!
 **  Performs the norm product of matrix V by matrix X: t(Y) * X * Y
 **
 ** \return Pointer to the newly created AMatrixC matrix
 **
 ** \param[in]  x          AMatrixCSquare matrix
 ** \param[in]  y          Second AMatrixC matrix
 **
 *****************************************************************************/
AMatrixCSquare* MatrixCFactory::matNorm(const AMatrixCSquare* x,
                                        const AMatrixC* y)
{
  if (x->getNCols() != y->getNRows())
  {
    my_throw("Incompatible dimensions when making norm product of two matrices");
  }

  const MatrixCSDiag*    mxsd  = dynamic_cast<const MatrixCSDiag*>(x);
  const MatrixCSDiagCst* mxsdc = dynamic_cast<const MatrixCSDiagCst*>(x);
  const MatrixCSSym*     mxsym = dynamic_cast<const MatrixCSSym*>(x);
  const MatrixCSDiag*    mysd  = dynamic_cast<const MatrixCSDiag*>(y);
  const MatrixCSDiagCst* mysdc = dynamic_cast<const MatrixCSDiagCst*>(y);
  const MatrixCSSym*     mysym = dynamic_cast<const MatrixCSSym*>(y);

  AMatrixCSquare* res = nullptr;

  if (mxsdc != nullptr && mysdc != nullptr)
  {
    res = new MatrixCSDiagCst();
  }
  else if ((mxsdc != nullptr || mxsd != nullptr)
      && (mysdc != nullptr || mysd != nullptr))
  {
    res = new MatrixCSDiag();
  }
  else if ((mxsdc != nullptr || mxsd != nullptr || mxsym != nullptr)
      && (mysdc != nullptr || mysd != nullptr || mysym != nullptr))
  {
    res = new MatrixCSSym();
  }
  else
  {
    res = new MatrixCSGeneral();
  }

  res->reset(y->getNCols(), y->getNCols());
  res->normMatrix(*x, *y);

  return res;

}

AMatrixC* MatrixCFactory::createIdentity(int nrow, bool sparse)
{
  return new MatrixCSDiagCst(nrow, sparse);
}
