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

/**
 * Create a submatrix from an input matrix
 * by specifying the list of rows and columns to be extracted or excluded
 * @param x             Input matrix
 * @param selRows       List of selected rows
 * @param selCols       List of selected columns
 * @param flagKeepRows  True if the selected rows must be kept (they must be excluded otherwise)
 * @param flagKeepCols  True if the selected columns must be kept (they must be excluded otherwise)
 * @return
 *
 * @remark If a list if not defined, the whole set of rows (resp. columns) is considered)
 * @remark (this assumes that flagKeep is True)
 */
AMatrix* MatrixFactory::createReduce(const AMatrix *x,
                                     const VectorInt &selRows,
                                     const VectorInt &selCols,
                                     bool flagKeepRows,
                                     bool flagKeepCols)
{
  int nrows = x->getNRows();
  VectorInt localSelRows;
  if (selRows.empty())
    localSelRows = VH::sequence(nrows);
  else
  {
    localSelRows = VH::filter(selRows, 0, x->getNRows());
    if (!flagKeepRows) localSelRows = VH::complement(VH::sequence(nrows), localSelRows);
  }
  int newNRows = (int) localSelRows.size();
  if (newNRows <= 0)
  {
    messerr("The new Matrix has no Row left");
    return nullptr;
  }

  int ncols = x->getNCols();
  VectorInt localSelCols;
  if (selCols.empty())
    localSelCols = VH::sequence(ncols);
  else
  {
    localSelCols = VH::filter(selCols, 0, x->getNCols());
    if (!flagKeepCols) localSelCols = VH::complement(VH::sequence(ncols), localSelCols);
  }
  int newNCols = (int) localSelCols.size();
  if (newNCols <= 0)
  {
    messerr("The new Matrix has no Column left");
    return nullptr;
  }
  bool flagSame = (localSelRows == localSelCols);

  /// TODO : use typeinfo
  AMatrix* res = nullptr;
  const MatrixRectangular*        mxrg  = dynamic_cast<const MatrixRectangular*>(x);
  const MatrixSquareGeneral*      mxsg  = dynamic_cast<const MatrixSquareGeneral*>(x);
  const MatrixSquareSymmetric*    mxsym = dynamic_cast<const MatrixSquareSymmetric*>(x);

  if (mxsym != nullptr)
  {
    // Case of a square symmetric input matrix
    if (flagSame)
      res = new MatrixSquareSymmetric(newNRows);
    else
      res = new MatrixRectangular(newNRows, newNCols);
  }
  else if (mxsg != nullptr)
  {
    // Case of a square general input matrix
    if (flagSame)
      res = new MatrixSquareGeneral(newNRows);
    else
      res = new MatrixRectangular(newNRows, newNCols);
  }
  else if (mxrg != nullptr)
  {
    // Case of a rectangular input matrix
    res = new MatrixRectangular(newNRows, newNCols);
  }
  else
    messageAbort("CreateReduce cannot be called for such matrix. This should never happen");

  res->copyReduce(x, localSelRows, localSelCols);

  return res;
}

/**
 * Create a submatrix from an input matrix by specifying the row and or column to be extracted or excluded
 * @param x             Input matrix
 * @param selRow        Rank of the selected row
 * @param selCol        Rank of the selected column
 * @param flagKeepRow   True if the selected row must be kept (it must be excluded otherwise)
 * @param flagKeepCol   True if the selected column must be kept (it must be excluded otherwise)
 * @return
 *
 * @remark If a rank is negative, the whole set of rows (resp. columns) is considered)
 * @remark (this assumes that flagKeep is True)
 */
AMatrix* MatrixFactory::createReduceOne(const AMatrix *x,
                                        int selRow,
                                        int selCol,
                                        bool flagKeepRow,
                                        bool flagKeepCol)
{
  VectorInt localSelRows;
  if (selRow >= 0)
  {
    localSelRows.resize(1);
    localSelRows[0] = selRow;
   }

  VectorInt localSelCols;
  if (selCol >= 0)
  {
    localSelCols.resize(1);
    localSelCols[0] = selCol;
  }
  return MatrixFactory::createReduce(x, localSelRows, localSelCols, flagKeepRow, flagKeepCol);
}

/*****************************************************************************/
/*!
 **  Concatenate two matrices
 **
 ** \return Pointer on the newly created concatenated matrix (or NULL)
 **
 ** \param[in]  a1   Pointer to the first matrix
 ** \param[in]  a2   Pointer to the second matrix
 ** \param[in]  flagShiftRow Concatenate by Row
 ** \param[in]  flagShiftCol Concatenate by Column
 **
 *****************************************************************************/
AMatrix* MatrixFactory::createGlue(const AMatrix* a1,
                                   const AMatrix* a2,
                                   bool flagShiftRow,
                                   bool flagShiftCol)
{
  AMatrix *a = nullptr;

  // Preliminary checks

  bool isSparse = a1->isSparse();
  if ((a1->isSparse() && ! a2->isSparse()) || (! a1->isSparse() && a2->isSparse()))
  {
    messerr("In 'createGlue()' both matrices should be sparse or not sparse");
    return a;
  }
  int n11 = a1->getNRows();
  int n12 = a1->getNCols();
  int n21 = a2->getNRows();
  int n22 = a2->getNCols();
  int nrows = n11;
  int ncols = n12;
  if (flagShiftRow && ! isSparse)
  {
    if (n12 != n22)
    {
      messerr("Binding by row: Input matrices must share same column number");
      return a;
    }
    nrows += n21;
  }
  if (flagShiftCol && ! isSparse)
  {
    if (n11 != n21)
    {
      messerr("Binding by column: Input matrices must share same row number");
      return a;
    }
    ncols += n22;
  }

  /* Core allocation */


  if (isSparse)
  {
    const MatrixSparse* aloc1 = dynamic_cast<const MatrixSparse*>(a1);
    const MatrixSparse* aloc2 = dynamic_cast<const MatrixSparse*>(a2);
    a = MatrixSparse::glue(aloc1, aloc2, flagShiftRow, flagShiftCol);
  }
  else
    a = MatrixRectangular::glue(a1, a2, flagShiftRow, flagShiftCol);

  return a;
}
