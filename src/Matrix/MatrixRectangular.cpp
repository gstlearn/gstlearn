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
#include "Matrix/MatrixRectangular.hpp"
#include "Basic/AStringable.hpp"
#include "Matrix/AMatrix.hpp"
#include "Basic/VectorHelper.hpp"

MatrixRectangular::MatrixRectangular(int nrows, int ncols)
  : AMatrixDense(nrows, ncols)
{
}

MatrixRectangular::MatrixRectangular(const MatrixRectangular &r)
  : AMatrixDense(r)
{
}

MatrixRectangular::MatrixRectangular(const AMatrix &m)
  : AMatrixDense(m)
{
}

MatrixRectangular& MatrixRectangular::operator= (const MatrixRectangular &r)
{
  if (this != &r)
  {
    AMatrixDense::operator=(r);
  }
  return *this;
}

MatrixRectangular::~MatrixRectangular()
{
}

MatrixRectangular* MatrixRectangular::create(const MatrixRectangular* mat)
{
  return new MatrixRectangular(*mat);
}
MatrixRectangular* MatrixRectangular::create(int nrow, int ncol)
{
  return new MatrixRectangular(nrow, ncol);
}

/**
 * Converts a VectorVectorDouble into a Matrix
 * Note: the input argument is stored by row (if coming from [] specification)
 * @param  X Input VectorVectorDouble argument
 * @return The returned rectangular matrix
 *
 * @remark: the matrix is transposed implicitly while reading
 */
MatrixRectangular* MatrixRectangular::createFromVVD(const VectorVectorDouble& X)
{
  int nrow = (int) X.size();
  int ncol = (int) X[0].size();

  MatrixRectangular* mat = new MatrixRectangular(nrow, ncol);
  mat->_fillFromVVD(X);
  return mat;
}

MatrixRectangular* MatrixRectangular::createFromVD(const VectorDouble &X,
                                                   int nrow,
                                                   int ncol,
                                                   bool byCol,
                                                   bool invertColumnOrder)
{
  if (nrow * ncol != (int) X.size())
  {
    messerr("Inconsistency between arguments 'nrow'(%d) and 'ncol'(%d)", nrow, ncol);
    messerr("and the dimension of the input Vector (%d)", (int) X.size());
  }
  MatrixRectangular* mat = new MatrixRectangular(nrow, ncol);

  int lec = 0;
  if (byCol)
  {
    for (int irow = 0; irow < nrow; irow++)
      for (int icol = 0; icol < ncol; icol++)
      {
        int jcol = (invertColumnOrder) ? ncol - icol - 1 : icol;
        mat->setValue(irow, jcol, X[lec++]);
      }
  }
  else
  {
    for (int icol = 0; icol < ncol; icol++)
      for (int irow = 0; irow < nrow; irow++)
      {
        int jcol = (invertColumnOrder) ? ncol - icol - 1 : icol;
        mat->setValue(irow, jcol, X[lec++]);
      }
  }
  return mat;
}

MatrixRectangular* MatrixRectangular::glue(const AMatrix *A1,
                                           const AMatrix *A2,
                                           bool flagShiftRow,
                                           bool flagShiftCol)
{
  // Create the new matrix
  int shiftRow = (flagShiftRow) ? A1->getNRows() : 0;
  int shiftCol = (flagShiftCol) ? A1->getNCols() : 0;

  int nrows = (flagShiftRow) ? A1->getNRows() + A2->getNRows() : MAX(A1->getNRows(), A2->getNRows());
  int ncols = (flagShiftCol) ? A1->getNCols() + A2->getNCols() : MAX(A1->getNCols(), A2->getNCols());

  MatrixRectangular* mat = new MatrixRectangular(nrows, ncols);
  mat->fill(0.);

  // Copy the first input matrix

  for (int irow = 0; irow < A1->getNRows(); irow++)
    for (int icol = 0; icol < A1->getNCols(); icol++)
      mat->setValue(irow, icol, A1->getValue(irow, icol));

  // Copy the second input matrix

  for (int irow = 0; irow < A2->getNRows(); irow++)
    for (int icol = 0; icol < A2->getNCols(); icol++)
      mat->setValue(irow + shiftRow, icol + shiftCol, A2->getValue(irow, icol));

  return mat;
}

void MatrixRectangular::addRow(int nrow_added)
{
  int nrows = getNRows();
  int ncols = getNCols();

  AMatrix* statsSave = this->clone();
  reset(nrows+nrow_added, ncols);
  for (int irow=0; irow< nrows; irow++)
    for (int icol=0; icol<ncols; icol++)
      setValue(irow, icol, statsSave->getValue(irow, icol));
  delete statsSave;
}

void MatrixRectangular::addColumn(int ncolumn_added)
{
  int nrows = getNRows();
  int ncols = getNCols();

  AMatrix* statsSave = this->clone();
  reset(nrows, ncols+ncolumn_added);
  for (int irow=0; irow< nrows; irow++)
    for (int icol=0; icol<ncols; icol++)
      setValue(irow, icol, statsSave->getValue(irow, icol));
  delete statsSave;
}

/**
 * @brief Create an output Rectangular Matrix by selecting some rows and columns
 *        of the Input matrix 'A'
 *
 * @param A        Input Rectangular Matrix
 * @param rowKeep  Set of Rows to be kept (all if not defined)
 * @param colKeep  Set of Columns to be kept (all if not defined)
 * @param flagInvertRow when True, transform 'rowKeep' into 'rowDrop'
 * @param flagInvertCol when True, transform 'colKeep' into 'colDrop'
 * @return Newly created Rectangular Matrix
 */
MatrixRectangular* MatrixRectangular::sample(const AMatrix* A,
                                             const VectorInt& rowKeep,
                                             const VectorInt& colKeep,
                                             bool flagInvertRow,
                                             bool flagInvertCol)
{
  int nrowtotal = A->getNRows();
  int ncoltotal = A->getNCols();
  VectorInt rows = rowKeep;
  if (rows.empty()) rows = VH::sequence(nrowtotal);
  if (flagInvertRow) rows = VH::complement(VH::sequence(nrowtotal), rows);
  VectorInt cols = colKeep;
  if (cols.empty()) cols = VH::sequence(ncoltotal);
  if (flagInvertCol) cols = VH::complement(VH::sequence(ncoltotal), cols);

  int nrows = (int)rows.size();
  int ncols = (int)cols.size();
  if (nrows <= 0 || ncols <= 0) return nullptr;

  for (int irow = 0; irow < nrows; irow++)
  {
    if (!checkArg("Selected Row index", rows[irow], nrowtotal)) return nullptr;
  }
  for (int icol = 0; icol < ncols; icol++)
  {
    if (!checkArg("Selected Column index", cols[icol], ncoltotal)) return nullptr;
  }

  MatrixRectangular* mat = new MatrixRectangular(nrows, ncols);
  for (int irow = 0; irow < nrows; irow++)
    for (int icol = 0; icol < ncols; icol++)
      mat->setValue(irow, icol, A->getValue(rows[irow], cols[icol]));
  return mat;
}

/**
 * @brief Set the values contained in 'A' into the current matrix
 *
 * @param A Input Matrix
 * @param rowFetch Set of row indices of 'this' where values of 'A' should be stored
 * @param colFetch Set of column indices of'this' where values of 'A' should be stored
 * @param flagInvertRow when True, transform 'rowFetch' into 'rowAvoid'
 * @param flagInvertCol when True, transform 'colFetch' into 'colAvoid'
 */
void MatrixRectangular::unsample(const AMatrix* A,
                                 const VectorInt& rowFetch,
                                 const VectorInt& colFetch,
                                 bool flagInvertRow,
                                 bool flagInvertCol)
{
  int nrowtotal = getNRows();
  int ncoltotal = getNCols();
  VectorInt rows = rowFetch;
  if (rows.empty()) rows = VH::sequence(nrowtotal);
  if (flagInvertRow) rows = VH::complement(VH::sequence(nrowtotal), rows);
  VectorInt cols = colFetch;
  if (cols.empty()) cols = VH::sequence(ncoltotal);
  if (flagInvertCol) cols = VH::complement(VH::sequence(ncoltotal), cols);

  int nrows = (int)rows.size();
  int ncols = (int)cols.size();
  if (nrows <= 0 || ncols <= 0) return;

  for (int irow = 0; irow < nrows; irow++)
  {
    if (!checkArg("Selected Row index", rows[irow], getNRows()))
      return;
  }
  for (int icol = 0; icol < ncols; icol++)
  {
    if (!checkArg("Selected Column index", cols[icol], getNCols()))
      return;
  }

  for (int irow = 0; irow < nrows; irow++)
    for (int icol = 0; icol < ncols; icol++)
      setValue(rows[irow], cols[icol], A->getValue(irow, icol));
}
