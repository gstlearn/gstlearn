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
#include "Matrix/NF_Triplet.hpp"
#include "Matrix/LinkMatrixSparse.hpp"

#include "csparse_d.h"

NF_Triplet::NF_Triplet()
    : _nrowmax(0),
      _ncolmax(0),
      _eigenT()
{
}

NF_Triplet::NF_Triplet(const NF_Triplet &r)
  : _nrowmax(r._nrowmax),
    _ncolmax(r._ncolmax),
    _eigenT(r._eigenT)
{
}

NF_Triplet& NF_Triplet::operator= (const NF_Triplet &r)
{
  if (this != &r)
  {
    _nrowmax = r._nrowmax;
    _ncolmax = r._ncolmax;
    _eigenT = r._eigenT;
  }
  return *this;
}

NF_Triplet::~NF_Triplet()
{
}

void NF_Triplet::add(int irow, int icol, double value)
{
  if (irow > _nrowmax) _nrowmax = irow;
  if (icol > _ncolmax) _ncolmax = icol;
  _eigenT.push_back(T(irow, icol, value));
}

/**
 * Force the dimension of the Sparse matrix
 * This is done by adding a fictitious sample at position 'nrow-1' and 'ncol-1' with value 0
 * @param nrow    Ultimate number of rows
 * @param ncol    Ultimate number of columns
 */
void NF_Triplet::force(int nrow, int ncol)
{
  // Check if the maximum positions have been reached
  int nrow_max = getNRows();
  int ncol_max = getNCols();
  if (nrow_max < nrow || ncol_max < ncol)
    add(nrow - 1, ncol - 1, 0.);
}

cs* NF_Triplet::buildCsFromTriplet() const
{
  cs* local = cs_spalloc2(0,0,1,1,1);
  for (int i = 0, n = getNumber(); i < n; i++)
    (void) cs_entry2(local, getRow(i), getCol(i), getValue(i));
  cs* Q = cs_triplet2(local);
  local = cs_spfree2(local);
  return Q;
}

Eigen::SparseMatrix<double> NF_Triplet::buildEigenFromTriplet() const
{
  Eigen::SparseMatrix<double> mat(_nrowmax+1, _ncolmax+1);
  mat.setFromTriplets(_eigenT.begin(), _eigenT.end());
  return mat;
}

NF_Triplet NF_Triplet::createFromCs(const cs* mat, int shiftRow, int shiftCol)
{
  return csToTriplet(mat, shiftRow, shiftCol);
}

NF_Triplet NF_Triplet::createFromEigen(const Eigen::SparseMatrix<double>& mat, int shiftRow, int shiftCol)
{
  NF_Triplet NF_T;
  std::vector<T> v;
  int row_max = 0;
  int col_max = 0;
  for(int i = 0; i < mat.outerSize(); i++)
    for(typename Eigen::SparseMatrix<double>::InnerIterator it(mat,i); it; ++it)
    {
      int irow = it.row() + shiftRow;
      int icol = it.col() + shiftCol;
      if (irow > row_max) row_max = irow;
      if (icol > col_max) col_max = icol;
      v.emplace_back(irow,icol,it.value());
    }
  NF_T._nrowmax = row_max;
  NF_T._ncolmax = col_max;
  NF_T._eigenT = v;
  return NF_T;
}

int NF_Triplet::getRow(int i) const
{
  if (i < 0 || i >= getNumber()) return ITEST;
  return _eigenT[i].row();
}

int NF_Triplet::getCol(int i) const
{
  if (i < 0 || i >= getNumber()) return ITEST;
  return _eigenT[i].col();
}

double NF_Triplet::getValue(int i) const
{
  if (i < 0 || i >= getNumber()) return TEST;
  return _eigenT[i].value();
}

VectorDouble NF_Triplet::getValues() const
{
  int n = getNumber();
  VectorDouble vec(n);
  for (int i = 0; i< n; i++)
    vec[i] = _eigenT[i].value();
  return vec;
}

VectorInt NF_Triplet::getRows(bool flag_from_1) const
{
  int n = getNumber();
  int shift = (flag_from_1) ? 1 : 0;
  VectorInt vec(n);
  for (int i = 0; i< n; i++)
    vec[i] = _eigenT[i].row() + shift;
  return vec;
}

VectorInt NF_Triplet::getCols(bool flag_from_1) const
{
  int n = getNumber();
  int shift = (flag_from_1) ? 1 : 0;
  VectorInt vec(n);
  for (int i = 0; i< n; i++)
    vec[i] = _eigenT[i].col() + shift;
  return vec;
}

/**
 * Append NF_Triplet 'T2' at the end of the current one
 * @param T2 NF_Triplet to be appended
 */
void NF_Triplet::appendInPlace(const NF_Triplet& T2)
{
  _eigenT.insert( _eigenT.end(), T2._eigenT.begin(), T2._eigenT.end() );

  // Update the maxima for rows and columns
  for (int i = 0, n = T2.getNumber(); i < n; i++)
  {
    int irow = T2.getRow(i);
    int icol = T2.getCol(i);
    if (irow > _nrowmax) _nrowmax = irow;
    if (icol > _ncolmax) _ncolmax = icol;
  }
}
