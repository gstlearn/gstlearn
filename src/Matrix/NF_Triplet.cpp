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

namespace gstlrn
{
NF_Triplet::NF_Triplet()
  : _nrowmax(0)
  , _ncolmax(0)
  , _eigenT()
{
}

NF_Triplet::NF_Triplet(const NF_Triplet& r)
  : _nrowmax(r._nrowmax)
  , _ncolmax(r._ncolmax)
  , _eigenT(r._eigenT)
{
}

NF_Triplet& NF_Triplet::operator=(const NF_Triplet& r)
{
  if (this != &r)
  {
    _nrowmax = r._nrowmax;
    _ncolmax = r._ncolmax;
    _eigenT  = r._eigenT;
  }
  return *this;
}

NF_Triplet::~NF_Triplet()
{
}

void NF_Triplet::add(Id irow, Id icol, double value)
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
void NF_Triplet::force(Id nrow, Id ncol)
{
  // Check if the maximum positions have been reached
  Id nrow_max = getNRows();
  Id ncol_max = getNCols();
  if (nrow_max < nrow || ncol_max < ncol)
    add(nrow - 1, ncol - 1, 0.);
}

Eigen::SparseMatrix<double> NF_Triplet::buildEigenFromTriplet() const
{
  Eigen::SparseMatrix<double> mat(_nrowmax + 1, _ncolmax + 1);
  mat.setFromTriplets(_eigenT.begin(), _eigenT.end());
  //  mat.prune(EPSILON10);
  return mat;
}

NF_Triplet NF_Triplet::createFromEigen(const Eigen::SparseMatrix<double>& mat, Id shiftRow, Id shiftCol)
{
  NF_Triplet NF_T;
  std::vector<T> v;
  Id row_max = 0;
  Id col_max = 0;
  for (Id i = 0; i < mat.outerSize(); i++)
    for (typename Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it)
    {
      Id irow = it.row() + shiftRow;
      Id icol = it.col() + shiftCol;
      if (irow > row_max) row_max = irow;
      if (icol > col_max) col_max = icol;
      v.emplace_back(irow, icol, it.value());
    }
  NF_T._nrowmax = row_max;
  NF_T._ncolmax = col_max;
  NF_T._eigenT  = v;
  return NF_T;
}

Id NF_Triplet::getRow(Id i) const
{
  if (i < 0 || i >= getNElements()) return ITEST;
  return _eigenT[i].row();
}

Id NF_Triplet::getCol(Id i) const
{
  if (i < 0 || i >= getNElements()) return ITEST;
  return _eigenT[i].col();
}

double NF_Triplet::getValue(Id i) const
{
  if (i < 0 || i >= getNElements()) return TEST;
  return _eigenT[i].value();
}

VectorDouble NF_Triplet::getValues() const
{
  auto n = getNElements();
  VectorDouble vec(n);
  for (Id i = 0; i < n; i++)
    vec[i] = _eigenT[i].value();
  return vec;
}

VectorInt NF_Triplet::getRows(bool flag_from_1) const
{
  auto n    = getNElements();
  Id shift = (flag_from_1) ? 1 : 0;
  VectorInt vec(n);
  for (Id i = 0; i < n; i++)
    vec[i] = _eigenT[i].row() + shift;
  return vec;
}

VectorInt NF_Triplet::getCols(bool flag_from_1) const
{
  auto n    = getNElements();
  Id shift = (flag_from_1) ? 1 : 0;
  VectorInt vec(n);
  for (Id i = 0; i < n; i++)
    vec[i] = _eigenT[i].col() + shift;
  return vec;
}

/**
 * Append NF_Triplet 'T2' at the end of the current one
 * @param T2 NF_Triplet to be appended
 */
void NF_Triplet::appendInPlace(const NF_Triplet& T2)
{
  _eigenT.insert(_eigenT.end(), T2._eigenT.begin(), T2._eigenT.end());

  // Update the maxima for rows and columns
  for (Id i = 0, n = T2.getNElements(); i < n; i++)
  {
    auto irow = T2.getRow(i);
    auto icol = T2.getCol(i);
    if (irow > _nrowmax) _nrowmax = irow;
    if (icol > _ncolmax) _ncolmax = icol;
  }
}
} // namespace gstlrn