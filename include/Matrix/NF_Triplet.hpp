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
#include "geoslib_define.h"

#include "Basic/ASerializable.hpp"
#include "Basic/AStringable.hpp"

#include <Eigen/Sparse>

class cs;

/**
 * Stores the contents of a sparse matrix in Triplet form
 * The format is adapter to Eigen or to csparse
 */

typedef Eigen::Triplet<double> T;

class GSTLEARN_EXPORT NF_Triplet
{
public:
  NF_Triplet();
  NF_Triplet(const NF_Triplet &m);
  NF_Triplet& operator= (const NF_Triplet &m);
  virtual ~NF_Triplet();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  void add(int irow, int icol, double value);
  int getNumber() const { return (int) _eigenT.size(); }
  int getNRows()  const { return _nrowmax; }
  int getNCols()  const { return _ncolmax; }
  void force(int nrow, int ncol);

  int getRow(int i) const;
  int getCol(int i) const;
  double getValue(int i) const;
  VectorDouble getValues() const;
  VectorInt getRows(bool flag_from_1 = false) const;
  VectorInt getCols(bool flag_from_1 = false) const;
  void appendInPlace(const NF_Triplet& T2);

  cs* buildCsFromTriplet() const;
  Eigen::SparseMatrix<double> buildEigenFromTriplet() const;

  static NF_Triplet createFromCs(const cs* mat, int shiftRow=0, int shiftCol=0);
  static NF_Triplet createFromEigen(const Eigen::SparseMatrix<double>& mat, int shiftRow=0, int shiftCol=0);

private:
  int  _nrowmax;
  int  _ncolmax;
  std::vector<T> _eigenT; // Triplet in Eigen format
};
