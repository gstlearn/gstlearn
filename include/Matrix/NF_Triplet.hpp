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

#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/WarningMacro.hpp"
#include "gstlearn_export.hpp"

#ifndef SWIG
DISABLE_WARNING_PUSH
DISABLE_WARNING_COND_EXPR_CONSTANT
DISABLE_WARNING_UNUSED_BUT_SET_VARIABLE
DISABLE_WARNING_DECLARATION_HIDE_GLOBAL
#  include <Eigen/Sparse>
DISABLE_WARNING_POP

typedef Eigen::Triplet<double, gstlrn::Id> T;
#endif

namespace gstlrn
{
/**
 * Stores the contents of a sparse matrix in Triplet form
 * The format is adapter to Eigen
 */
class GSTLEARN_EXPORT NF_Triplet
{
public:
  NF_Triplet();
  NF_Triplet(const NF_Triplet& r);
  NF_Triplet& operator=(const NF_Triplet&);
  virtual ~NF_Triplet();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  void add(Id irow, Id icol, double value);
  Id getNElements() const { return static_cast<Id>(_eigenT.size()); }
  Id getNRows() const { return _nrowmax; }
  Id getNCols() const { return _ncolmax; }
  void force(Id nrow, Id ncol);

  Id getRow(Id i) const;
  Id getCol(Id i) const;
  double getValue(Id i) const;
  VectorDouble getValues() const;
  VectorInt getRows(bool flag_from_1 = false) const;
  VectorInt getCols(bool flag_from_1 = false) const;
  void appendInPlace(const NF_Triplet& T2);

#ifndef SWIG
  ::Eigen::SparseMatrix<double> buildEigenFromTriplet() const;

  static NF_Triplet createFromEigen(const Eigen::SparseMatrix<double>& mat, Id shiftRow = 0, Id shiftCol = 0);
#endif

private:
  Id _nrowmax;
  Id _ncolmax;
#ifndef SWIG
  std::vector<T> _eigenT; // Triplet in Eigen format
#endif
};
} // namespace gstlrn
