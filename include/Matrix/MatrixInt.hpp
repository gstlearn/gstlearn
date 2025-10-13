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

#include "Basic/VectorNumT.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/AStringable.hpp"

namespace gstlrn
{
/**
 * Rectangular matrices are stored by columns
 */
class GSTLEARN_EXPORT MatrixInt : public AStringable, public ICloneable {

public:
  MatrixInt(Id nrow = 0, Id ncol = 0);
  MatrixInt(const MatrixInt &r);
  MatrixInt& operator= (const MatrixInt &r);
	virtual ~MatrixInt();

  /// ICloneable interface
  IMPLEMENT_CLONING(MatrixInt)

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  static MatrixInt* createFromVI(const VectorInt &X,
                                 Id nrow,
                                 Id ncol,
                                 bool byCol = false);
  static MatrixInt* createFromVVI(const VectorVectorInt &X);

  void   reset(Id nrows, Id ncols);
  void   resetFromArray(Id nrows, Id ncols, const Id* tab, bool byCol = true);
  Id    getValue(Id irow, Id icol) const;
  Id    getValue(Id irank) const;

  void   setValueByRank(Id rank, Id value);
  void   setValue(Id irow, Id icol, Id value);
  Id    getMatrixSize() const;
  Id    size() const { return getMatrixSize(); }
  VectorInt getValues() const;
  VectorInt getValuesPerRow(Id irow) const;
  VectorInt getValuesPerColumn(Id icol) const;
  VectorVectorInt getMatrix() const;
  void   setValues(const VectorInt& values, bool byCol = true);
  void   setValuesOldStyle(const Id* values, bool byCol = true);
  void   transposeInPlace();
  bool   empty() const { return _nRows <= 0 || _nCols <= 0; }
  void   fill(Id value);

  Id  getNCols() const { return _nCols; }
  void setNCols(Id cols) { _nCols = cols; }
  Id  getNRows() const { return _nRows; }
  void setNRows(Id rows) { _nRows = rows; }

  /*! Get value operator override */
  Id  operator()(Id irow, Id icol) const { return getValue(irow, icol); }
  /*! Set value operator override */
  Id &operator()(Id irow, Id icol)       { return _getValueRef(irow, icol); }

private:
  void _allocate();
  void _deallocate();
  bool _isIndexValid(Id irow, Id icol) const;
  bool _isRankValid(Id rank) const;
  Id  _getIndexToRank(Id irow, Id icol) const;
  Id& _getValueRef(Id irow, Id icol);
  static bool _isNumbersValid(Id nrows, Id ncols);
  static void _transposeInPlace(Id n1, Id n2, const Id* v1, Id* w1);

private:
  Id _nRows;
  Id _nCols;
  VectorInt _rectMatrix;
};
}