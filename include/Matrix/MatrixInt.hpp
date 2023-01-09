/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 9 avr. 2019 by N. Desassis                                     */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "Basic/VectorNumT.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/AStringable.hpp"
#include "Matrix/AMatrix.hpp"

/**
 * Rectangular matrices are stored by columns
 */
class GSTLEARN_EXPORT MatrixInt : public AStringable, public ICloneable {

public:
  MatrixInt(int nrow = 0, int ncol = 0);
  MatrixInt(const MatrixInt &m);
  MatrixInt& operator= (const MatrixInt &r);
	virtual ~MatrixInt();

  /// ICloneable interface
  IMPLEMENT_CLONING(MatrixInt)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void   reset(int nrows, int ncols);
  int    getValue(int irow, int icol) const;
  int    getValue(int irank) const;
  void   setValue(int rank, int value);
  void   setValue(int irow, int icol, int value);
  int    getMatrixSize() const;
  int    size() const { return getMatrixSize(); }
  VectorInt getValues() const;
  VectorInt getValuesPerRow(int irow) const;
  VectorInt getValuesPerColumn(int icol) const;
  VectorVectorInt getMatrix() const;
  void   setValues(const VectorInt& values, bool byCol = true);
  void   setValuesOldStyle(const int* values, bool byCol = true);
  void   transposeInPlace();
  bool   empty() const { return _nRows <= 0 || _nCols <= 0; }
  void   fill(int value);

  int getNCols() const { return _nCols; }
  void setNCols(int cols) { _nCols = cols; }
  int getNRows() const { return _nRows; }
  void setNRows(int rows) { _nRows = rows; }

  static MatrixInt* createFromVVD(const VectorVectorInt& X);

private:
  int    _getIndexToRank(int irow,int icol) const;
  void   _allocate();
  void   _deallocate();
  bool   _isIndexValid(int irow, int icol) const;
  bool   _isRankValid(int rank) const;
  bool   _isNumbersValid(int nrows, int ncols) const;

private:
  int _nRows;
  int _nCols;
  VectorInt _rectMatrix;
};
