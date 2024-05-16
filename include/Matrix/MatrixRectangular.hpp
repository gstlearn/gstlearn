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

#include "Basic/VectorNumT.hpp"
#include "Matrix/AMatrixDense.hpp"

/**
 * Rectangular matrices are stored by columns
 */
class GSTLEARN_EXPORT MatrixRectangular : public AMatrixDense {

public:
  MatrixRectangular(int nrow = 0, int ncol = 0, int opt_eigen=-1);
  MatrixRectangular(const MatrixRectangular &m);
  MatrixRectangular(const AMatrix &m);
  MatrixRectangular& operator= (const MatrixRectangular &r);
	virtual ~MatrixRectangular();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

	/// Cloneable interface
  IMPLEMENT_CLONING(MatrixRectangular)

  /// Interface for AMatrix
  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const override { return false; }
  /*! Say if the matrix must be diagonal */
  bool mustBeDiagonal() const override { return false; }
  /*! Say if the matrix must be diagonal constant */
  bool mustBeDiagCst() const override { return false; }

  static MatrixRectangular* createFromVVD(const VectorVectorDouble& X, int opt_eigen = -1);
  static MatrixRectangular* createFromVD(const VectorDouble &X,
                                         int nrow,
                                         int ncol,
                                         bool byCol = false,
                                         int opt_eigen = -1,
                                         bool invertColumnOrder = false);
  static MatrixRectangular* glue(const AMatrix *A1,
                                 const AMatrix *A2,
                                 bool flagShiftRow,
                                 bool flagShiftCol);

  /*! Adding a Row or a Column (at the bottom or right of Rectangular Matrix) */
  void addRow(int nrow_added=1);
  void addColumn(int ncolumn_added = 1);

protected:
  virtual int _getIndexToRank(int irow,int icol) const override;

private:
  /// Interface for AMatrix
  virtual int     _getMatrixPhysicalSize() const override;
  virtual double& _getValueRef(int irow, int icol) override;
  virtual void    _allocate() override;
  virtual double  _getValue(int irow, int icol) const override;
  virtual double  _getValueByRank(int irank) const override;
  virtual void    _setValueByRank(int rank, double value) override;
  virtual void    _setValue(int irow, int icol, double value) override;
  virtual void    _updValue(int irow, int icol, const EOperator& oper, double value) override;

  virtual void    _transposeInPlace() override;
  virtual void    _prodMatVecInPlacePtr(const double *x,double *y, bool transpose = false) const override;
  virtual void    _prodVecMatInPlacePtr(const double *x,double *y, bool transpose = false) const override;
  virtual int     _invert() override;
  virtual int     _solve(const VectorDouble& b, VectorDouble& x) const override;

private:
  /// ========================================================================
  /// The subsequent methods rely on the specific local storage ('rectMatrix')
  /// ========================================================================
  void    _allocateLocal();
  void    _recopyLocal(const MatrixRectangular& r);
  double  _getValueLocal(int irow, int icol) const;
  double  _getValueLocal(int irank) const;
  double& _getValueRefLocal(int irow, int icol);
  void    _setValueLocal(int irank, double value);
  void    _setValueLocal(int irow, int icol, double value);
  void    _updValueLocal(int irow, int icol, const EOperator& oper, double value);
  void    _prodMatVecInPlacePtrLocal(const double *x, double *y, bool transpose = false) const;
  void    _prodVecMatInPlacePtrLocal(const double *x, double *y, bool transpose = false) const;
  void    _transposeInPlaceLocal();
  int     _getIndexToRankLocal(int irow, int icol) const;

private:
  VectorDouble _rectMatrix; // Classical storage
};
