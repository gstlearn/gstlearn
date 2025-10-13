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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Basic/WarningMacro.hpp"
#include "Matrix/AMatrix.hpp"

#ifndef SWIG
DISABLE_WARNING_PUSH
DISABLE_WARNING_COND_EXPR_CONSTANT
DISABLE_WARNING_UNUSED_BUT_SET_VARIABLE
DISABLE_WARNING_MAYBE_UNINITIALIZED
DISABLE_WARNING_DECLARATION_HIDE_GLOBAL
#  include <Eigen/Dense>
#  include <Eigen/Eigenvalues>
DISABLE_WARNING_POP
#endif

namespace gstlrn
{
class MatrixSquare;
class MatrixSymmetric;
class EOperator;

/**
 * Dense Matrix
 * This class provides all the functions that can be performed using a Matrix stored
 * in "Dense" format (in opposition to the "Sparse" format).
 * This class can be derived in the case the matrix is Square, and even more if it is
 * Square and Symmetric.
 */

class GSTLEARN_EXPORT MatrixDense: public AMatrix
{

public:
  MatrixDense(Id nrow = 0, Id ncol = 0);
  MatrixDense(const MatrixDense& r);
  MatrixDense(const AMatrix& r);
  MatrixDense& operator=(const MatrixDense& r);
  virtual ~MatrixDense();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// Cloneable interface
  IMPLEMENT_CLONING(MatrixDense)

  /// Interface for AMatrix

  /*! Returns if the current matrix is Dense */
  bool isDense() const override { return true; }
  /*! Returns if the current matrix is Sparse */
  bool isSparse() const override { return false; }
  /*! Get the value from a matrix cell */
  double getValue(Id irow, Id icol) const override;
  /*! Set the value for in a matrix cell */
  void setValue(Id irow, Id icol, double value) override;
  /*! Update the contents of a matrix cell */
  void updValue(Id irow,
                Id icol,
                const EOperator& oper,
                double value) override;
  /*! Set the contents of a Column */
  void setColumn(Id icol, const VectorDouble& tab) override;
  /*! Set the contents of a Column to a constant value */
  void setColumnToConstant(Id icol, double value) override;
  /*! Set the contents of a Row */
  void setRow(Id irow, const VectorDouble& tab) override;
  /*! Set the contents of a Row to a constant value*/
  void setRowToConstant(Id irow, double value) override;
  /*! Set the contents of the (main) Diagonal */
  void setDiagonal(const VectorDouble& tab) override;
  /*! Set the contents of the (main) Diagonal to a constant value */
  void setDiagonalToConstant(double value = 1.) override;
  /*! Add a value to each matrix component */
  void addScalar(double v) override;
  /*! Add value to matrix diagonal */
  void addScalarDiag(double v) override;
  /*! Multiply each matrix component by a value */
  void prodScalar(double v) override;
  /*! Set all the values of the Matrix at once */
  void fill(double value) override;
  /*! Multiply a Matrix row-wise */
  void multiplyRow(const VectorDouble& vec) override;
  /*! Multiply a Matrix column-wise */
  void multiplyColumn(const VectorDouble& vec) override;
  /*! Divide a Matrix row-wise */
  void divideRow(const VectorDouble& vec) override;
  /*! Divide a Matrix column-wise */
  void divideColumn(const VectorDouble& vec) override;

  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const override { return false; }

  /*! Extract a Row */
  VectorDouble getRow(Id irow) const override;
  /*! Extract a Column */
  VectorDouble getColumn(Id icol) const override;
  /*! Multiply matrix 'x' by matrix 'y' and store the result in 'this' */
  void prodMatMatInPlace(const AMatrix* x,
                         const AMatrix* y,
                         bool transposeX = false,
                         bool transposeY = false) override;
  /*! Product 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' stored in 'this'*/
  void prodNormMatMatInPlace(const AMatrix* a,
                             const AMatrix* m,
                             bool transpose = false) override;
  /*! Perform 'this' = 't(A)' %*% 'vec' %*% 'A' or 'A' %*% 'vec' %*% 't(A)' */
  void prodNormMatVecInPlace(const AMatrix* a,
                             const VectorDouble& vec,
                             bool transpose = false) override;
  /*! Perform 'this' = 't(A)' %*% 'A' or 'A' %*% 't(A)' */
  void prodNormMatInPlace(const AMatrix* a,
                          bool transpose = false) override;
  /*! Perform 'this' = 'val1' * 'mat1' + 'val2' * 'mat2' + 'val3' * 'mat3' */
  void linearCombination(double val1,
                         const AMatrix* mat1,
                         double val2         = 1.,
                         const AMatrix* mat2 = nullptr,
                         double val3         = 1.,
                         const AMatrix* mat3 = nullptr) override;
  /*! Add a matrix (multiplied by a constant) */
  void addMat(const AMatrix& y, double cx = 1., double cy = 1.) override;

  const VectorDouble& getEigenValues() const { return _eigenValues; }
  const MatrixSquare* getEigenVectors() const { return _eigenVectors; }
  Id invert2(MatrixDense& res) const;
  void unsample(const AMatrix* A,
                const VectorInt& rowFetch,
                const VectorInt& colFetch,
                bool flagInvertRow = false,
                bool flagInvertCol = false);
  MatrixDense compressMatLC(const MatrixDense& matLC, bool transpose = false);

  // Adding a Row or a Column (at the bottom or right of Rectangular Matrix)
  void addRow(Id nrow_added = 1);
  void addColumn(Id ncolumn_added = 1);
  constvect getColumnPtr(Id icol) const;

  static MatrixDense* create(const MatrixDense* mat);
  static MatrixDense* create(Id nrow, Id ncol);
  static double traceProd(const MatrixDense& a, MatrixDense& b); // Warning: b is modified to reduce memory allocation
  static MatrixDense* createFromVVD(const VectorVectorDouble& X);
  static MatrixDense* createFromVD(const VectorDouble& X,
                                   Id nrow,
                                   Id ncol,
                                   bool byCol             = false,
                                   bool invertColumnOrder = false);
  static MatrixDense* glue(const AMatrix* A1,
                           const AMatrix* A2,
                           bool flagShiftRow,
                           bool flagShiftCol);
  static bool sample(MatrixDense& res,
                     const AMatrix& A,
                     const VectorInt& rowKeep = VectorInt(),
                     const VectorInt& colKeep = VectorInt(),
                     bool flagInvertRow       = false,
                     bool flagInvertCol       = false);
#ifndef SWIG
  static void sum(const MatrixDense* mat1,
                  const MatrixDense* mat2,
                  MatrixDense* mat3);
#endif // !SWIG

protected:
  void _allocate() override;
  void _deallocate() override;
  Id _getIndexToRank(Id irow, Id icol) const override;
  void _setValueByRank(Id rank, double value) override;
  void _transposeInPlace() override;
  Id _invert() override;
  Id _solve(const VectorDouble& b, VectorDouble& x) const override;
  Id _getMatrixPhysicalSize() const override;
  double _getValueByRank(Id rank) const override;
  double& _getValueRef(Id irow, Id icol) override;

#ifndef SWIG
  void _addProdVecMatInPlacePtr(constvect x, vect y, bool transpose = false) const override;
  void _addProdMatVecInPlacePtr(constvect x, vect y, bool transpose = false) const override;
#endif

  Id _computeEigen(bool optionPositive = true);
  Id _computeGeneralizedEigen(const MatrixSymmetric& b, bool optionPositive = true);

private:
  void _recopy(const MatrixDense& r);
  Id _terminateEigen(const Eigen::VectorXd& eigenValues,
                     const Eigen::MatrixXd& eigenVectors,
                     bool optionPositive = true,
                     bool changeOrder    = false);
  bool _needToReset(Id nrows, Id ncols) override;

#ifndef SWIG

public:
  constvect getViewOnColumn(Id icol) const;
  vect getViewOnColumnModify(Id icol);
#endif

#ifndef SWIG

public:
  Eigen::Map<const Eigen::MatrixXd> eigenMat() const
  {
    return Eigen::Map<const Eigen::MatrixXd>(_eigenMatrix.data(), getNRows(), getNCols());
  }
  Eigen::Map<Eigen::MatrixXd> eigenMat()
  {
    return Eigen::Map<Eigen::MatrixXd>(_eigenMatrix.data(), getNRows(), getNCols());
  }
#endif

protected:
  VectorDouble _eigenValues;
  MatrixSquare* _eigenVectors;
  VectorDouble _eigenMatrix;

private:
  Id _maxSize;
};
} // namespace gstlrn
