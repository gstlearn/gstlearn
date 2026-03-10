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

#include "Basic/Message.hpp"
#include "gstlearn_export.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"

namespace gstlrn
{
class NF_Triplet;
class EOperator;

/// TODO : Transform into template for storing something else than double

/**
 * This class is the root of the Matrix organization in gstlearn
 * A matrix is a 2-D organization: it is characterized by its number of rows
 * and its number of columns.
 * Although the user should not bother with this remark, the elements of a matrix
 * processed in 'gstlearn' are stored in a Row-major format.
 * This is to say that the internal rank of an element characterized by its row and column numbers is:
 *  (icol * getNRows() + irow)
 *
 * Since gstlearn version v1.3:
 * - Dense Matrices storage and algebra rely on Eigen3 library
 * - Sparse Matrices storage and algebra rely on Eigen3 library
 */
class GSTLEARN_EXPORT AMatrix: public AStringable, public ICloneable
{
public:
  AMatrix(Id nrow = 0, Id ncol = 0);
  AMatrix(const AMatrix& m);
  AMatrix& operator=(const AMatrix& m);
  virtual ~AMatrix();

  /// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /*! Returns if the current matrix is Sparse */
  virtual bool isSparse() const = 0;
  /*! Returns if the matrix belongs to the MatrixDense class */
  virtual bool isDense() const = 0;
  /*! Check if the matrix is (non empty) square */
  virtual bool isSquare(bool printWhyNot = false) const;
  /*! Check if the input matrix is (non empty and square) symmetric */
  virtual bool isSymmetric(double eps = EPSILON10, bool printWhyNot = false) const;
  /*! Say if the matrix must be symmetric */
  virtual bool mustBeSymmetric() const { return false; }

  /*! Sets the value at row 'irow' and column 'icol' */
  virtual void setValue(Id irow, Id icol, double value) = 0;
  /*! Update the value at row 'irow' and column 'icol' */
  virtual void updValue(Id irow, Id icol, const EOperator& oper, double value) = 0;
  /*! Set the contents of a Column */
  virtual void setColumn(Id icol, const VectorDouble& tab) = 0;
  /*! Set the contents of a Column to a constant value*/
  virtual void setColumnToConstant(Id icol, double value) = 0;
  /*! Set the contents of a Row */
  virtual void setRow(Id irow, const VectorDouble& tab) = 0;
  /*! Set the contents of a Row to a constant value*/
  virtual void setRowToConstant(Id irow, double value) = 0;
  /*! Set the contents of the (main) Diagonal */
  virtual void setDiagonal(const VectorDouble& tab) = 0;
  /*! Set the contents of the (main) Diagonal to a constant value */
  virtual void setDiagonalToConstant(double value = 1.) = 0;

  /*! Multiply a Matrix row-wise */
  virtual void multiplyRow(const VectorDouble& vec) = 0;
  /*! Multiply a Matrix column-wise */
  virtual void multiplyColumn(const VectorDouble& vec) = 0;
  /*! Divide a Matrix row-wise */
  virtual void divideRow(const VectorDouble& vec) = 0;
  /*! Divide a Matrix column-wise */
  virtual void divideColumn(const VectorDouble& vec) = 0;

  /*! Gets the value at row 'irow' and column 'icol' */
  virtual double getValue(Id irow, Id icol) const = 0;
  /*! Extract a Row */
  virtual VectorDouble getRow(Id irow) const;
  /*! Extract a Column */
  virtual VectorDouble getColumn(Id icol) const;

  /*! Set all the values of the Matrix at once */
  virtual void fill(double value) = 0;

  /**
   * @brief Operate addition a constant to all terms of a matrix and store the result in res = other1 + cst
   *
   * @tparam T Type of the matrix
   * @param res Output Matrix
   * @param other Input matrix
   * @param cst Constant value to be added
   */
  template<typename T>
  static void add(T& res, const T& other, double cst)
  {
    static_assert(std::is_base_of_v<AMatrix, T>,
                  "Invalid type for Matrix addition (template method). "
                  "Only MatrixDense and MatrixSparse are allowed");
    if (res.getNRows() != other.getNRows() || res.getNCols() != other.getNCols())
    {
      res.resize(other.getNRows(), other.getNCols());
    }
    res.addGeneral(res, other, cst);
  }

  /**
   * @brief Operate addition of two matrices and store the result in a third one (res = other1 + other2)
   *
   * @tparam T Type of the matrix
   * @param res Output Matrix
   * @param other1 First input matrix
   * @param other2 Second input matrix
   * @note The dimensions of the matrices must be the same. If 'res' has different dimensions, it will be resized to fit the dimensions of 'other1' and 'other2
   */
  template<typename T>
  static void add(T& res, const T& other1, const T& other2)
  {
    static_assert(std::is_base_of_v<AMatrix, T>,
                  "Invalid type for Matrix addition (template method). "
                  "Only MatrixDense and MatrixSparse are allowed");
    // Check dimensions
    if (other1.getNRows() != other2.getNRows() ||
        other1.getNCols() != other2.getNCols())
    {
      messerr("Matrix addition error: matrices have different dimensions");
      return;
    }
    if (res.getNRows() != other1.getNRows() || res.getNCols() != other1.getNCols())
    {
      res.resize(other1.getNRows(), other1.getNCols());
    }
    res.addGeneral(res, other1, other2);
  }

  /**
   * @brief Operate product of a constant to all terms of a matrix and store the result in res = other1 * cst
   *
   * @tparam T Type of the matrix
   * @param res Output Matrix
   * @param other Input matrix
   * @param cst Constant value to be multiplied
   */
  template<typename T>
  static void prod(T& res, const T& other, double cst)
  {
    static_assert(std::is_base_of_v<AMatrix, T>,
                  "Invalid type for Matrix product (template method). "
                  "Only MatrixDense and MatrixSparse are allowed");
    if (res.getNRows() != other.getNRows() || res.getNCols() != other.getNCols())
    {
      res.resize(other.getNRows(), other.getNCols());
    }
    res.prodGeneral(res, other, cst);
  }

  /**
   * @brief Operate element-wise product of two matrices and store the result in a third one (res = other1 .* other2)
   *
   * @tparam T Type of the matrix
   * @param res Output Matrix
   * @param other1 First input matrix
   * @param other2 Second input matrix
   * @note The dimensions of the matrices must be the same. If 'res' has different dimensions, it will be resized to fit the dimensions of 'other1' and 'other2
   */
  template<typename T>
  static void prodHadamard(T& res, const T& other1, const T& other2)
  {
    static_assert(std::is_base_of_v<AMatrix, T>,
                  "Invalid type for Matrix Hadamard product (template method). "
                  "Only MatrixDense and MatrixSparse are allowed");
    // Check dimensions
    if (other1.getNRows() != other2.getNRows() ||
        other1.getNCols() != other2.getNCols())
    {
      messerr("Matrix Hadamard product error: matrices have different dimensions");
      return;
    }
    if (res.getNRows() != other1.getNRows() || res.getNCols() != other1.getNCols())
    {
      res.resize(other1.getNRows(), other1.getNCols());
    }
    res.prodGeneral(res, other1, other2);
  }

  /**
   * @brief Linear combination of up to 3 matrices and a constant
   *
   * @tparam T
   * @param res Resulting matrix of the linear combination
   * @param addition Constant to be added to the linear combination
   * @param val1 Weight of the first matrix
   * @param other1 First matrix
   * @param val2  Weight of the second matrix
   * @param other2 Second matrix
   * @param val3  Weight of the third matrix
   * @param other3 Third matrix
   */
  template<typename T>
  static void linearCombination(T& res,
                                double addition,
                                double val1,
                                const T& other1,
                                double val2     = 0.,
                                const T& other2 = T(),
                                double val3     = 0.,
                                const T& other3 = T())
  {
    static_assert(std::is_base_of_v<AMatrix, T>,
                  "Invalid type for Matrix addition (template method). "
                  "Only MatrixDense and MatrixSparse are allowed");

    // Check dimensions
    auto checkDims = [](const T& a, const T& b, const char* nameA, const char* nameB)
    {
      if (a.getNRows() != b.getNRows() || a.getNCols() != b.getNCols())
      {
        messerr("Matrix linear combination error: matrices '%s' and '%s' have different dimensions",
                nameA, nameB);
        return false;
      }
      return true;
    };

    if (val1 != 0. && val2 != 0. && !checkDims(other1, other2, "other1", "other2")) return;
    if (val1 != 0. && val3 != 0. && !checkDims(other1, other3, "other1", "other3")) return;
    if (val2 != 0. && val3 != 0. && !checkDims(other2, other3, "other2", "other3")) return;

    if (res.getNRows() != other1.getNRows() || res.getNCols() != other1.getNCols())
    {
      res.resize(other1.getNRows(), other1.getNCols());
    }

    T temp(res.getNRows(), res.getNCols());
    T cumul(res.getNRows(), res.getNCols());

    if (val1 != 0.)
    {
      T::prodGeneral(temp, other1, val1);
      T::addGeneral(cumul, cumul, temp);
    }

    if (val2 != 0.)
    {
      T::prodGeneral(temp, other2, val2);
      T::addGeneral(cumul, cumul, temp);
    }

    if (val3 != 0.)
    {
      T::prodGeneral(temp, other3, val3);
      T::addGeneral(cumul, cumul, temp);
    }
    if (addition != 0.)
      T::addGeneral(cumul, cumul, addition);

    res = cumul;
  }

  virtual void reset(Id nrows, Id ncols);
  virtual void resetFromValue(Id nrows, Id ncols, double value);
  virtual void resetFromArray(Id nrows, Id ncols, const double* tab, bool byCol = true);
  virtual void resetFromVD(Id nrows, Id ncols, const VectorDouble& tab, bool byCol = true);
  virtual void resetFromVVD(const VectorVectorDouble& tab, bool byCol = true);

  /*! Transpose the matrix in place*/
  virtual void transposeInPlace();
  /*! Transpose the matrix and return it as a copy*/
  virtual AMatrix* transpose() const;

  /*! Perform 'this' = 'x' * 'y' */
  virtual void prodMatMatInPlace(const AMatrix* x,
                                 const AMatrix* y,
                                 bool transposeX = false,
                                 bool transposeY = false);
  /*! Perform 'this' = 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' */
  virtual void prodNormMatMatInPlace(const AMatrix* a,
                                     const AMatrix* m,
                                     bool transpose = false);
  /*! Perform 'this' = 't(A)' %*% 'vec' %*% 'A' or 'A' %*% 'vec' %*% 't(A)' */
  virtual void prodNormMatVecInPlace(const AMatrix* a,
                                     const VectorDouble& vec,
                                     bool transpose = false);
  /*! Perform 'this' = 't(A)' %*% 'A' or 'A' %*% 't(A)' */
  virtual void prodNormMatInPlace(const AMatrix* a,
                                  bool transpose = false);
  /*! Add a matrix (multiplied by a constant) */
  virtual void addMat(const AMatrix& y, double cx = 1., double cy = 1.);

  /*! Extract the contents of the matrix */
  virtual NF_Triplet getMatrixToTriplet(Id shiftRow = 0, Id shiftCol = 0) const;

  void clear();
  /*! Modify the dimension of the matrix (if needed) */
  void resize(Id nrows, Id ncols);
  /*! Indicate if the given indices are valid for the current matrix size */
  bool isValid(Id irow, Id icol, bool printWhyNot = false) const;
  /*! Check if the matrix is square and Identity */
  bool isIdentity(bool printWhyNot = false) const;

  /*! Add a value to a matrix term */
  void addValue(Id irow, Id icol, double value);
  /*! Check if a matrix is the same as me (norm L1) */
  bool isSame(const AMatrix& m, double eps = EPSILON4, bool printWhyNot = false);
  /*! Check that 'm' has the same dimensions as 'this' */
  bool isSameSize(const AMatrix& m) const;
  /*! Returns if the current matrix is Empty */
  bool empty() const
  {
    return (_nRows == 0 || _nCols == 0);
  }
  /*! Returns the sum of absolute difference between argument and this */
  double compare(const AMatrix& mat) const;
  /*! Returns the number of rows */
  Id getNRows() const
  {
    return _nRows;
  }
  /*! Returns the number of columns */
  Id getNCols() const
  {
    return _nCols;
  }
  /*! Get the total number of elements of the (full) matrix (as for VectorT) */
  Id size() const
  {
    return _nRows * _nCols;
  }
  /*! Returns the contents of the whole matrix as a VectorDouble */
  VectorDouble getValues(bool byCol = true) const;
  /*! Extract a Diagonal (main or secondary) of this */
  const VectorDouble& getDiagonal(Id shift = 0) const;
  /*! Checks if a Column is valid (contains a non TEST value) */
  bool isColumnDefined(Id icol) const;
  /*! Checks if a Row is valid (contains a non TEST value) */
  bool isRowDefined(Id irow) const;
  /*! Define the number of defined columns */
  Id getNColDefined() const;
  /*! Define the number of defined rows */
  Id getNRowDefined() const;
  /*! Extract a portion of a Column */
  VectorDouble getColumnByRowRange(Id icol, Id rowFrom, Id rowTo) const;
  /*! Check if the matrix does not contain any negative element */
  bool isNonNegative(bool verbose = false) const;

  /*! Perform 'y' = 'this' * 'x' */
  VectorDouble prodMatVec(const VectorDouble& x, bool transpose = false) const;
  void prodMatVecInPlace(const VectorDouble& x, VectorDouble& y, bool transpose = false) const;
#ifndef SWIG
  void prodMatVecInPlaceC(const constvect x, vect y, bool transpose = false) const;
  void addProdMatVecInPlaceC(const constvect x, vect y, bool transpose = false) const;
#endif

  /*! Perform 'y' = 'x' * 'this' */
  VectorDouble prodVecMat(const VectorDouble& x, bool transpose = false) const;
  void prodVecMatInPlace(const VectorDouble& x, VectorDouble& y, bool transpose = false) const;
#ifndef SWIG
  void prodVecMatInPlaceC(const constvect x, vect y, bool transpose = false) const;
  void addProdVecMatInPlaceC(const constvect x, vect y, bool transpose = false) const;
#endif

  /*! Perform x %*% 'this' %*% y */
  double prodVecMatVec(const VectorDouble& x, const VectorDouble& y) const;
  /*! Perform 'this' = 'y' %*% 'this' or 'this' %*% 'y' */
  void prodMat(const AMatrix* matY, bool transposeY = false);

  /*! Matrix inversion in place */
  Id invert();
  /*! Solving the Matrix Linear system */
  Id solve(const VectorDouble& b, VectorDouble& x) const;
  /*! Dump a specific range of samples from the internal storage */
  void dumpElements(const String& title, Id ifrom, Id ito) const;
  /*! Dump statistics on the Matrix */
  void dumpStatistics(const String& title) const;
  /*! Sets the matrix as Identity */
  void setIdentity(double value = 1.);
  void fillRandom(double zeroPercent = 0, Id seed = 432432);
  void setValues(const VectorDouble& values, bool byCol = true);
  double getMeanByColumn(Id icol) const;
  double getMinimum() const;
  double getMaximum() const;
  double getNormInf() const;
  void copyReduce(const AMatrix* x,
                  const VectorInt& validRows,
                  const VectorInt& validCols);
  void copyElements(const AMatrix& m, double factor = 1.);

  void makePositiveColumn();

  void dumpRange(const char* title);

#ifndef SWIG
  /*! Get value operator */
  double operator()(Id row, Id col) const
  {
    return getValue(row, col);
  }
  /*! Set value operator */
  double& operator()(Id row, Id col)
  {
    return _getValueRef(row, col);
  }
#endif

protected:
  virtual void _allocate()                                        = 0;
  virtual void _deallocate()                                      = 0;
  virtual Id _getIndexToRank(Id irow, Id icol) const              = 0;
  virtual void _setValueByRank(Id rank, double value)             = 0;
  virtual void _transposeInPlace()                                = 0;
  virtual Id _invert()                                            = 0;
  virtual Id _solve(const VectorDouble& b, VectorDouble& x) const = 0;
  virtual Id _getMatrixPhysicalSize() const                       = 0;
  virtual double _getValueByRank(Id rank) const                   = 0;
  virtual double& _getValueRef(Id irow, Id icol)                  = 0;

#ifndef SWIG
  virtual void _addProdMatVecInPlacePtr(constvect x,
                                        vect y,
                                        bool transpose = false) const = 0;
  virtual void _addProdVecMatInPlacePtr(constvect x,
                                        vect y,
                                        bool transpose = false) const = 0;
#endif

  virtual bool _needToReset(Id nrows, Id ncols);
  virtual bool _isPhysicallyPresent(Id /*irow*/, Id /*icol*/) const
  {
    return true;
  }
  virtual void _setValues(const double* values, bool byCol);

  virtual void _clear();
  virtual bool _isNumbersValid(Id nrows, Id ncols) const;

  void _setNCols(Id ncols)
  {
    _nCols = ncols;
  }
  void _setNRows(Id nrows)
  {
    _nRows = nrows;
  }
  bool _isColumnValid(Id icol) const;
  bool _isRowValid(Id irow) const;
  bool _isIndexValid(Id irow, Id icol) const;
  bool _isRowVectorConsistent(const VectorDouble& tab) const;
  bool _isColVectorConsistent(const VectorDouble& tab) const;
  bool _isVectorSizeConsistent(const VectorDouble& tab) const;
  bool _isColumnSizeConsistent(const VectorDouble& tab) const;
  bool _isRowSizeConsistent(const VectorDouble& tab) const;
  bool _isRankValid(Id rank) const;
  void _fillFromVVD(const VectorVectorDouble& X);

  // Static functions
  static bool _isMatrixCompatible(const String& name,
                                  const AMatrix* mat1 = nullptr,
                                  Id vsize1           = 0,
                                  bool transpose1     = false,
                                  const AMatrix* mat2 = nullptr,
                                  Id vsize2           = 0,
                                  bool transpose2     = false,
                                  const AMatrix* mat3 = nullptr,
                                  Id vsize3           = 0,
                                  bool transpose3     = false);
  static bool _identifyRowAndCol(const AMatrix* mat,
                                 Id vsize,
                                 bool transpose,
                                 Id* nrow,
                                 Id* ncol);

private:
  bool _matrixNeedToReset(Id nrows, Id ncols);

private:
  mutable VectorDouble _diagonal;
  Id _nRows;
  Id _nCols;
  double _nullTerm; // Used for returning a null constant address
};

/* Shortcut functions for C style aficionados */
GSTLEARN_EXPORT I32 getMultiThread();
GSTLEARN_EXPORT void setMultiThread(I32 nthreads);
GSTLEARN_EXPORT bool isMultiThread();
GSTLEARN_EXPORT bool getFlagMatrixCheck();
GSTLEARN_EXPORT void setFlagMatrixCheck(bool flag);
} // namespace gstlrn
