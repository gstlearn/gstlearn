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
#include "Matrix/AMatrix.hpp"
#include "Matrix/LinkMatrixSparse.hpp"
#include <Eigen/Sparse>

/**
 * Sparse Matrix
 */
class GSTLEARN_EXPORT MatrixSparse : public AMatrix {

public:
  MatrixSparse(int nrow = 0, int ncol = 0);
#ifndef SWIG
  MatrixSparse(const cs* A);
#endif
  MatrixSparse(const MatrixSparse &m);
  MatrixSparse& operator= (const MatrixSparse &m);
  virtual ~MatrixSparse();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// Cloneable interface
  IMPLEMENT_CLONING(MatrixSparse)

  /// Interface for AMatrix
  /*! Returns if the current matrix is Sparse */
  bool isSparse() const { return true; }
  /*! Returns if the matrix belongs to the MatrixSparse class (avoids dynamic_cast) */
  virtual bool isMatrixSparse() const { return true; }

  /*! Set the contents of a Column */
  virtual void setColumn(int icol, const VectorDouble& tab) override;
  /*! Set the contents of a Row */
  virtual void setRow(int irow, const VectorDouble& tab) override;
  /*! Set the contents of the (main) Diagonal */
  virtual void setDiagonal(const VectorDouble& tab) override;
  /*! Set the contents of the (main) Diagonal to a constant value */
  virtual void setDiagonalToConstant(double value = 1.) override;
  /*! Transpose the matrix and return it as a copy*/
  virtual MatrixSparse* transpose() const override;
  /*! Add a value to each matrix component */
  virtual void addScalar(double v) override;
  /*! Add value to matrix diagonal */
  virtual void addScalarDiag(double v) override;
  /*! Multiply each matrix component by a value */
  virtual void prodScalar(double v) override;
  /*! Set a set of values simultaneously from an input array */
  void setValuesFromTriplet(const Triplet& T) override;
  /*! Set all the values of the Matrix at once */
  virtual void fill(double value) override;
  /*! Multiply a Matrix row-wise */
  virtual void multiplyRow(const VectorDouble& vec) override;
  /*! Multiply a Matrix column-wise */
  virtual void multiplyColumn(const VectorDouble& vec) override;
  /*! Divide a Matrix row-wise */
  virtual void divideRow(const VectorDouble& vec) override;
  /*! Divide a Matrix column-wise */
  virtual void divideColumn(const VectorDouble& vec) override;
  /*! Perform M * 'vec' */
  virtual VectorDouble prodVector(const VectorDouble& vec) const override;
  /*! Perform 'vec'^T * M */
  virtual VectorDouble prodTVector(const VectorDouble& vec) const override;

#ifndef SWIG
  /*! Extract the contents of the matrix */
  virtual Triplet getValuesAsTriplets() const override;
#endif

  //// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  void init(int nrows, int ncols);

  /// The next functions use specific definition of matrix (to avoid dynamic_cast)
  /// rather than manipulating AMatrix. They are no more generic of AMatrix
  /*! Add a matrix (multiplied by a constant) */
  virtual void addMatrix(const MatrixSparse& y, double value = 1.);
  /*! Multiply a matrix by another and store the result in the current matrix */
  virtual void prodMatrix(const MatrixSparse& x, const MatrixSparse& y);
  /*! Linear combination of matrices */
  virtual void linearCombination(double cx, double cy, const MatrixSparse& y);

#ifndef SWIG
  /*! Returns a pointer to the Sparse storage */
  const cs* getCs() const { return _csMatrix; }
#endif
  Triplet getSparseToTriplet(bool flag_from_1 = false) const;


  void reset(int nrows, int ncols);
  void reset(int nrows, int ncols, double value);
  void reset(int nrows, int ncols, const double* tab, bool byCol = true);
  void reset(int nrows,
             int ncols,
             const VectorDouble &tab,
             bool byCol = true);
  void reset(const VectorVectorDouble& tab, bool byCol = true);


  /*! Dump a specific range of samples from the internal storage */
  void dumpElements(const String& title, int ifrom, int ito) const;

  /*! Set all the values of the Matrix with random values */
  void fillRandom(int seed = 432432, double zeroPercent = 0.1);


protected:
  /// Interface for AMatrix
  bool    _isPhysicallyPresent(int irow, int icol) const { DECLARE_UNUSED(irow, icol); return true; }
  bool    _isCompatible(const AMatrix& m) const override { DECLARE_UNUSED(m); return (isSparse()); }
  void    _allocate() override;
  void    _deallocate() override;

  virtual double& _getValueRef(int irow, int icol) override;
  virtual int     _getMatrixPhysicalSize() const override;
  virtual void    _setValueByRank(int rank, double value) override;
  virtual void    _setValue(int irow, int icol, double value) override;
  virtual void    _setValues(const double* values, bool byCol) override;
  virtual double  _getValueByRank(int rank) const override;
  virtual double  _getValue(int irow, int icol) const override;
  virtual int     _getIndexToRank(int irow,int icol) const override;
  virtual void    _transposeInPlace() override;

  virtual void    _prodVectorInPlace(const double *inv,double *outv) const override;
  virtual int     _invert() override;
  virtual int     _solve(const VectorDouble& b, VectorDouble& x) const override;

  void _clear();
  bool _isElementPresent(int irow, int icol) const;

private:
  void _forbiddenForSparse(const String& func) const;

private:
  cs*  _csMatrix; // Classical storage for Sparse matrix
  Eigen::SparseMatrix<double> _eigenMatrix; // Eigen storage in Eigen Library (always stored Eigen::ColMajor)
};

/*! Transform any matrix in a Sparse format */
GSTLEARN_EXPORT MatrixSparse *createFromAnyMatrix(const AMatrix* mat);
GSTLEARN_EXPORT void setUpdateNonZeroValue(int status = 2);
GSTLEARN_EXPORT int getUpdateNonZeroValue();
