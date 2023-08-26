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

  void init(int nrows, int ncols);

#ifndef SWIG
  /*! Extract the contents of the matrix */
  void getValuesAsTriplets(VectorInt &irows,
                           VectorInt &icols,
                           VectorDouble &values) const;
#endif

#ifndef SWIG
  /*! Returns a pointer to the Sparse storage */
  const cs* getCs() const { return _csMatrix; }
#endif
  Triplet getCsToTriplet(bool flag_from_1 = false) const;

  /*! Returns if the current matrix is Sparse */
  bool isSparse() const { return true; }

  void reset(int nrows, int ncols);
  void reset(int nrows, int ncols, double value);
  void reset(int nrows, int ncols, const double* tab, bool byCol = true);
  void reset(int nrows,
             int ncols,
             const VectorDouble &tab,
             bool byCol = true);
  void reset(const VectorVectorDouble& tab, bool byCol = true);

  /*! Transpose the matrix and return it as a copy*/
  virtual MatrixSparse* transpose() const override;
  /*! Add a value to each matrix component */
  virtual void addScalar(double v) override;
  /*! Add value to matrix diagonal */
  virtual void addScalarDiag(double v) override;
  /*! Multiply each matrix component by a value */
  virtual void prodScalar(double v) override;

  /*! Add a matrix to this component by component */
  void addMatrix(const MatrixSparse& y);
  /*! Multiply a matrix by another and store the result in the current matrix */
  void prodMatrix(const MatrixSparse& x, const MatrixSparse& y);
  /*! Linear combination of matrices */
  void linearCombination(double cx, double cy, const MatrixSparse& y);

  /*! Dump a specific range of samples from the internal storage */
  void dumpElements(const String& title, int ifrom, int ito) const;
  /*! Conversion to a string */
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /*! Set all the values of the Matrix at once */
  void fill(double value);
  void fillRandom(int seed = 432432, double zeroPercent = 0.1);
  void setValuesByArrays(const VectorInt &irows,
                         const VectorInt &icols,
                         const VectorDouble &values) override;

protected:
#ifndef SWIG
  virtual double& _getValueRef(int irow, int icol) override;

protected:
  bool    _isPhysicallyPresent(int /*irow*/, int /*icol*/) const { return true; }
  bool    _isCompatible(const AMatrix& m) const override { return (isSparse()); }
  void    _allocate() override;
  void    _deallocate() override;

  /*! Returns the number of elements actually stored as members in subsequent classes */
  virtual int     _getMatrixSize() const override;
  virtual void    _setValue(int rank, double value) override;
  virtual void    _setValue(int irow, int icol, double value) override;
  virtual void    _setValues(const double* values, bool byCol) override;
  virtual double  _getValue(int rank) const override;
  virtual double  _getValue(int irow, int icol) const override;
  virtual void    _transposeInPlace() override;
  virtual void    _prodVector(const double *inv,double *outv) const override;
  virtual int     _invert() override;
  virtual int     _solve(const VectorDouble& b, VectorDouble& x) const override;
  virtual void    _clearContents() {};
  virtual int     _getIndexToRank(int irow,int icol) const override;

  void _clear();

private:
  void _initiateSparse();
  void _recopySparse(const cs* cs);
  void _forbiddenForSparse(const String& func) const;

private:
  cs*  _csMatrix;
#endif
};

/*! Transform any matrix in a Sparse format */
GSTLEARN_EXPORT MatrixSparse* toSparse(const AMatrix* mat);

