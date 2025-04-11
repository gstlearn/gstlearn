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
#include "Matrix/MatrixDense.hpp"

/**
 * Square Matrix
 */
class GSTLEARN_EXPORT MatrixSquare : public MatrixDense {

public:
  MatrixSquare(int nrow = 0);
  MatrixSquare(const MatrixSquare &r);
  MatrixSquare(const AMatrix &m);
  MatrixSquare& operator= (const MatrixSquare &r);
  virtual ~MatrixSquare();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// ICloneable interface
  IMPLEMENT_CLONING(MatrixSquare)

  /// Interface for AMatrix
  virtual double determinant(void) const;
  /*! Check if the matrix is (non empty) square */
  bool isSquare(bool printWhyNot = false) const override { DECLARE_UNUSED(printWhyNot); return true; }
  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const override { return false; }

  /*! Returns the size of the matrix (nrows=ncols) */
  int getNSize() const { return getNRows(); }
  void resetFromVVD(const VectorVectorDouble& tab, bool byCol = true) override;

  static MatrixSquare* createFromVVD(const VectorVectorDouble& X);
  static MatrixSquare* createFromVD(const VectorDouble& X,
                                     int nrow,
                                     bool byCol             = false,
                                     bool invertColumnOrder = false);
  double trace() const;

  /*! Perform inner product */
  void innerMatrix(const MatrixSquare& x,
                   const AMatrix& r1,
                   const AMatrix& r2);
  /*! Multiply the diagonal by a vector */
  void prodDiagByVector(const VectorDouble& diag);
  /*! Divide the diagonal by a vector */
  void divideDiagByVector(const VectorDouble& diag);
  /*! Multiply by a Diagonal matrix provided as VectorDouble (in place) */
  void prodByDiagInPlace(int mode, const VectorDouble& c);

  double normVec(const VectorDouble& vec);
  int decomposeLU(MatrixSquare& tls,
                  MatrixSquare& tus,
                  double eps = EPSILON20);

protected:
  bool _isNumbersValid(int nrows, int ncols) const override;
  void _setNSize(int nval);

private:
  int _invertLU();
  int _solveLU(const MatrixSquare& tus,
               const MatrixSquare& tls,
               const double* b,
               double* x);
  int _forwardLU(const MatrixSquare& tls, const double* b, double* x, double eps = EPSILON20);
  int _backwardLU(const MatrixSquare& tus, const double* b, double* x, double eps = EPSILON20);
};

/*! Product 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' */
GSTLEARN_EXPORT MatrixSquare* prodNormMatMat(const MatrixDense* a,
                                              const MatrixDense* m,
                                              bool transpose = false);
/*! Product 't(A)' %*% 'A' or 'A' %*% 't(A)' */
GSTLEARN_EXPORT MatrixSquare* prodNormMat(const MatrixDense& a,
                                           const VectorDouble& vec = VectorDouble(),
                                           bool transpose          = false);
