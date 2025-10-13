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

#include "Matrix/MatrixDense.hpp"
#include "gstlearn_export.hpp"

/**
 * Square Matrix
 */
namespace gstlrn
{
class GSTLEARN_EXPORT MatrixSquare: public MatrixDense
{

public:
  MatrixSquare(Id nrow = 0);
  MatrixSquare(const MatrixSquare& r);
  MatrixSquare(const AMatrix& m);
  MatrixSquare& operator=(const MatrixSquare& r);
  virtual ~MatrixSquare();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// ICloneable interface
  IMPLEMENT_CLONING(MatrixSquare)

  /// Interface for AMatrix
  virtual double determinant(void) const;
  /*! Check if the matrix is (non empty) square */
  bool isSquare(bool printWhyNot = false) const override
  {
    DECLARE_UNUSED(printWhyNot);
    return true;
  }
  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const override { return false; }

  /*! Returns the size of the matrix (nrows=ncols) */
  Id getNSize() const { return getNRows(); }
  void resetFromVVD(const VectorVectorDouble& tab, bool byCol = true) override;

  static MatrixSquare* createFromVVD(const VectorVectorDouble& X);
  static MatrixSquare* createFromVD(const VectorDouble& X,
                                    Id nrow,
                                    bool byCol             = false,
                                    bool invertColumnOrder = false);
  static MatrixSquare* createFromTridiagonal(const VectorDouble& vecdiag,
                                             const VectorDouble& vecinf,
                                             const VectorDouble& vecsup);

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
  void prodByDiagInPlace(Id mode, const VectorDouble& c);

  double normVec(const VectorDouble& vec);
  Id decomposeLU(MatrixSquare& tls,
                 MatrixSquare& tus,
                 double eps = EPSILON20);

  Id computeEigen(bool optionPositive = true);

protected:
  bool _isNumbersValid(Id nrows, Id ncols) const override;
  void _setNSize(Id nval);

private:
  Id _invertLU();
  Id _solveLU(const MatrixSquare& tus,
              const MatrixSquare& tls,
              const double* b,
              double* x);
  Id _forwardLU(const MatrixSquare& tls, const double* b, double* x, double eps = EPSILON20);
  Id _backwardLU(const MatrixSquare& tus, const double* b, double* x, double eps = EPSILON20);
};

/*! Product 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' */
GSTLEARN_EXPORT MatrixSquare* prodNormMatMat(const MatrixDense* a,
                                             const MatrixDense* m,
                                             bool transpose = false);
/*! Product 't(A)' %*% 'A' or 'A' %*% 't(A)' */
GSTLEARN_EXPORT MatrixSquare* prodNormMat(const MatrixDense& a,
                                          bool transpose = false);
/*! Product 't(A)' %*% 'vec' %*% 'A' or 'A' %*% 'vec' %*% 't(A)' */
GSTLEARN_EXPORT MatrixSquare* prodNormMatVec(const MatrixDense& a,
                                             const VectorDouble& vec,
                                             bool transpose = false);
} // namespace gstlrn