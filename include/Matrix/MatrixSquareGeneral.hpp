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
#include "Matrix/AMatrixSquare.hpp"

/**
 * Square Matrix General
 */
class GSTLEARN_EXPORT MatrixSquareGeneral : public AMatrixSquare {

public:
  MatrixSquareGeneral(int nrow = 0);
  MatrixSquareGeneral(const MatrixSquareGeneral &r);
  MatrixSquareGeneral(const AMatrix &m);
  MatrixSquareGeneral& operator= (const MatrixSquareGeneral &r);
	virtual ~MatrixSquareGeneral();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// ICloneable interface
  IMPLEMENT_CLONING(MatrixSquareGeneral)

  /// Interface for AMatrix
  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const override { return false; }

  static MatrixSquareGeneral* createFromVVD(const VectorVectorDouble& X);
  static MatrixSquareGeneral* createFromVD(const VectorDouble &X,
                                           int nrow,
                                           bool byCol = false,
                                           bool invertColumnOrder = false);

  int decomposeLU(MatrixSquareGeneral& tls,
                  MatrixSquareGeneral& tus,
                  double eps = EPSILON20);

private:
  int     _invertLU();
  int     _solveLU(const MatrixSquareGeneral& tus,
                   const MatrixSquareGeneral& tls,
                   const double *b,
                   double *x);
  int     _forwardLU(const MatrixSquareGeneral& tls, const double *b, double *x, double eps = EPSILON20);
  int     _backwardLU(const MatrixSquareGeneral& tus, const double *b, double *x, double eps = EPSILON20);
};

/*! Product 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' */
GSTLEARN_EXPORT MatrixSquareGeneral* prodNormMatMat(const AMatrixDense* a,
                                                    const AMatrixDense* m,
                                                    bool transpose = false);
/*! Product 't(A)' %*% 'A' or 'A' %*% 't(A)' */
GSTLEARN_EXPORT MatrixSquareGeneral* prodNormMat(const AMatrixDense &a,
                                                 const VectorDouble& vec = VectorDouble(),
                                                 bool transpose = false);
