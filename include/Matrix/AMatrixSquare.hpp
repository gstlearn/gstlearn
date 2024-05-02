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
#include "Matrix/MatrixRectangular.hpp"

/**
 * Square Matrix
 */
class GSTLEARN_EXPORT AMatrixSquare : public MatrixRectangular {

public:
  AMatrixSquare(int nrow = 0, int opt_eigen = -1);
  AMatrixSquare(const AMatrixSquare &m);
  AMatrixSquare(const AMatrix &m);
  AMatrixSquare& operator= (const AMatrixSquare &r);
	virtual ~AMatrixSquare();

	/// Interface for AMatrix
  virtual double determinant(void) const;
  /*! Check if the matrix is (non empty) square */
  bool isSquare(bool printWhyNot = false) const override { DECLARE_UNUSED(printWhyNot); return true; }

  /*! Returns the size of the matrix (nrows=ncols) */
  int getNSize() const { return getNRows(); }

  double trace() const;

  /*! Perform inner product */
  void innerMatrix(const AMatrixSquare& x,
                   const AMatrix& r1,
                   const AMatrix& r2);
  /*! Multiply the diagonal by a vector */
  void prodDiagByVector(const VectorDouble& diag);
  /*! Divide the diagonal by a vector */
  void divideDiagByVector(const VectorDouble& diag);
  /*! Multiply by a Diagonal matrix provided as VectorDouble (in place) */
  void prodByDiagInPlace(int mode, const VectorDouble& c);

  double normVec(const VectorDouble& vec);

protected:
  void   _setNSize(int nval);
};
