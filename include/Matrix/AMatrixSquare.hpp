/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Matrix/AMatrix.hpp"

/**
 * Square Matrix
 */
class GSTLEARN_EXPORT AMatrixSquare : public AMatrix {

protected:
  AMatrixSquare(int nrow = 0);
  AMatrixSquare(const AMatrixSquare &m);
  AMatrixSquare& operator= (const AMatrixSquare &r);

public:
	virtual ~AMatrixSquare();

  /*! Returns the size of the matrix (nrows=ncols) */
  int  getNSize() const { return getNRows(); }
  /*! Perform Norm matrix */
  void normMatrix(const AMatrixSquare& x, const AMatrix& y);
  void normTMatrix(const AMatrixSquare& x, const AMatrix& y);

  double trace() const;

  /*! Perform inner product */
  void innerMatrix(const AMatrixSquare& x,
                   const AMatrix& r1,
                   const AMatrix& r2);
  /*! Multiply the diagonal by a vector */
  void prodDiagByVector(const VectorDouble& diag);
  /*! Divide the diagonal by a vector */
  void divideDiagByVector(const VectorDouble& diag);
  /*! Returns the Determinant value */
  virtual double determinant(void) const;

protected:
  void   _setNSize(int nval);
  bool   _isNumberValid(int nrows,int ncols) const;
};
