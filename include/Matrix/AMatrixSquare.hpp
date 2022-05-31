/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 9 avr. 2019 by N. Desassis                                     */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Matrix/AMatrix.hpp"

/**
 * Square Matrix
 */
class GSTLEARN_EXPORT AMatrixSquare : public AMatrix {

protected:
  AMatrixSquare(int nrow = 0, bool sparse = false);
  AMatrixSquare(const AMatrixSquare &m);
  AMatrixSquare& operator= (const AMatrixSquare &r);
public:
	virtual ~AMatrixSquare();

  /*! Returns the size of the matrix (nrows=ncols) */
  int  getNSize() const { return getNRows(); }
  /*! Perform Norm matrix */
  void normMatrix(const AMatrixSquare& x, const AMatrix& y);
  void normTMatrix(const AMatrixSquare& x, const AMatrix& y);

  /*! Perform inner product */
  void innerMatrix(const AMatrixSquare& x,
                   const AMatrix& r1,
                   const AMatrix& r2);
  /*! Multiply the diagonal by a vector */
  void prodDiagByVector(const VectorDouble& diag);
  /*! Divide the diagonal by a vector */
  void divideDiagByVector(const VectorDouble& diag);
  /*! Returns the Determinant value */
  virtual double _determinant(void) const override;

protected:
  void   _setNSize(int nval);
  bool   _isNumberValid(int nrows,int ncols) const;
};
