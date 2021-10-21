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

#include "MatrixC/AMatrixC.hpp"

/**
 * Square Matrix
 */
class AMatrixCSquare : public AMatrixC {

protected:
  AMatrixCSquare(int nrow = 0, bool sparse = false);
  AMatrixCSquare(const AMatrixCSquare &m);
  AMatrixCSquare& operator= (const AMatrixCSquare &r);
public:
	virtual ~AMatrixCSquare();

  /*! Returns the size of the matrix (nrows=ncols) */
  int  getNSize() const { return getNRows(); }
  /*! Perform Norm matrix */
  void normMatrix(const AMatrixCSquare& x, const AMatrixC& y);
  /*! Perform inner product */
  void innerMatrix(const AMatrixCSquare& x,
                   const AMatrixC& r1,
                   const AMatrixC& r2);
  /*! Multiply the diagonal by a vector */
  void prodDiagByVector(const VectorDouble& diag);
  /*! Divide the diagonal by a vector */
  void divideDiagByVector(const VectorDouble& diag);

protected:
  void   _setNSize(int nval);
  bool   _isNumberValid(int nrows,int ncols) const;
};
