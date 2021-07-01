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

class AMatrixC;
class AMatrixCSquare;

class MatrixCFactory {

public:
  MatrixCFactory();
  virtual ~MatrixCFactory();

  /// TODO : Use smartpointer
  static AMatrixC* matProduct(const AMatrixC* x, const AMatrixC* y);
  static AMatrixCSquare* matNorm(const AMatrixCSquare* x, const AMatrixC* y);
  static AMatrixC* createIdentity(int nrow, bool sparse);
};
