/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

class AMatrix;
class AMatrixSquare;

class GSTLEARN_EXPORT MatrixFactory {

public:
  MatrixFactory();
  virtual ~MatrixFactory();

  /// TODO : Use smartpointer
  static AMatrix* matProduct(const AMatrix* x, const AMatrix* y);
  static AMatrixSquare* matNorm(const AMatrixSquare* x, const AMatrix* y);
  static AMatrix* createIdentity(int nrow, bool sparse);
  static AMatrixSquare* createMatrixSquare(const AMatrixSquare* x,int nrow);
};
