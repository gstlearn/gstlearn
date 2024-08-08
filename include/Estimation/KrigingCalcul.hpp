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

#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/AMatrix.hpp"

class GSTLEARN_EXPORT KrigingCalcul
{
public:
  KrigingCalcul(const MatrixSquareSymmetric* sigCma = nullptr,
                const MatrixRectangular* X = nullptr);
  KrigingCalcul(const KrigingCalcul &r) = delete;
  KrigingCalcul& operator=(const KrigingCalcul &r) = delete;
  virtual ~KrigingCalcul();

  int setC(const MatrixSquareSymmetric* C);
  int setX(const MatrixRectangular* X);
  int setC0(const MatrixRectangular* C0);
  int setX0(const MatrixRectangular* X0);
  int setZ(const VectorDouble& Z);
  int setBeta(const VectorDouble& beta);

private:
  bool _matchDimensions(const AMatrix* mat, int nrowsRef, int ncolsRef);
  int _invertC();

private:
  const MatrixSquareSymmetric* _C; // Pointer not to be deleted
  const MatrixRectangular* _X;     // Pointer not to be deleted
  const MatrixRectangular* _C0;    // Pointer not to be deleted
  const MatrixRectangular* _X0;    // Pointer not to be deleted
  VectorDouble _Z;
  VectorDouble _beta;
  MatrixRectangular _lambda;
  MatrixSquareSymmetric* _Cm1;
  int _neq;
  int _nbfl;
  int _nrhs;
};
