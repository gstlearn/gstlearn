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

#include "LinearOp/ALinearOp.hpp"
#include "LinearOp/Cholesky.hpp"
#include "gstlearn_export.hpp"
#include "LinearOp/SPDEOp.hpp"
#include "Matrix/MatrixSparse.hpp"

class PrecisionOpMultiMatrix;
class ProjMultiMatrix;
class MatrixSparse;


class GSTLEARN_EXPORT SPDEOpMatrix : public SPDEOp
{
public:
  SPDEOpMatrix(const PrecisionOpMultiMatrix* pop = nullptr, const ProjMultiMatrix* A = nullptr, const MatrixSparse* invNoise = nullptr);
  virtual ~SPDEOpMatrix();

#ifndef SWIG
private:
  int _addToDestImpl(const constvect inv, vect outv) const override;
  int _solve(const constvect inv, vect outv) const override;
#endif

private:
  mutable MatrixSparse _QpAinvNoiseAt; //mutable is required to perform the Cholesky decompositions
  mutable Cholesky*    _chol;         // when needed, e.g in a const method.
};
                               

