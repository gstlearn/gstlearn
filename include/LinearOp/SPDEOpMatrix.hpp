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
  int _addToDestImpl(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const override;
  int _solve(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const override;
#endif

private:
  mutable MatrixSparse _QpAinvNoiseAt; //mutable is required to perform the Cholesky decompositions
                                       // when needed, e.g in a const method.
};
                               

