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

#include "LinearOp/SPDEOp.hpp"
#include "gstlearn_export.hpp"

class PrecisionOpMultiMatrix;
class ProjMultiMatrix;
class MatrixSparse;


class GSTLEARN_EXPORT SPDEOpMatrix : public SPDEOp
{
public:
  SPDEOpMatrix(const PrecisionOpMultiMatrix* pop = nullptr, const ProjMultiMatrix* A = nullptr, const MatrixSparse* invNoise = nullptr);
  virtual ~SPDEOpMatrix();

  int getSize() const override;

#ifndef SWIG
protected:
  int _addToDest(const Eigen::VectorXd& inv, Eigen::VectorXd& outv) const override;
#endif

private:
  MatrixSparse* _QpAinvNoiseAt;
};

