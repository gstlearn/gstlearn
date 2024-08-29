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

#include "Matrix/MatrixSparse.hpp"
#include "LinearOp/ProjMulti.hpp"

#include "gstlearn_export.hpp"
#ifndef SWIG
  #include <Eigen/Core>
  #include <Eigen/Dense>
  #include <Eigen/src/Core/Matrix.h>
#endif

#include <algorithm>
class ProjMatrix;


class GSTLEARN_EXPORT ProjMultiMatrix : public ProjMulti, public MatrixSparse
{
public:
  ProjMultiMatrix(const std::vector<std::vector<const ProjMatrix*>> &proj);
  static std::vector<std::vector<const ProjMatrix*>> create(std::vector<const ProjMatrix*> &vectproj, int nvariable);
  virtual ~ProjMultiMatrix(){}

#ifndef SWIG           
  protected:
  virtual int _addPoint2mesh(const Eigen::VectorXd& inv,
                        Eigen::VectorXd& outv) const override;
  virtual int _addMesh2point(const Eigen::VectorXd& inv,
                        Eigen::VectorXd& outv) const override;
#endif

};