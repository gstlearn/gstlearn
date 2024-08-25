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
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/Matrix.h>
#include "LinearOp/IProjMatrix.hpp"
#include "Basic/VectorNumT.hpp"

class ProjMatrix;

class GSTLEARN_EXPORT ProjMatrixMulti : public IProjMatrix
{
public:
  ProjMatrixMulti(const std::vector<ProjMatrix*> &proj,
                  int nvar = 1);
  //int point2mesh(const VectorDouble& inv, VectorDouble& outv) const override;
  //int mesh2point(const VectorDouble& inv, VectorDouble& outv) const override;
  int getApexNumber() const override;
  int getPointNumber() const override;
  virtual ~ProjMatrixMulti(){}
private:
  const std::vector<ProjMatrix*> _projs;
  int _apicesNumber;
  int _pointsNumber;
  const int _nvar;

#ifndef SWIG           
  protected:
  int _point2mesh(const Eigen::VectorXd& inv,
                        Eigen::VectorXd& outv) const override;
  int _mesh2point(const Eigen::VectorXd& inv,
                        Eigen::VectorXd& outv) const override;
#endif

};