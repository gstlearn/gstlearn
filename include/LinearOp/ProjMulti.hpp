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
#ifndef SWIG
  #include <Eigen/Core>
  #include <Eigen/Dense>
  #include <Eigen/src/Core/Matrix.h>
#endif
#include "LinearOp/IProjMatrix.hpp"

class GSTLEARN_EXPORT ProjMulti : public IProjMatrix
{
public:
  ProjMulti(const std::vector<std::vector<const IProjMatrix*>> &proj);
  int getApexNumber() const override;
  int getPointNumber() const override;
  int getNVariable() const { return _nvariable; }
  int getNLatent() const { return _nlatent; }
  virtual ~ProjMulti(){}

#ifndef SWIG           
  protected:
  virtual int _addPoint2mesh(const Eigen::VectorXd& inv,
                        Eigen::VectorXd& outv) const override;
  virtual int _addMesh2point(const Eigen::VectorXd& inv,
                        Eigen::VectorXd& outv) const override;
#endif

protected:
int findFirstNoNullOnRow(int j) const;
int findFirstNoNullOnCol(int j) const;
const std::vector<int>& getPointNumbers() const {return _pointNumbers;}
const std::vector<int>& getApexNumbers()  const {return _apexNumbers;}

protected:
std::vector<std::vector<const IProjMatrix*> >_projs;

private: 
  void _makeEmpty();
private:
int _pointNumber;
int _apexNumber;
int _nlatent;
int _nvariable;
std::vector<int> _pointNumbers;
std::vector<int> _apexNumbers;
mutable Eigen::VectorXd _work;
mutable Eigen::VectorXd _workmesh;


};