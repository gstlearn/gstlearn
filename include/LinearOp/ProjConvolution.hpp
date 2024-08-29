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

#include "Db/DbGrid.hpp"
#include "IProjMatrix.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Basic/VectorNumT.hpp"
#include "Matrix/MatrixSparse.hpp"
#include <Eigen/src/Core/Matrix.h>

#ifndef SWIG
  #include <Eigen/Core>
  #include <Eigen/Dense>
#endif
/**
 * Projection matrix for vertical convolution
 */
class GSTLEARN_EXPORT ProjConvolution: public IProjMatrix
{

public:
  ProjConvolution(const VectorDouble &convolution = VectorDouble(),
                  const DbGrid *grid_point = nullptr,
                  const VectorInt& nodeRes2D = VectorInt(),
                  const VectorDouble& gext = VectorDouble());
  ProjConvolution(const ProjConvolution &m)= delete;
  ProjConvolution& operator= (const ProjConvolution &m)= delete;
  virtual ~ProjConvolution();

  int getApexNumber() const override;
  int getPointNumber() const override;

  DbGrid* getResolutionGrid() const;

  /// TODO : return a shared pointer ?
  const MatrixSparse* getAProjHoriz() const { return _AProjHoriz; }

  const VectorDouble& getConvolution() const { return _convolution; }
  const VectorInt& getShiftVector() const { return _shiftVector; }

private:
  int  _getConvSize() const { return (int) _convolution.size(); }
  int  _getHalfSize() const { return (_getConvSize() - 1) / 2; }
  void _buildGridSeis2D();
  void _buildGridRes2D();
  void _buildShiftVector();
  int  _buildAprojHoriz();
  int  _getNDim() const { return _gridSeismic->getNDim(); }
  Grid _getGridCharacteristicsRR(bool delLastDim = false) const;
  Grid _getGridCharacteristicsRS() const;

  

  #ifndef SWIG        
  private:
  void _convolve(const Eigen::VectorXd &valonvertex,
                 Eigen::VectorXd &valonseismic) const;
  void _convolveT(const Eigen::VectorXd &valonseismic,
                   Eigen::VectorXd &valonvertex) const;
  bool _isVecDimCorrect(const Eigen::VectorXd &valonseismic,
                        const Eigen::VectorXd &valonvertex) const;   
  protected:
  int _addPoint2mesh(const Eigen::VectorXd& valonseismic,
                        Eigen::VectorXd& valonvertex) const override;
  int _addMesh2point(const Eigen::VectorXd& valonvertex,
                        Eigen::VectorXd& valonseismic) const override;
  #endif

private:
  VectorDouble            _convolution;
  const DbGrid*           _gridSeismic;
  VectorInt               _nodeRes2D;
  VectorDouble            _gext;
  VectorInt               _shiftVector;
  DbGrid*                 _gridSeis2D;
  DbGrid*                 _gridRes2D;
  MatrixSparse*           _AProjHoriz;
  mutable Eigen::VectorXd _work;
};

