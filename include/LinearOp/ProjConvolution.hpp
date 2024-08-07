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

  int point2mesh(const VectorDouble& valonseismic, VectorDouble& valonvertex) const override;
  int mesh2point(const VectorDouble& valonvertex, VectorDouble& valonseismic) const override;

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
  bool _isVecDimCorrect(const VectorDouble &valonseismic,
                        const VectorDouble &valonvertex) const;

  void _convolve(const VectorDouble &valonvertex,
                 VectorDouble &valonseismic) const;
  void _convolveT(const VectorDouble &valonseismic,
                   VectorDouble &valonvertex) const;

private:
  VectorDouble         _convolution;
  const DbGrid*        _gridSeismic;
  VectorInt            _nodeRes2D;
  VectorDouble         _gext;
  VectorInt            _shiftVector;
  DbGrid*              _gridSeis2D;
  DbGrid*              _gridRes2D;
  MatrixSparse*        _AProjHoriz;
  mutable VectorDouble _work;
};

