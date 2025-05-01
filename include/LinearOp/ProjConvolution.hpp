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
#include "IProj.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Basic/VectorNumT.hpp"
#include "Matrix/MatrixSparse.hpp"

/**
 * Projection matrix for vertical convolution
 */
class GSTLEARN_EXPORT ProjConvolution: public IProj
{

public:
  ProjConvolution(const VectorDouble &convolution = VectorDouble(),
                  const DbGrid *grid_point = nullptr,
                  const VectorInt& nodeRes2D = VectorInt(),
                  const VectorDouble& gext = VectorDouble());
  ProjConvolution(const ProjConvolution &m)= delete;
  ProjConvolution& operator= (const ProjConvolution &m)= delete;
  virtual ~ProjConvolution();

  int getNApex() const override;
  int getNPoint() const override;

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
    void _convolve(const constvect valonvertex, vect valonseismic) const;
    void _convolveT(const constvect valonseismic, vect valonvertex) const;
    bool _isVecDimCorrect(const constvect valonseismic,
                          const constvect valonvertex) const;

  protected:
    int _addPoint2mesh(const constvect valonseismic,
                       vect valonvertex) const override;
    int _addMesh2point(const constvect valonvertex,
                       vect valonseismic) const override;
#endif

private:
  VectorDouble                _convolution;
  const DbGrid*               _gridSeismic;
  VectorInt                   _nodeRes2D;
  VectorDouble                _gext;
  VectorInt                   _shiftVector;
  DbGrid*                     _gridSeis2D;
  DbGrid*                     _gridRes2D;
  MatrixSparse*               _AProjHoriz;
  mutable std::vector<double> _work;
};

