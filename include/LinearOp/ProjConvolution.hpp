/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Db/DbGrid.hpp"
#include "IProjMatrix.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Basic/Vector.hpp"


/**
 * Projection matrix for vertical convolution
 */
class GSTLEARN_EXPORT ProjConvolution: public IProjMatrix
{

public:
  ProjConvolution(const VectorDouble &convolution = VectorDouble(),
                  const DbGrid *grid_point = nullptr);
  ProjConvolution(const ProjConvolution &m)= delete;
  ProjConvolution& operator= (const ProjConvolution &m)= delete;
  virtual ~ProjConvolution();

  int point2mesh(const VectorDouble& in, VectorDouble& out) const override;
  int mesh2point(const VectorDouble& valonvertex, VectorDouble& valonseismic) const override;

  int getApexNumber() const override;
  int getPointNumber() const override;

private:
  int _getConvSize() const { return (int) _convolution.size(); }
  int _getHalfSize() const { return (_getConvSize() - 1) / 2; }
  VectorInt _getShiftVector() const;

private:
  VectorDouble  _convolution;
  const DbGrid* _gridPoint;
  mutable cs* _Aproj; // Stockage temporaire de la matrice creuse de Projection
};

