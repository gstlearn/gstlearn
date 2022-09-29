/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 9 avr. 2019 by N. Desassis                                     */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/

#include "LinearOp/ProjConvolution.hpp"

ProjConvolution::ProjConvolution(const VectorDouble& convolution, const DbGrid* seismicOut)
 : _convolution(convolution)
 , _meshOut()
{
 _meshOut = MeshETurbo(seismicOut);
 /* TODO : tester 3D */
}

ProjConvolution::~ProjConvolution() { }

int ProjConvolution::point2mesh(const VectorDouble& in, VectorDouble& out) const
{
  return 0;
}

int ProjConvolution::ProjConvolution::mesh2point(const VectorDouble& valonvertex, VectorDouble& valonseismic) const
{
  return 0;
}


int ProjConvolution::ProjConvolution::getApexNumber() const
{
  Grid grid = _meshOut.getGrid();
  VectorInt nxs = grid.getNXs();
  nxs[grid.getNDim()-1]-= _convolution.size();
  return ut_vector_prod(nxs);
}

int ProjConvolution::ProjConvolution::getPointNumber() const
{
  return _meshOut.getNApices();
}
