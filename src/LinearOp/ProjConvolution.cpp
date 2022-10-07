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
#include "csparse_f.h"

ProjConvolution::ProjConvolution(const VectorDouble &convolution,
                                 const DbGrid *grid_point)
    : _convolution(convolution),
      _gridPoint(grid_point),
      _shiftVector(),
      _Aproj(nullptr)
{
  _constructShiftVector();

  _constructAprojCS();
}

ProjConvolution::~ProjConvolution() { }

/**
 * Calculate the Aproj sparse matrix.
 * This method is kept for establishing time bench marks.
 * It emulates mesh2point algorithm.
 * Note that this algorithm does not handle the presence of undefined values
 */
int ProjConvolution::_constructAprojCS()
{
  cs* Atriplet;
  Atriplet = cs_spalloc(0, 0, 1, 1, 1);
  if (_Aproj != nullptr) _Aproj = cs_spfree(_Aproj);

  for (int is = 0; is < getPointNumber(); is++)
  {
    for (int i = -_getHalfSize(); i <= _getHalfSize(); i++)
    {
      int j = i + _getHalfSize();
      int id = is + _shiftVector[j];
      if (id < 0) return 1;

      (void) cs_entry(Atriplet,is,id,_convolution[j]);
    }
  }

  if (_Aproj != nullptr) _Aproj = cs_spfree(_Aproj);
  _Aproj = cs_triplet(Atriplet);
  Atriplet  = cs_spfree(Atriplet);
}

/**
 * Calculate the vector of grid index shifts
 * This vector is calculated for the cell located in the center of the grid
 */
void ProjConvolution::_constructShiftVector()
{
  int ndim = _gridPoint->getNDim();
  int center = 1;
  for (int idim = 0; idim < ndim; idim++)
    center *= _gridPoint->getNX(idim);
  center /= 2;

  VectorInt indp(ndim);
  VectorInt indm(ndim);
  _shiftVector.resize(_getConvSize());

  _gridPoint->rankToIndice(center, indp);
  for (int idim = 0; idim < ndim; idim++) indm[idim] = indp[idim];

  // Shift the index of last coordinate by the shift of the grid
  indp[ndim - 1] += _getHalfSize();

  for (int i = -_getHalfSize(); i <= _getHalfSize(); i++)
  {
    indm[ndim - 1] = indp[ndim - 1] + i;
    int id = _gridPoint->indiceToRank(indm);
    _shiftVector[i + _getHalfSize()] = id - center;
  }
}

int ProjConvolution::point2mesh(const VectorDouble &valonseismic,
                                VectorDouble &valonvertex) const
{

  if ((int) valonvertex.size() != getApexNumber())
  {
    messerr("Dimension of 'valonvertex'(%d) incorrect. If should be %d",
            (int) valonvertex.size(), getApexNumber());
    return 1;
  }
  if ((int) valonseismic.size() != getPointNumber())
  {
    messerr("Dimension of 'valonseismic'(%d) incorrect. If should be %d",
            (int) valonseismic.size(), getPointNumber());
    return 1;
  }

  for (auto &e : valonvertex)
     e = 0.;

  for (int is = 0; is < (int) valonseismic.size(); is++)
  {
    for (int i = -_getHalfSize(); i <= _getHalfSize(); i++)
    {
      int j = i + _getHalfSize();
      int id = is + _shiftVector[j];
      if (id < 0) return 1;

      double valm1 = valonseismic[is];
      if (FFFF(valm1))
      {
        valonvertex[id] = TEST;
        break;
      }
      double valm = valm1 * _convolution[j];
      valonvertex[id] += valm;
    }
  }

  // Comparing with the Aproj method (if initiated)

  if (_Aproj != nullptr)
  {
    VectorDouble valcheck = valonvertex;
    cs_tmulvec(_Aproj,(int) valcheck.size(),valonseismic.data(),valcheck.data());
    valcheck.subtract(valonvertex);
    message("Point2Mesh: norme de la difference = %lf\n",valcheck.norm());
  }

  return 0;
}

int ProjConvolution::mesh2point(const VectorDouble &valonvertex,
                                VectorDouble &valonseismic) const
{
  if ((int) valonvertex.size() != getApexNumber())
  {
    messerr("Dimension of 'valonvertex'(%d) incorrect. If should be %d",
            (int) valonvertex.size(), getApexNumber());
    return 1;
  }
  if ((int) valonseismic.size() != getPointNumber())
  {
    messerr("Dimension of 'valonseismic'(%d) incorrect. If should be %d",
            (int) valonseismic.size(), getPointNumber());
    return 1;
  }

  for (int is = 0; is < (int) valonseismic.size(); is++)
  {
    double valp = 0;
    for (int i = -_getHalfSize(); i <= _getHalfSize(); i++)
    {
      int j = i + _getHalfSize();
      int id = is + _shiftVector[j];
      if (id < 0) return 1;

      double valm1 = valonvertex[id];
      if( FFFF(valm1))
      {
        valp = TEST;
        break;
      }
      double valm = valm1  * _convolution[j];
      valp += valm;
    }
    valonseismic[is] = valp;
  }

  // Comparing with the Aproj method (if initiated)

  if (_Aproj != nullptr)
  {
    VectorDouble valcheck = valonseismic;
    cs_mulvec(_Aproj,(int) valcheck.size(),valonvertex.data(),valcheck.data());
    valcheck.subtract(valonseismic);
    message("Mesh2point: norme de la difference = %lf\n",valcheck.norm());
  }
  return 0;
}

int ProjConvolution::getApexNumber() const
{
  VectorInt nxs = _gridPoint->getNXs();
  nxs[_gridPoint->getNDim() - 1] += (_getConvSize() - 1);
  return ut_vector_prod(nxs);
}

int ProjConvolution::getPointNumber() const
{
  VectorInt nxs = _gridPoint->getNXs();
  return ut_vector_prod(nxs);
}
