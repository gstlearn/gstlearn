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
#include "Basic/AStringable.hpp"
#include "csparse_f.h"

ProjConvolution::ProjConvolution(const VectorDouble &convolution,
                                 const DbGrid *grid_point,
                                 const VectorInt& nmult,
                                 bool useAProj)
    : _convolution(convolution),
      _gridPoint(grid_point),
      _nmult(nmult),
      _shiftVector(),
      _weights(),
      _Aproj(nullptr)
{
  int ndim = grid_point->getNDim();
  _nmult.resize(ndim, 1);
  _nmult[ndim-1] = 1;

  _constructShiftVector();

  _constructWeights();

  if (useAProj) _constructAprojCS();
}

ProjConvolution::~ProjConvolution()
{
  if (_Aproj != nullptr) _Aproj = cs_spfree(_Aproj);
}

void ProjConvolution::_constructWeights()
{
  int ndim = _getNDim();

  if (ndim == 2)
  {
    _weights = MatrixRectangular(_nmult[0],1);
    for (int ix = 0; ix < _nmult[0]; ix++)
    {
      double w1 = (double) ix / (double) _nmult[0];
      _weights.setValue(ix, 0, w1);
    }
  }
  else if (ndim == 3)
  {
    _weights = MatrixRectangular(_nmult[0], _nmult[1]);
    for (int ix = 0; ix < _nmult[0]; ix++)
    {
      double w1 = (double) ix / (double) _nmult[0];
      for (int iy = 0; iy < _nmult[1]; iy++)
      {
        double w2 = (double) ix / (double) _nmult[0];
        _weights.setValue(ix, iy, w1 * w2);
      }
    }
  }
  else
  {
    messerr("Not Programmed");
  }
}


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
  _Aproj = cs_triplet(Atriplet);
  Atriplet  = cs_spfree(Atriplet);
  return 0;
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

bool ProjConvolution::_isVecDimCorrect(const VectorDouble &valonseismic,
                                       const VectorDouble &valonvertex) const
{
  if ((int) valonvertex.size() != getApexNumber())
  {
    messerr("Dimension of 'valonvertex'(%d) incorrect. If should be %d",
            (int) valonvertex.size(), getApexNumber());
    return false;
  }
  if ((int) valonseismic.size() != getPointNumber())
  {
    messerr("Dimension of 'valonseismic'(%d) incorrect. If should be %d",
            (int) valonseismic.size(), getPointNumber());
    return false;
  }
  return true;
}

int ProjConvolution::point2mesh(const VectorDouble &valonseismic,
                                VectorDouble &valonvertex) const
{
  if (! _isVecDimCorrect(valonseismic, valonvertex)) return 1;

  for (auto &e : valonvertex)
     e = 0.;

  int count = (int) valonseismic.size();
  int size  = _getConvSize();
  double valm = 0.;
  int id = 0;
  for (int is = 0; is < count; is++)
  {
    for (int j = 0; j < size; j++)
    {
      id = is + _shiftVector[j];
      valm = valonseismic[is];
      if (FFFF(valm))
      {
        valonvertex[id] = TEST;
        break;
      }
      valonvertex[id] += valm * _convolution[j];
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
  if (! _isVecDimCorrect(valonseismic, valonvertex)) return 1;

  int count = (int) valonseismic.size();
  int size  = _getConvSize();
  double valp  = 0.;
  double valm = 0.;
  int id = 0;
  for (int is = 0; is < count; is++)
  {
    valp = 0;
    for (int j = 0; j < size; j++)
    {
      id = is + _shiftVector[j];
      if (id < 0) return 1;

      valm = valonvertex[id];
      if( FFFF(valm))
      {
        valp = TEST;
        break;
      }
      valp += valm * _convolution[j];
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

DbGrid* ProjConvolution::getResolutionGrid() const
{
  int ndim = _gridPoint->getNDim();

  VectorInt nxs(ndim);
  VectorDouble dx(ndim);
  VectorDouble x0(ndim);
  _gridPoint->getGrid().divider(_nmult, 1, nxs, dx, x0);

  // Create the new grid
  DbGrid* dbgrid = DbGrid::create(nxs, dx, x0, _gridPoint->getAngles());
  return dbgrid;
}

VectorInt ProjConvolution::_getNXResolutionGrid() const
{
  int ndim = _gridPoint->getNDim();

  VectorInt nxs(ndim);
  VectorDouble dx(ndim);
  VectorDouble x0(ndim);
  _gridPoint->getGrid().divider(_nmult, 1, nxs, dx, x0);

  // Correct the last dimension
  nxs[_gridPoint->getNDim() - 1] += (_getConvSize() - 1);
  return nxs;
}

int ProjConvolution::getApexNumber() const
{
  VectorInt nxs = _getNXResolutionGrid();
  return ut_vector_prod(nxs);
}

int ProjConvolution::getPointNumber() const
{
  VectorInt nxs = _gridPoint->getNXs();
  return ut_vector_prod(nxs);
}
