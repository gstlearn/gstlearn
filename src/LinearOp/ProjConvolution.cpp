/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 9 avr. 2019 by N. Desassis                                     */
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
      _gridSeismic(grid_point),
      _nmult(nmult),
      _shiftVector(),
      _weightx(),
      _weighty(),
      _Aproj(nullptr)
{
  int ndim = grid_point->getNDim();
  if (ndim != 2 && ndim != 3)
  {
    messerr("ProjConvolution is limited to 2-D or 3-D case");
    _nmult.clear();
    return;
  }
  _nmult.resize(ndim, 1);
  _nmult[ndim-1] = 1;

  _buildShiftVector();

  _buildWeights();

  if (useAProj) _buildAprojCS();
}

ProjConvolution::~ProjConvolution()
{
  if (_Aproj != nullptr) _Aproj = cs_spfree(_Aproj);
}

void ProjConvolution::_buildWeights()
{
  _weightx.resize(_nmult[0]);
  for (int ix = 0; ix < _nmult[0]; ix++)
    _weightx[ix] = 1. - (double) ix / (double) _nmult[0];

  if (_getNDim() <= 2) return;

  _weighty.resize(_nmult[1]);
  for (int iy = 0; iy < _nmult[1]; iy++)
    _weighty[iy] = 1. - (double) iy / (double) _nmult[1];
}


/**
 * Calculate the Aproj sparse matrix.
 * This method is kept for establishing time bench marks.
 * It emulates mesh2point algorithm.
 * Note that this algorithm does not handle the presence of undefined values
 */
int ProjConvolution::_buildAprojCS()
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
 * Calculate the vector of grid index shifts (in Point Grid)
 * This vector is calculated for the cell located in the center of the grid
 */
void ProjConvolution::_buildShiftVector()
{
  // Creating the characteristics of the Point Grid

  Grid grid = _getResolutionGridCharacteristics();

  int ndim = _gridSeismic->getNDim();
  int center = 1;
  for (int idim = 0; idim < ndim; idim++)
    center *= grid.getNX(idim);
  center /= 2;

  VectorInt indp(ndim);
  VectorInt indm(ndim);
  _shiftVector.resize(_getConvSize());

  grid.rankToIndice(center, indp);
  for (int idim = 0; idim < ndim; idim++) indm[idim] = indp[idim];

  // Shift the index of last coordinate by the shift of the grid
  indp[ndim - 1] += _getHalfSize();

  for (int i = -_getHalfSize(); i <= _getHalfSize(); i++)
  {
    indm[ndim - 1] = indp[ndim - 1] + i;
    int id = grid.indiceToRank(indm);
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
  if (_shiftVector.empty())
  {
    messerr("The ProjConvolution object has not been built correctly");
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

  if (_getNMultProd() == 1 && true)  // DR: modif pour provoquer le passage dans nouveau code
  {
    if (_mesh2pointRef(valonvertex, valonseismic)) return 1;
  }
  else if (_getNDim() == 2)
  {
    if (_mesh2point2D(valonvertex, valonseismic)) return 1;
  }
  else
  {
    if (_mesh2point3D(valonvertex, valonseismic)) return 1;
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

int ProjConvolution::_mesh2point3D(const VectorDouble &valonvertex,
                                   VectorDouble &valonseismic) const
{
  int size  = _getConvSize();
  VectorInt    indp(3);
  VectorInt    ranks(4);
  VectorDouble wgt(4);

  // Get the characteristics of the coarse grid

  Grid grid = _getResolutionGridCharacteristics();

  // Loop on the nodes of the seismic grid

  int ecrs    = 0;
  int lecp    = 0;
  int decale  = 0;
  double valp = 0.;
  double valm = 0.;
  double wloc = 0.;
  double vloc = 0.;
  double wx   = 0.;
  double wy   = 0.;

  for (int iz = 0; iz < _gridSeismic->getNX(2); iz++)
  {
    indp[2] = iz;
    for (int iy = 0; iy < grid.getNX(1); iy++)
    {
      indp[1] = iy;
      int nymax = (iy < grid.getNX(1) - 1) ? _nmult[1] : 1;
      for (int iym = 0; iym < nymax; iym++)
      {
        wy = _weighty[iym];
        for (int ix = 0; ix < grid.getNX(0); ix++)
        {
          indp[0] = ix;
          int nxmax = (ix < grid.getNX(0) - 1) ? _nmult[0] : 1;
          for (int ixm = 0; ixm < nxmax; ixm++)
          {
            wx = _weightx[ixm];

            lecp = grid.indiceToRank(indp);
            ranks[0] = lecp;
            ranks[1] = lecp + 1;
            ranks[2] = lecp + 1 + grid.getNX(0);
            ranks[3] = lecp +     grid.getNX(0);

            wgt[0] = wx        * wy;
            wgt[1] = (1. - wx) * wy;
            wgt[2] = (1. - wx) * (1. - wy);
            wgt[3] = wx        * (1. - wy);

            // Loop on the convolution
            valp = 0;
            for (int j = 0; j < size; j++)
            {

              // Derive the value of the 2-D center of gravity
              valm = 0.;
              decale = _shiftVector[j];
              for (int i = 0; i < 4; i++)
              {
                wloc = wgt[i];
                if (wloc > 0.)
                {
                  vloc = valonvertex[ranks[i] + decale];
                  if (FFFF(vloc))
                  {
                    valp = TEST;
                    break;
                  }
                  valm += vloc * wloc;
                }
              }

              // Add to the convolution
              valp += valm * _convolution[j];
            }

            VectorInt indloc(3);
            indloc[0] = ix * _nmult[0] + ixm;
            indloc[1] = iy * _nmult[1] + iym;
            indloc[2] = iz;
            int lecs = _gridSeismic->getGrid().indiceToRank(indloc);
//            if (lecs != ecrs)
//              message("Erreur lecs=%d ecrs=%d\n",lecs,ecrs);
//            valonseismic[ecrs++] = valp;
            valonseismic[lecs] = valp;
          }
        }
      }
    }
  }
  return 0;
}

int ProjConvolution::_mesh2point2D(const VectorDouble &valonvertex,
                                   VectorDouble &valonseismic) const
{
  int size  = _getConvSize();
  VectorInt    indp(2);
  VectorInt    ranks(2);
  VectorDouble wgt(2);

  // Get the characteristics of the coarse grid

  Grid grid = _getResolutionGridCharacteristics();

  // Loop on the nodes of the seismic grid

  int ecrs    = 0;
  int lecp    = 0;
  int decale  = 0;
  double valp = 0.;
  double valm = 0.;
  double wloc = 0.;
  double vloc = 0.;
  double wx   = 0.;

  for (int iy = 0; iy < _gridSeismic->getNX(1); iy++)
  {
    indp[1] = iy;
    for (int ix = 0; ix < grid.getNX(0); ix++)
    {
      indp[0] = ix;
      int nxmax = (ix < grid.getNX(0) - 1) ? _nmult[0] : 1;
      for (int ixm = 0; ixm < nxmax; ixm++)
      {
        wx = _weightx[ixm];

        lecp = grid.indiceToRank(indp);
        ranks[0] = lecp;
        ranks[1] = lecp + 1;

        wgt[0] = wx;
        wgt[1] = 1. - wx;

        // Loop on the convolution

        valp = 0;
        for (int j = 0; j < size; j++)
        {

          // Derive the value of the 2-D center of gravity
          decale = _shiftVector[j];
          valm = 0.;
          for (int i = 0; i < 2; i++)
          {
            wloc = wgt[i];
            if (wloc > 0.)
            {
              vloc = valonvertex[ranks[i] + decale];
              if (FFFF(vloc))
              {
                valp = TEST;
                break;
              }
              valm += vloc * wloc;
            }
          }

          // Add to the convolution
          valp += valm * _convolution[j];
        }
        valonseismic[ecrs++] = valp;
      }
    }
  }
  return 0;
}

int ProjConvolution::_mesh2pointRef(const VectorDouble &valonvertex,
                                    VectorDouble &valonseismic) const
{
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
  return 0;
}

Grid ProjConvolution::_getResolutionGridCharacteristics() const
{
  int ndim = _gridSeismic->getNDim();

  VectorInt nxs(ndim);
  VectorDouble dx(ndim);
  VectorDouble x0(ndim);
  _gridSeismic->getGrid().multiple(_nmult, 0, nxs, dx, x0);

  // Correct the last dimension
  nxs[ndim - 1] += (_getConvSize() - 1);
  x0[ndim - 1]  -= (_getConvSize() - 1) * dx[ndim - 1];
  Grid grid(ndim, nxs, x0, dx);

  return grid;
}

DbGrid* ProjConvolution::getResolutionGrid() const
{
  // Get the characteristics of the Point Grid
  Grid grid = _getResolutionGridCharacteristics();

  // Create the new Point grid
  DbGrid* dbgrid = DbGrid::create(grid.getNXs(),
                                  grid.getDXs(),
                                  grid.getX0s(),
                                  _gridSeismic->getAngles());
  return dbgrid;
}

int ProjConvolution::getApexNumber() const
{
  Grid grid = _getResolutionGridCharacteristics();
  return ut_vector_prod(grid.getNXs());
}

int ProjConvolution::getPointNumber() const
{
  VectorInt nxs = _gridSeismic->getNXs();
  return ut_vector_prod(nxs);
}
