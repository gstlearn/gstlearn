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

ProjConvolution::ProjConvolution(const VectorDouble &convolution,
                                 const DbGrid *grid_point)
 : _convolution(convolution)
 , _gridPoint(grid_point)
{
}

ProjConvolution::~ProjConvolution() { }



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

  int ndim = _gridPoint->getNDim();
  VectorInt indp(ndim);
  VectorInt indm(ndim);
  for (int iechp = 0; iechp < (int) valonseismic.size(); iechp++)
  {
    _gridPoint->rankToIndice(iechp, indp);

    for (int idim = 0; idim < ndim; idim++)
      indm[idim] = indp[idim];

    indp[ndim-1] += _getHalfSize();

    for (int i = -_getHalfSize(); i <= _getHalfSize(); i++)
    {
      int j = i + _getHalfSize();

      // First value of the difference
      indm[ndim-1] = indp[ndim-1] + i;
      int rankm = _gridPoint->indiceToRank(indm);
      if (rankm < 0)
      {
        //messerr("Error indexing for target iechp=%d\n",iechp);
        //ut_ivector_display("Indices in convolution (out of grid)", indm);
        continue; // This is in order to avoid crash on next line
      }
      double valm1 = valonseismic[iechp];
      if( FFFF(valm1))
      {
        valonvertex[rankm] = TEST;
        break;
      }
      // Second value of the difference
//      indm[ndim-1] += 1;
//      rankm = _gridPoint->indiceToRank(indm);
//      if (rankm < 0)
//      {
//        //messerr("Error indexing for target iechp=%d\n",iechp);
//        //ut_ivector_display("Indices in convolution (out of grid)", indm);
//        continue; // This is in order to avoid crash on next line
//      }
//      double valm2 = valonvertex[rankm];

      //double valm = (valm1 - valm2) * _convolution[j];
      double valm = valm1  * _convolution[j];
      valonvertex[rankm] = valm;
    }
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

  int ndim = _gridPoint->getNDim();
  VectorInt indp(ndim);
  VectorInt indm(ndim);
  for (int iechp = 0; iechp < (int) valonseismic.size(); iechp++)
  {
    _gridPoint->rankToIndice(iechp, indp);

    for (int idim = 0; idim < ndim; idim++)
      indm[idim] = indp[idim];

    indp[ndim-1] += _getHalfSize();

    double valp = 0;
    for (int i = -_getHalfSize(); i <= _getHalfSize(); i++)
    {
      int j = i + _getHalfSize();

      // First value of the difference
      indm[ndim-1] = indp[ndim-1] + i;
      int rankm = _gridPoint->indiceToRank(indm);
      if (rankm < 0)
      {
        //messerr("Error indexing for target iechp=%d\n",iechp);
        //ut_ivector_display("Indices in convolution (out of grid)", indm);
        continue; // This is in order to avoid crash on next line
      }
      double valm1 = valonvertex[rankm];
      if( FFFF(valm1))
      {
        valp = TEST;
        break;
      }
      // Second value of the difference
//      indm[ndim-1] += 1;
//      rankm = _gridPoint->indiceToRank(indm);
//      if (rankm < 0)
//      {
//        //messerr("Error indexing for target iechp=%d\n",iechp);
//        //ut_ivector_display("Indices in convolution (out of grid)", indm);
//        continue; // This is in order to avoid crash on next line
//      }
//      double valm2 = valonvertex[rankm];

      //double valm = (valm1 - valm2) * _convolution[j];
      double valm = valm1  * _convolution[j];
      valp += valm;
    }
    valonseismic[iechp] = valp;
  }
  return 0;
}

int ProjConvolution::getApexNumber() const
{
  VectorInt nxs = _gridPoint->getNXs();
  nxs[_gridPoint->getNDim() - 1]+= _getConvSize();
  return ut_vector_prod(nxs);
}

int ProjConvolution::getPointNumber() const
{
  VectorInt nxs = _gridPoint->getNXs();
  return ut_vector_prod(nxs);
}
