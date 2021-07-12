/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 18 dec. 2020 by N. Desassis                                    */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "LinearOp/PrecisionOpMultiConditional.hpp"
#include "geoslib_e.h"

PrecisionOpMultiConditional::PrecisionOpMultiConditional()
  :_multiPrecisionOp(std::vector<PrecisionOp*>())
  ,_multiProjData(std::vector<IProjMatrix*>())
  ,_nugget(0.)
  ,_ndat(0)
  ,_work1(VectorDouble())
  ,_work2(VectorVectorDouble())
{

}


void PrecisionOpMultiConditional::push_back(PrecisionOp* pmatElem,
                                            IProjMatrix* projDataElem)
{
  _multiPrecisionOp.push_back(pmatElem);
  _work2.push_back(VectorDouble(pmatElem->getSize()));
  if (_multiProjData.size() == 0)
  {
    _ndat = projDataElem->getPointNumber(); //TODO Vérifier la cohérence. _ndat doit coïncider pour tous les projDataElem.
    _work1.resize(_ndat);
  }
  _multiProjData.push_back(projDataElem);
  _updated();
}

/*****************************************************************************/
/*!
** Compute diag(Q1,...,Qncova) x + 1/nugget [A1,...,Ancova]^t [A1,...,Ancova] x
** in a block form where ncova is the number of basic structures excluding the
** nugget effect. Qi are the precision matrices associated to each structure and Ai
** are the projection matrices from the meshing vertices to the data locations.
**
** \param[in]  in     Array of input values
**
** \param[out] out    Array of output values
**
*******************************************************************************/
void PrecisionOpMultiConditional::_evalDirect(const VectorVectorDouble& in,
                                              VectorVectorDouble& out) const
{
  _init();
  double invNugg = 1. / _nugget;
  for (int imod = 0; imod < size(); imod++)
  {
    _multiPrecisionOp[imod]->eval(in[imod], out[imod]);

    for (int jmod = 0; jmod < size(); jmod++)
    {
      _multiProjData[jmod]->mesh2point(in[jmod], _work1);
      for (int idat = 0; idat < _ndat; idat++)
      {
         _work1[idat] *= invNugg;
      }
      _multiProjData[imod]->point2mesh(_work1, _work2[imod]);
      _linearComb(1., _work2, 1., out, out);

    }
  }
}

PrecisionOpMultiConditional::~PrecisionOpMultiConditional(){}
