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
#include <LinearOp/HessianOp.hpp>
#include "geoslib_e.h"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"

HessianOp::HessianOp() 
  : ALinearOp()
  , _isInitialized(false)
  , _flagSeismic(false)
  , _pMat(nullptr)
  , _projData(nullptr)
  , _projSeis(nullptr)
  , _indic()
  , _propSeis()
  , _varSeis()
  , _lambda()
  , _workp()
  , _workx()
  , _workv()
  , _works()
{
}

HessianOp::~HessianOp() 
{
}

/*****************************************************************************/
/*!
**  Initialize the Hessian Operator
**
** \param[in]  pmat     The precision matrix to be optimized
** \param[in]  projdata The Projection operator between Data and Meshing
** \param[in]  projseis The Projection operator between Seismic and Meshing
** \param[in]  indic    Array of facies values
** \param[in]  propseis Array of facies proportions
** \param[in]  varseis  Array of variance attached to the seismic
**
*****************************************************************************/
int HessianOp::init(const PrecisionOp*  pmat,
                    const ProjMatrix*   projdata,
                    const ProjMatrix*   projseis,
                    const VectorDouble& indic,
                    const VectorDouble& propseis,
                    const VectorDouble& varseis)
{
  // Initialization

  int error = 0;

  try
  {
    _pMat     = pmat;
    _projData = projdata;
    _projSeis = projseis;
    _indic    = indic;
    _propSeis = propseis;
    _varSeis  = varseis;
    
    int nvertex = _projData->getApexNumber();
    int npoint  = _projData->getPointNumber();
    
    // Particular case of the Seismic
    _flagSeismic = (projseis != (ProjMatrix *) NULL && 
                    projseis->getPointNumber() > 0);
    if (_flagSeismic)
    {
      int nseis = _projSeis->getPointNumber();
      _works.resize(nseis);
    }
    
    // Auxiliary working arrays
    _workp.resize(npoint);
    _workx.resize(npoint);
    _workv.resize(nvertex);
    
    // Set the initialization flag
    _isInitialized = true;
  }

  catch(const char * str)
  {
    error = 1;
    std::cout << str << std::endl;
  }
  return error;
}

/*****************************************************************************/
/*!
**  Operate the operation: 'out' = HESS * 'in'
**
** \param[in]  in       Array of input values
**
** \param[out] out      Array of output values
**
*****************************************************************************/
void HessianOp::_evalDirect(const VectorDouble& in,
                            VectorDouble& out) const
{
  if (! _isInitialized)
    my_throw("'HessianOp' must be initialized beforehand");

  // Contribution of the spatial structure

  _pMat->getShiftOp()->evalDirect(in,out);

  // Contribution of the Data

  _projData->mesh2point(_lambda,_workp);
  _projData->mesh2point(in,_workx);

  double denom,dl;
  for (int i=0; i<_projData->getPointNumber(); i++)
  {
    double ratio = 0.;
    if (! FFFF(_indic[i]))
    {
      denom = _indic[i] - law_cdf_gaussian(_workp[i]);
      dl    = law_df_gaussian(_workp[i]);
      ratio = dl / denom;
    }
    _workp[i] = (- _workp[i] * ratio + pow(ratio,2)) * _workx[i];
  }
  _projData->point2mesh(_workp, _workv);
  for (int i=0; i<_projData->getApexNumber(); i++) out[i] += _workv[i];

  // Contribution of Seismic (optional)

  if (_flagSeismic)
  {
    for (int i=0; i<_projSeis->getApexNumber(); i++) 
      _workv[i] = law_cdf_gaussian(_lambda[i]);
    _projSeis->mesh2point(_workv, _works);
    for (int i=0; i<_projSeis->getPointNumber(); i++) 
    {
      _works[i] -= _propSeis[i]; 
      _works[i] *= _varSeis[i];
    }
    _projSeis->point2mesh(_works, _workv);

    for (int i=0; i<_projData->getApexNumber(); i++) 
      out[i] -= _lambda[i] * law_df_gaussian(_lambda[i]) * _workv[i] * in[i];

    for (int i=0; i<_projSeis->getApexNumber(); i++)
      _workv[i] = in[i] * law_df_gaussian(_lambda[i]);
    _projSeis->mesh2point(_workv, _works);
    for (int i=0; i<_projSeis->getPointNumber(); i++)
      _works[i] *= _varSeis[i];
    _projSeis->point2mesh(_works, _workv);
    for (int i=0; i<_projSeis->getApexNumber(); i++)
      _workv[i] *= law_df_gaussian(_lambda[i]); // d2u

    for (int i=0; i<_projData->getApexNumber(); i++) 
      out[i] += _workv[i];
  }
}
