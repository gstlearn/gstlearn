/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "LinearOp/HessianOp.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/Law.hpp"

#include <math.h>

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
int HessianOp::init(PrecisionOp*  pmat,
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
    messerr("%s", str);
  }
  return error;
}

/*****************************************************************************/
/*!
**  Operate the operation: 'outv' = HESS * 'inv'
**
** \param[in]  inv       Array of input values
**
** \param[out] outv      Array of output values
**
*****************************************************************************/
void HessianOp::_evalDirect(const VectorDouble& inv,
                            VectorDouble& outv) const
{
  if (! _isInitialized)
    my_throw("'HessianOp' must be initialized beforehand");

  // Contribution of the spatial structure

  _pMat->eval(inv,outv);

  // Contribution of the Data

  _projData->mesh2point(_lambda,_workp);
  _projData->mesh2point(inv,_workx);

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
  for (int i=0; i<_projData->getApexNumber(); i++) outv[i] += _workv[i];

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
      outv[i] -= _lambda[i] * law_df_gaussian(_lambda[i]) * _workv[i] * inv[i];

    for (int i=0; i<_projSeis->getApexNumber(); i++)
      _workv[i] = inv[i] * law_df_gaussian(_lambda[i]);
    _projSeis->mesh2point(_workv, _works);
    for (int i=0; i<_projSeis->getPointNumber(); i++)
      _works[i] *= _varSeis[i];
    _projSeis->point2mesh(_works, _workv);
    for (int i=0; i<_projSeis->getApexNumber(); i++)
      _workv[i] *= law_df_gaussian(_lambda[i]); // d2u

    for (int i=0; i<_projData->getApexNumber(); i++) 
      outv[i] += _workv[i];
  }
}
