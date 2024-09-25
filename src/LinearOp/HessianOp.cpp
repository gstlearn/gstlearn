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
#include "LinearOp/HessianOp.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/AException.hpp"
#include "Basic/Law.hpp"

#include "LinearOp/ALinearOp.hpp"
#include "LinearOp/PrecisionOp.hpp"

#include <math.h>

HessianOp::HessianOp()
  : _isInitialized(false)
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

/**
 * @brief Return the size of the operator
 * 
 * @return int 
 */
int HessianOp::getSize() const
{
  return _pMat->getSize();
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
    _lambda.resize(nvertex);
    // Set the initialization flag
    _isInitialized = true;
  }
  catch(const std::string& str)
  {
    // TODO : Check if std::exception can be used
    error = 1;
    messerr("%s", str.c_str());
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
int HessianOp::_addToDest(const constvect inv, vect outv) const
{
  if (!_isInitialized) my_throw("'HessianOp' must be initialized beforehand");
  
  // Map Eigen Vector to VectorDouble arguments
  // TODO : VectorXd => VectorDouble = Memory copy !!
 // VectorDouble einv(inv.data(), inv.data() + inv.size());
 // VectorDouble eoutv(outv.size());

  // Contribution of the spatial structure

  _pMat->addToDest(inv,outv);

  // Contribution of the Data
  constvect lambdas(_lambda);
  vect wps(_workp);
  vect wxs(_workx);
  _projData->mesh2point(lambdas,wps);
  _projData->mesh2point(inv,wxs);

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
    _workp[i] = _workx[i] * (- _workp[i] * ratio + pow(ratio,2) ) ;
  }
  vect wvs(_workv);
  _projData->point2mesh(wps, wvs);
  for (int i=0; i<_projData->getApexNumber(); i++) 
  {
    outv[i]+= _workv[i];
  }
  // Contribution of Seismic (optional)

  if (_flagSeismic)
  {
    for (int i=0; i<_projSeis->getApexNumber(); i++) 
      _workv[i] = law_cdf_gaussian(_lambda[i]);
    vect wvs(_workv);
    vect wss(_works);
    _projSeis->mesh2point(wvs, wss);
    for (int i=0; i<_projSeis->getPointNumber(); i++) 
    {
      _works[i]-= _propSeis[i];
      _works[i]*= _varSeis[i];
    }
    _projSeis->point2mesh(wss, wvs);

    for (int i=0; i<_projData->getApexNumber(); i++) 
    { 
      double val = _lambda[i];
      outv[i] -= val * law_df_gaussian(val) * _workv[i] * inv[i];
    }
    for (int i=0; i<_projSeis->getApexNumber(); i++)
      _workv[i] = inv[i] * law_df_gaussian(_lambda[i]);
    _projSeis->mesh2point(wvs, wss);
    for (int i=0; i<_projSeis->getPointNumber(); i++)
      _works[i] *= _varSeis[i];
    _projSeis->point2mesh(wss, wvs);
    for (int i=0; i<_projSeis->getApexNumber(); i++)
      _workv[i] *= law_df_gaussian(_lambda[i]);

    for (int i=0; i<_projData->getApexNumber(); i++) 
      outv[i] += _workv[i];
  }
  return 0;
}
