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
#include "Model/ModelOptim.hpp"

#include "Enum/EConsElem.hpp"
#include "geoslib_define.h"

#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Model/Model.hpp"
#include "Variogram/Vario.hpp"

#include <nlopt.h>

#define IJDIR(ijvar, ipadir) ((ijvar)*npadir + (ipadir))
#define WT(ijvar, ipadir)    wt[IJDIR(ijvar, ipadir)]

ModelOptim::ModelOptim()
  : _modelPart()
{
}

ModelOptim::ModelOptim(const ModelOptim& m)
  : _modelPart()
{
   _copyModelPart(m._modelPart);
}

void ModelOptim::_copyModelPart(const Model_Part& modelPart)
{
  _modelPart._model   = modelPart._model;
  _modelPart._params  = modelPart._params;
  _modelPart._tabval  = modelPart._tabval;
  _modelPart._tablow  = modelPart._tablow;
  _modelPart._tabupp  = modelPart._tabupp;
  _modelPart._verbose = modelPart._verbose;
}

ModelOptim& ModelOptim::operator=(const ModelOptim& m)
{
  if (this != &m)
  {
    _copyModelPart(m._modelPart);
  }
  return (*this);
}

ModelOptim::~ModelOptim()
{
}

int ModelOptim::_buildModelParamList()
{
  if (_modelPart._model == nullptr)
  {
    messerr("Argument '_model' must be defined beforehand");
    return 1;
  }

  // Clear previous contents
  _modelPart._tabval.clear();
  _modelPart._tablow.clear();
  _modelPart._tabupp.clear();

  // Loop on the covariances
  const Model* model = _modelPart._model;
  int nvar = model->getVariableNumber();
  for (int icov = 0, ncov = model->getCovaNumber(); icov < ncov; icov++)
  {
    const CovAniso* cova = model->getCova(icov);

    if (cova->hasRange())
    {
      // Add the 'Range' (scalar) attribute
      _addOneModelParam(icov, EConsElem::RANGE, 0, EPSILON2, TEST);
    }
    if (cova->hasParam())
    {
      // Add the 'Param' (scalar) attribute
      _addOneModelParam(icov, EConsElem::PARAM, 0, 0., 2.);
    }

    // Add the 'Sill' (vectorial)attribute
    int ijvar = 0;
    for (int ivar = ijvar = 0; ivar < nvar; ivar++)
      for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
        _addOneModelParam(icov, EConsElem::SILL, ijvar, TEST, TEST);
  }
  return 0;
}

void ModelOptim::_addOneModelParam(int icov,
                                   const EConsElem& type,
                                   int rank,
                                   double lbound,
                                   double ubound)
{
  OneParam param;
  param._icov = icov;
  param._type = type;
  param._rank = rank;
  _modelPart._params.push_back(param);
  _modelPart._tabval.push_back(1.);
  if (FFFF(lbound)) lbound = -1.e30;
  _modelPart._tablow.push_back(lbound);
  if (FFFF(ubound)) ubound = 1.e30;
  _modelPart._tabupp.push_back(ubound);
}

void ModelOptim::_patchModel(Model_Part& modelPart, const double* current)
{
  // Initializations
  int nvar    = modelPart._model->getVariableNumber();
  int nvs2    = nvar * (nvar + 1) / 2;
  int nparams = (int)modelPart._params.size();

  // Loop on the parameters
  int iparam = 0;
  while (iparam < nparams)
  {
    const OneParam& param = modelPart._params[iparam];
    int icov              = param._icov;
    CovAniso* cova        = modelPart._model->getCova(icov);

    if (param._type == EConsElem::RANGE)
    {
      // Range
      cova->setRangeIsotropic(current[iparam++]);
    }

    if (param._type == EConsElem::SILL)
    {
      // Sill (through the AIC matrix)
      MatrixSquareGeneral aic(nvar);
      int ijvar = 0;
      for (int ivar = ijvar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          aic.setValue(ivar, jvar, current[iparam++]);
      if (ijvar == nvs2)
      {
        MatrixSquareSymmetric sills(nvar);
        sills.prodNormMatVecInPlace(aic, VectorDouble());
        cova->setSill(sills);
      }
    }

    if (param._type == EConsElem::PARAM)
    {
      // Set the 'param' attribute
      cova->setParam(current[iparam++]);
    }
  }
}
