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
#include "LinearOp/CholeskyDense.hpp"

#include <cmath>
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
    // TODO: inference of PARAM is deactivated (DR on 2024/11/14)
    // if (cova->hasParam())
    // {
    //   // Add the 'Param' (scalar) attribute
    //   _addOneModelParam(icov, EConsElem::PARAM, 0, 0., 2.);
    // }

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
  param._scale = 1.;
  _modelPart._params.push_back(param);
  _modelPart._tabval.push_back(1.);
  if (FFFF(lbound)) lbound = -HUGE_VAL;
  _modelPart._tablow.push_back(lbound);
  if (FFFF(ubound)) ubound = +HUGE_VAL;
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

  // Add 1 to the internal counter
  modelPart._niter++;
}

void ModelOptim::updateModelParamList(double hmax,
                                      const MatrixSquareSymmetric& vars)
{
  double value = TEST;
  double scale = 1.;
  int nparams  = _getParamNumber();
  Model* model = _modelPart._model;
  int ncov     = model->getCovaNumber(true);

  // Cholesky decomposition of the matrix of variances
  CholeskyDense cholesky(&vars);
  VectorDouble varchol = cholesky.getLowerTriangle();

    for (int iparam = 0; iparam < nparams; iparam++)
  {
    OneParam& param      = _modelPart._params[iparam];
    int icov             = param._icov;
    const CovAniso* cova = model->getCova(icov);

    switch (param._type.toEnum())
    {
      case EConsElem::E_SILL:
      {
        value = varchol[param._rank] / sqrt(ncov);
        scale = ABS(value);
        break;
      }

      case EConsElem::E_RANGE:
      {
        double dunit = hmax / ncov / 2.;
        value        = dunit * (icov + 1);
        scale        = dunit;

        break;
      }

      case EConsElem::E_PARAM:
      {
        value = 1.;
        if (cova->getType() == ECov::COSEXP) value = hmax / 3.;
        scale = 1.;
        break;
      }

      case EConsElem::E_ANGLE:
      {
        value = 0.;
        scale = 1800.;
        break;
      }

        default: break;
    }

    _modelPart._tabval[iparam] = value;
    param._scale = scale;
  }
}

void ModelOptim::dumpParamList() const
{
  mestitle(1, "Model Optimization parameters");
  int nparams = _getParamNumber();
  for (int iparam = 0; iparam < nparams; iparam++)
  {
    _dumpOneModelParam(_modelPart._params[iparam], _modelPart._tabval[iparam]);
  }
}

void ModelOptim::_dumpOneModelParam(const OneParam& param, double value)
{
  message("Covariance %d - Type = %s - Rank = %d - Scale = %lf - Current = %lf\n",
    param._icov, param._type.getDescr().data(), param._rank, param._scale, value);
}

void ModelOptim::_printResult(const String& title,
                              const Model_Part& modelPart,
                              double result)
{
  if (!modelPart._verbose) return;
  message("Iteration %3d - %s = %lf\n", modelPart._niter, title.c_str(), result);
}