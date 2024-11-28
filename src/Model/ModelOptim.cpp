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
#include "Variogram/Vario.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Model/Model.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Model/Constraints.hpp"

#include <cmath>
#include <nlopt.h>

#define IJDIR(ijvar, ipadir) ((ijvar)*npadir + (ipadir))
#define WT(ijvar, ipadir)    wt[IJDIR(ijvar, ipadir)]

ModelOptim::ModelOptim(Model* model,
                       Constraints* constraints,
                       const Option_AutoFit& mauto,
                       const Option_VarioFit& optvar)
  : _modelPart()
  , _constraints(constraints)
  , _mauto(mauto)
  , _optvar(optvar)
{
  _modelPart._model = model;
}

ModelOptim::ModelOptim(const ModelOptim& m)
  : _modelPart()
  , _constraints(m._constraints)
  , _mauto(m._mauto)
  , _optvar(m._optvar)
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
    _constraints = m._constraints;
    _mauto       = m._mauto;
    _optvar      = m._optvar;
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

  _modelPart._calcmode.setAsVario(true);
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
    double scale          = param._scale;
    CovAniso* cova        = modelPart._model->getCova(icov);

    if (param._type == EConsElem::RANGE)
    {
      // Range
      cova->setRangeIsotropic(scale * current[iparam++]);
    }

    if (param._type == EConsElem::SILL)
    {
      // Sill (through the AIC matrix)
      MatrixSquareGeneral aic(nvar);
      int ijvar = 0;
      for (int ivar = ijvar = 0; ivar < nvar; ivar++)
        for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
          aic.setValue(ivar, jvar, scale * current[iparam++]);
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
      cova->setParam(scale * current[iparam++]);
    }
  }

  // Add 1 to the internal counter
  modelPart._niter++;
}

void ModelOptim::updateModelParamList(double distmax_def,
                                      const MatrixSquareSymmetric& vars_def)
{
  double value = TEST;
  double scale = 1.;
  int nparams  = _getParamNumber();
  Model* model = _modelPart._model;
  int ncov     = model->getCovaNumber(true);

  // Cholesky decomposition of the matrix of variances
  VectorDouble varchol;
  if (!vars_def.empty())
  {
    CholeskyDense cholesky(&vars_def);
    varchol = cholesky.getLowerTriangle();
  }

    for (int iparam = 0; iparam < nparams; iparam++)
  {
    OneParam& param      = _modelPart._params[iparam];
    int icov             = param._icov;
    const CovAniso* cova = model->getCova(icov);

    value = 1.;
    scale = 1.;
    switch (param._type.toEnum())
    {
      case EConsElem::E_SILL:
      {
        if (!vars_def.empty())
        {
          value = varchol[param._rank] / sqrt(ncov);
          scale = ABS(value);
        }
        break;
      }

      case EConsElem::E_RANGE:
      {
        if (!FFFF(distmax_def))
        {
          double dunit = distmax_def / ncov / 2.;
          value        = dunit * (icov + 1);
          scale        = dunit;
        }
        break;
      }

      case EConsElem::E_PARAM:
      {
        if (!FFFF(distmax_def))
        {
          if (cova->getType() == ECov::COSEXP) value = distmax_def / 3.;
        }
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

    _modelPart._tabval[iparam] = value / scale;
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
  int nparams = (int)modelPart._params.size();

  message("Iteration %3d - %s (", modelPart._niter, title.c_str());
  for (int iparam = 0; iparam < nparams; iparam++)
    message("%lf ", modelPart._params[iparam]._scale * modelPart._tabval[iparam]);
  message(") = %lf\n", result);
}

void ModelOptim::_setSill(int icov, int ivar, int jvar, double value) const
{
  _modelPart._model->setSill(icov, ivar, jvar, value);
}

void ModelOptim::_performOptimization(nlopt_func f,
                                      void* f_data,
                                      double distmax_def,
                                      const MatrixSquareSymmetric& vars_def)
{
  // Define the optimization criterion
  int npar      = _getParamNumber();
  nlopt_opt opt = nlopt_create(NLOPT_LN_NELDERMEAD, npar);
  nlopt_set_lower_bounds(opt, _modelPart._tablow.data());
  nlopt_set_upper_bounds(opt, _modelPart._tabupp.data());
  nlopt_srand(12345);
  nlopt_set_ftol_rel(opt, 1e-6);

  // Update the initial optimization values (due to variogram)
  updateModelParamList(distmax_def, vars_def);

  // Define the cost function
  nlopt_set_min_objective(opt, f, f_data);

  // Perform the optimization (store the minimized value in 'minf')
  double minf;
  if (_modelPart._verbose) mestitle(1, "Model Fitting from Variogram");
  nlopt_optimize(opt, _modelPart._tabval.data(), &minf);
  nlopt_destroy(opt);
}