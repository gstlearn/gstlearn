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
#include "Model/AModelOptim.hpp"

#include "Basic/AStringable.hpp"
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
#define TAKE_ROT                                                               \
  ((_optvar.getLockSamerot() && first_covrot < 0) || !_optvar.getLockSamerot())
#define DEFINE_THIRD   (cova->hasParam())
#define DEFINE_RANGE   (cova->hasRange() > 0)
#define DEFINE_ANICOEF (cova->hasRange() != 0 && _optvar.getAuthAniso())
#define DEFINE_ANIROT                                                          \
  (cova->hasRange() != 0 && _optvar.getAuthAniso() && _optvar.getAuthRotation())
#define UNDEFINE_ANIROT                                                        \
  (cova->hasRange() == 0 || !_optvar.getAuthAniso() ||                         \
   !_optvar.getAuthRotation())

AModelOptim::AModelOptim(Model* model,
                         Constraints* constraints,
                         const Option_AutoFit& mauto,
                         const Option_VarioFit& optvar)
  : _modelPart()
  , _constraints(constraints)
  , _mauto(mauto)
{
  _modelPart._model = model;
  _modelPart._optvar = optvar;
}

AModelOptim::AModelOptim(const AModelOptim& m)
  : _modelPart()
  , _constraints(m._constraints)
  , _mauto(m._mauto)
{
   _copyModelPart(m._modelPart);
}

void AModelOptim::_copyModelPart(const Model_Part& modelPart)
{
  _modelPart._model   = modelPart._model;
  _modelPart._optvar  = modelPart._optvar;
  _modelPart._params  = modelPart._params;
  _modelPart._tabval  = modelPart._tabval;
  _modelPart._tablow  = modelPart._tablow;
  _modelPart._tabupp  = modelPart._tabupp;
  _modelPart._verbose = modelPart._verbose;
}

AModelOptim& AModelOptim::operator=(const AModelOptim& m)
{
  if (this != &m)
  {
    _copyModelPart(m._modelPart);
    _constraints = m._constraints;
    _mauto       = m._mauto;
  }
  return (*this);
}

AModelOptim::~AModelOptim()
{
}

void AModelOptim::_addOneModelParam(int icov,
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

void AModelOptim::_updateModelParamList(double distmax_def,
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

void AModelOptim::_dumpParamList() const
{
  mestitle(1, "List of the Model parameters to be infered");
  int nparams = _getParamNumber();
  for (int iparam = 0; iparam < nparams; iparam++)
  {
    _dumpOneModelParam(_modelPart._params[iparam], _modelPart._tabval[iparam]);
  }
  message("\n");
}

void AModelOptim::_dumpOneModelParam(const OneParam& param, double value)
{
  message("Covariance %d - %s(%d) - Scale = %lf - Current = %lf\n",
    param._icov, param._type.getDescr().data(), param._rank, param._scale, value);
}

void AModelOptim::_printResult(const String& title,
                               const Model_Part& modelPart,
                               double result)
{
  if (!modelPart._verbose) return;
  int nparams = (int)modelPart._params.size();

  message("%s (", title.c_str());
  for (int iparam = 0; iparam < nparams; iparam++)
    message("%lf ", modelPart._params[iparam]._scale * modelPart._tabval[iparam]);
  message(") : %lf (%d)\n", result, modelPart._niter);
}

void AModelOptim::_setSill(int icov, int ivar, int jvar, double value) const
{
  _modelPart._model->setSill(icov, ivar, jvar, value);
}

void AModelOptim::_performOptimization(double (*optim_func)(unsigned n,
                                                            const double* x,
                                                            double* gradient,
                                                            void* func_data),
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
  nlopt_set_ftol_rel(opt, EPSILON6);

  // Update the initial optimization values (due to variogram)
  _updateModelParamList(distmax_def, vars_def);

  // Define the cost function
  nlopt_set_min_objective(opt, optim_func, f_data);

  // Perform the optimization (store the minimized value in 'minf')
  double minf;
  nlopt_optimize(opt, _modelPart._tabval.data(), &minf);
  nlopt_destroy(opt);
}

int AModelOptim::_buildModelParamList()
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
  const Model* model       = _modelPart._model;
  int nvar                 = model->getVariableNumber();
  int ndim                 = model->getDimensionNumber();
  bool flagRotationDefined = false;

  for (int icov = 0, ncov = model->getCovaNumber(); icov < ncov; icov++)
  {
    const CovAniso* cova = model->getCova(icov);
    bool flagSill        = true;
    bool flagRange       = cova->hasRange() > 0;
    bool flagAniso       = cova->hasRange() != 0 && _modelPart._optvar.getAuthAniso();
    bool flagAnisoRot    = cova->hasRange() != 0 && _modelPart._optvar.getAuthAniso() &&
                        _modelPart._optvar.getAuthRotation();
    bool flagMultiRanges = false;

    // Distance parameters
    if (flagRange)
    {

      // Add the 'Range(s)' attribute
      if (!flagAniso)
      {
        // Add the 'Isotropic Range' (scalar)
        _addOneModelParam(icov, EConsElem::RANGE, 0, EPSILON2, TEST);
      }
      else
      {
        // Add the 'Anisotropic Range' (vectorial)
        if (ndim == 2)
        {
          // For anisotropy in 2-D, one value is sufficient
          _addOneModelParam(icov, EConsElem::RANGE, 0, EPSILON2, TEST);
        }
        else if (ndim == 3)
        {
          // For anisotropy in 3-D, vectorial definition is needed
          if (_modelPart._optvar.getLockIso2d())
          {
            // If Isotropy is forced in 2-D, back to a single angle
            _addOneModelParam(icov, EConsElem::RANGE, 0, EPSILON2, TEST);
          }
          else
          {
            for (int idim = 0; idim < ndim; idim++)
              _addOneModelParam(icov, EConsElem::RANGE, idim, EPSILON2, TEST);
            flagMultiRanges = true;
          }
        }
        else
        {
          // For other space dimension, consider Isotropic Range
          _addOneModelParam(icov, EConsElem::RANGE, 0, EPSILON2, TEST);
        }

        // Define possible Anisotropy rotation by angle(s)

        if (flagMultiRanges && flagAnisoRot)
        {
          // Skip rotation definition is already defined and should be shared
          if (flagRotationDefined && _modelPart._optvar.getLockSamerot()) continue;

          if (ndim == 2)
          {
            if (!_modelPart._optvar.getLockIso2d())
            {
              // For anisotropy in 2-D, one value is sufficient
              _addOneModelParam(icov, EConsElem::ANGLE, 0);
              flagRotationDefined = true;
            }
          }
          else if (ndim == 3)
          {
            for (int idim = 0; idim < ndim; idim++)
            {
              if (idim == 0 && _modelPart._optvar.getLockIso2d()) continue;
              _addOneModelParam(icov, EConsElem::ANGLE, idim);
            }
            flagRotationDefined = true;
          }
        }
      }
    }

    // Add the 'Sill' (vectorial) attribute
    if (flagSill)
    {
      if (!_modelPart._optvar.getFlagGoulardUsed())
      {
        int ijvar = 0;
        for (int ivar = ijvar = 0; ivar < nvar; ivar++)
          for (int jvar = 0; jvar <= ivar; jvar++, ijvar++)
            _addOneModelParam(icov, EConsElem::SILL, ijvar, TEST, TEST);
      }
    }

    // TODO: suppressed by DR
    // if (cova->hasParam())
    // {
    //   // Add the 'Param' (scalar) attribute
    //   _addOneModelParam(icov, EConsElem::PARAM, 0, 0., 2.);
    // }
  }

  // Final assignment
  _modelPart._calcmode.setAsVario(true);

  // Optional printout
  if (_modelPart._verbose) _dumpParamList();
  return 0;
}

  void AModelOptim::_patchModel(Model_Part & modelPart, const double* current)
  {
    // Initializations
    int ncov    = modelPart._model->getCovaNumber();
    int nvar    = modelPart._model->getVariableNumber();
    int nvs2    = nvar * (nvar + 1) / 2;
    int nparams = (int)modelPart._params.size();
    bool samerot = modelPart._optvar.getLockSamerot();

    // Loop on the parameters
    int iparam = 0;
    while (iparam < nparams)
    {
      const OneParam& param = modelPart._params[iparam];
      int icov              = param._icov;
      int rank              = param._rank;
      double scale          = param._scale;
      CovAniso* cova        = modelPart._model->getCova(icov);

      if (param._type == EConsElem::RANGE)
      {
        // Ranges

        // Use first direction range as isotropic
        if (rank == 0)
          cova->setRangeIsotropic(scale * current[iparam++]);
        else
          cova->setRange(rank, scale * current[iparam++]);
      }

      else if (param._type == EConsElem::ANGLE)
      {
        double angle = scale * current[iparam++];
        cova->setAnisoAngle(rank, angle);
        if (samerot)
        {
          // Export the Anisotropy Rotation information to all covariances
          for (int jcov = 0; jcov < ncov; jcov++)
          {
            CovAniso* mcova = modelPart._model->getCova(jcov);
            if (mcova->hasRange() > 0) mcova->setAnisoAngle(rank, angle);
          }
        }
      }

      else if (param._type == EConsElem::SILL)
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

      else if (param._type == EConsElem::PARAM)
      {
        // Set the 'param' attribute
        cova->setParam(scale * current[iparam++]);
      }

      else {
        messageAbort("AModelOptim: This should never happen");
      }
    }

    // Add 1 to the internal counter
    modelPart._niter++;
  }
