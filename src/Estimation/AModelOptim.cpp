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
#include "Estimation/AModelOptim.hpp"
#include "Basic/OptCustom.hpp"
#include "Basic/Optim.hpp"
#include "Basic/VectorHelper.hpp"
#include "Model/ModelCovList.hpp"

namespace gstlrn
{
AModelOptim::AModelOptim(ModelGeneric* model, bool verbose)
  : _model(model)
  , _verbose(verbose)
  , _trace(false)
{
  if (_model == nullptr)
    throw std::invalid_argument("Model cannot be null");

  bool useGradient = static_cast<bool>(OptCustom::query("UseGradient", 1));
  _params          = _model->generateListParams();

  // Id nvar                 = _model->getNVar();
  // MatrixSymmetric varsUnit = MatrixSymmetric(nvar);
  // for (Id ivar = 0; ivar < nvar; ivar++) varsUnit.setValue(ivar, ivar, 1.);
  // _model->initParams(varsUnit, 1.);
  _x = _params->getOptimizableValues();
  if (useGradient)
    _opt = new Optim(LBFGS, static_cast<Id>(_x.size()));
  else
    _opt = new Optim(NELDERMEAD, static_cast<Id>(_x.size()));

  _opt->setXtolRel(EPSILON6);
  _opt->setObjective([this](const std::vector<double>& x)
                     { return this->eval(x); });
  _opt->setGradient([this](vect grad)
                    { this->evalGrad(grad); },
                    _params->getDispatch(),
                    _params->getDispatchIndex());
  resetIter();
};

void AModelOptim::setEnvironment(const MatrixSymmetric& vars, double href, double epsilon, double min, double max)
{
  _model->initParams(vars, href, min, max);
  _opt->setLowerBounds(_params->getMinValues(epsilon), _params->getDispatch());
  _opt->setUpperBounds(_params->getMaxValues(epsilon), _params->getDispatch());
  _x = _params->getOptimizableValues();
}

AModelOptim& AModelOptim::operator=(const AModelOptim& r)
{
  if (this != &r)
  {
    DECLARE_UNUSED(r)
    messerr("Assignment operator not implemented for AModelOptim");
  }
  return *this;
}

void AModelOptim::setAuthorizedAnalyticalGradients(bool authorized)
{
  if (_opt == nullptr)
  {
    messerr("Optimizer is not initialized");
    return;
  }
  _opt->setAuthorizedAnalyticalGradients(authorized);
}

bool AModelOptim::getAuthorizedAnalyticalGradients() const
{
  if (_opt == nullptr)
  {
    messerr("Optimizer is not initialized");
    return false;
  }
  return _opt->getAuthorizedAnalyticalGradients();
}

AModelOptim::~AModelOptim()
{
  delete _opt;
}

void AModelOptim::setGradients(std::vector<std::function<double(const std::vector<double>&)>>& gradients)
{
  if (_opt == nullptr)
  {
    messerr("Optimizer is not initialized");
    return;
  }
  _opt->setGradientComponents(gradients);
}
void AModelOptim::setVerbose(bool verbose, bool trace)
{
  _verbose = verbose;
  _trace   = trace;
  if (trace) _verbose = true;

  // Export 'verbose' and 'trace' flags down to FitSill (if defined)
  auto* mcv = dynamic_cast<ModelCovList*>(_model);
  if (mcv != nullptr)
  {
    AModelFitSills* amf = mcv->getFitSills();
    if (amf != nullptr)
    {
      // In this internal Fitting Sills procedure, the verbose flag is switched OFF
      // in order to avoid intermediate printouts
      amf->setVerbose(false);
      amf->setTrace(trace);
    }
  }

  // In the verbose case, first print the list of parameters
  if (verbose || trace)
    _params->display();
}

double AModelOptim::eval(const std::vector<double>& x)
{
  _iter++;

  // Set the current parameters inside the Model
  _params->setValues(x);

  // Update the different parameters of the Model
  _model->updateModel();

  // Calculate the cost
  double result = computeCost(false, false);

  if (_trace)
  {
    message("Iteration %4d - Cost = %lf", _iter, result);
    VH::dump(" - Current parameters", x, false);
  }

  return result;
};

void AModelOptim::evalGradInEffectiveDimension(vect res)
{
  if (_opt == nullptr)
  {
    messerr("Optimizer is not initialized");
    return;
  }
  _opt->evalGrad(res);
}

void AModelOptim::evalGrad(vect res) {
  DECLARE_UNUSED(res)
};
void AModelOptim::_printSummary(double minf, const std::vector<double>& x) const
{
  message("Summary of Optimization procedure:\n");
  message("Count of Iterations = %4d - Final Cost = %lf\n",
          _iter, minf);
  VH::dump("- Final parameters", x, false);
  auto* mcv           = dynamic_cast<ModelCovList*>(_model);
  AModelFitSills* amf = mcv->getFitSills();
  if (amf != nullptr)
  {
    Id nitergCum = mcv->getCovList()->getNitergCum();
    amf->printFitSillSummary(nitergCum);
  }
}

double AModelOptim::run()
{
  double minf = _opt->minimize(_x);

  if (_verbose) _printSummary(minf, _x);
  resetIter();
  return minf;
}

void AModelOptim::resetIter()
{
  _iter = 0;
}
} // namespace gstlrn