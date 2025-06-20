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
#pragma once

#include "Basic/AStringable.hpp"
#include "Covariances/ACov.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Model/AModelFitSills.hpp"
#include "Model/ModelGeneric.hpp"
#include "Model/ModelCovList.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
#include "Basic/Optim.hpp"
#include "Basic/VectorHelper.hpp"

class ModelGeneric;

class GSTLEARN_EXPORT AModelOptim
{
public:
  AModelOptim(ModelGeneric* model = nullptr,
              bool verbose        = false)
    : _model(model)
    , _verbose(verbose)
    , _trace(false)
  {
    if (_model == nullptr)
      throw std::invalid_argument("Model cannot be null");

    bool useGradient = true;
    _params          = _model->generateListParams();

    int nvar                 = _model->getNVar();
    MatrixSymmetric varsUnit = MatrixSymmetric(nvar);
    for (int ivar = 0; ivar < nvar; ivar++) varsUnit.setValue(ivar, ivar, 1.);
    _model->initParams(varsUnit, 1.);
    _x    = _params->getOptimizableValues();
    _xmin = _params->getMinValues();
    _xmax = _params->getMaxValues();
    if (useGradient)
      _opt = new Optim(LBFGS, (int)_x.size());
    else
      _opt = new Optim(NELDERMEAD, (int)_x.size());
    _opt->setLowerBounds(_xmin);
    _opt->setUpperBounds(_xmax);
    _opt->setXtolRel(EPSILON6);
    _opt->setObjective([this](const std::vector<double>& x)
                       { return this->eval(x); });
    _opt->setGradient([this](vect grad)
                      { this->evalGrad(grad); },
                      _params->getDispatch(),
                      _params->getDispatchIndex());
    resetIter();
  };

  void setEnvironment(const MatrixSymmetric& vars, double href)
  {
    _model->initParams(vars, href);
    _x = _params->getOptimizableValues();
  }

  AModelOptim& operator=(const AModelOptim& r)
  {
    if (this != &r)
    {
      DECLARE_UNUSED(r)
      messerr("Assignment operator not implemented for AModelOptim");
    }
    return *this;
  }

  void setAuthorizedAnalyticalGradients(bool authorized)
  {
    if (_opt == nullptr)
    {
      messerr("Optimizer is not initialized");
      return;
    }
    _opt->setAuthorizedAnalyticalGradients(authorized);
  }

  bool getAuthorizedAnalyticalGradients() const
  {
    if (_opt == nullptr)
    {
      messerr("Optimizer is not initialized");
      return false;
    }
    return _opt->getAuthorizedAnalyticalGradients();
  }

  virtual ~AModelOptim()
  {
    delete _opt;
  }

  void setGradients(std::vector<std::function<double(const std::vector<double>&)>>& gradients)
  {
    if (_opt == nullptr)
    {
      messerr("Optimizer is not initialized");
      return;
    }
    _opt->setGradientComponents(gradients);
  }
  void setVerbose(bool verbose = false, bool trace = false)
  {
    _verbose = verbose;
    _trace   = trace;
    if (trace) _verbose = true;

    // Export 'verbose' and 'trace' flags down to FitSill (if defined)
    ModelCovList* mcv = dynamic_cast<ModelCovList*>(_model);
    if (mcv != nullptr)
    {
      AModelFitSills* amf = mcv->getFitSills();
      if (amf != nullptr)
      {
        amf->setVerbose(verbose);
        amf->setTrace(trace);
      }
    }

    // In the verbose case, first print the list of parameters
    if (verbose || trace)
      _params->display();
  }

  double eval(const std::vector<double>& x)
  {
    _iter++;

    // Set the current parameters inside the Model
    _params->setValues(x);

    // Update the different parameters of the Model
    _model->updateModel();

    // Calculate the cost
    double result = computeCost(false);

    if (_trace)
    {
      message("Iteration %4d - Cost = %lf", _iter, result);
      VH::dump(" - Current parameters", x, false);
    }

    return result;
  };

  virtual void evalGrad(vect res) {
    DECLARE_UNUSED(res)
  };
  void _printSummary(double minf, const std::vector<double>& x) const
  {
    message("Summary of Optimization procedure:\n");
    message("Count of Iterations = %4d - Final Cost = %lf\n",
            _iter, minf);
    VH::dump("- Final parameters", x, false);
    ModelCovList* mcv   = dynamic_cast<ModelCovList*>(_model);
    AModelFitSills* amf = mcv->getFitSills();
    if (amf != nullptr)
    {
      int nitergCum = mcv->getCovList()->getNitergCum();
      amf->printFitSillSummary(nitergCum);
    }
  }

  void run()
  {
    double minf = _opt->minimize(_x);

    if (_verbose) _printSummary(minf, _x);
  }

  void resetIter()
  {
    _iter = 0;
  }
  virtual double computeCost(bool verbose = false) = 0;

protected:
  ModelGeneric* _model; // Pointer to the model being optimized

private:
  std::shared_ptr<ListParams> _params; // Parameters of the model to be optimized
  Optim* _opt;
  bool _verbose;
  bool _trace;
  std::vector<double> _x;
  std::vector<double> _xmin;
  std::vector<double> _xmax;
  int _iter;
};
