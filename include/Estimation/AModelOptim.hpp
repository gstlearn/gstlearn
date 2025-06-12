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
  AModelOptim(ModelGeneric* model = nullptr, bool verbose = false)
    : _model(model)
    , _verbose(verbose)
  {
    if (_model == nullptr)
    {
      throw std::invalid_argument("Model cannot be null");
    }
    _params = _model->generateListParams();
    _model->initParams();
    _x    = _params->getValues();
    _xmin = _params->getMinValues();
    _xmax = _params->getMaxValues();
    _opt  = new Optim(opt_algorithm::NELDERMEAD, (int) _x.size());
    _opt->setLowerBounds(_xmin);
    _opt->setUpperBounds(_xmax);
    _opt->setXtolRel(EPSILON6);
    _opt->setObjective([this](const std::vector<double>& x)
                       { return this->eval(x); });
    _iter = 0;
  };

  AModelOptim(const AModelOptim& r) 
  {
    DECLARE_UNUSED(r)
    throw std::runtime_error("You shoudln't arrive here!");
  };

  AModelOptim& operator=(const AModelOptim& r)
  {
    if (this != &r)
    {
      DECLARE_UNUSED(r)
      messerr("Assignment operator not implemented for AModelOptim");
    }
    return *this;
  }

  virtual ~AModelOptim()
  {
    delete _opt;
  }

  void setVerbose(bool verbose)
  {
    _verbose = verbose;

    // In the verbose case, first print the list of parameters
    if (verbose) _params->display();
  }

  double eval(const std::vector<double>& x)
  {
    _iter++;
    _params->setValues(x);
    _model->updateModel();

    // Calculate the cost
    double result = computeCost(false);

    if (_verbose)
    {
      message("Iteration %4d - Cost = %lf", _iter, result);
      VH::dump(" - Current parameters", x, false);
    }

    // Check if Goulard must be applied
    ModelCovList* mcv = dynamic_cast<ModelCovList*>(_model);
    if (mcv != nullptr)
    {
      AModelFitSills* amf = mcv->getFitSills();
      if (amf != nullptr) amf->fitSills(_verbose);
    }

    return -result;
  };

  void run()
  {
    _opt->optimize(_x);
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
  std::vector<double> _x;
  std::vector<double> _xmin;
  std::vector<double> _xmax;
  int _iter;
};
