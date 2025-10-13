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

#include "Model/GaussianProcess.hpp"

#include "Model/Model.hpp"
#include "Model/ModelGeneric.hpp"
#include "Db/Db.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Basic/Law.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "geoslib_define.h"
#include <memory>

namespace gstlrn
{
GaussianProcess::GaussianProcess(ModelGeneric* model, Db* data)
  : AStringable()
  , _model(nullptr)
  , _data(nullptr)
{
  init(model, data);
}

GaussianProcess::GaussianProcess(const GaussianProcess& gp)
  : AStringable()
  , _model(gp._model)
  , _data(gp._data)
{
}

GaussianProcess& GaussianProcess::operator=(const GaussianProcess& gp)
{
  if (this != &gp)
  {
    _model = gp._model;
    _data  = gp._data;
  }
  return *this;
}

GaussianProcess::~GaussianProcess()
{
}

String GaussianProcess::toString(const AStringFormat* strfmt) const
{
  String sstr              = "Gaussian Process Model:\n";
  const auto* modelderived = dynamic_cast<const Model*>(_model.get());
  if (modelderived == nullptr)
  {
    sstr.append("No model defined.\n");
    return sstr;
  }
  sstr.append("Model:\n");
  sstr.append(modelderived->toString(strfmt));
  sstr.append("Data:\n");
  if (_data == nullptr)
  {
    sstr.append("No data defined.\n");
    return sstr;
  }
  sstr.append(_data->toString(strfmt));
  return sstr;
}

void GaussianProcess::init(const ModelGeneric* model, const Db* data)
{
  if (model != nullptr)
    _model = std::shared_ptr<ModelGeneric>(model->clone());
  if (data != nullptr)
    _data = std::shared_ptr<Db>(data->clone());
}

void GaussianProcess::fit(Id nb_neighVecchia,
                          bool reml,
                          bool verbose,
                          bool trace)
{
  if (_model && _data)
  {
    _model->fitNew(_data.get(), nullptr, nullptr, nullptr,
                   ModelOptimParam(), nb_neighVecchia, verbose, trace, reml);
  }
}

void GaussianProcess::predict(Db* out)
{
  NeighUnique neigh;
  if (_model && _data && out)
  {
    kriging(getData().get(), out, getModel().get(), &neigh);
  }
}

void GaussianProcess::simulate(Db* out, Id nbsimus)
{
  NeighUnique neigh;
  auto seed = law_get_random_seed();
  if (_model && _data)
  {
    simtub(getData().get(), out, dynamic_cast<Model*>(getModel().get()), &neigh, nbsimus, seed);
  }
}
}