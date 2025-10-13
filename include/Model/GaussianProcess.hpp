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

#include "gstlearn_export.hpp"
#include "Basic/AStringable.hpp"
#include "Model/Model.hpp"
#include <memory>

namespace gstlrn
{
class Db;
class ModelGeneric;
class Model;

class GSTLEARN_EXPORT GaussianProcess: public AStringable
{
public:
  GaussianProcess(ModelGeneric* model = nullptr, Db* data = nullptr);
  GaussianProcess(const GaussianProcess& gp);
  GaussianProcess& operator=(const GaussianProcess& gp);
  virtual ~GaussianProcess();

  void init(const ModelGeneric* model, const Db* data);
  String toString(const AStringFormat* strfmt = nullptr) const override;
  void fit(Id nb_neighVecchia = ITEST,
           bool reml           = false,
           bool verbose        = false,
           bool trace          = false);
  void predict(Db* out);
  void simulate(Db* out, Id nbsimus = 1);
  Model* getRawModel() { return _getRawModel(); }
#ifndef SWIG
  auto getModel() { return _model; }
  auto getData() { return _data; }

#endif

private:
#ifndef SWIG
  Model* _getRawModel() const { return dynamic_cast<Model*>(_model.get()); }
#endif

private:
  std::shared_ptr<ModelGeneric> _model; // The model used for the Gaussian process
  std::shared_ptr<Db> _data;            // The database containing the data for the Gaussian process
};
}