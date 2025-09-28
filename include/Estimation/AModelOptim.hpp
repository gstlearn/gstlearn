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

#include "Covariances/ACov.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Model/AModelFitSills.hpp"
#include "Model/ModelGeneric.hpp"
#include "Basic/Optim.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class ModelGeneric;

class GSTLEARN_EXPORT AModelOptim
{
public:
  AModelOptim(ModelGeneric* model = nullptr,
              bool verbose        = false);
  void setEnvironment(const MatrixSymmetric& vars, double href, double epsilon = EPSILON6);

  AModelOptim& operator=(const AModelOptim& r);

  void setAuthorizedAnalyticalGradients(bool authorized);

  bool getAuthorizedAnalyticalGradients() const;

  virtual ~AModelOptim();

  void setGradients(std::vector<std::function<double(const std::vector<double>&)>>& gradients);

  void setVerbose(bool verbose = false, bool trace = false);

  double eval(const std::vector<double>& x);

  virtual void evalGrad(vect res);
  double run();

  void resetIter();

  virtual double computeCost(bool verbose = false) = 0;
  std::shared_ptr<ListParams> getParams() const { return _params; }
  void evalGradInEffectiveDimension(vect res);
private:
  void _printSummary(double minf, const std::vector<double>& x) const;

protected:
  ModelGeneric* _model; // Pointer to the model being optimized

private:
  std::shared_ptr<ListParams> _params; // Parameters of the model to be optimized
  Optim* _opt;
  bool _verbose;
  bool _trace;
  std::vector<double> _x;
  Id _iter;
};
}