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

#include "Basic/Optim.hpp"
#include "Covariances/ACov.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Model/AModelFitSills.hpp"
#include "Model/ModelGeneric.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
  class ModelGeneric;

  class GSTLEARN_EXPORT AModelOptim
  {
  public:
    AModelOptim(ModelGeneric* model = nullptr, bool verbose = false);
#ifndef SWIG
    void setEnvironment(const MatrixSymmetric& vars,
                        double href,
                        double epsilon = EPSILON6,
                        double min = 0.,
                        double max = INF);
#endif
    AModelOptim& operator=(const AModelOptim& r);
    AModelOptim(const AModelOptim& r);
    void setAuthorizedAnalyticalGradients(bool authorized);

    bool getAuthorizedAnalyticalGradients() const;

    virtual ~AModelOptim();

    void setGradients(
      std::vector<std::function<double(const VectorDouble&)>>& gradients);

    void setVerbose(bool verbose = false, bool trace = false);

    double eval(const VectorDouble& x);

    virtual void evalGrad(vect res);
    double run();

    void resetIter();

    virtual double
      computeCost(bool flagPrint = false, bool verbose = false) = 0;

    std::shared_ptr<ListParams> getParams() const { return _params; }

    void evalGradInEffectiveDimension(vect res);

  private:
    void _printSummary(double minf, const VectorDouble& x) const;

  protected:
    ModelGeneric* _model; // Pointer to the model being optimized

  private:
    std::shared_ptr<ListParams>
      _params; // Parameters of the model to be optimized
    Optim* _opt;
    bool _verbose;
    bool _trace;
    VectorDouble _x;
    Id _iter;
  };
} // namespace gstlrn
