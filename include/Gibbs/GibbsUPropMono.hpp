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

#include "Gibbs/GibbsMultiMono.hpp"

class Db;
class Model;

/**
 * This class is designated for Gibbs with the following properties
 * - Unique (absent) Neighborhood
 * - Monovariate case only
 * - Propagation algorithm (no need to establish and invert Covariance matrix)
 * - No bound provided
 */
class GSTLEARN_EXPORT GibbsUPropMono : public GibbsMultiMono
{
public:
  GibbsUPropMono();
  GibbsUPropMono(Db* db, std::vector<Model *> models, double rho);
  GibbsUPropMono(const GibbsUPropMono &r);
  GibbsUPropMono& operator=(const GibbsUPropMono &r);
  virtual ~GibbsUPropMono();

  void update(VectorVectorDouble& y,
              int isimu,
              int ipgs,
              int iter) override;
  int covmatAlloc(bool verbose, bool verboseTimer = false) override;

  double getEps() const { return _eps; }
  void setEps(double eps) { _eps = eps; }
  double getRval() const { return _rval; }
  void setRval(double rval) { _rval = rval; }

private:
  double _rval;
  double _eps;
};
