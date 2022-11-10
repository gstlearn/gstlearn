/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
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
