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

#include "GibbsMultiMono.hpp"
#include "Gibbs/AGibbs.hpp"

class Db;
class Model;

/**
 * This class is designated for Gibbs with the following properties
 * - Unique (absent) Neighborhood
 * - Multivariate case: Multiple Monovariate systems
 * (even if the model is provided as multivariate)
 */
class GSTLEARN_EXPORT GibbsUMultiMono : public GibbsMultiMono
{
public:
  GibbsUMultiMono();
  GibbsUMultiMono(Db* db, std::vector<Model *> models, double rho);
  GibbsUMultiMono(const GibbsUMultiMono &r);
  GibbsUMultiMono& operator=(const GibbsUMultiMono &r);
  virtual ~GibbsUMultiMono();

  void update(VectorVectorDouble &y, int isimu, int ipgs, int iter) override;
  int covmatAlloc(bool verbose, bool verboseTimer = false) override;

private:
  double _getVariance(int ivar, int iact) const;
  double _getEstimate(int icase, int ivar, int iact, VectorVectorDouble& y) const;

private:
  VectorVectorDouble _covmat; // One matrix per variable
};
