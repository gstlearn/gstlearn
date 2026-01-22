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

#include "Gibbs/AGibbs.hpp"
#include "GibbsMultiMono.hpp"

namespace gstlrn
{
class Db;
class Model;

/**
 * This class is designated for Gibbs with the following properties
 * - Unique (absent) Neighborhood
 * - Multivariate case: Multiple Monovariate systems
 * (even if the model is provided as multivariate)
 */
class GSTLEARN_EXPORT GibbsUMultiMono: public GibbsMultiMono
{
public:
  GibbsUMultiMono();
  GibbsUMultiMono(Db* db, const std::vector<Model*>& models, double rho);
  GibbsUMultiMono(const GibbsUMultiMono& r);
  GibbsUMultiMono& operator=(const GibbsUMultiMono& r);
  virtual ~GibbsUMultiMono();

  void update(VectorVectorDouble& y, Id isimu, Id ipgs, Id iter) override;
  Id covmatAlloc(bool verbose, bool verboseTimer = false) override;

private:
  double _getVariance(Id ivar, Id iact) const;
  double _getEstimate(Id icase, Id ivar, Id iact, VectorVectorDouble& y) const;

private:
  std::vector<MatrixSquare> _covmat; // One matrix per variable
};
} // namespace gstlrn