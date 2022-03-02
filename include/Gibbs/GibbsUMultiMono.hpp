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
#include "GibbsMultiMono.hpp"
#include "Gibbs/AGibbs.hpp"
#include "Basic/Vector.hpp"

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

  void update(VectorVectorDouble& y,
              int isimu,
              int ipgs,
              int iter) override;
  int covmatAlloc(bool verbose, bool verboseTimer = false) override;

private:
  VectorVectorDouble _covmat; // One matrix per variable
};
