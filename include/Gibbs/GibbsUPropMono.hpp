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

#include "Gibbs/AGibbs.hpp"
#include "Basic/Vector.hpp"

class Db;
class Model;

/**
 * This class is designated for Gibbs with the following properties
 * - Unique (absent) Neighborhood
 * - Monovariate case only
 * - Propagation algorithm (no need to establish and invert Covariance matrix)
 * - No bound provided
 */
class GibbsUPropMono : public AGibbs
{
public:
  GibbsUPropMono();
  GibbsUPropMono(Db* db, Model* model);
  GibbsUPropMono(const GibbsUPropMono &r);
  GibbsUPropMono& operator=(const GibbsUPropMono &r);
  virtual ~GibbsUPropMono();

  void update(VectorVectorDouble& y,
              int isimu,
              int ipgs,
              int iter) override;
  int covmatAlloc(bool verbose) override;
};
