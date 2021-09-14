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
class Neigh;

class GibbsMoving : public AGibbs
{
private:
  struct GibbsWeights {
    VectorInt _ivars;
    VectorInt _iechs;
    VectorInt _pivot;
    VectorDouble _sigma;
    VectorVectorDouble _ll;
  };

public:
  GibbsMoving();
  GibbsMoving(const GibbsMoving &r);
  GibbsMoving& operator=(const GibbsMoving &r);
  virtual ~GibbsMoving();

  int calculInitialize(int flag_category,
                       int flag_order,
                       Db *dbin,
                       Model *model,
                       int isimu,
                       int igrf,
                       int ipgs,
                       bool verbose);
  int calculIteration(Db *dbin,
                      Model *model,
                      int isimu,
                      int ipgs,
                      int igrf,
                      int verbose);
  int covmatAlloc(Db *dbin, Model *model, Neigh *neigh, bool verbose);

private:
  std::vector<GibbsWeights> _wgt;
};
