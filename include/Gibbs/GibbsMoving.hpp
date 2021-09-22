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
    int _neq;
    int _pivot;
    VectorInt _ranks;
    VectorVectorDouble _ll;
  }; // Per sample

public:
  GibbsMoving();
  GibbsMoving(Db* db, Model* model, Neigh* neigh);
  GibbsMoving(const GibbsMoving &r);
  GibbsMoving& operator=(const GibbsMoving &r);
  virtual ~GibbsMoving();

  void update(VectorVectorDouble& y,
              int isimu,
              int ipgs,
              int iter) override;
  int covmatAlloc(bool verbose) override;

  Neigh* getNeigh() const { return _neigh; }

private:
  Neigh* _neigh;
  std::vector<GibbsWeights> _wgt; // For each sample
};
