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

#include "Gibbs/GibbsMulti.hpp"
#include "Basic/Vector.hpp"

class Db;
class Model;

class GibbsMMulti : public GibbsMulti
{
private:
  struct GibbsWeights {
    int _pivot;
    VectorInt _ranks;
    VectorVectorDouble _ll;
  }; // Per sample

public:
  GibbsMMulti();
  GibbsMMulti(Db* db, Model* model);
  GibbsMMulti(const GibbsMMulti &r);
  GibbsMMulti& operator=(const GibbsMMulti &r);
  virtual ~GibbsMMulti();

  void update(VectorVectorDouble& y,
              int isimu,
              int ipgs,
              int iter) override;
  int covmatAlloc(bool verbose) override;

  const cs* getQ() const { return _Q; }
  void setEpsilon1(double epsilon) { _epsilon1 = epsilon; }
  void setEpsilon2(double epsilon) { _epsilon2 = epsilon; }

private:
  int _testConditioning(bool verbose);
  void _display(int iact) const;
  void _display() const;
  int _getVariableNumber() const;
  void _makeQSymmetric(cs* Q) const;
  int  _buildQ();
  void _extractWeightFromQ();
  void _tableStore(const Db* db, const cs* Cmat, const csn* N, bool verbose);

private:
  std::vector<GibbsWeights> _wgt; // For each sample
  cs*  _Q;
  double _epsilon1;
  double _epsilon2;
};
