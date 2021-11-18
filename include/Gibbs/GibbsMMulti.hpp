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

  void setEpsilon1(double epsilon) { _epsilon1 = epsilon; }
  void setEpsilon2(double epsilon) { _epsilon2 = epsilon; }
  void setFlagCheckCovariance(bool flagCheckCovariance) { _flagCheckCovariance = flagCheckCovariance; }

private:
  void _display(int iact) const;
  void _display() const;
  int  _getVariableNumber() const;
  void _tableStore(const Db* db, const cs* Cmat, bool verbose);
  double _checkForSampleIdentity(int iact, const cs* Cmat) const;
  void   _checkForIdentity(const cs* Cmat, bool verbose = false) const;
  VectorVectorDouble _getWeights(int iech,
                                 int *nbgh_arg,
                                 int *pivot_arg) const;

private:
  cs*       _Ln;
  VectorInt _Pn;
  double    _epsilon1;
  double    _epsilon2;
  bool      _flagCheckCovariance;

  // Mutable arrays (declared to speed up the process)
  mutable VectorInt _ranks;
  mutable VectorDouble _b;
  mutable VectorDouble _x;
};
