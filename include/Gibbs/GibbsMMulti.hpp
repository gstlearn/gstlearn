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
class Neigh;

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
  GibbsMMulti(Db* db, Model* model, Neigh* neigh);
  GibbsMMulti(const GibbsMMulti &r);
  GibbsMMulti& operator=(const GibbsMMulti &r);
  virtual ~GibbsMMulti();

  void update(VectorVectorDouble& y,
              int isimu,
              int ipgs,
              int iter) override;
  int covmatAlloc(bool verbose) override;

  Neigh* getNeigh() const { return _neigh; }

private:
  int _isValidConditioning() const;
  void _print(int iact) const;
  int _getVariableNumber() const;
  void _setQCont(VectorBool& QCont,
                 int nech,
                 int iech,
                 const VectorInt& ranks,
                 bool flagSym) const;
  VectorInt _getQCont(VectorBool& Qcont, int nech, int iech);

private:
  Neigh* _neigh;
  std::vector<GibbsWeights> _wgt; // For each sample
};
