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
  const cs* getQ() const { return _Q; }
  bool getFlagSymNeigh() const { return _flagSymNeigh; }
  void setFlagSymNeigh(bool flagSymNeigh) { _flagSymNeigh = flagSymNeigh; }
  bool getFlagSymQ() const { return _flagSymQ; }
  void setFlagSymQ(bool flagSymQ) { _flagSymQ = flagSymQ; }
  bool getFlagPrintQ() const { return _flagPrintQ; }
  void setFlagPrintQ(bool flagPrintQ) { _flagPrintQ = flagPrintQ; }

private:
  int _improveConditioning(bool verbose);
  void _print(int iact) const;
  int _getVariableNumber() const;
  void _setQFlag(VectorBool& QFlag,
                 int nech,
                 int iech,
                 const VectorInt& ranks) const;
  VectorInt _getQFlag(VectorBool& QFlag, int nech, int iech);
  void _makeQSymmetric(cs* Q) const;
  int  _buildQ();
  void _extractWeightFromQ();

private:
  Neigh* _neigh;
  std::vector<GibbsWeights> _wgt; // For each sample
  bool _flagSymNeigh;
  bool _flagSymQ;
  bool _flagPrintQ;
  cs*  _Q;
};
