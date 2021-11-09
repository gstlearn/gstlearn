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
#include "Basic/AStringable.hpp"

class Db;
class Model;

class GibbsMultiMono : public AGibbs
{
public:
  GibbsMultiMono();
  GibbsMultiMono(Db* db, std::vector<Model*> models, double rho);
  GibbsMultiMono(const GibbsMultiMono &r);
  GibbsMultiMono& operator=(const GibbsMultiMono &r);
  virtual ~GibbsMultiMono();

  Model* getModels(int ivar) const { return _models[ivar]; } // TODO: protect by const asap
  double getRho() const { return _rho; }
  int getVariableNumber() const { return static_cast<int>(_models.size()); }

  int calculInitialize(VectorVectorDouble& y,
                       int isimu,
                       int ipgs);
  double getSimulate(VectorVectorDouble& y,
                     double yk,
                     double sk,
                     int ipgs,
                     int ivar,
                     int iact,
                     int iter);
  int checkGibbs(const VectorVectorDouble& y, int isimu, int ipgs);

private:
  std::vector<Model *> _models;
  double _rho;
};
