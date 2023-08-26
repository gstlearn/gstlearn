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
#include "Basic/AStringable.hpp"

class Db;
class Model;

class GSTLEARN_EXPORT GibbsMultiMono : public AGibbs
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
