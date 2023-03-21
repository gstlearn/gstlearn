/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "geoslib_define.h"
#include "Gibbs/AGibbs.hpp"

class Db;
class Model;

class GSTLEARN_EXPORT GibbsMulti: public AGibbs
{
public:
  GibbsMulti();
  GibbsMulti(Db* db, Model* model);
  GibbsMulti(const GibbsMulti &r);
  GibbsMulti& operator=(const GibbsMulti &r);
  virtual ~GibbsMulti();

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

  Model* getModel() const { return _model; } // protect using const asap

private:
  Model* _model;
};
