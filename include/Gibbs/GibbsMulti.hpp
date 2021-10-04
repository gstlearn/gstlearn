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

class GibbsMulti : public AGibbs
{
public:
  GibbsMulti();
  GibbsMulti(Db* db, Model* model);
  GibbsMulti(const GibbsMulti &r);
  GibbsMulti& operator=(const GibbsMulti &r);
  virtual ~GibbsMulti();

  int calculInitialize(VectorVectorDouble& y,
                       int isimu,
                       int ipgs,
                       bool verbose);
  double getSimulate(VectorVectorDouble& y,
                     double yk,
                     double sk,
                     int iact,
                     int ipgs,
                     int ivar,
                     int iter);
  int checkGibbs(const VectorVectorDouble& y, int isimu, int ipgs);

  Model* getModel() const { return _model; } // protect using const asap

private:
  Model* _model;
};
