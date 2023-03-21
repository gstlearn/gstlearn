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

#include "Gibbs/GibbsMulti.hpp"

class Db;
class Model;

class GSTLEARN_EXPORT GibbsUMulti : public GibbsMulti
{
public:
  GibbsUMulti();
  GibbsUMulti(Db* db, Model* model);
  GibbsUMulti(const GibbsUMulti &r);
  GibbsUMulti& operator=(const GibbsUMulti &r);
  virtual ~GibbsUMulti();

  void update(VectorVectorDouble& y,
              int isimu,
              int ipgs,
              int iter) override;
  int covmatAlloc(bool verbose, bool verboseTimer = false) override;

private:
  VectorDouble _covmat;
};
