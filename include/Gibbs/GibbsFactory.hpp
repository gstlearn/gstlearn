/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Gibbs/AGibbs.hpp"

class Db;
class Model;

class GSTLEARN_EXPORT GibbsFactory
{
public:
  GibbsFactory();
  virtual ~GibbsFactory();

  static AGibbs *createGibbs(Db* db,
                             Model* model,
                             bool flagMoving);
  static AGibbs *createGibbs(Db* db,
                             std::vector<Model *> models,
                             double rho,
                             bool flag_propagation);
};
