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
                             const std::vector<Model *>& models,
                             double rho,
                             bool flag_propagation);
};
