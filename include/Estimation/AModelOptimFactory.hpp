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

#include "Estimation/Vecchia.hpp"
#include "Estimation/Likelihood.hpp"
#include "Estimation/AModelOptimNew.hpp"
#include "Model/ModelOptimVMap.hpp"
#include "Model/ModelOptimVario.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

class AModelOptimNew;

class GSTLEARN_EXPORT AModelOptimFactory
{
  public:
  static AModelOptimNew* create(
    ModelGeneric* model,
    const Db* db,
    Vario* vario,
    const DbGrid* dbmap,
    Constraints* constraints,
    const Option_AutoFit& mauto,
    const Option_VarioFit& optvar,
    int nb_neighVecchia = ITEST)
  {
 
    if (db != nullptr && nb_neighVecchia != ITEST)
    {
      return Vecchia::createForOptim(model, db, nb_neighVecchia);
    }
    if (dbmap != nullptr)
    {
      return ModelOptimVMap::createForOptim(model, dbmap, constraints, mauto, optvar);
    }
    if (vario != nullptr)
    {
      return ModelOptimVario::createForOptim(model, vario, constraints, mauto, optvar);
    }
    return Likelihood::createForOptim(model, db);
  }  
};