/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/NamingConvention.hpp"

#include "Enum/ESPDECalcMode.hpp"
#include "Enum/ECalcMember.hpp"

#include <vector>

class Db;
class DbGrid;
class Model;
class SPDE;
class RuleProp;

class GSTLEARN_EXPORT PGSSPDE
{
public:
  PGSSPDE(std::vector<Model*> models,
          const Db* field,
          const RuleProp* ruleprop,
          const Db* data = nullptr);
  PGSSPDE(const PGSSPDE& r) = delete;
  PGSSPDE& operator=(const PGSSPDE& r) = delete;
  virtual ~PGSSPDE();

  void compute(Db *dbout,
               int seed = 131351,
               int nitergibbs = 0,
               const NamingConvention &namconv = NamingConvention("Facies"));

private:
  Db*                _data;
  std::vector<SPDE*> _spdeTab;
  const RuleProp*    _ruleProp;
  ESPDECalcMode      _calcul;
};
