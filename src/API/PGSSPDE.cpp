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
#include "API/PGSSPDE.hpp"
#include "API/SPDE.hpp"

#include "Drifts/DriftList.hpp"
#include "Db/Db.hpp"
#include "Model/Model.hpp"
#include "Drifts/ADrift.hpp"
#include "Basic/String.hpp"
#include "LithoRule/RuleProp.hpp"
#include "Basic/Law.hpp"

PGSSPDE::PGSSPDE(const std::vector<Model*>& models,
                 const Db* field,
                 const RuleProp* ruleprop,
                 const Db* data)
    : _data(),
      _spdeTab(),
      _ruleProp(ruleprop),
      _calcul()
{
  _calcul = (data == nullptr) ? ESPDECalcMode::SIMUNONCOND :
                                ESPDECalcMode::SIMUCOND;
  for (const auto &e : models)
  {
    _spdeTab.push_back(new SPDE(e, field, data, _calcul));
  }
}

void PGSSPDE::compute(Db *dbout,
                      int seed,
                      int /*nitergibbs*/,
                      const NamingConvention &namconv)
{
  // Set the seed for random number generator
  law_set_random_seed(seed);

  int ngrf = (int) _spdeTab.size();
  VectorInt iuids(ngrf);

  VectorString namesG = generateMultipleNames("simuGauss", ngrf);
  int seed_loc = seed;
  for(int igrf = 0; igrf < ngrf; igrf++)
  {
    iuids[igrf] = _spdeTab[igrf]->compute(dbout, 1, seed_loc);
    seed_loc = 0;
  }
  dbout->setLocatorsByUID(iuids, ELoc::SIMU);

  if (_calcul == ESPDECalcMode::SIMUCOND)
    _ruleProp->categoryToThresh(_data);

  _ruleProp->gaussToCategory(dbout, namconv);

  dbout->deleteColumnsByUID(iuids);
}

PGSSPDE::~PGSSPDE()
{
  for(auto &e:_spdeTab)
  {
    delete e;
  }
}
