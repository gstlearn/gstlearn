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

#include "API/newAPIs.hpp"
#include "Basic/VectorT.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovContext.hpp"
#include "Db/Db.hpp"
#include "Model/GaussianProcess.hpp"
#include "Model/Model.hpp"

namespace gstlrn
{
Db* createDbFromDataFrame(const DataFrame* dat,
                          const VectorString& coordinates)
{
  auto* db = new Db();
  for (const auto& [name, values]: *dat)
  {

    db->setColumn(values, name);
  }
  db->setLocators(coordinates, ELoc::X);
  return db;
}

GaussianProcess* createModelFromData(const Db* dat,
                                     const VectorString& variables,
                                     const std::vector<ECov>& structs,
                                     bool addMeasurementError)
{

  auto* gp = new GaussianProcess();
  CovContext ctxt(static_cast<Id>(variables.size()), dat->getNDim());
  Model model(ctxt);
  if (structs.empty())
  {
    messerr("No covariance structures provided.");
    return nullptr;
  }

  if (addMeasurementError)
  {
    CovAniso nugget(ECov::NUGGET, ctxt);
    model.addCov(nugget);
  }

  for (size_t i = 0; i < structs.size(); ++i)
  {
    CovAniso covi(structs[i], ctxt);
    model.addCov(covi);
  }
  model.setDriftIRF(0);
  gp->init(&model, dat);
  auto data = gp->getData();
  data->setLocators(variables, ELoc::Z);
  return gp;
}
} // namespace gstlrn