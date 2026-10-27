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
#include "Basic/ASerializable.hpp"
#include "Db/Db.hpp"
#include "Enum/ECov.hpp"
#include "Model/Model.hpp"
#include "Simulation/Simulations.hpp"
#include "geoslib_define.h"
#include "geoslib_f.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Unless you test a specific feature, bring it to its minimal expansion
  // e.g. the few lines below
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-");

  auto np = 100;
  auto seed = 14313;

  auto* data = Db::createFillRandom(np, 2, 0);
  auto* model = Model::createFromParam(ECov::EXPONENTIAL, 0.5, 1.);
  model->setMeans(0);
  model->setDriftIRF(-1);
  (void)simtub(NULL, data, model, nullptr, 1, seed);

  auto nburn = 10;
  auto niter = 30; // Was 300 but the time consumption was too long for tests
  auto nbsimu = 5;

  (void)gibbs_sampler(
    data, model, nbsimu, seed, nburn, niter, true, true, false, false, false, 0,
    5., false, false, true);

  data->getContentsAsTable().display();
  return 0;
}
