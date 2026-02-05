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
#include "Simulation/SimuSpectralRN.hpp"
#include "geoslib_define.h"

#include "Basic/ASerializable.hpp"
#include "Model/Model.hpp"
#include "Simulation/CalcSimuSpectral.hpp"
#include "geoslib_f.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  // type_cov <- c("GAUSSIAN", "EXPONENTIAL", "MATERN")[idx_type]
  // nu <- c(10, 1/2, 3/4)[idx_type]
  // ranges = c(5, 10, 15)
  // angles = c(30., 0, 0)

  Id ndim = 1;
  Id nb   = 1000;
  Id seed = 13112;

  defineDefaultSpace(ESpaceType::RN, ndim);
  VectorInt nx    = {10000};
  VectorDouble dx = {1.};
  auto* db        = DbGrid::create(nx, dx);

  auto* mod = Model::createFromParam(ECov::MATERN, 5., 1, 0.75);
  auto sim  = SimuSpectralRN(1, nb, 100, seed, NULL, true);
  (void)sim.setModelGeneric(mod);
  (void)sim.simulate(mod->getCovAniso(0));
  auto gamma = sim.getGamma();
  auto omega = sim.getOmega();
  (void)sim.compute(db);

  db->getStatsAsTable().display();
  return 0;
}
