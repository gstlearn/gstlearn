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
#include "Covariances/CorGneiting.hpp"
#include "Db/DbGrid.hpp"
#include "Model/ModelGeneric.hpp"
#include "Simulation/Simulations.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  auto ctxt = CovContext(1, 1);

  VectorDouble scales  = {1};
  double alpha         = 1.5;
  double beta          = 1;
  double nu            = 3. / 4.;
  double timerange     = 2.5;
  auto* mg_in_gstlearn = CorGneiting::create(ctxt, ECov::GNEITING_G, alpha, beta, timerange,
                                             {nu}, {1.0}, scales, VectorDouble(), 1.0, false);
  auto modelGeneric    = ModelGeneric(ctxt);
  (void)modelGeneric.setCov(mg_in_gstlearn);

  VectorInt nx    = {100, 50};
  VectorDouble dx = {0.15, 0.3};
  auto* grd       = DbGrid::create(nx, dx);

  Id nsim = 20;
  Id ns   = 5000;
  (void)simuSpectral(nullptr, grd, &modelGeneric, nullptr, nsim, 0, ns, 100, nullptr, true,
                     NamingConvention("test_a_template"));

  return 0;
}
