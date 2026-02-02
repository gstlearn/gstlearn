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
#include "geoslib_define.h"

#include "Basic/ASerializable.hpp"
#include "Covariances/CorAniso.hpp"
#include "Db/DbGrid.hpp"
#include "Model/ModelGeneric.hpp"
#include "Simulation/CalcSimuSpectral.hpp"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  Id ndim                 = 2;
  Id nvar                 = 1;
  auto gs_type            = ECov::EXPONENTIAL;
  auto gs_ctxt            = CovContext(nvar, ndim);
  VectorDouble ani_scales = {3. / 4., 3. / 2.};
  VectorDouble ani_angles = {30., 0.};
  double nus              = .3 / 4.;
  auto* gs_cov            = CorAniso::createAnisotropic(gs_ctxt, gs_type,
                                                        ani_scales,
                                                        nus,
                                                        ani_angles,
                                                        false);
  auto* modelGeneric      = new ModelGeneric(gs_ctxt);
  modelGeneric->setCov(gs_cov);

  VectorDouble x0 = {0., 0.};
  VectorInt nx    = {200, 200};
  VectorDouble dx = {0.1, 0.1};
  auto* grid      = DbGrid::create(nx, dx, x0);
  grid->display();

  Id nbsimu    = 3;
  Id seed      = 234567;
  Id ns        = 1000;
  Id nd        = 100;
  bool verbose = true;
  (void)simuSpectral(nullptr, grid, modelGeneric, nullptr, nbsimu, seed, ns, nd, nullptr, verbose,
                     NamingConvention("mySimu"));
  grid->display();

  return 0;
}
