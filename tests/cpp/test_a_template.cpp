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
#include "Model/Model.hpp"
#include "geoslib_define.h"
#include "geoslib_f.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  auto* data = Db::createFromNF("/home/drenard/DataBefore.NF");
  data->deleteColumn("g*");
  data->deleteColumn("es.y.*");
  data->deleteColumn("*Bound");
  data->deleteColumn("code");
  data->deleteColumn("Simu");
  data->deleteColumn("Ind");

  auto* model = Model::createFromParam(ECov::EXPONENTIAL, 0.5);
  model->display();

  auto nbsimu          = 4;
  auto nburn           = 10;
  auto niter           = 30;
  auto gibbsOptStats   = 0;
  auto percent         = 0.;
  auto flagMoving      = true;
  auto flagNorm        = false;
  auto flagMultiMono   = false;
  auto flagPropagation = false;
  auto flagSymNeigh    = false;
  auto flagCE          = false;
  auto flagSTD         = false;
  auto verbose         = false;
  VectorInt seed_gb    = {2113, 7555, 6002, 4124};

  data->display();
  for (Id s = 0; s < nbsimu; s++)
  {
    message(">>>> processing simulation #%d\n", s + 1);
    std::string s_nm = "cs.y." + std::to_string(s);
    data->clearLocators(ELoc::Z);
    (void)gibbs_sampler(data, model, 1, seed_gb[s], nburn, niter, flagMoving, flagNorm, flagMultiMono,
                        flagPropagation, flagSymNeigh, gibbsOptStats, percent, flagCE, flagSTD, verbose,
                        NamingConvention(s_nm));
  }
  data->getStatsAsTable({"cs.y.*"}).display();

  return 0;
}
