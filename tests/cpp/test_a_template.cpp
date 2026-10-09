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
  // Unless you test a specific feature, bring it to its minimal expansion
  // e.g. the few lines below
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-");

  auto* db = Db::createFillRandom(100, 3, 1);
  auto* model = Model::createFromParam(ECov::CUBIC);

  (void)gibbs_sampler(db, model, 10, 426831, 10, 30, true);

  // GSTLEARN_EXPORT Id gibbs_sampler(
  //   Db * dbin, Model * model, Id nbsimu = 1, Id seed = 14341,
  //   Id gibbs_nburn = 10, Id gibbs_niter = 30, bool flag_moving = false,
  //   bool flag_norm = true, bool flag_multi_mono = false,
  //   bool flag_propagation = false, bool flag_sym_neigh = false,
  //   Id gibbs_optstats = 0, double percent = 5., bool flag_ce = true,
  //   bool flag_cstd = false, bool verbose = false,
  //   const NamingConvention& namconv = NamingConvention("Gibbs"));

  return 0;
}
