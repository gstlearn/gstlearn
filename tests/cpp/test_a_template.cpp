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
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Estimation/Vecchia.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  DECLARE_UNUSED(argc);
  DECLARE_UNUSED(argv);

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-");

  auto* model = Model::createFromParam(ECov::EXPONENTIAL);

  Id ndat  = 10;
  Id ndim  = 2;
  Id nvar  = 0;
  auto* db = Db::createFillRandom(ndat, ndim, nvar);

  auto* grid = DbGrid::createFillRandom({4, 4});

  auto neigh = NeighUnique(ndim);

  Id err = kriging(db, grid, model, &neigh);
  messerr("Error = %ld", err);
  return 0;
}
