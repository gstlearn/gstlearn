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
#include "Basic/OptCst.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorT.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
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
  OptCst::define(ECst::NTCOL, -1);
  OptCst::define(ECst::NTROW, -1);

  Id nvar    = 3;
  Id ndim    = 2;
  Id ndat    = 15;
  auto dbfmt = DbStringFormat(FLAG_ARRAY, VectorString(), VectorInt(), false);
  VectorDouble hetero(nvar, 0.1);
  auto* db           = Db::createFillRandom(ndat, ndim, nvar, 0, 0, 0., 0.1, hetero);
  VectorString names = {"z-1", "z-3"};
  VectorInt ranks    = VH::sequence(5., 3, 1);
  Db* dbaux;

  // Complete file
  mestitle(1, "Initial data set");
  db->display(&dbfmt);

  message("\n---> Reducing by:\n");
  message(" - selecting some variables (z-1 and z-3)\n");
  message(" - suppressing masked samples\n");
  dbaux = Db::createReduce(db, names);
  dbaux->display(&dbfmt);
  delete dbaux;

  message("\n---> Reducing by:\n");
  message(" - selecting some variables (z-1 and z-3)\n");
  message(" - selecting some sample ranks (5 samples starting from rank 3)\n");
  dbaux = Db::createReduce(db, names, ranks);
  dbaux->display(&dbfmt);
  delete dbaux;

  message("\n---> Reducing by:\n");
  message(" - selecting some variables (z-1 and z-3)\n");
  message(" - selecting only samples where all variables are defined\n");
  dbaux = Db::createReduce(db, names, VectorInt(), true);
  dbaux->display(&dbfmt);
  delete dbaux;

  // Checking operators
  db->display(&dbfmt);
  double value = (*db)(3, "z-1");
  message("Valeur initiale = %f\n", value);
  (*db)(3, "z-1") = 12.;
  db->display(&dbfmt);
  return 0;
}
