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
#include "Basic/VectorT.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Estimation/Vecchia.hpp"
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

  auto* model = Model::createFromParam(ECov::EXPONENTIAL);
  model->setDriftIRF(0);

  Id nbvecchia = 4;
  double like;
  bool verbose = true;

  Id nvar  = 1;
  Id dim   = 2;
  Id ndat  = 6;
  Id indNa = 3;
  Db* db   = Db::createFillRandom(ndat, dim, nvar, 0, 0, 0);
  db->setValue("z", indNa, TEST);
  auto dbfmt = DbStringFormat(FLAG_ARRAY, VectorString(), VectorInt(), false);
  db->display(&dbfmt);

  like = logLikelihoodVecchia(db, model, nbvecchia, verbose);
  message("Likelihood without filtering %lf\n", like);

  VectorDouble sel(ndat, 1.);
  sel[indNa] = 0.;
  db->addSelection(sel, "sel");
  db->display(&dbfmt);
  like = logLikelihoodVecchia(db, model, nbvecchia, verbose);
  message("Likelihood with selection %lf\n", like);

  db->deleteSample(indNa);
  db->display(&dbfmt);
  like = logLikelihoodVecchia(db, model, nbvecchia, verbose);
  message("Likelihood after deleting the NA sample %lf\n", like);

  return 0;
}
