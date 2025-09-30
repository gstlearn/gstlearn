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

  Id nbVecchia = 8;
  double like;
  bool verbose = true;

  Id nvar    = 1;
  Id dim     = 2;
  Id ndat    = 6;
  Id indNa   = 3;
  auto dbfmt = DbStringFormat(FLAG_VARS | FLAG_ARRAY, VectorString(), VectorInt(), false);
  Db* dbref  = Db::createFillRandom(ndat, dim, nvar, 0, 0, 0);
  dbref->setValue("z", indNa, TEST);
  VectorDouble sel(ndat, 1.);
  sel[indNa] = 0.;
  dbref->addSelection(sel, "sel");
  Db* db;

  mestitle(0, "Removing the bad sample");
  db = dbref->clone();
  db->deleteSample(indNa);
  db->display(&dbfmt);
  like = logLikelihoodVecchia(db, model, nbVecchia, verbose);
  message("Likelihood = %lf\n", like);
  delete db;

  mestitle(0, "With Selection");
  db = dbref->clone();
  db->display(&dbfmt);
  like = logLikelihoodVecchia(db, model, nbVecchia, verbose);
  message("Likelihood = %lf\n", like);
  delete db;

  mestitle(0, "Without filtering");
  db = dbref->clone();
  db->clearSelection();
  db->display(&dbfmt);
  like = logLikelihoodVecchia(db, model, nbVecchia, verbose);
  message("Likelihood = %lf\n", like);
  delete db;

  return 0;
}
