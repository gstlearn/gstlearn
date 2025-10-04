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
#include "Basic/VectorT.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Estimation/Vecchia.hpp"
#include "Matrix/MatrixSymmetric.hpp"
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

  auto* model1 = Model::createFromParam(ECov::EXPONENTIAL);
  model1->setDriftIRF(0);
  MatrixSymmetric* sills = MatrixSymmetric::createRandomDefinitePositive(2);
  auto* model2           = Model::createFromParam(ECov::EXPONENTIAL, 1, 0., 1., VectorDouble(), *sills);
  model2->setDriftIRF(0);

  Id nbVecchia = 8;
  double like;
  bool verbose = false;
  int mode     = 0;

  Id nvar  = 1;
  Id dim   = 2;
  Id ndat  = 6;
  Id indNa = 3;
  VectorDouble sel(ndat, 1.);
  sel[indNa] = 0.;
  auto dbfmt = DbStringFormat(FLAG_VARS | FLAG_ARRAY, VectorString(), VectorInt(), false);

  // Constructing a Monovariate Db
  nvar       = 1;
  Db* dbref1 = Db::createFillRandom(ndat, dim, nvar, 0, 0, 0);
  dbref1->setValue("z", indNa, TEST);
  dbref1->addSelection(sel, "sel");

  // Constructing a bivariate Db
  nvar       = 2;
  Db* dbref2 = Db::createFillRandom(ndat, dim, nvar, 0, 0, 0);
  dbref2->setValue("z-1", indNa, TEST);
  dbref2->setValue("z-2", indNa, TEST);
  dbref2->addSelection(sel, "sel");

  if (mode == 0 || mode == 10 || mode == 11)
  {
    message("\n>>>>> Db1 when removing the bad sample\n");
    Db* db = dbref1->clone();
    db->deleteSample(indNa);
    if (verbose) db->display(&dbfmt);
    like = logLikelihoodVecchia(db, model1, nbVecchia, verbose);
    message("Likelihood = %lf\n", like);
    delete db;
  }

  if (mode == 0 || mode == 10 || mode == 12)
  {
    message("\n>>>>> Db1 with Selection\n");
    Db* db = dbref1->clone();
    if (verbose) db->display(&dbfmt);
    like = logLikelihoodVecchia(db, model1, nbVecchia, verbose);
    message("Likelihood = %lf\n", like);
    delete db;
  }

  if (mode == 0 || mode == 10 || mode == 13)
  {
    message("\n>>>>> Db1 with TEST values\n");
    Db* db = dbref1->clone();
    db->clearSelection();
    if (verbose) db->display(&dbfmt);
    like = logLikelihoodVecchia(db, model1, nbVecchia, verbose);
    message("Likelihood = %lf\n", like);
    delete db;
  }

  if (mode == 0 || mode == 10 || mode == 14)
  {
    message("\n>>>>> Db1 with TEST values using traditional LogLikelihood\n");
    Db* db = dbref1->clone();
    db->clearSelection();
    if (verbose) db->display(&dbfmt);
    like = -model1->computeLogLikelihood(db, verbose);
    message("Likelihood = %lf\n", like);
    delete db;
  }

  if (mode == 0 || mode == 20 || mode == 21)
  {
    message("\n>>>>> Db2 when removing the bad sample\n");
    Db* db = dbref2->clone();
    db->deleteSample(indNa);
    if (verbose) db->display(&dbfmt);
    like = logLikelihoodVecchia(db, model2, nbVecchia, verbose);
    message("Likelihood = %lf\n", like);
    delete db;
  }

  if (mode == 0 || mode == 20 || mode == 22)
  {
    message("\n>>>>> Db2 with Selection\n");
    Db* db = dbref2->clone();
    if (verbose) db->display(&dbfmt);
    like = logLikelihoodVecchia(db, model2, nbVecchia, verbose);
    message("Likelihood = %lf\n", like);
    delete db;
  }

  if (mode == 0 || mode == 20 || mode == 23)
  {
    message("\n>>>>> Db2 with TEST values\n");
    Db* db = dbref2->clone();
    db->clearSelection();
    if (verbose) db->display(&dbfmt);
    like = logLikelihoodVecchia(db, model2, nbVecchia, verbose);
    message("Likelihood = %lf\n", like);
    delete db;
  }

  if (mode == 0 || mode == 20 || mode == 24)
  {
    message("\n>>>>> Db2 with TEST values using traditional LogLikelihood\n");
    Db* db = dbref2->clone();
    db->clearSelection();
    if (verbose) db->display(&dbfmt);
    like = -model2->computeLogLikelihood(db, verbose);
    message("Likelihood = %lf\n", like);
    delete db;
  }
  return 0;
}
