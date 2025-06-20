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

/**
 * This test aims to test the various possibilities of Model Fitting
 *
 * The example contains several subsets:
 * - Fitting sills only using old technology (Goulard)
 * - Fitting sills only using new technology (AModelFitSills)
 * - Fitting a Model using the old technology (Foxleg based)
 * - Fitting using the new technology (AModelOptim*)
 *
 * As options, you can define:
 * - nvar: to run a Monovariate or a Multivariate (independent variables) case
 * - calcul: tu run the fitting in variogram or in covariance
 */

#include "Basic/VectorNumT.hpp"
#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VMap.hpp"
#include "Variogram/VarioParam.hpp"
#include "Model/ModelOptimVario.hpp"
#include "Model/ModelOptimVMap.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

#include "geoslib_define.h"
#include "geoslib_old_f.h"

static Vario* _computeVariogram(Db* db2D, const ECalcVario& calcul)
{
  double hmax            = db2D->getExtensionDiagonal();
  int nlag               = 10;
  double dlag            = hmax / 2. / nlag;
  VarioParam* varioparam = VarioParam::createOmniDirection(nlag, dlag);
  Vario* vario           = Vario::computeFromDb(*varioparam, db2D, calcul);
  (void)vario->dumpToNF("Vario.ascii");
  delete varioparam;
  return vario;
}

static void _firstTest(Db* db2D, Model* model, const ECalcVario& calcul, bool verbose)
{
  DECLARE_UNUSED(verbose);
  mestitle(0, "Sill fitting from Variogram (old version)");

  Vario* vario = _computeVariogram(db2D, calcul);

  (void) model_fitting_sills(vario, model);
  (void)model->dumpToNF("Model_Sills.ascii");
  model->display();

  delete vario;
}

static void _secondTest(Db* db2D,
                        Model* model,
                        const ECalcVario& calcul,
                        bool verbose,
                        bool trace)
{
  mestitle(0, "Model fitting from Variogram (new version)");

  Vario* vario = _computeVariogram(db2D, calcul);
  ModelOptimParam mop = ModelOptimParam();
  mop.setWmode(2);
  mop.setFlagGoulard(false);
  model->fitNew(nullptr, vario, nullptr, nullptr, mop, ITEST,
                verbose, trace);
  (void)model->dumpToNF("Model_Vario.ascii");
  model->display();

  delete vario;
}

static void _thirdTest(DbGrid* dbgrid,
                       Model* model,
                       const ECalcVario& calcul,
                       bool verbose,
                       bool trace)
{
  mestitle(0, "Sill Fitting from Variogram Map (new version)");

  DbGrid* dbmap = db_vmap(dbgrid, calcul, {50,50});
  (void) dbmap->dumpToNF("VMap.ascii");

  ModelOptimParam mop = ModelOptimParam();
  mop.setWmode(2);
  mop.setFlagGoulard(false);
  model->fitNew(nullptr, nullptr, dbmap, nullptr, mop, ITEST,
                verbose, trace);
  (void)model->dumpToNF("Model_VMap.ascii");
  model->display();

  delete dbmap;
}

static MatrixSymmetric _buildSillMatrix(int nvar, double value)
{
  MatrixSymmetric mat(nvar);
  for (int ivar = 0; ivar < nvar; ivar++) mat.setValue(ivar, ivar, value);
  return mat;
}

int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("ModelFit-");

  // Global parameters
  int nvar                = 1;
  const ECalcVario calcul = ECalcVario::VARIOGRAM;

  // Creating the Model used to simulate the Data
  Model* model_simu           = new Model(nvar);
  MatrixSymmetric sill_nugget = _buildSillMatrix(nvar, 2.);
  model_simu->addCovFromParam(ECov::NUGGET, 0., 0., 0., VectorDouble(),
                              sill_nugget);

  double range1         = 0.25;
  double param1         = 1.;
  MatrixSymmetric sill1 = _buildSillMatrix(nvar, 3.);
  model_simu->addCovFromParam(ECov::SPHERICAL, range1, 0., param1,
                              VectorDouble(), sill1);
  message("Model used for simulating the Data\n");
  model_simu->display();
  model_simu->dumpToNF("Model_Ref.ascii");

  // Data set
  int nech = 100;
  Db* db2D = Db::createFromBox(nech, {0., 0.}, {1., 1.});
  (void)simtub(nullptr, db2D, model_simu);
  (void)db2D->dumpToNF("db.ascii");

  // Grid Data set
  int nx         = 100;
  double dx      = 1. / nx;
  DbGrid* dbgrid = DbGrid::create({nx, nx}, {dx, dx});
  (void)simtub(nullptr, dbgrid, model_simu);
  (void)dbgrid->dumpToNF("dbgrid.ascii");

  // Creating the testing Model
  message("Initial Model used for Test\n");
  model_simu->display();

  // Optimization tests
  int mode     = 2;
  bool verbose = true;
  bool trace = true;
  Model* model_test;

  if (mode == 0 || mode == 1)
  {
    model_test = model_simu->clone();
    _firstTest(db2D, model_test, calcul, verbose);
  }
  if (mode == 0 || mode == 2)
  {
    model_test = model_simu->clone();
    _secondTest(db2D, model_test, calcul, verbose, trace);
  }
  if (mode == 0 || mode == 3)
  {
    model_test = model_simu->clone();
    _thirdTest(dbgrid, model_test, calcul, verbose, trace);
  }

  delete db2D;
  delete dbgrid;
  delete model_simu;
  return 0;
}
