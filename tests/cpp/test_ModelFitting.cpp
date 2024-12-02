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
 * - Fitting sills only using new technology (AModelOptimSills)
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
#include "Model/ModelOptimSillsVario.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

#include "geoslib_define.h"
#include "geoslib_old_f.h"

static Vario* _computeVariogram(Db* db2D, const ECalcVario& calcul)
{
  double hmax            = db2D->getExtensionDiagonal();
  int npas               = 10;
  double dpas            = hmax / 2. / npas;
  VarioParam* varioparam = VarioParam::createOmniDirection(npas, dpas);
  Vario* vario           = Vario::computeFromDb(*varioparam, db2D, calcul);
  (void)vario->dumpToNF("vario.ascii");
  delete varioparam;
  return vario;
}

// Old style Sill fitting (Goulard) from Variogram
static void _firstTest(Db* db2D, Model* model, const ECalcVario& calcul)
{
  mestitle(0, "Using Goulard");

  // Calculating the experimental variogram
  Vario* vario = _computeVariogram(db2D, calcul);

  Constraints constraints;

  // Fit the Model
  (void) model_fitting_sills(vario, model, constraints);
  (void)model->dumpToNF("Goulard.ascii");
  model->display();

  delete vario;
}

// New style Sill fitting (Goulard) from Variogram
static void _secondTest(Db* db2D, Model* model, const ECalcVario& calcul)
{
  mestitle(0, "Using ModelOptimSillsVario");

  // Calculating the experimental variogram
  Vario* vario = _computeVariogram(db2D, calcul);

  // Fit the Model
  ModelOptimSillsVario model_opt(model);
  model_opt.fit(vario);
  (void) model->dumpToNF("ModelOptimSillsVario.ascii");
  model->display();

  delete vario;
}

// New style Sill fitting (Goulard) from VMap
static void _thirdTest(DbGrid* dbgrid, Model* model, const ECalcVario& calcul)
{
  mestitle(0, "Using Variogram Map");

  // Calculating the experimental variogram
  DbGrid* dbmap = db_vmap(dbgrid, calcul, {50,50});
  (void) dbmap->dumpToNF("VMap.ascii");

  // Fit the Model
  Option_VarioFit* optvar  = new Option_VarioFit();
  optvar->setFlagGoulardUsed(false);
  ModelOptimVMap model_opt(model);
  model_opt.fit(dbmap, true);
  (void)model->dumpToNF("ModelVMap.ascii");
  model->display();

  delete dbmap;
  delete optvar;
}

static VectorDouble _buildSillMatrix(int nvar, double value)
{
  VectorDouble mat(nvar * nvar);
  int ecr = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++, ecr++)
      mat[ecr] = (ivar == jvar) ? value : 0.;
  return mat;
}
  
int main(int argc, char *argv[])
{
    std::stringstream sfn;
    sfn << gslBaseName(__FILE__) << ".out";
    StdoutRedirect sr(sfn.str(), argc, argv);
    ASerializable::setContainerName(true);
    ASerializable::setPrefixName("ModelFit-");

    // Global parameters
    int nvar = 1;
    const ECalcVario calcul = ECalcVario::VARIOGRAM;
    bool verbose = false;

    // Creating the Model used to simulate the Data
    VectorDouble sill_nugget = _buildSillMatrix(nvar, 2.);
    Model* model_simu  = new Model(nvar);
    model_simu->addCovFromParam(ECov::NUGGET, 0., 0., 0.,
                                 VectorDouble(), sill_nugget);

    double range1 = 0.25;
    double param1 = 1.;
    VectorDouble sill1 = _buildSillMatrix(nvar, 3.);
    model_simu->addCovFromParam(ECov::SPHERICAL, range1, 0., param1,
                                VectorDouble(), sill1);
    if (verbose)
    {
      message("Model used for simulating the Data\n");
      model_simu->display();
    }

    // Data set
    int nech = 100;
    Db* db2D   = Db::createFromBox(nech, {0., 0.}, {1., 1.});
    (void)simtub(nullptr, db2D, model_simu);
    (void)db2D->dumpToNF("db.ascii");

    // Grid Data set
    int nx = 100;
    double dx = 1. / nx;
    DbGrid* dbgrid = DbGrid::create({nx, nx}, {dx, dx});
    (void)simtub(nullptr, dbgrid, model_simu);
    (void)dbgrid->dumpToNF("dbgrid.ascii");

    // Creating the testing Model
    Model* model_test  = new Model(nvar);
    model_test->addCovFromParam(ECov::NUGGET);
    model_test->addCovFromParam(ECov::SPHERICAL);
    if (verbose)
    {
      message("Model used for Test\n");
      model_test->display();
    }

    // Optimization tests
    int mode = 0;

    if (mode == 0 || mode == 1) _firstTest(db2D, model_test, calcul);

    if (mode == 0 || mode == 2) _secondTest(db2D, model_test, calcul);

    if (mode == 0 || mode == 3) _thirdTest(dbgrid, model_test, calcul);

    delete db2D;
    delete dbgrid;
    delete model_simu;
    return 0;
}
