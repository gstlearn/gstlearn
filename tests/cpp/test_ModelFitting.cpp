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
 * - Fitting sills only using new technology (ModelOptimSills)
 * - Fitting a Model using the old technology (Foxleg based)
 * - Fitting using the new technology (ModelOptim*)
 *
 * As options, you can define:
 * - nvar: to run a Monovariate or a Multivariate (independent variables) case
 * - calcul: tu run the fitting in variogram or in covariance
 */

#include "Basic/VectorNumT.hpp"
#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"
#include "Model/ModelOptimVario.hpp"
#include "Model/ModelOptimSills.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

#include "geoslib_define.h"
#include "geoslib_old_f.h"

static Vario* _computeVariogram(Db* db, const ECalcVario& calcul)
{
  double hmax            = db->getExtensionDiagonal();
  int npas               = 10;
  double dpas            = hmax / 2. / npas;
  VarioParam* varioparam = VarioParam::createOmniDirection(npas, dpas);
  Vario* vario           = Vario::computeFromDb(*varioparam, db, calcul);
  (void)vario->dumpToNF("vario.ascii");
  delete varioparam;
  return vario;
}

static void _firstTest(Db* db, Model* model, const ECalcVario& calcul)
{
  mestitle(0, "Using Goulard");

  // Calculating the experimental variogram
  Vario* vario = _computeVariogram(db, calcul);

  Constraints* constraints = new Constraints();

  // Fit the Model
  (void) model_fitting_sills(vario, model, *constraints);
  (void)model->dumpToNF("Goulard.ascii");
  model->display();

  delete vario;
  delete constraints;
}

static void _secondTest(Db* db, Model* model, const ECalcVario& calcul)
{
  mestitle(0, "Using ModelOptimSills");

  // Calculating the experimental variogram
  Vario* vario = _computeVariogram(db, calcul);

  // Fit the Model
  Constraints* constraints = new Constraints();
  Option_AutoFit* mauto = new Option_AutoFit();
  Option_VarioFit* optvar = new Option_VarioFit();
  ModelOptimSills model_opt;
  model_opt.fit(vario, model, constraints, mauto, optvar);
  (void) model->dumpToNF("ModelOptimSills.ascii");
  model->display();

  delete vario;
  delete constraints;
  delete mauto;
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
    int nvar = 2;
    const ECalcVario calcul = ECalcVario::VARIOGRAM;

    // Creating the Model used to simulate the Data
    VectorDouble sill_nugget = _buildSillMatrix(nvar, 2.);
    Model* model_simu  = new Model(nvar);
    model_simu->addCovFromParam(ECov::NUGGET, 0., 0., 0.,
                                 VectorDouble(), sill_nugget);

    double range1 = 0.25;
    double param1 = 1.;
    VectorDouble sill1 = _buildSillMatrix(nvar, 3.);
    model_simu->addCovFromParam(ECov::MATERN, range1, 0., param1,
                                 VectorDouble(), sill1);
    message("Model used for simulating the Data\n");
    model_simu->display();

    // Data set
    int nech = 100;
    Db* db1   = Db::createFromBox(nech, {0., 0.}, {1., 1.});
    (void)simtub(nullptr, db1, model_simu);
    (void)db1->dumpToNF("db.ascii");

    // Creating the testing Model
    Model* model_test  = new Model(nvar);
    model_test->addCovFromParam(ECov::NUGGET);
    model_test->addCovFromParam(ECov::MATERN);
    // model_test->setDriftIRF(0);
    message("Model used for Test\n");
    model_test->display();

    // Optimization tests
    int mode = 0;

    if (mode == 0 || mode == 1) _firstTest(db1, model_test, calcul);

    if (mode == 0 || mode == 2) _secondTest(db1, model_test, calcul);

    delete db1;
    delete model_simu;
    return 0;
}
