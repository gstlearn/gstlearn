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
 * This test aims to test the use of the nlopt library for optimization.
 * Only a quadratic function is minimized in this example.
 *
 * The example contains two subsets:
 * 1) Find the minimum of a simple quadratic function
 * 2) On a numerical field simulated on a Data set, calculate an experimental
 *    variogram and retrieve the parameters of the Model
 */

#include "Basic/VectorNumT.hpp"
#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"
#include "Model/ModelOptimVario.hpp"
#include "Model/ModelOptimLikelihood.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

#include "geoslib_define.h"

#include <iostream>
#include <nlopt.h>

// Function to minimize
double myfunc(unsigned n, const double *x, double *grad, void *my_func_data = nullptr)
{
    DECLARE_UNUSED(n);
    DECLARE_UNUSED(my_func_data);
    if (grad) {
        grad[0] = 2 * (x[0] - 3);
    }
    double value = (x[0] - 3) * (x[0] - 3);
    // std::cout << "current value = " << x[0] << " -> Minimum = " << value
    //           << std::endl;
    return value;
}

static void _firstTest()
{
  mestitle(0,"Minimization of a simple function");
  int npar = 1;
  VectorDouble x = {1.};
  nlopt_opt opt = nlopt_create(NLOPT_LN_NELDERMEAD, npar);

  // Bounds for each parameter
  VectorDouble lb = {1., 10.};
  nlopt_set_lower_bounds(opt, lb.data());
  VectorDouble ub = {5., 10.};
  nlopt_set_upper_bounds(opt, ub.data());

  nlopt_set_min_objective(opt, myfunc, nullptr);

  // Set the tolerance for the stopping criteria
  nlopt_set_ftol_rel(opt, EPSILON4);

  // Minimization
  double minf    = 1.e30;
  try
  {
    nlopt_optimize(opt, x.data(), &minf);
    std::cout << "Optimum: x = " << x[0] << " -> Minimum value = " << minf << std::endl;
  }
  catch (std::exception& e)
  {
    std::cerr << "Error in the optimization: " << e.what() << std::endl;
  }
}

static void _secondTest(Db* db, Model* model, bool converge)
{
  mestitle(0, "Fitting a Model from a Variogram");

  // Calculating the experimental variogram
  double hmax = db->getExtensionDiagonal();
  int npas = 10;
  double dpas = hmax / 2. / npas;
  VarioParam* varioparam = VarioParam::createOmniDirection(npas, dpas);
  Vario* vario           = Vario::computeFromDb(*varioparam, db);
  (void)vario->dumpToNF("vario2.ascii");

  // Fit the Model
  ModelOptimVario model_opt(model);
  model_opt.fit(vario, false, 2, converge);
  (void) model->dumpToNF("model2.ascii");
  model->display();

  delete varioparam;
  delete vario;
}

static void _thirdTest(Db* db, Model* model, bool flagSPDE, bool converge)
{
  if (flagSPDE)
    mestitle(0, "Fitting a Model using Loglikelihood (SPDE)");
  else
    mestitle(0, "Fitting a Model using Loglikelihood (Covariance)");

  // Fit the Model
  ModelOptimLikelihood model_opt(model);
  model_opt.fit(db, flagSPDE, converge);
  (void)model->dumpToNF("model3.ascii");
  model->display();
}

int main(int argc, char *argv[])
{
    std::stringstream sfn;
    sfn << gslBaseName(__FILE__) << ".out";
    StdoutRedirect sr(sfn.str(), argc, argv);
    ASerializable::setContainerName(true);
    ASerializable::setPrefixName("NlOpt-");

       // Creating the Model used to simulate the Data
    double sill_nugget = 2.;
    Model* model_simu  = new Model();
    model_simu->addCovFromParam(ECov::NUGGET, 0., sill_nugget);
    double range1 = 0.25;
    double sill1  = 3.;
    double param1 = 1.;
    model_simu->addCovFromParam(ECov::MATERN, range1, sill1, param1);
    message("Model used for simulating the Data\n");
    model_simu->display();

    // Data set
    int nech = 100;
    Db* db   = Db::createFromBox(nech, {0., 0.}, {1., 1.});
    (void)simtub(nullptr, db, model_simu);
    (void)db->dumpToNF("db.ascii");

    // Creating the testing Model
    Model* model_test  = new Model();
    model_test->addCovFromParam(ECov::NUGGET);
    model_test->addCovFromParam(ECov::MATERN);
    // model_test->setDriftIRF(0);
    message("Model used for Test\n");
    model_test->display();

    // Optimization tests
    int mode      = 0;
    bool converge = true;
    bool flagSPDE = false;

    if (mode == 0 || mode == 1) _firstTest();

    if (mode == 0 || mode == 2) _secondTest(db, model_test, converge);

    if (mode == 0 || mode == 3) _thirdTest(db, model_test, flagSPDE, converge);

    delete db;
    delete model_simu;
    return 0;
}
