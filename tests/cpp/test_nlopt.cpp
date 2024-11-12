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
#include "Db/DbStringFormat.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"
#include "Model/Model_Optim.hpp"
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
    return (x[0] - 3) * (x[0] - 3);
}

static void _firstTest()
{
  int npar      = 2;
  nlopt_opt opt = nlopt_create(NLOPT_LD_LBFGS, npar);

  // Bounds for each parameter
  VectorDouble lb = {1., 10.};
  nlopt_set_lower_bounds(opt, lb.data());
  VectorDouble ub = {5., 10.};
  nlopt_set_upper_bounds(opt, ub.data());

  nlopt_set_min_objective(opt, myfunc, nullptr);

  // Set the tolerance for the stopping criteria
  nlopt_set_xtol_rel(opt, EPSILON4);

  // Starting point
  double x[2] = {1., 5.0};

  double minf; // minimum value of the objective function
  try
  {
    nlopt_optimize(opt, x, &minf);
    std::cout << "End optimization" << std::endl;
    std::cout << "x = " << x[0] << " " << x[1] << std::endl;
    std::cout << "Minimum value of the objective " << minf << std::endl;
  }
  catch (std::exception& e)
  {
    std::cerr << "Error in the optimization: " << e.what() << std::endl;
  }
}

static void _secondTest()
{
  bool verbose = true;

  // Creating the Model
  double sill_nugget = 2.;
  Model* model = new Model();;
  model->addCovFromParam(ECov::NUGGET, 0., sill_nugget);
  double range_spherical = 0.2;
  double sill_spherical  = 3.;
  model->addCovFromParam(ECov::SPHERICAL, range_spherical, sill_spherical);
  if (verbose) model->display();

  // Creating the Data set
  int nech = 500.;
  Db* db   = Db::createFromBox(nech, {0., 0.}, {1., 1.});
  (void) simtub(nullptr, db, model);
  if (verbose)
  {
    DbStringFormat* dbfmt = DbStringFormat::createFromFlags(true, true, true);
    db->display(dbfmt);
  }

  // Calculating the experimental variogram
  int npas = 10;
  double dpas = 0.3 / npas;
  VarioParam* varioparam = VarioParam::createOmniDirection(npas, dpas);
  Vario* vario           = Vario::computeFromDb(*varioparam, db);
  (void)vario->dumpToNF("vario.ascii");
  if (verbose) vario->display();

  // Fit the Model
  Model_Optim model_opt(vario);
  model_opt.fit(model, 2);
  (void) model->dumpToNF("model.ascii");
  if (verbose) model->display();

  delete model;
  delete db;
  delete varioparam;
  delete vario;
}

int main(int argc, char *argv[])
{
    std::stringstream sfn;
    sfn << gslBaseName(__FILE__) << ".out";
    StdoutRedirect sr(sfn.str(), argc, argv);
    ASerializable::setContainerName(true);
    ASerializable::setPrefixName("NlOpt-");

    int mode = 2;

    if (mode == 0 || mode == 1) _firstTest();

    if (mode == 0 || mode == 2) _secondTest();
    
    return 0;
}
