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

#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/Optim.hpp"
#include "Basic/VectorNumT.hpp"

#include "geoslib_define.h"

#include <iostream>

// Function to minimize
double myfunc2(const std::vector<double>& x)
{
  // if (grad)
  //   grad[0] = 2 * (x[0] - 3);

  double value = (x[0] - 3) * (x[0] - 3);
  // std::cout << "current value = " << x[0] << " -> Minimum = " << value
  //           << std::endl;
  return value;
}

static void _firstTest()
{
  mestitle(0,"Minimization of a Function");
  int npar = 1;
  std::vector<double> x = {1.};
  Optim* opt     = new Optim(NELDERMEAD, npar);

  // Bounds for each parameter
  VectorDouble lb = {1., 10.};
  opt->setLowerBounds(lb);
  VectorDouble ub = {5., 10.};
  opt->setUpperBounds(ub);
  auto func = [](const std::vector<double>& x) { return myfunc2(x);};
  opt -> setObjective(func);
  opt->setXtolRel(EPSILON4);
  double minf = opt->optimize(x);
  std::cout << "Optimum: x = " << x[0] << " -> Minimum value = " << minf << std::endl;
}

int main(int argc, char *argv[])
{
    std::stringstream sfn;
    sfn << gslBaseName(__FILE__) << ".out";
    StdoutRedirect sr(sfn.str(), argc, argv);
    ASerializable::setPrefixName("NlOpt-");

    // Optimization tests

   _firstTest();

    return 0;
}
