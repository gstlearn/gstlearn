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
#include "Basic/NamingConvention.hpp"
#include "Db/DbGrid.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

/**
 * Test the NamingConvention class, particularly the new simulation naming functionality
 */
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Create a simple grid for testing
  VectorInt nx = {10, 10};
  VectorDouble dx = {1., 1.};
  DbGrid* grid = DbGrid::create(nx, dx);

  mestitle(0, "Testing NamingConvention for Multivariate Simulations");

  // Test 1: Non-conditional multivariate simulations (2 variables, 2 simulations)
  // Storage order: simulation-first (default)
  message(
    "\nTest 1: Non-conditional multivariate (nvar=2, nsim=2, simu-first)\n");
  {
    Id nvar = 2;
    Id nsim = 2;
    Id iatt = grid->addColumnsByConstant(nvar * nsim, 0.);

    NamingConvention namconv("Test1", false, false, false);
    namconv.setOutputForSimulations(
      VectorString(), nvar, grid, iatt, nsim, true, false);

    message("Variable names: ");
    for (Id i = 0; i < nvar * nsim; i++)
    {
      message("%s", grid->getNameByUID(iatt + i).c_str());
      if (i < nvar * nsim - 1) message(", ");
    }
    message("\n");
  }

  // Test 2: Conditional multivariate simulations (2 variables with names, 2 simulations)
  message("\nTest 2: Conditional multivariate (nvar=2, nsim=2, simu-first)\n");
  {
    Id nvar = 2;
    Id nsim = 2;
    Id iatt = grid->addColumnsByConstant(nvar * nsim, 0.);

    VectorString varnames = {"Fe", "Al"};
    NamingConvention namconv("Test2", false, false, false);
    namconv.setOutputForSimulations(
      varnames, nvar, grid, iatt, nsim, true, false);

    message("Variable names: ");
    for (Id i = 0; i < nvar * nsim; i++)
    {
      message("%s", grid->getNameByUID(iatt + i).c_str());
      if (i < nvar * nsim - 1) message(", ");
    }
    message("\n");
  }

  // Test 3: Variable-first storage order
  message("\nTest 3: Non-conditional (nvar=2, nsim=2, variable-first)\n");
  {
    Id nvar = 2;
    Id nsim = 2;
    Id iatt = grid->addColumnsByConstant(nvar * nsim, 0.);

    NamingConvention namconv("Test3", false, false, false);
    namconv.setOutputForSimulations(
      VectorString(), nvar, grid, iatt, nsim, false, false);

    message("Variable names: ");
    for (Id i = 0; i < nvar * nsim; i++)
    {
      message("%s", grid->getNameByUID(iatt + i).c_str());
      if (i < nvar * nsim - 1) message(", ");
    }
    message("\n");
  }

  // Test 4: Single simulation consistency (nsim=1 should still show S1)
  message("\nTest 4: Single simulation (nvar=2, nsim=1)\n");
  {
    Id nvar = 2;
    Id nsim = 1;
    Id iatt = grid->addColumnsByConstant(nvar * nsim, 0.);

    NamingConvention namconv("Test4", false, false, false);
    namconv.setOutputForSimulations(
      VectorString(), nvar, grid, iatt, nsim, true, false);

    message("Variable names: ");
    for (Id i = 0; i < nvar * nsim; i++)
    {
      message("%s", grid->getNameByUID(iatt + i).c_str());
      if (i < nvar * nsim - 1) message(", ");
    }
    message("\n");
  }

  // Test 5: Larger multivariate case (3 variables, 3 simulations)
  message("\nTest 5: Larger case (nvar=3, nsim=3, simu-first)\n");
  {
    Id nvar = 3;
    Id nsim = 3;
    Id iatt = grid->addColumnsByConstant(nvar * nsim, 0.);

    VectorString varnames = {"Cu", "Pb", "Zn"};
    NamingConvention namconv("Test5", false, false, false);
    namconv.setOutputForSimulations(
      varnames, nvar, grid, iatt, nsim, true, false);

    message("Variable names:\n");
    for (Id i = 0; i < nvar * nsim; i++)
    {
      message("  [%2d] %s\n", i, grid->getNameByUID(iatt + i).c_str());
    }
  }

  delete grid;
  return 0;
}
