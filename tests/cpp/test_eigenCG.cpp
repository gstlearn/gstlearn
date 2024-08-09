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
#include "Basic/File.hpp"
#include "LinearOp/ScaleOp.hpp"
#include "LinearOp/LinearOpCGSolver.hpp"
#include "Matrix/VectorEigen.hpp"

#include <Eigen/src/Core/Matrix.h>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  int n = 5;

  VectorDouble b({2, 4, 6, 8, 10});
  VectorEigen B(b);

  VectorDouble o(n);
  VectorEigen O(n);

  VectorDouble c(n);
  VectorEigen C(n);

  ScaleOp I(n, 2.0);
  LinearOpCGSolver<ScaleOp> s(&I);

  std::cout << "b = " << b.toString() << std::endl;

  s.solve(b, o);
  std::cout << "o = " << o.toString() << std::endl;

  s.solve(B, O);
  std::cout << "O = " << O << std::endl;

  I.evalDirect(o, c);
  std::cout << "c = " << c.toString() << std::endl;

  I.evalDirect(O, C);
  std::cout << "C = " << C << std::endl;
}
