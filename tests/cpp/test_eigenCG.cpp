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
#include "LinearOp/MatrixReplacement.hpp"
#include "LinearOp/IdentityEigenCG.hpp"
#include "LinearOp/LinearOpEigenCGSolver.hpp"
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
  VectorDouble x(n, 0.);
  VectorEigen X(x);

  IdentityEigenCG I(n);

  std::cout << "b = " << b.toString() << std::endl;

  I.evalInverse(b, x);
  std::cout << "x = " << x << std::endl;

  I.evalInverse(B, X);
  std::cout << "X = " << X << std::endl;

  LinearOpEigenCGSolver<IdentityEigenCG> s(&I);
  VectorEigen out(n);
  s.solve(B, out);
  std::cout << "X' = " << out<< std::endl;
}
