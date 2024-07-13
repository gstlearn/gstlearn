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
#include "LinearOp/ALinearOpEigenCG.hpp"
#include "Matrix/MatrixSparse.hpp"

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

  int n = 10;

/*
  Eigen::SparseMatrix<double> S = Eigen::MatrixXd::Random(n,n).sparseView(0.5,1);
  S = S.transpose()*S;
 
  ALinearOpEigenCG A;
  A.attachMyMatrix(S);
*/
///*
  MatrixSparse s(n, n, 0); // use flagEigen = 0 while issue https://github.com/gstlearn/gstlearn/issues/227
  s.fillRandom(123456, 0.5);
  MatrixSparse* S = createFromAnyMatrix(&s, 1);
  S->prodMatInPlace(S, true); // Make S definite positive
 
  ALinearOpEigenCG A;
  A.attachMyMatrix(S);
 //*/

  Eigen::VectorXd b(n), x;
  b.setRandom();

  std::cout << "S = " << S->toString() << std::endl;
  std::cout << "b = " << b << std::endl;
 
  // Solve Ax = b using various iterative solver with matrix-free version:
  {
    Eigen::ConjugateGradient<ALinearOpEigenCG, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> cg;
    cg.compute(A);
    x = cg.solve(b);
    std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
    std::cout << "x = " << x << std::endl;
  }
 
  {
    Eigen::BiCGSTAB<ALinearOpEigenCG, Eigen::IdentityPreconditioner> bicg;
    bicg.compute(A);
    x = bicg.solve(b);
    std::cout << "BiCGSTAB: #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << std::endl;
    std::cout << "x = " << x << std::endl;
  }
 
  {
    Eigen::GMRES<ALinearOpEigenCG, Eigen::IdentityPreconditioner> gmres;
    gmres.compute(A);
    x = gmres.solve(b);
    std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
    std::cout << "x = " << x << std::endl;
  }
 
  {
    Eigen::DGMRES<ALinearOpEigenCG, Eigen::IdentityPreconditioner> gmres;
    gmres.compute(A);
    x = gmres.solve(b);
    std::cout << "DGMRES:   #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
    std::cout << "x = " << x << std::endl;
  }
 
  {
    Eigen::MINRES<ALinearOpEigenCG, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> minres;
    minres.compute(A);
    x = minres.solve(b);
    std::cout << "MINRES:   #iterations: " << minres.iterations() << ", estimated error: " << minres.error() << std::endl;
    std::cout << "x = " << x << std::endl;
  }

  delete S;
}
