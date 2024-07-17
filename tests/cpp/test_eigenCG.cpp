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
#include "LinearOp/Cholesky.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Basic/VectorHelper.hpp"

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

  MatrixSparse s(n, n, 0); // use flagEigen = 0 while issue https://github.com/gstlearn/gstlearn/issues/227
  s.fillRandom(123456, 0.5);
  s.prodMatInPlace(&s, true); // Make s definite positive
  MatrixSparse* S = createFromAnyMatrix(&s, 1); // Create eigen sparse (big 'S')
 
  VectorDouble b({2, 4, 6, 8, 10});
  Eigen::Map<const Eigen::VectorXd> B(b.data(), b.size());

  VectorDouble x(n);
  Eigen::VectorXd X;

  MatrixReplacement M;
  M.attachMyMatrix(S); // Use Eigen sparse matrix
 
  IdentityEigenCG I(n);

  std::cout << "S = " << S->toString() << std::endl;
  std::cout << "B = " << B << std::endl;
/*
  // Solve Ax = B using various methods and iterative solver with matrix-free version:
  {
    S->solve(b, x);
    std::cout << "solve:   #iterations: N/A, estimated error: N/A" << std::endl;
    std::cout << "x = ";
    for (auto v : x) std::cout << v << std::endl;
  }

  {
    Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> cg;
    cg.compute(M);
    X = cg.solve(B);
    std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
    std::cout << "X = " << X << std::endl;
  }
*/
/*
  {
    Eigen::ConjugateGradient<Cholesky, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> cg2;
    cg2.compute(C);
    X = cg2.solve(B);
    std::cout << "CG2:      #iterations: " << cg2.iterations() << ", estimated error: " << cg2.error() << std::endl;
    std::cout << "X = " << X << std::endl;
  }

  {
    Eigen::ConjugateGradient<CholeskyEigenCG, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> cg2b;
    cg2b.compute(C);
    X = cg2b.solve(B);
    std::cout << "CG2b:      #iterations: " << cg2b.iterations() << ", estimated error: " << cg2b.error() << std::endl;
    std::cout << "X = " << X << std::endl;
  }
*/
  {
    x.fill(0);
    I.evalInverse(b, x);
    std::cout << "x = " << x << std::endl;
  }
 /*
  {
    Eigen::BiCGSTAB<MatrixReplacement, Eigen::IdentityPreconditioner> bicg;
    bicg.compute(M);
    X = bicg.solve(B);
    std::cout << "BiCGSTAB: #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << std::endl;
    std::cout << "X = " << X << std::endl;
  }

  {
    Eigen::GMRES<MatrixReplacement, Eigen::IdentityPreconditioner> gmres;
    gmres.compute(M);
    X = gmres.solve(B);
    std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
    std::cout << "X = " << X << std::endl;
  }
 
  {
    Eigen::DGMRES<MatrixReplacement, Eigen::IdentityPreconditioner> gmres;
    gmres.compute(M);
    X = gmres.solve(B);
    std::cout << "DGMRES:   #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
    std::cout << "X = " << X << std::endl;
  }
 
  {
    Eigen::MINRES<MatrixReplacement, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> minres;
    minres.compute(M);
    X = minres.solve(B);
    std::cout << "MINRES:   #iterations: " << minres.iterations() << ", estimated error: " << minres.error() << std::endl;
    std::cout << "X = " << X << std::endl;
  }
*/
  delete S;
}
