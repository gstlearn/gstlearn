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
#include "Covariances/CovAniso.hpp"
#include "Matrix/MatrixRectangular.hpp"
#include "Basic/Timer.hpp"
#include <vector>
/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 * It will be compiled but not run nor diff'ed.
 */
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  MatrixRectangular mat(4,4);
  
  Timer timer;
  mat.reset(5, 5);
  std::vector<int> a = {200,400,800};
  timer.reset();
  for (int i = 0; i < 10000; i++)
  {
    mat.reset(a[i%3], a[i%3]);
  }
  timer.displayIntervalMilliseconds("Matrix resizing");

  timer.reset();
  for (int i = 0; i < 10000; i++)
  {
    mat.resize(a[i%3], a[i%3]);
  }
  timer.displayIntervalMilliseconds("Matrix resizing");

  return(0);
}
