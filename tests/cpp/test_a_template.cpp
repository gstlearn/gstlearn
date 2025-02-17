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
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"

#include "Basic/File.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Timer.hpp"

/**
 * This file is meant to perform any test that needs to be coded for a quick trial
 * It will be compiled but not run nor diff'ed.
 */
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Testing the optimal way to program
  int n1 = 10000;
  int n2 = 10000;
  VectorDouble v1 = VH::simulateGaussian(n1);
  VectorDouble v2 = VH::simulateGaussian(n2);

  double vv1;
  double vv2;
  double total;

  Timer timer;

  // Avec les iterateurs
  timer.reset();
  VectorDouble::const_iterator iv1(v1.begin());

  total = 0.;
  while (iv1 < v1.end())
  {
    vv1 = *iv1;
    iv1++;

    VectorDouble::const_iterator iv2(v2.begin());
    while (iv2 < v2.end())
    {
      vv2 = *iv2;
      iv2++;

      total += vv1 * vv2;
    }
  }
  message("Valeur obtenue = %lf\n", total);
  timer.displayIntervalMilliseconds("Avec les iterateurs");

  // Avec les crochets
  int nn1 = (int)v1.size();
  int nn2 = (int)v2.size();
  total   = 0.;
  for (int i1 = 0; i1 < nn1; i1++)
  {
    vv1 = v1[i1];
    for (int i2 = 0; i2 < nn2; i2++)
    {
      vv2 = v2[i2];
      total += vv1 * vv2;
    }
  }
  message("Valeur obtenue = %lf\n", total);
  timer.displayIntervalMilliseconds("Avec les crochets");

  // Avec les pointeurs
  int np1            = (int)v1.size();
  int np2            = (int)v2.size();
  const double* ptr1_s = v1.data();
  const double* ptr2_s = v2.data();
  const double* ptr1;
  const double* ptr2;
  total                = 0.;
  ptr1 = ptr1_s;
  for (int i1 = 0; i1 < np1; i1++)
  {
    vv1 = (*ptr1);
    ptr2 = ptr2_s;
    for (int i2 = 0; i2 < np2; i2++)
    {
      vv2 = (*ptr2);
      total += vv1 * vv2;
      ptr2++;
    }
    ptr1++;
  }
  message("Valeur otenue = %lf\n", total);
  timer.displayIntervalMilliseconds("Avec les pointeurs");

  // Avec les 'auto'
  total = 0.;
  for (const auto& vv1: v1)
  {
    for (const auto& vv2: v2)
    {
      total += vv1 * vv2;
    }
  }
  message("Valeur otenue = %lf\n", total);
  timer.displayIntervalMilliseconds("Avec les auto");

  return (0);
}
