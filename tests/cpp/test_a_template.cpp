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
  timer.displayIntervalMilliseconds("Avec les iterateurs begin() et end()");

  // Avec les crochets
  total   = 0.;
  for (int i1 = 0; i1 < n1; i1++)
  {
    vv1 = v1[i1];
    for (int i2 = 0; i2 < n2; i2++)
    {
      vv2 = v2[i2];
      total += vv1 * vv2;
    }
  }
  message("Valeur obtenue = %lf\n", total);
  timer.displayIntervalMilliseconds("Avec les VectorDouble et les operateurs crochets");

  // Avec les pointeurs
  const double* ptr1_s = v1.data();
  const double* ptr2_s = v2.data();
  const double* ptr1;
  const double* ptr2;
  total                = 0.;
  ptr1 = ptr1_s;
  for (int i1 = 0; i1 < n1; i1++)
  {
    vv1 = *ptr1++;
    ptr2 = ptr2_s;
    for (int i2 = 0; i2 < n2; i2++)
    {
      vv2 = *ptr2++;
      total += vv1 * vv2;
    }
  }
  message("Valeur otenue = %lf\n", total);
  timer.displayIntervalMilliseconds("Avec les pointeurs sur les data() des vecteur doubles");

  // Avec les 'auto'
  total = 0.;
  for (const auto& vv1: v1) {
    for (const auto& vv2: v2) {
      total += vv1 * vv2;
    }
  }
  message("Valeur otenue = %lf\n", total);
  timer.displayIntervalMilliseconds("Avec les auto");

  // Avec span
  timer.reset();
  total           = 0.;
  const auto v1sp = std::span(v1);
  const auto v2sp = std::span(v2);
  for (const auto vv1: v1sp) {
    for (const auto vv2: v2sp) {
      total += vv1 * vv2;
    }
  }
  message("Valeur otenue = %lf\n", total);
  timer.displayIntervalMilliseconds("Avec span");

  // auto + getVector
  total = 0.;
  for (const auto vv1: v1.getVector()) {
    for (const auto vv2: v2.getVector()) {
      total += vv1 * vv2;
    }
  }
  message("Valeur otenue = %lf\n", total);
  timer.displayIntervalMilliseconds("Avec auto + getVector");

  // Auto + iterateur
  total          = 0.;
  const auto beg = v2.cbegin();
  const auto end = v2.cend();
  for (const auto vv1: v1) {
    for (auto it = beg; it != end; ++it) {
      total += vv1 * *it;
    }
  }
  message("Valeur otenue = %lf\n", total);
  timer.displayIntervalMilliseconds("Avec auto + const iterator");

  // Avec data() + size()
  total           = 0.;
  const auto v1s  = v1.size();
  const auto* v1d = v1.data();
  const auto v2s  = v2.size();
  const auto* v2d = v2.data();
  for (size_t i = 0; i < v1s; ++i) {
    for (size_t j = 0; j < v2s; ++j) {
      total += v1d[i] * v2d[j];
    }
  }
  message("Valeur otenue = %lf\n", total);
  timer.displayIntervalMilliseconds("Avec data() + size()");

  return (0);
}
