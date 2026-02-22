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
#include "Basic/ASerializable.hpp"
#include "Basic/VectorHelper.hpp"
#include "geoslib_define.h"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  auto nech = 5;
  auto V1   = VH::simulateGaussian(nech, 0., 1.);
  V1.dump("Vector V1");
  auto V2 = VH::simulateGaussian(nech, 0., 1.);
  V2.dump("Vector V2");
  VectorDouble V4;

  auto V3 = VH::add(V1, V2);
  V3.dump("Checking VH::add(V1,V2)");

  VH::add(V4, V1, V2);
  V4.dump("Checking VH::add(v4,V1,V2)");

  auto V5 = V1 + V2;
  V5.dump("Checking V = V1 + V2");

  auto V6 = V1;
  V6 += V2;
  V6.dump("Checking V6 += V2");

  auto V7 = V1 + V2 + V1;
  V7.dump("Checking V = V1 + V2 + V1");

  auto V8 = V1 + 3.1;
  V8.dump("Checking V = V1 + 3.1");

  auto V9 = 2.1 + V1;
  V9.dump("Checking V = 2.1 + V1");

  auto V10 = 2.1 + V1 + 2.3;
  V10.dump("Checking V = 2.1 + V1 + 2.3");
  return 0;
}
