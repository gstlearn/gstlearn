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
  V1.dump("Vector V1", false);
  auto V2 = VH::simulateGaussian(nech, 0., 1.);
  V2.dump("Vector V2", false);

  VectorDouble Vres;
  message("\n");

  Vres = VH::add(V1, V2);
  Vres.dump("Checking VH::add(V1,V2)", false);

  VH::add(Vres, V1, V2);
  Vres.dump("Checking VH::add(Vres,V1,V2)", false);

  Vres = V1 + V2;
  Vres.dump("Checking Vres = V1 + V2", false);

  Vres = V1;
  Vres += V2;
  Vres.dump("Checking Vres(V1) += V2", false);

  Vres = V1 + V2 + V1;
  Vres.dump("Checking Vres = V1 + V2 + V1", false);

  Vres = V1 + 3.1;
  Vres.dump("Checking Vres = V1 + 3.1", false);

  Vres = 2.1 + V1;
  Vres.dump("Checking Vres = 2.1 + V1", false);

  Vres = 2.1 + V1 + 2.3;
  Vres.dump("Checking Vres = 2.1 + V1 + 2.3", false);

  Vres = V1 - V2;
  Vres.dump("Checking Vres = V1 - V2", false);

  Vres = V1 - 3.1;
  Vres.dump("Checking Vres = V1 - 3.1", false);

  Vres = -3.1 + V1;
  Vres.dump("Checking Vres = -3.1 + V1", false);

  Vres = V1 * V2;
  Vres.dump("Checking Vres = V1 * V2", false);

  Vres = V1 * 3.1;
  Vres.dump("Checking Vres = V1 * 3.1", false);

  Vres = 2.1 * V1;
  Vres.dump("Checking Vres = 2.1 * V1", false);

  Vres = V1 / V2;
  Vres.dump("Checking Vres = V1 / V2", false);

  Vres = V1 / 3.1;
  Vres.dump("Checking Vres = V1 / 3.1", false);

  Vres = 2.1 / V1;
  Vres.dump("Checking Vres = 2.1 / V1", false);

  Vres = V1;
  Vres /= V2;
  Vres.dump("Checking Vres (V1) /= V2", false);

  Vres = V1 + 4. * V2;
  Vres.dump("Checking Vres = V1 + 4. * V2", false);

  Vres = (V1 + 4.) * V2;
  Vres.dump("Checking Vres = (V1 + 4.) * V2 -> NO", false);

  Vres = V1 + (4. * V2);
  Vres.dump("Checking Vres = V1 + (4. * V2)", false);

  Vres = V2 * 4. + V1;
  Vres.dump("Checking Vres = V2 * 4. + V1", false);
  return 0;
}
