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
#include "Basic/Law.hpp"
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

  mestitle(1, "Testing Double");
  auto V1 = VH::simulateGaussian(nech, 0., 1.);
  V1.dump("Vector V1", false);
  auto V2 = VH::simulateGaussian(nech, 0., 1.);
  V2.dump("Vector V2", false);
  VectorDouble Vres;

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

  mestitle(1, "Testing Integers");
  VectorInt IV1(nech);
  for (Id i = 0; i < nech; i++)
    IV1[i] = static_cast<Id>(10. * law_uniform());
  IV1.dump("Vector IV1", false);
  VectorInt IV2(nech);
  for (Id i = 0; i < nech; i++)
    IV2[i] = static_cast<Id>(10. * law_uniform());
  IV2.dump("Vector IV2", false);
  VectorInt IVres;

  IVres = VH::add(IV1, IV2);
  IVres.dump("Checking VH::add(IV1,IV2)", false);

  VH::add(IVres, IV1, IV2);
  IVres.dump("Checking VH::add(IVres,IV1,IV2)", false);

  IVres = IV1 + IV2;
  IVres.dump("Checking IVres = IV1 + IV2", false);

  IVres = IV1;
  IVres += IV2;
  IVres.dump("Checking IVres(IV1) += IV2", false);

  IVres = IV1 + IV2 + IV1;
  IVres.dump("Checking IVres = IV1 + IV2 + IV1", false);

  IVres = IV1 + 3;
  IVres.dump("Checking IVres = IV1 + 3", false);

  IVres = 2 + IV1;
  IVres.dump("Checking IVres = 2 + IV1", false);

  IVres = 2 + IV1 + 5;
  IVres.dump("Checking IVres = 2 + IV1 + 5", false);

  IVres = IV1 - IV2;
  IVres.dump("Checking IVres = IV1 - IV2", false);

  IVres = IV1 - 3;
  IVres.dump("Checking IVres = IV1 - 3", false);

  IVres = -3 + IV1;
  IVres.dump("Checking IVres = -3 + IV1", false);

  IVres = IV1 * IV2;
  IVres.dump("Checking IVres = IV1 * IV2", false);

  IVres = IV1 * 3;
  IVres.dump("Checking IVres = IV1 * 3", false);

  IVres = 2 * IV1;
  IVres.dump("Checking IVres = 2 * IV1", false);

  IVres = IV1 / IV2;
  IVres.dump("Checking IVres = IV1 / IV2", false);

  IVres = IV1 / 3;
  IVres.dump("Checking IVres = IV1 / 3", false);

  IVres = 2 / IV1;
  IVres.dump("Checking IVres = 2 / IV1", false);

  IVres = IV1;
  IVres /= IV2;
  IVres.dump("Checking IVres (IV1) /= IV2", false);

  IVres = IV1 + 4 * IV2;
  IVres.dump("Checking IVres = IV1 + 4 * IV2", false);

  return 0;
}
