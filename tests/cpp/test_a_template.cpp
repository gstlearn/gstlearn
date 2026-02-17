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
#include "geoslib_define.h"
#include "Space/SpaceRN.hpp"
#include "Space/ASpaceObject.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/CorAniso.hpp"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  // Do not remove
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-"); // Here set the test name

  setDefaultSpace(SpaceRN::create(3));
  std::cout << "Default space dimension: " << getDefaultSpace()->getNDim() << std::endl;

  CovContext ctxt(1, 1);
  CorAniso* cov1 = CorAniso::create(ctxt, ECov::MATERN, {0.5}, {1.});
  std::cout << "CorAniso dimension: " << cov1->getNDim() << std::endl;
  std::cout << "eval0 = " << cov1->eval0() << std::endl;
}
