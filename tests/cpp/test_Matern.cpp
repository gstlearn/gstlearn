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
#include "Covariances/CorMatern.hpp"
#include "Covariances/CovContext.hpp"
#include "Db/DbGrid.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Model/ModelGeneric.hpp"
#include "Simulation/ASimuSpectral.hpp"
#include "geoslib_define.h"
#include "Matrix/MatrixSymmetric.hpp"

using namespace gstlrn;

int main(int argc, char* argv[])
{
  DECLARE_UNUSED(argc);
  DECLARE_UNUSED(argv);

  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);
  ASerializable::setPrefixName("test_a_template-");

  VectorDouble ranges = {1.0, 2.0};
  VectorDouble angles = {30.0, 0.0};     // en degrés
  VectorDouble rr     = {1.0, 1.0, 2.0}; // on va multiplier les ranges par 1. et 2.
  VectorDouble params = {0.5, 1.0, 2.0};
  MatrixSymmetric gsSigma(3);
  gsSigma.setIdentity(1.);
  
  bool flagRange = true;
  CorMatern cor_tri(ranges,
                    angles,
                    rr,
                    params,
                    gsSigma,
                    flagRange);

  Id nvar            = 3;
  CovContext ctxt    = CovContext(nvar);
  ModelGeneric model = ModelGeneric(ctxt);
  model.setCov(&cor_tri);
  DbGrid* grid = DbGrid::create({100, 100});
  simuSpectral(nullptr, grid, &cor_tri, 1, 43431, 100, 100, nullptr, true, 
               NamingConvention("Simu"));

  return 0;
}
