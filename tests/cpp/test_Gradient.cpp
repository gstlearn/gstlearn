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
#include "API/SPDE.hpp"
#include "Basic/File.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/OptDbg.hpp"
#include "Covariances/CovContext.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Enum/ECst.hpp"
#include "Enum/ESpaceType.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Estimation/CalcKrigingGradient.hpp"
#include "LinearOp/PrecisionOpMatrix.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Space/ASpaceObject.hpp"
#include "geoslib_old_f.h"

using namespace gstlrn;

/****************************************************************************/
/*!
** This test is meant to check the elaboration of the CovGradientGenric class
** and its use in Kriging (or Depth adn Gradient)
**
*****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  /***********************/
  /* 1 - Initializations */
  /***********************/
  Id ndim = 2;
  defineDefaultSpace(ESpaceType::RN, ndim);
  ASerializable::setPrefixName("test_Gradient-");
  bool flagForceNumeric;

  // Setup constants
  OptDbg::reset();
  OptCst::define(ECst::NTCAR, 10.);
  OptCst::define(ECst::NTDEC, 6.);

  // Define the data set
  Db* data = Db::create();
  data->addColumns({0., 1., 1., 2., 3.}, "x1", ELoc::X, 0);
  data->addColumns({0., 3., 0., 0., 0.}, "x2", ELoc::X, 1);
  data->addColumns({0., 0., 1., 2., 3.}, "z1", ELoc::Z, 0);
  data->addColumns({1., -1., 0., 1., 2.}, "g1", ELoc::G, 0);
  data->addColumns({1., -1., 0., 2., 1.}, "g2", ELoc::G, 1);
  DbStringFormat* dbfmt = DbStringFormat::create(FLAG_ARRAY);
  data->display(dbfmt);
  delete dbfmt;

  // Define the output Grid
  DbGrid* grid = DbGrid::create({4, 4});
  grid->display();

  // Define the (monovariate) Model
  Model* model = Model::createFromParam(ECov::GAUSSIAN, 2.);
  model->setDriftIRF(1);
  model->display();

  // Define the (unique) Neighborhood
  NeighUnique* neigh = NeighUnique::create();
  neigh->display();

  // Options
  OptDbg::setReference(1);
  double ballradius = 0.01;

  // Gradient with Analytic solution
  mestitle(0, "With Analytical Solution");
  flagForceNumeric = false;
  (void)krigingGradient(data, grid, model, neigh, true, true, ballradius, flagForceNumeric,
                        NamingConvention("KrigGradAnalytic"));

  // Gradient with Numeric solution
  mestitle(0, "With Numerical Solution");
  flagForceNumeric = true;
  (void)krigingGradient(data, grid, model, neigh, true, true, ballradius, flagForceNumeric,
                        NamingConvention("KrigGradNumeric"));

  // Free memory
  delete data;
  delete grid;
  delete model;
  delete neigh;

  return (0);
}
