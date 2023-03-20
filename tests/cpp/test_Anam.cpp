/******************************************************************************/
/* gstlearn C++ Library                                                       */
/*                                                                            */
/* Authors: <authors>                                                         */
/*                                                                            */
/* License: BSD 3 Clause                                                      */
/******************************************************************************/
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

#include "Simulation/CalcSimuTurningBands.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "Basic/File.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/Law.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Anam-");
  int seed = 10355;
  law_set_random_seed(seed);

  OptCst::defineByKey("ASP",0);

  ///////////////////////
  // Creating the Db
  VectorInt nx_S = {100,100};
  VectorDouble dx_S = {0.01, 0.01};
  double dx_P = 0.25;
  int nx_B = 5;
  double dx_B = dx_P / nx_B;

  // Generate initial grid
  DbGrid* grid = DbGrid::create(nx_S, dx_S);
  grid->display();
  // Create grid of Panels covering the simulated area
  DbGrid* panel = DbGrid::createCoveringDb(grid, VectorInt(), {dx_P,dx_P});
  panel->display();
  // Discretization with a grid of blocks which cover the simulated area
  DbGrid* blocs = DbGrid::createCoveringDb(grid, VectorInt(), {dx_B,dx_B});
  blocs->display();

  // Simulation of the Gaussian variable
  Model* model_init = Model::createFromParam(ECov::EXPONENTIAL, 0.1, 1.);
  (void) simtub(NULL, grid, model_init);
  grid->setName("Simu", "Y");

  // Nonlinear transform (lognormal)
  double m_Z = 1.5;
  double s_Z = 0.5;
  VectorDouble Y = grid->getColumn("Y");
  VectorDouble Z(Y.size());
  for (int i = 0; i < (int) Y.size(); i++)
    Z[i] = m_Z * exp(s_Z * Y[i] - 0.5*s_Z*s_Z);
  grid->addColumns(Z, "Z");

  // Extracting a subset
  int np = 500;
  Db* data = Db::createSamplingDb(grid, 0., np, {"x1","x2","Y","Z"});
  data->setLocator("Z", ELoc::Z);
  data->display();

  // Gaussian Anamorphosis
  AAnam* anam = AnamHermite::create(20);
  anam->fit(data, "Z");
  anam->display();

  return 0;
}

