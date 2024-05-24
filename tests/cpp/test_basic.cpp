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
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Model/Model.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Variogram/VarioParam.hpp"
#include "Variogram/Vario.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Basic/Law.hpp"

/**
 * This file is used as a demonstration showing the versatility of 'gstlearn'.
 * It also exists in RMarkdown and Jupyter-notebook formats.
 */
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  defineDefaultSpace(ESpaceType::RN, 2);
  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Basic-");
  law_set_old_style(true); // Added to ensure the similarity of non-regression tests per platform

  // We create a grid of 150 by 100 square cells of 1m edge, called mygrid.

  VectorInt nx = { 100, 150 };
  VectorDouble dx = { 1., 1. };
  DbGrid* mygrid = DbGrid::create(nx, dx);

  // We create a Geostatistical Model constituted of a single Spherical anisotropic structure
  // with a sill of 1, a shortest range of 30m and a longest one of 50m.
  // The orientation of the long range is in direction 30 degrees counted counter-clockwise from East.

  double sill = 1.;
  VectorDouble ranges = { 50., 30.};
  VectorDouble angles = { 30., 0.};
  Model *mymodel = Model::createFromParam(ECov::SPHERICAL, 1., sill, 1., ranges,
                                          VectorDouble(), angles);
  mymodel->display();

  // We perform one non-conditional simulation using the Turning Band method
  // (the number of bands is set to 1000).
  // This results in adding a field (commonly referred to as a variable) in mygrid which is called Data.
  // Typing the name of the grid data base is an easy way to get a summary of its contents.

  int nbtuba = 1000;
  int seed = 14341;
  (void) simtub(nullptr, mygrid, mymodel, nullptr, 1, seed, nbtuba, false, false,
                NamingConvention("Data"));
  mygrid->display();

  // We sample the grid on a set of 100 samples randomly located within the area covered by the grid.
  // Typing the name of the newly created data base gives the list of variables stored in this data base.

  Db* mypoints = Db::createSamplingDb(mygrid, 0., 100);
  mypoints->display();

  // We calculate a variogram along the two main directions of the Model, i.e. 30 and 120 degrees.

  int npas = 25;
  double dpas = 2.;
  VarioParam* varioparam = VarioParam::createSeveral2D(angles, npas, dpas);
  Vario* myvario = Vario::computeFromDb(*varioparam, mypoints);
  myvario->display();

  // In the next step, we consider that the initial model (mymodel) has been correctly
  // rendered through the simulated information at the 100 samples,
  // Therefore we will use this initial model for estimation and conditional simulations.

  // For this subsequent steps, we need to define the neighborhood which will serve in designating
  // the set of samples, close to the target, which will be used for the estimation or simulation
  // at the target location (called myneigh).
  // Due to the small number of points, we decide to use the whole set of available samples,
  // wherever the target is located: this is known as a Unique neighborhood.

  auto myneigh = std::shared_ptr<ANeigh>(NeighUnique::create());

  // We have all the ingredients to perform the estimation using Kriging.
  // When typing the name of the (output) grid data base, we can check that 2 variables have been added:
  // - the one corresponding to the estimation result (called Kriging.Data.estim)
  // - the one corresponding to the standard deviation of the estimation error (called Kriging.Data.stdev)

  (void) kriging(mypoints, mygrid, mymodel, myneigh);
  mygrid->display();

  // We can now construct 2 conditional simulations, using the Turning Bands algorithm again,
  // but conditioned by the variable informed at the samples of mypoints.
  // The newly created variables in mygrid are called *Simu.Data.1" and "Simu.Data.2".
  // Each simulation outcome reproduces the spatial characteristics provided by the Model and
  // honors the information provided at sample points.

  (void) simtub(mypoints, mygrid, mymodel, myneigh, 2, seed, nbtuba);
  mygrid->display();

  // We can save the final version of the data base in Neutral File format

  mypoints->dumpToNF("Data");
  mygrid->dumpToNF("Grid");

  // Free the pointers
  delete mygrid;
  delete mymodel;
  delete mypoints;

  return 0;
}
