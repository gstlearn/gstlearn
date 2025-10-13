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
#include "Basic/File.hpp"
#include "Basic/Law.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "LinearOp/ProjConvolution.hpp"
#include "Model/Model.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

using namespace gstlrn;

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  ASerializable::setPrefixName("test_ProjConvolution-");
  Id seed = 10355;
  Id ndim = 2;
  law_set_random_seed(seed);
  defineDefaultSpace(ESpaceType::RN, ndim);

  ///////////////////////
  // Creating the Db
  Id ngrid          = 101;
  Id nxval          = 101;
  VectorInt nx      = {nxval, ngrid};
  DbGrid* grid_data = DbGrid::create(nx);

  Id conv_dim  = 11;
  double range = 3.;
  double total = 0.;
  VectorDouble convolution(conv_dim);
  for (Id i = 0; i < conv_dim; i++)
  {
    double dist    = (i - conv_dim / 2.) / range;
    convolution[i] = exp(-dist * dist);
    total += convolution[i];
  }
  for (Id i = 0; i < conv_dim; i++)
    convolution[i] /= total;
  VH::dump("Convolution", convolution);

  Id ngrid_seismic     = ngrid - (conv_dim - 1);
  nx                   = VectorInt({nxval, ngrid_seismic});
  DbGrid* grid_seismic = DbGrid::create(nx);

  // Creating the Model
  Model* model = Model::createFromParam(ECov::SPHERICAL, 10.);

  // Simulating on the initial grid
  (void)simtub(nullptr, grid_data, model);
  auto uid_in = grid_data->getLastUID();

  // Save the initial grid in a NF
  (void)grid_data->dumpToNF("Initial.NF");

  // Operate the convolution
  ProjConvolution Aproj(convolution, grid_seismic);

  // Extract the simulation from the input Grid
  VectorDouble data = grid_data->getColumnByUID(uid_in);

  // Perform the convolution
  VectorDouble seismic(grid_seismic->getNSample());
  if (Aproj.mesh2point(data, seismic)) return 1;
  grid_seismic->addColumns(seismic, "Seismic", ELoc::Z);

  // Save the final grid in a NF
  (void)grid_seismic->dumpToNF("Seismic.NF");

  // Perform convolution back-transform
  VectorDouble data2(grid_data->getNSample());
  if (Aproj.point2mesh(seismic, data2)) return 1;
  grid_data->addColumns(data2, "Data", ELoc::Z);

  // Save the final grid in a NF
  (void)grid_data->dumpToNF("Point.NF");

  delete grid_data;
  delete grid_seismic;
  delete model;
  return 0;
}
