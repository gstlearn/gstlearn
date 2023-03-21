/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/Law.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Basic/File.hpp"
#include "Basic/VectorHelper.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "LinearOp/ProjConvolution.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Convolution-");
  int seed = 10355;
  int ndim = 2;
  law_set_random_seed(seed);
  defineDefaultSpace(ESpaceType::RN, ndim);

  ///////////////////////
  // Creating the Db
  int ngrid = 101;
  int nxval = 101;
  VectorInt nx = { nxval, ngrid };
  DbGrid* grid_data = DbGrid::create(nx);

  int conv_dim = 11;
  double range = 3.;
  double total = 0.;
  VectorDouble convolution(conv_dim);
  for (int i = 0; i < conv_dim; i++)
  {
    double dist = (i - conv_dim / 2) / range;
    convolution[i] = exp(- dist * dist);
    total += convolution[i];
  }
  for (int i = 0; i < conv_dim; i++)
    convolution[i] /= total;
  VH::display("Convolution", convolution);

  int ngrid_seismic = ngrid - (conv_dim - 1);
  nx = VectorInt({nxval, ngrid_seismic});
  DbGrid* grid_seismic = DbGrid::create(nx);

  // Creating the Model
  Model* model = Model::createFromParam(ECov::SPHERICAL, 10.);

  // Simulating on the initial grid
  (void) simtub(nullptr, grid_data, model);
  int uid_in = grid_data->getLastUID();

  // Save the initial grid in a NF
  (void) grid_data->dumpToNF("Initial.ascii");

  // Operate the convolution
  ProjConvolution Aproj(convolution, grid_seismic);

  // Extract the simulation from the input Grid
  VectorDouble data = grid_data->getColumnByUID(uid_in);

  // Perform the convolution
  VectorDouble seismic = VectorDouble(grid_seismic->getSampleNumber());
  if (Aproj.mesh2point(data, seismic)) return 1;
  grid_seismic->addColumns(seismic,"Seismic",ELoc::Z);

  // Save the final grid in a NF
  (void) grid_seismic->dumpToNF("Seismic.ascii");

  // Perform convolution back-transform
  VectorDouble data2 = VectorDouble(grid_data->getSampleNumber());
  if (Aproj.point2mesh(seismic, data2)) return 1;
  grid_data->addColumns(data2,"Data",ELoc::Z);

  // Save the final grid in a NF
  (void) grid_data->dumpToNF("Point.ascii");

  delete grid_data;
  delete grid_seismic;
  delete model;
  return 0;
}

