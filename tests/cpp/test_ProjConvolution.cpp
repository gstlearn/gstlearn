#include "Basic/Vector.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Db/DbGrid.hpp"
#include "Basic/Law.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Basic/File.hpp"
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
//  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("Convolution-");
  int seed = 10355;
  int ndim = 2;
  law_set_random_seed(seed);
  ASpaceObject::defineDefaultSpace(ESpaceType::SPACE_RN, ndim);

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
  ut_vector_display("Convolution", convolution);

  int ngrid_seismic = ngrid - (conv_dim - 1);
  nx = VectorInt({nxval, ngrid_seismic});
  DbGrid* grid_seismic = DbGrid::create(nx);

  ///////////////////////
  // Creating the Model
  Model* model = Model::createFromParam(ECov::SPHERICAL, 10.);

  // Simulating on the initial grid
  (void) simtub(nullptr, grid_data, model);
  int uid_in = grid_data->getLastUID();

  // Save the initial grid in a NF
  (void) grid_data->dumpToNF("Initial_Grid.ascii");

  // Operate the convolution
  ProjConvolution Aproj(convolution, grid_seismic);

  // Extract the simulation from the input Grid
  VectorDouble data = grid_data->getColumnByUID(uid_in);

  // Perform the convolution
  VectorDouble seismic = VectorDouble(grid_seismic->getSampleNumber());
  if (Aproj.mesh2point(data, seismic)) return 1;
  grid_seismic->addColumns(seismic);

  // Save the final grid in a NF
  (void) grid_seismic->dumpToNF("Final_Grid.ascii");

  // Perform convolution back-transform
  VectorDouble data2 = VectorDouble(grid_data->getSampleNumber());
  if (Aproj.point2mesh(seismic, data2)) return 1;
  grid_data->addColumns(data2);

  // Save the final grid in a NF
  (void) grid_data->dumpToNF("Final_Data.ascii");

  delete grid_data;
  delete grid_seismic;
  delete model;
  return 0;
}

