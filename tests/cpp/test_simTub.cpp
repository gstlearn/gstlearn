/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_d.h"

#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Drifts/Drift1.hpp"
#include "Drifts/DriftX.hpp"
#include "Drifts/DriftY.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/OptCustom.hpp"
#include "Basic/VectorHelper.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"

static Db* createLocalDb(int nech, int ndim, int nvar)
{
  // Coordinates
  VectorDouble tab = VH::simulateUniform(ndim * nech, 0., 50.);
  // Variable
  for (int ivar=0; ivar<nvar; ivar++)
  {
    VectorDouble tabvar = VH::simulateGaussian(nech);
    tab.insert(tab.end(), tabvar.begin(), tabvar.end());
  }

  Db* data = Db::createFromSamples(nech,ELoadBy::COLUMN,tab);
  data->setNameByUID(1,"x1");
  data->setNameByUID(2,"x2");

  data->setLocatorByUID(1,ELoc::X,0);
  data->setLocatorByUID(2,ELoc::X,1);

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    data->setNameByUID(3+ivar,"Var");
    data->setLocatorByUID(3+ivar,ELoc::Z,ivar);
  }
  return data;
}

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
  ASerializable::setPrefixName("Simtub-");

  // Global parameters
  law_set_random_seed(32131);
  bool verbose = 1;
  int ndim = 2;
  int nvar = 1;
  int nbsimu = 3;
  DbGrid* grid_res;
  ASpaceObject::defineDefaultSpace(ESpaceType::SPACE_RN, ndim);
  DbStringFormat dbfmt(FLAG_STATS,{"Simu*"});

  // Generate the output grid
  VectorInt nx = {50,50};
  DbGrid* grid = DbGrid::create(nx);
  grid->display();

  // Generate the data base
  int nech = 100;
  Db* data = createLocalDb(nech, ndim, nvar);
  data->display(&dbfmt);

  // Create the Model
  Model* model = Model::createFromParam(ECov::SPHERICAL, 10., 1.);
  model->display();

  // Creating a Moving Neighborhood
  NeighMoving* neighM = NeighMoving::create(false, 25);
  neighM->display();

  // Unique Neighborhood
  NeighUnique* neighU = NeighUnique::create();
  neighU->display();

  // ====================== Simulation (turning bands) ====================
  message("\n<----- Simulation (Moving Neighborhood) ----->\n");
  grid_res = grid->clone();
  simtub(data, grid_res, model, neighM, nbsimu);
  grid_res->display(&dbfmt);
  (void) grid_res->dumpToNF("Moving.ascii",verbose);

  message("\n<----- Simulation (Unique Neighborhood) ----->\n");
  grid_res = grid->clone();
  simtub(data, grid_res, model, neighU, nbsimu);
  grid_res->display(&dbfmt);
  (void) grid_res->dumpToNF("Unique.ascii",verbose);

  // ====================== Free pointers ==================================
  if (neighM    != nullptr) delete neighM;
  if (neighU    != nullptr) delete neighU;
  if (data      != nullptr) delete data;
  if (grid      != nullptr) delete grid;
  if (grid_res  != nullptr) delete grid_res;
  if (model     != nullptr) delete model;

  return (0);
}
