#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Basic/OptCst.hpp"
#include "Basic/ECst.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "API/SPDE.hpp"
#include "Model/Model.hpp"
#include "Model/NoStatArray.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Basic/FunctionalSpirale.hpp"
#include "geoslib_f.h"
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This program is meant to check the manipulation of the Db
 **
 *****************************************************************************/
int main(int /*argc*/, char */*argv*/[])
{
  // Standard output redirection to file
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  ASerializable::setContainerName(true);
  ASerializable::setPrefixName("TestDb-");
  int seed = 10355;
  law_set_random_seed(seed);

  // Creating the Grid Rotated Db
  DbGrid* grid = DbGrid::create({6,4}, {1.,2.}, {10.,20.}, {10.,0.});
  grid->display();

  // Creating the Model
  Model* model = Model::createFromDb(grid);
  model->display();
  CovContext ctxt = model->getContext();
  CovLMC covs(ctxt.getSpace());
  CovAniso cova = CovAniso(ECov::CUBIC,ctxt);
  cova.setRanges({10,45});
  cova.setAnisoAngles({30.,0.});
  covs.addCov(&cova);
  model->setCovList(&covs);
  model->display();

  // Creating the MeshTurbo which contains the Db
  MeshETurbo mesh;
  mesh.initFromCova(cova,grid,10,2,true,true);

  /////////////////////////
  // Testing the selections
  /////////////////////////

  int nech = grid->getSampleNumber();

  // First selection generated with Bernoulli (proba=0.6)
  VectorDouble sel1 = ut_vector_simulate_bernoulli(nech, 0.6);
  ut_vector_display("sel1", sel1);
  grid->addSelection(sel1, "Sel1");

  // Second selection generated with Bernoulli (proba=0.4) combined with previous one
  VectorDouble sel2 = ut_vector_simulate_bernoulli(nech, 0.4);
  ut_vector_display("sel2", sel2);
  grid->addSelection(sel2, "Sel2","and");

  // Retrieve resulting selection for check
  VectorDouble sel3 = grid->getSelection();
  ut_vector_display("sel1 && sel2",sel3);

  // Testing Filters on Db printout (only Statistics on the variables "Sel*")
  DbStringFormat dbfmt(FLAG_VARS | FLAG_STATS,{"Sel*"});
  grid->display(&dbfmt);

  // Creating a Selection by setting individual values
  OptCst::define(ECst::NTROW,-1);
  DbStringFormat dbfmt2(FLAG_VARS | FLAG_ARRAY);
  grid->addSelection(VectorDouble(), "mySel");
  grid->setValue("mySel", 12, 0.);
  grid->setValue("mySel", 19, 0.);
  grid->display(&dbfmt2);

  delete grid;
  delete model;
  return 0;
}

