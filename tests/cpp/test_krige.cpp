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
#include "geoslib_f.h"
#include "Space/Space.hpp"
#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighMoving.hpp"

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

  // Global parameters
  int ndim = 2;
  int nvar = 1;
  law_set_random_seed(32131);

  setup_license("Demonstration");
  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  DbStringFormat dbfmt(FLAG_STATS);

  // Generate the output grid
  VectorInt nx = {100,100};
  DbGrid* grid = DbGrid::create(nx);
  grid->display();

  // Generate the data base
  int nech = 100;
  // Coordinates
  VectorDouble tab = ut_vector_simulate_uniform(ndim * nech, 0., 100.);
  // Variable
  for (int ivar=0; ivar<nvar; ivar++)
  {
    VectorDouble tabvar = ut_vector_simulate_gaussian(nech, 0., 5);
    tab.insert(tab.end(), tabvar.begin(), tabvar.end());
  }

  Db* data = Db::createFromSamples(nech,ELoadBy::COLUMN,tab);
  data->setNameByUID(1,"x1");
  data->setNameByUID(2,"x2");
  data->setNameByUID(3,"Var");
  data->setLocatorByUID(1,ELoc::X,0);
  data->setLocatorByUID(2,ELoc::X,1);
  data->setLocatorByUID(3,ELoc::Z);
  data->display(&dbfmt);

  // Create the Model
  CovContext ctxt(nvar); // use default space
  Model* model = Model::create(ctxt);
  CovLMC covs(ctxt.getSpace());
  CovAniso cova1(ECov::SPHERICAL, 40., 0., 45., ctxt);
  covs.addCov(&cova1);
  CovAniso cova2(ECov::NUGGET, 0., 0., 12., ctxt);
  covs.addCov(&cova2);
  model->setCovList(&covs);
  model->setMean(0,123.);
  model->display();

  // ====================== Moving Neighborhood case ===========================
  mestitle(0,"Test in Moving Neighborhood");

  // Creating a Moving Neighborhood
  NeighMoving* neighM = NeighMoving::create(ndim, false, 100);
  neighM->display();

  // Launch Kriging
  data->setLocatorByUID(3,ELoc::Z);
  kriging(data, grid, model, neighM);
  grid->display(&dbfmt);

  // Launch Cross-Validation
  data->setLocatorByUID(3,ELoc::Z);
  xvalid(data, model, neighM, 0, -1, -1);
  data->display(&dbfmt);

  // ====================== Unique Neighborhood case ===========================
  mestitle(0,"Test in Unique Neighborhood");

  // Unique Neighborhrood
  NeighUnique* neighU = NeighUnique::create(ndim,false);

  // Launch Cross-Validation
  data->setLocatorByUID(3,ELoc::Z);
  xvalid(data, model, neighU, 0, -1, -1);
  data->display(&dbfmt);
  OptDbg::setReference(-1);

  // Launch Kriging
  data->setLocatorByUID(3,ELoc::Z);
  kriging(data, grid, model, neighU);
  grid->display(&dbfmt);

  // ====================== Image Neighborhood case ===========================
  mestitle(0,"Test in Image Neighborhood");

  // Generate the Image
  DbGrid* image = DbGrid::create(nx);
  (void) simtub(nullptr,image,model);
  image->display();

  // Image Neighborhood
  NeighImage* neighI = NeighImage::create(ndim, {2,2}, 2);
  neighI->display();

  // Modify the Model (for filtering)
  model->setCovaFiltered(1, true);
  model->display();

  // Perform Image filtering
  image->setLocator("Simu", ELoc::Z);
  krimage(image, model, neighI);
  image->display(&dbfmt);

  delete neighM;
  delete neighU;
  delete neighI;
  delete data;
  delete grid;
  delete image;
  delete model;

  return (0);
}
