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
#include "Drifts/Drift1.hpp"
#include "Drifts/DriftX.hpp"
#include "Drifts/DriftY.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/OptCustom.hpp"
#include "Neigh/ANeighParam.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighMoving.hpp"

static Db* createLocalDb(int nech, int ndim, int nvar)
{
  // Coordinates
  VectorDouble tab = ut_vector_simulate_uniform(ndim * nech, 0., 50.);
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
  return data;
}

/**
 * Creating internal Model
 * @param nvar    Number of variables
 * @param typecov 1 for Spherical + Nugget; 2 for Linear
 * @param typedrift 0 None, 1 Linear drift
 * @param typemean 1 for defining a constant Mean
 * @return
 */
static Model* createModel(int nvar, int typecov, int typedrift, int typemean)
{
  CovContext ctxt(nvar); // use default space
  Model* model = Model::create(ctxt);
  CovLMC covs(ctxt.getSpace());

  if (typecov == 1)
  {
    CovAniso cova1(ECov::SPHERICAL, 40., 0., 45., ctxt);
    covs.addCov(&cova1);
    CovAniso cova2(ECov::NUGGET, 0., 0., 12., ctxt);
    covs.addCov(&cova2);
    model->setCovList(&covs);

  }
  else if (typecov == 2)
  {
    CovAniso covaL(ECov::LINEAR, 1., 0., 1., ctxt);
    covs.addCov(&covaL);
    model->setCovList(&covs);
  }

  if (typedrift == 1)
  {
    Drift1 drift1 = Drift1(ctxt);
    model->addDrift(&drift1);
  }
  else if (typedrift == 2)
  {
    Drift1 drift1 = Drift1(ctxt);
    DriftX driftx = DriftX(ctxt);
    DriftY drifty = DriftY(ctxt);
    model->addDrift(&drift1);
    model->addDrift(&driftx);
    model->addDrift(&drifty);
  }

  if (typemean == 1)
  {
    model->setMean(0, 123.);
  }
  return model;
}

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

  DbGrid* grid_res;
  DbGrid* image_res;
  Db* data_res;
  Db* data_short_res;

  // Global parameters
  int ndim = 2;
  int nvar = 1;
  law_set_random_seed(32131);

  setup_license("Demonstration");
  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  DbStringFormat dbfmt(FLAG_STATS);
  DbStringFormat dbfmtXvalid(FLAG_STATS,{"Xvalid*"});
  DbStringFormat dbfmtKriging(FLAG_STATS,{"Krig*"});
  DbStringFormat dbfmtImage(FLAG_STATS,{"Filter*"});
  DbStringFormat dbfmtBayes(FLAG_STATS,{"Bayes*"});
  DbStringFormat dbfmtSimu(FLAG_STATS,{"Sim*"});

  // Generate the output grid
  VectorInt nx = {50,50};
  DbGrid* grid = DbGrid::create(nx);
  grid->addColumnsByConstant(2, 1., "Extend", ELoc::BLEX);
  grid->display();

  // Generate the data base
  int nech = 100;
  Db* data = createLocalDb(nech, ndim, nvar);
  data->display(&dbfmt);

  // Create the Model
  Model* model1 = createModel(nvar, 1, 1, 1);
  model1->display();

  // ====================== Moving Neighborhood case ===========================
  // Creating a Moving Neighborhood
  NeighMoving* neighM = NeighMoving::create(ndim, false, 25);
  neighM->display();

  message("\n<----- Cross-Validation in Moving Neighborhood ----->\n");
  data_res = dynamic_cast<Db*>(data->clone());
  xvalid(data_res, model1, neighM, 0, -1, -1);
  data_res->display(&dbfmtXvalid);

  message("\n<----- Kriging in Moving Neighborhood ----->\n");
  grid_res = dynamic_cast<DbGrid*>(grid->clone());
  kriging(data, grid_res, model1, neighM);
  grid_res->display(&dbfmtKriging);

  message("\n<----- Kriging Test in Moving Neighborhood ----->\n");
  grid_res = dynamic_cast<DbGrid*>(grid->clone());
  Krigtest_Res ktest = krigtest(data, grid_res, model1, neighM, 0);
  message("\nTesting KrigTest facility\n");
  message("- Space Dimension = %d\n",ktest.ndim);
  message("- Number of Neighbors = %d\n",ktest.nech);
  message("- Number of Kriging System equations = %d\n",ktest.neq);
  ut_ivector_display("- Neighboring Sample Indices", ktest.nbgh);

  // ====================== Unique Neighborhood case ===========================
  // Unique Neighborhood
  NeighUnique* neighU = NeighUnique::create(ndim,false);
  neighU->display();

  message("\n<----- Cross-Validation in Unique Neighborhood ----->\n");
  data_res = dynamic_cast<Db*>(data->clone());
  xvalid(data_res, model1, neighU, 0, -1, -1);
  data->display(&dbfmtXvalid);

  message("\n<----- Kriging in Unique Neighborhood ----->\n");
  grid_res = dynamic_cast<DbGrid*>(grid->clone());
  kriging(data, grid_res, model1, neighU);
  grid_res->display(&dbfmtKriging);

  message("\n<----- Simulations in Unique Neighborhood ----->\n");
  grid_res = dynamic_cast<DbGrid*>(grid->clone());
  simtub(data, grid_res, model1, neighU, 3, 12345);
  grid_res->display(&dbfmtSimu);

  message("\n<----- Bayesian Simulations in Unique Neighborhood ----->\n");
  grid_res = dynamic_cast<DbGrid*>(grid->clone());
  simbayes(data, grid_res, model1, neighU, 3, 12345);
  grid_res->display(&dbfmtSimu);

  message("\n<----- Global Estimate (average) ----->\n");
  grid_res = dynamic_cast<DbGrid*>(grid->clone());
  Global_Res gres = global_arithmetic(data, grid_res, model1, 0, true);

  message("\n<----- Global Estimate (Kriging) ----->\n");
  gres = global_kriging(data, grid_res, model1, 0, true);

  // ====================== Block Kriging case ===========================
  message("\n<----- Block Kriging (fixed size) ----->\n");
  grid_res = dynamic_cast<DbGrid*>(grid->clone());
  kriging(data, grid_res, model1, neighU, EKrigOpt::BLOCK, 1, 1, 0, {3,3});
  grid_res->display(&dbfmtKriging);

  message("\n<----- Block Kriging (variable size) ----->\n");
  grid_res = dynamic_cast<DbGrid*>(grid->clone());
  krigcell(data, grid_res, model1, neighU, 1, 1, {3,3});
  grid_res->display(&dbfmtKriging);

  // ====================== Image Neighborhood case ===========================
  // Generate the Image
  DbGrid* image = DbGrid::create(nx);
  (void) simtub(nullptr,image,model1);
  image->setLocator("Simu", ELoc::Z);
  image->display();

  // Image Neighborhood
  NeighImage* neighI = NeighImage::create(ndim, {2,2}, 2);
  neighI->display();

  // Modify the Model (for filtering)
  Model* model4 = dynamic_cast<Model*>(model1->clone());
  model4->setCovaFiltered(1, true);
  model4->display();

  message("\n<----- Image Filtering ----->\n");
  image_res = dynamic_cast<DbGrid*>(image->clone());
  krimage(image_res, model4, neighI);
  image_res->display(&dbfmtImage);

  // ====================== Testing Baysian Kriging ===========================
  // Create the Model
  Model* model2 = createModel(nvar, 2, 1, 2);
  model2->display();

  message("\n<----- Bayesian Kriging in Unique Neighborhood ----->\n");
  grid_res = dynamic_cast<DbGrid*>(grid->clone());
  OptDbg::define(EDbg::BAYES);
  kribayes(data, grid_res, model2, neighU);
  OptDbg::undefine(EDbg::BAYES);
  grid_res->display(&dbfmtBayes);

  // ====================== Testing Printout ==================================
  message("\n<----- Test Kriging Printout ----->\n");

  // Create the Local Data Base
  Db* data_short = createLocalDb(4, 2, 1);

  // Create the Local Model
  Model* model3 = createModel(nvar, 2, 1, 0);
  model3->display();

  message("\n---> Kriging in Place (checking Exact Interpolator)\n");
  OptDbg::setReference(1);
  data_short_res = dynamic_cast<Db*>(data_short->clone());
  kriging(data_short_res, data_short_res, model3, neighU);
  OptDbg::setReference(0);

  message("\n---> Kriging in general\n");
  OptDbg::setReference(1);
  kriging(data_short, grid_res, model3, neighU);
  OptDbg::setReference(0);

  delete neighM;
  delete neighU;
  delete neighI;
  delete data;
  delete data_res;
  delete data_short;
  delete grid;
  delete grid_res;
  delete image;
  delete image_res;
  delete model1;
  delete model2;
  delete model3;
  delete model4;

  return (0);
}
