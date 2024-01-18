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
#include "geoslib_d.h"
#include "geoslib_f.h"

#include "Enum/ESpaceType.hpp"

#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Model/Model.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovLMC.hpp"
#include "Drifts/DriftM.hpp"
#include "Basic/Law.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/OptCustom.hpp"
#include "Basic/VectorHelper.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighImage.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamContinuous.hpp"
#include "Anamorphosis/CalcAnamTransform.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Estimation/CalcImage.hpp"
#include "Estimation/CalcGlobal.hpp"

static Db* createLocalDb(int nech, int ndim, int nvar, int seed)
{
  // Define seed
  law_set_random_seed(seed);

  // Coordinates
  VectorDouble tab = VH::simulateGaussian(ndim * nech, 0., 50.);
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
  else if (typecov == 3)
  {
    CovAniso cova1(ECov::SPHERICAL, 40., 0., 1., ctxt);
    covs.addCov(&cova1);
    model->setCovList(&covs);
  }

  if (typedrift == 1)
  {
    DriftM drift1 = DriftM();
    model->addDrift(&drift1);
  }
  else if (typedrift == 2)
  {
    DriftM drift1 = DriftM();
    DriftM driftx = DriftM(VectorInt({1}));
    DriftM drifty = DriftM({0,1});
    model->addDrift(&drift1);
    model->addDrift(&driftx);
    model->addDrift(&drifty);
  }

  if (typemean == 1)
  {
    model->setMean(123.);
  }
  return model;
}

/****************************************************************************/
/*!
 ** Main Program
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  DbGrid* grid_res  = nullptr;
  DbGrid* image_res = nullptr;
  Db* data_res      = nullptr;
  Model* model_res  = nullptr;
  AnamHermite* anam = nullptr;
  Global_Result gres;
  Krigtest_Res ktest;
  VectorDouble tab;

  // Global parameters
  int ndim = 2;
  int nvar = 1;
  law_set_random_seed(32131);

  defineDefaultSpace(ESpaceType::RN, ndim);
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
  Db* data = createLocalDb(nech, ndim, nvar, 342673);
  data->display(&dbfmt);

  // Generate an image file
  DbGrid* image = DbGrid::create(nx);

  // Create the Model
  Model* model = createModel(nvar, 1, 1, 1);
  model->display();

  // Image Neighborhood
  NeighImage* neighI = NeighImage::create({2,2}, 2);
  neighI->display();

  // Creating a Moving Neighborhood
  NeighMoving* neighM = NeighMoving::create(false, 25);
  neighM->display();

  // Unique Neighborhood
  NeighUnique* neighU = NeighUnique::create();
  neighU->display();

  // Block discretization
  VectorInt ndiscs = {3,3};

  // ====================== Testing Neighborhood Storage ===========================
  message("\n---> Testing Neighborhood storage\n");
  grid_res = grid->clone();
  test_neigh(data, grid_res, model, neighM);
  grid_res->display(&dbfmtKriging);

  // ====================== Moving Neighborhood case ===========================
  message("\n<----- Cross-Validation in Moving Neighborhood ----->\n");
  data_res = data->clone();
  xvalid(data_res, model, neighM, 0, -1, -1, 0);
  data_res->display(&dbfmtXvalid);

  message("\n<----- Kriging in Moving Neighborhood ----->\n");
  grid_res = grid->clone();
  kriging(data, grid_res, model, neighM);
  grid_res->display(&dbfmtKriging);

  message("\n<----- Declustering in Moving Neighborhood ----->\n");
  declustering(data_res, model, 3, neighM, grid, VectorDouble(), ndiscs, false, true);

  message("\n<----- Kriging Test in Moving Neighborhood ----->\n");
  grid_res = grid->clone();
  ktest = krigtest(data, grid_res, model, neighM, 0);
  message("\nTesting KrigTest facility\n");
  message("- Space Dimension = %d\n",ktest.ndim);
  message("- Number of Neighbors = %d\n",ktest.nech);
  message("- Number of Kriging System equations = %d\n",ktest.neq);
  VH::display("- Neighboring Sample Indices", ktest.nbgh);

  // ====================== Unique Neighborhood case ===========================
  message("\n<----- Cross-Validation in Unique Neighborhood ----->\n");
  data_res = data->clone();
  xvalid(data_res, model, neighU, 0, -1, -1, 0);
  data_res->display(&dbfmtXvalid);

  message("\n<----- Kriging in Unique Neighborhood ----->\n");
  grid_res = grid->clone();
  kriging(data, grid_res, model, neighU);
  grid_res->display(&dbfmtKriging);

  message("\n<----- Simulations in Unique Neighborhood ----->\n");
  grid_res = grid->clone();
  simtub(data, grid_res, model, neighU, 3, 12345);
  grid_res->display(&dbfmtSimu);

  message("\n<----- Bayesian Simulations in Unique Neighborhood ----->\n");
  grid_res = grid->clone();
  simbayes(data, grid_res, model, neighU, 3, 12345);
  grid_res->display(&dbfmtSimu);

  message("\n<----- Declustering in Unique Neighborhood ----->\n");
  data_res = data->clone();
  declustering(data_res, model, 2, neighU, nullptr, VectorDouble(), VectorInt(), false, true);

  message("\n<----- Global Estimate (Average) ----->\n");
  grid_res = grid->clone();
  gres = global_arithmetic(data, grid_res, model, 0, true);

  message("\n<----- Global Estimate (Kriging) ----->\n");
  gres = global_kriging(data, grid_res, model, 0, true);

  // ====================== Block Kriging case ===========================
  message("\n<----- Block Kriging (fixed size) ----->\n");
  grid_res = grid->clone();

  kriging(data, grid_res, model, neighU, EKrigOpt::BLOCK, 1, 1, 0, ndiscs);
  grid_res->display(&dbfmtKriging);

  message("\n<----- Block Kriging (variable size) ----->\n");
  grid_res = grid->clone();
  krigcell(data, grid_res, model, neighU, 1, 1, ndiscs);
  grid_res->display(&dbfmtKriging);

  // ====================== Image Neighborhood case ===========================
  // Generate the Image
  (void) simtub(nullptr,image,model);
  image->setLocator("Simu", ELoc::Z);
  image->display();

  // Modify the Model (for filtering)
  model_res = model->clone();
  model_res->setCovaFiltered(1, true);
  model_res->display();

  message("\n<----- Image Filtering ----->\n");
  image_res = image->clone();
  krimage(image_res, model_res, neighI);
  image_res->display(&dbfmtImage);

  // ====================== Testing Bayesian Kriging ===========================
  // Create the Model
  model = createModel(nvar, 2, 1, 2);
  model->display();

  message("\n<----- Bayesian Kriging in Unique Neighborhood ----->\n");
  grid_res = grid->clone();
  OptDbg::define(EDbg::BAYES);
  kribayes(data, grid_res, model, neighU);
  OptDbg::undefine(EDbg::BAYES);
  grid_res->display(&dbfmtBayes);

  // ====================== Testing Printout ==================================
  message("\n<----- Test Kriging Printout ----->\n");
  // Create the Local Data Base
  data = createLocalDb(4, 2, 1, 4291);

  // Create the Local Model
  model = createModel(nvar, 2, 1, 0);
  model->display();

  message("\n---> Kriging in Place (checking Exact Interpolator)\n");
  OptDbg::setReference(1);
  data_res = data->clone();
  kriging(data_res, data_res, model, neighU);
  OptDbg::setReference(0);

  message("\n---> Kriging in general\n");
  OptDbg::setReference(1);
  kriging(data, grid_res, model, neighU);
  OptDbg::setReference(0);

  // ====================== Testing Specials ==================================
  // Create the Local Data Base
  data = createLocalDb(10, 2, 3, 4901);

  message("\n<----- Test Kriging Multiple Variables under Constraints ----->\n");
  grid_res = grid->clone();
  tab = VH::simulateUniform(grid->getSampleNumber(), 10., 20.);
  grid_res->addColumns(tab, "Constraints", ELoc::SUM);
  krigsum(data, grid_res, model, neighU, true);
  grid_res->display(&dbfmtKriging);

  // Create the Local Data Base
  data = createLocalDb(100, 2, 1, 42781);

  // Create the Model
  model = createModel(1, 3, 0, 0);

  // Create the Gaussian
  anam = AnamHermite::create(20);
  anam->fitFromArray(data->getColumn("Var"));
  (void) anam->rawToGaussianByLocator(data);

  message("\n<----- Test Kriging Anamorphosed Gaussian ----->\n");
  grid_res = grid->clone();
  kriggam(data, grid_res, model, neighU, anam);
  grid_res->display(&dbfmtKriging);

  // ====================== Testing Multivariate===============================
  // Create the Local Data Base
  nvar = 3;
  data = createLocalDb(10, 2, 3, 4901);
  model = createModel(nvar, 1, 0, 0);

  message("\n<----- Test Kriging Multiple Variables with matCL ----->\n");
  grid_res = grid->clone();
  MatrixRectangular* matCL = MatrixRectangular::createFromVD({2., 2., 1., 1., 0., 1.}, 2, 3);
  kriging(data, grid_res, model, neighU, EKrigOpt::POINT, true, true, false,
          VectorInt(), VectorInt(), matCL);
  grid_res->display(&dbfmtKriging);

  // ====================== Free pointers ==================================
  if (neighM    != nullptr) delete neighM;
  if (neighU    != nullptr) delete neighU;
  if (neighI    != nullptr) delete neighI;
  if (data      != nullptr) delete data;
  if (data_res  != nullptr) delete data_res;
  if (grid      != nullptr) delete grid;
  if (grid_res  != nullptr) delete grid_res;
  if (image     != nullptr) delete image;
  if (image_res != nullptr) delete image_res;
  if (model     != nullptr) delete model;
  if (model_res != nullptr) delete model_res;
  if (anam      != nullptr) delete anam;

  return (0);
}
