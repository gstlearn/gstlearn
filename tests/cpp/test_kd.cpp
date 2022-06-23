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
#include "Neigh/NeighMoving.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamContinuous.hpp"

static Db* createLocalDb(int nech, int ndim, int nvar)
{
  // Coordinates
  VectorDouble tab = ut_vector_simulate_gaussian(ndim * nech, 0., 50.);
  // Variable
  for (int ivar=0; ivar<nvar; ivar++)
  {
    VectorDouble tabvar = ut_vector_simulate_gaussian(nech);
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
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str());

  DbGrid* grid_res  = nullptr;
  Db* data_res      = nullptr;
  Model* model_res  = nullptr;
  AnamHermite* anam = nullptr;
  Global_Res gres;
  Krigtest_Res ktest;
  VectorDouble tab;

  // Global parameters
  int ndim = 2;
  int nvar = 1;
  law_set_random_seed(32131);

  setup_license("Demonstration");
  ASpaceObject::defineDefaultSpace(SPACE_RN, ndim);
  DbStringFormat dbfmt(FLAG_STATS);
  DbStringFormat dbfmtXvalid(FLAG_STATS,{"Xvalid*"});
  DbStringFormat dbfmtKriging(FLAG_STATS,{"Krig*"});
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
  Model* model = createModel(nvar, 1, 1, 1);
  model->display();

  // Creating a Moving Neighborhood
  NeighMoving* neighM = NeighMoving::create(ndim, false, 25);
  neighM->display();

  // Create the Anamorphosis (and Gaussian Data)
  anam = AnamHermite::create(20);
  anam->fit(data->getColumn("Var"));
  anam->RawToGaussian(data, ELoc::Z);
  anam->display();
  data->display(&dbfmt);

  // ====================== KD =============================================

  message("\n<----- Test KD ----->\n");
  grid_res = dynamic_cast<DbGrid*>(grid->clone());
  int nfactor = 5;
  // Estimate Hermite polynomials at Data locations
  (void) calculateHermiteFactors(data, nfactor);
  model->addAnam(anam);
  dk(data, grid_res, model, neighM);
  grid_res->display(&dbfmtKriging);

  // ====================== Free pointers ==================================
  if (neighM    != nullptr) delete neighM;
  if (data      != nullptr) delete data;
  if (data_res  != nullptr) delete data_res;
  if (grid      != nullptr) delete grid;
  if (grid_res  != nullptr) delete grid_res;
  if (model     != nullptr) delete model;
  if (model_res != nullptr) delete model_res;
  if (anam      != nullptr) delete anam;

  return (0);
}
