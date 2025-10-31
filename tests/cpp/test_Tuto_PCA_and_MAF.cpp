/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/* This file is meant to test PCA and MAF features in a way similar to       */
/* the Tuto_PCA_and_MAF.Rmd script to diagnose memory issues with ASAN       */
/*                                                                            */
/******************************************************************************/
#include "Basic/AStringable.hpp"
#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Enum/ECov.hpp"
#include "Enum/EDbg.hpp"
#include "Model/Constraints.hpp"
#include "Model/Model.hpp"
#include "Model/Option_AutoFit.hpp"
#include "Model/Option_VarioFit.hpp"
#include "Simulation/CalcSimuTurningBands.hpp"
#include "Space/ASpaceObject.hpp"
#include "Stats/PCA.hpp"
#include "Variogram/Vario.hpp"
#include "Variogram/VarioParam.hpp"
#include "geoslib_define.h"
#include "utils.hpp"

using namespace gstlrn;

/****************************************************************************/
/*!
** Main Program for testing PCA and MAF similar to Tuto_PCA_and_MAF.Rmd
**
*****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  DECLARE_UNUSED(argc, argv)

  Id error = 0;
  Id ndim  = 2;
  Id nvar  = 3;

  defineDefaultSpace(ESpaceType::RN, ndim);

  // Grid of samples
  VectorInt nx_S    = {100, 100};
  VectorDouble dx_S = {0.01, 0.01};
  DbGrid* grid      = DbGrid::create(nx_S, dx_S);
  if (grid == nullptr)
  {
    messerr("Error creating grid");
    return 1;
  }

  mestitle(0, "Created Grid");
  grid->display();

  Id np = grid->getNSample();

  // Simulation of the Gaussian factors for structure #1 (Nugget Effect)
  mestitle(1, "Simulating Structure #1 (Nugget)");
  Model* m1 = Model::createFromParam(ECov::NUGGET, 0., 1.0);
  if (m1 == nullptr)
  {
    messerr("Error creating model 1");
    delete grid;
    return 1;
  }
  error = simtub(nullptr, grid, m1, nullptr, nvar, 432423, 100, false, false,
                 NamingConvention("U1"));
  if (error != 0)
  {
    messerr("Error in simtub for structure 1");
    delete grid;
    delete m1;
    return 1;
  }

  // Simulation of the Gaussian factors for structure #2 (Exponential)
  mestitle(1, "Simulating Structure #2 (Exponential)");
  Model* m2 = Model::createFromParam(ECov::EXPONENTIAL, 0.1, 1.0);
  if (m2 == nullptr)
  {
    messerr("Error creating model 2");
    delete grid;
    delete m1;
    return 1;
  }
  error = simtub(nullptr, grid, m2, nullptr, nvar, 432424, 100, false, false,
                 NamingConvention("U2"));
  if (error != 0)
  {
    messerr("Error in simtub for structure 2");
    delete grid;
    delete m1;
    delete m2;
    return 1;
  }

  // Simulation of the Gaussian factors for structure #3 (Cubic)
  mestitle(1, "Simulating Structure #3 (Cubic)");
  Model* m3 = Model::createFromParam(ECov::CUBIC, 0.25, 1.0);
  if (m3 == nullptr)
  {
    messerr("Error creating model 3");
    delete grid;
    delete m1;
    delete m2;
    return 1;
  }
  error = simtub(nullptr, grid, m3, nullptr, nvar, 432425, 100, false, false,
                 NamingConvention("U3"));
  if (error != 0)
  {
    messerr("Error in simtub for structure 3");
    delete grid;
    delete m1;
    delete m2;
    delete m3;
    return 1;
  }

  // Create correlated variables from the simulated factors
  // This is a simplified version - in reality, we'd need to apply correlation matrices
  mestitle(1, "Creating correlated variables");

  // For simplicity, we'll just copy and scale the simulated data
  // Z1, Z2, Z3 are created from linear combinations
  VectorDouble z1 = grid->getColumnByLocator(ELoc::Z, 0);
  VectorDouble z2 = grid->getColumnByLocator(ELoc::Z, 1);
  VectorDouble z3 = grid->getColumnByLocator(ELoc::Z, 2);

  // Simple linear combinations (simplified from the Rmd script)
  for (Id i = 0; i < np; i++)
  {
    z1[i] = 1.0 + 0.25 * z1[i];
    z2[i] = 2.0 + 3.0 * z2[i];
    z3[i] = 3.0 + 1.5 * z3[i];
  }

  grid->setColumn(z1, "Z1");
  grid->setColumn(z2, "Z2");
  grid->setColumn(z3, "Z3");

  // Data extraction - create sampling from grid
  mestitle(1, "Extracting data samples");
  Id npSamples       = 500;
  VectorString names = {"x1", "x2", "Z1", "Z2", "Z3"};
  Db data;
  data.resetSamplingDb(grid, 0., npSamples, names, 432426);

  data.setLocator("Z*", ELoc::Z, 0);
  mestitle(0, "Extracted Data");
  data.display();

  // Computing the experimental variogram
  mestitle(1, "Computing experimental variogram on raw data");
  Id nlag                = 10;
  double dlag            = 0.025;
  VarioParam* varioparam = VarioParam::createOmniDirection(nlag, dlag);
  if (varioparam == nullptr)
  {
    messerr("Error creating variogram parameters");
    delete grid;
    delete m1;
    delete m2;
    delete m3;
    return 1;
  }

  Vario vario_raw(*varioparam);
  vario_raw.compute(&data);

  // Fitting the variogram model on the experimental variogram
  mestitle(1, "Fitting variogram model on raw data");
  Model model_raw {};

  Constraints ctr;
  Option_VarioFit ovf;
  Option_AutoFit oaf;
  oaf.setVerbose(false);
  auto types = ECov::fromKeys({"NUGGET", "EXPONENTIAL", "CUBIC"});
  error      = model_raw.fit(&vario_raw, types, ctr, ovf, oaf, false);
  if (error != 0)
  {
    messerr("Error fitting model");
  }

  mestitle(0, "Fitted Model for Raw Data");
  model_raw.display();

  // ============
  // Evaluate PCA
  // ============
  mestitle(0, "Testing PCA");
  data.setLocator("Z*", ELoc::Z, 0);
  PCA pca(nvar);

  error = pca.pca_compute(&data, true);
  if (error != 0)
  {
    messerr("Error computing PCA");
  }
  pca.display();

  // Store the transformed variables
  error = pca.dbZ2F(&data, true, NamingConvention("U", false));
  if (error != 0)
  {
    messerr("Error transforming Z to PCA factors");
  }

  // Set locators for PCA factors
  data.setLocator("U*", ELoc::Z, 0);

  // Fitting the variogram model on PCA factors
  mestitle(1, "Computing and fitting variogram on PCA factors");
  Vario* vario_PCA = Vario::computeFromDb(*varioparam, &data);
  if (vario_PCA == nullptr)
  {
    messerr("Error computing PCA variogram");
  }
  else
  {
    vario_PCA->display();
    ctr.display();
    ovf.display();
    oaf.display();
    Model model_PCA {};
    error = model_PCA.fit(vario_PCA, types, ctr, ovf, oaf, true);
    if (error == 0)
    {
      mestitle(0, "Fitted Model for PCA");
      model_PCA.display();
    }

    delete vario_PCA;
  }

  // ============
  // Evaluate MAF
  // ============
  mestitle(0, "Testing MAF");
  data.setLocator("Z*", ELoc::Z, 0);
  PCA maf(nvar);

  // MAF computation using variogram at a specific lag
  Id ilag = 3; // lag index 3 (corresponding to ilag-1 in 0-based indexing)
  error   = maf.maf_compute(&data, *varioparam, ilag - 1, 0, true);
  if (error != 0)
  {
    messerr("Error computing MAF");
  }
  maf.display();

  // Store the transformed variables
  error = maf.dbZ2F(&data, true, NamingConvention("F", false));
  if (error != 0)
  {
    messerr("Error transforming Z to MAF factors");
  }

  // Set locators for MAF factors
  data.setLocator("F*", ELoc::Z, 0);

  // Fitting the variogram model on MAF factors
  mestitle(1, "Computing and fitting variogram on MAF factors");
  Vario* vario_MAF = Vario::computeFromDb(*varioparam, &data);
  if (vario_MAF == nullptr)
  {
    messerr("Error computing MAF variogram");
  }
  else
  {
    vario_MAF->display();
    ctr.display();
    ovf.display();
    oaf.display();
    Model model_MAF {};
    auto types_maf = ECov::fromKeys({"NUGGET", "EXPONENTIAL", "SPHERICAL"});
    error          = model_MAF.fit(vario_MAF, types_maf, ctr, ovf, oaf, true);
    if (error == 0)
    {
      mestitle(0, "Fitted Model for MAF");
      model_MAF.display();
    }
    delete vario_MAF;
  }

  // Cleanup
  delete grid;
  delete m1;
  delete m2;
  delete m3;
  delete varioparam;

  return static_cast<int>(error);
}
