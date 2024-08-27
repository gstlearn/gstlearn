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
#include "Basic/NamingConvention.hpp"
#include "Enum/ESpaceType.hpp"

#include "Matrix/MatrixRectangular.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Space/ASpaceObject.hpp"
#include "Db/Db.hpp"
#include "Db/DbStringFormat.hpp"
#include "Basic/Law.hpp"
#include "Model/Model.hpp"
#include "Basic/File.hpp"
#include "Basic/OptDbg.hpp"
#include "Basic/OptCst.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Estimation/KrigingCalcul.hpp"

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This test is composed of two parts:
 ** 1) Comparing the results of traditional Kriging with the results
 **    provided by the Algebraic calculations provided within 'KrigingCalcul'
 ** 2) Comparing the results for Collocated CoKriging or KFold Cross-validation
 **    in 'Unique' Neighborhood, whether they are programmed in the plain
 **    manner, or if they benefit from the inversion of the permanent part
**     of the Kriging System (performed a single time)
 **
 *****************************************************************************/
int main(int argc, char *argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Global parameters
  int ndim = 2;
  law_set_random_seed(32131);

  defineDefaultSpace(ESpaceType::RN, ndim);
  OptCst::define(ECst::NTCOL, -1);
  OptCst::define(ECst::NTROW, -1);
  OptDbg::define(EDbg::BAYES);

  // Parameters
  bool debugPrint = false;
  int nech        = 10;
  int nvar        = 2;
  int nfex        = 0;
  int nbfl        = (nfex + 1) * nvar;
  bool flagSK     = false;
  bool flagBayes  = false;
  if (flagBayes) flagSK = false;

  // Generate the data base
  Db* data = Db::createFillRandom(nech, ndim, nvar, nfex);
  DbStringFormat* dbfmt = DbStringFormat::createFromFlags(true, true, false, false, true);
  if (debugPrint) data->display(dbfmt);

  // Generate the target file
  Db* target = Db::createFillRandom(1, ndim, 0, nfex);
  if (debugPrint) target->display();

  // Create the constant Mean vector
  VectorDouble means(nvar, 0.);
  if (flagSK) means = VH::simulateGaussian(nvar);

  // Create the Model
  Model* model;
  double range = 0.2;
  MatrixSquareSymmetric sills(nvar);
  sills.setDiagonal(VH::simulateUniform(nvar, 1., 2.));
  model = Model::createFromParam(ECov::SPHERICAL, range, 0., 0., VectorDouble(), sills.getValues());
  model->setMeans(means);
  if (!flagSK) model->setDriftIRF(0, nfex);
  if (debugPrint) model->display();

  // Create the Bayesian Priors
  VectorDouble PriorMean = VH::simulateGaussian(nbfl);
  MatrixSquareSymmetric PriorCov(nbfl);
  PriorCov.setDiagonal(VH::simulateUniform(nbfl, 0.1, 0.5));

  // Unique Neighborhood
  NeighUnique* neigh = NeighUnique::create();
  if (debugPrint) neigh->display();

  // Define the verbose option
  int iech0 = 0;
  if (debugPrint) OptDbg::setReference(iech0 + 1);

  // ==============================
  // === First part of the test ===
  // ==============================

  // ---------------------- Using Standard Kriging procedure ---------------
  mestitle(0, "Using Standard Kriging procedure");
  Table table;
  if (!flagBayes)
  {
    kriging(data, target, model, neigh, EKrigOpt::POINT, true, true, true);
    table = target->printOneSample(iech0, {"Kriging*"}, true, true);
  }
  else
  {
    kribayes(data, target, model, neigh, PriorMean, PriorCov, true, true);
    table = target->printOneSample(iech0, {"Bayes*"}, true, true);
  }
  table.display();

  // ---------------------- Using Schur Class ------------------------------
  mestitle(0, "Using Schur class");
  MatrixSquareSymmetric C00   = model->getSillValues(0);
  MatrixSquareSymmetric Sigma = model->evalCovMatrixSymmetric(data);
  MatrixRectangular X         = model->evalDriftMatrix(data);
  MatrixRectangular Sigma0    = model->evalCovMatrix(data, target);
  MatrixRectangular X0        = model->evalDriftMatrix(target);
  VectorDouble Z = data->getMultipleValuesActive(VectorInt(), VectorInt(), means);

  KrigingCalcul Kcalc(Z, &Sigma, &X, &C00, means);
  Kcalc.setSigma0(&Sigma0);
  Kcalc.setX0(&X0);
  if (flagBayes) Kcalc.setBayes(PriorMean, &PriorCov);

  if (flagBayes)
  {
    VH::display("Prior Mean", PriorMean);
    message("Prior Variance-Covariance Matrix\n");
    PriorCov.display();
    VectorDouble beta = Kcalc.getPostMean();
    if (!beta.empty()) VH::display("Posterior Mean", beta);
    message("Posterior Variance-Covariance Matrix\n");
    Kcalc.getPostCov()->display();
  }

  VH::display("Kriging Value(s)", Kcalc.getEstimation());
  VH::display("Standard Deviation of Estimation Error", Kcalc.getStdv());
  VH::display("Variance of Estimator", Kcalc.getVarianceZstar());

  // ===============================
  // === Second part of the test ===
  // ===============================

  if (nvar > 1)
  {
    // Set of ranks of collocated variables
    VectorInt varsCol = {0};

    // Create the secondary Data base with the Target sample added as Data
    int iech = data->addSamples(1);
    for (int idim = 0; idim < ndim; idim++)
      data->setCoordinate(iech, idim, target->getCoordinate(0, idim));

    VectorDouble values = VH::simulateGaussian(nvar);
    for (int ivar = 0; ivar < (int)varsCol.size(); ivar++)
      values[varsCol[ivar]] = TEST;
    for (int ivar = 0; ivar < nvar; ivar++)
      data->setLocVariable(ELoc::Z, iech, ivar, values[ivar]);
    data->display(dbfmt);

    // ---------------------- With complemented Data Base ---------------------
    mestitle(0, "With complemented input Data Base");
    MatrixSquareSymmetric SigmaP = model->evalCovMatrixSymmetric(data);
    MatrixRectangular XP         = model->evalDriftMatrix(data);
    MatrixRectangular Sigma0P    = model->evalCovMatrix(data, target);
    VectorDouble ZP = data->getMultipleValuesActive(VectorInt(), VectorInt(), means);

    KrigingCalcul KcalcP(ZP, &SigmaP, &XP, &C00, means);
    KcalcP.setSigma0(&Sigma0P);
    KcalcP.setX0(&X0);

    VH::display("Kriging Value(s)", KcalcP.getEstimation());
    VH::display("Standard Deviation of Estimation Error", KcalcP.getStdv());
    VH::display("Variance of Estimator", KcalcP.getVarianceZstar());
  }
  
  // ====================== Free pointers ==================================
  delete neigh;
  delete data;
  delete target;
  delete model;

  return (0);
}
