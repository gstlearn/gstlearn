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

  // Parameters
  bool debugPrint = false;
  int mode = -1; // 0: Standard; 1: Schur; -1 Both

  bool flagSK         = false;
  bool flagBayes      = true;
  if (flagBayes) flagSK = false;

  int nech    = 10;
  int nvar    = 2;
  int nfex    = 0;
  int nbfl    = (nfex + 1) * nvar;

  VectorDouble means(nvar, 0.);
  if (flagSK) means = VH::simulateGaussian(nvar);
  VectorDouble PriorMean = VH::simulateGaussian(nbfl);
  MatrixSquareSymmetric PriorCov(nbfl);
  PriorCov.setDiagonal(VH::simulateUniform(nbfl, 0.1, 0.5));

  // Generate the data base
  Db* data = Db::createFillRandom(nech, ndim, nvar, nfex);
  DbStringFormat* dbfmt =
    DbStringFormat::createFromFlags(true, true, false, false, true);
  if (debugPrint) data->display(dbfmt);

  // Generate the target file
  Db* target = Db::createFillRandom(1, ndim, 0, nfex);
  if (debugPrint) target->display();

  // Create the Model
  Model* model;
  double range = 0.2;
  if (nvar == 1)
  {
    double sill = 2.;
    model = Model::createFromParam(ECov::SPHERICAL, range, sill);
  }
  else
  {
    MatrixSquareSymmetric sills(nvar);
    sills.setDiagonal(VH::simulateUniform(nvar, 1., 2.));
    model = Model::createFromParam(ECov::SPHERICAL, range, 0., 0., VectorDouble(), sills.getValues());
  }
  model->setMeans(means);
  if (! flagSK) model->setDriftIRF(0, nfex);
  if (debugPrint) model->display();

  // Unique Neighborhood
  NeighUnique* neigh = NeighUnique::create();
  if (debugPrint) neigh->display();

  // Define the verbose option
  int iech0 = 0;
  if (debugPrint) OptDbg::setReference(iech0 + 1);

  // ====================== Using Standard Kriging procedure ===============
  if (mode < 0 || mode == 0)
  {
    mestitle(0, "Using Standard Kriging procedure");
    Table table;
    if (!flagBayes)
    {
      kriging(data, target, model, neigh, EKrigOpt::POINT, true, true, true);
      table = target->printOneSample(iech0, {"Kriging*"});
    }
    else
    {
      OptDbg::define(EDbg::BAYES);
      kribayes(data, target, model, neigh, PriorMean, PriorCov, true, true);
      table = target->printOneSample(iech0, {"Bayes*"});
    }
    table.display();
  }

  // ====================== Using Schur Class ==============================
  if (mode < 0 || mode == 1)
  {
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

    VectorDouble Zstar = Kcalc.getEstimation();
    if (!Zstar.empty())
      VH::display("Kriging Value(s)", Zstar);

    VectorDouble stdv = Kcalc.getStdv();
    if (!stdv.empty())
      VH::display("Standard Deviation of Estimation Error", stdv);

    VectorDouble varianceZstar = Kcalc.getVarianceZstar();
    if (!varianceZstar.empty())
      VH::display("Variance of Estimator", varianceZstar);

  //  if (debugPrint) Kcalc.printStatus();
  }

  // ====================== Free pointers ==================================
  delete neigh;
  delete data;
  delete target;
  delete model;

  return (0);
}
