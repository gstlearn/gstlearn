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
#include "Basic/AStringFormat.hpp"
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

static Db* _dataComplement(Db* data, Db* target, const VectorDouble& valuesTarget)
{
  // Complement the Data Set

  Db* datap = data->clone();
  int iech = datap->addSamples(1);
  datap->setSampleCoordinates(iech, target->getSampleCoordinates(0));
  datap->setLocVariables(ELoc::Z, iech, valuesTarget);
  return datap;
}

static Db* _dataTargetDeplement(Db* data, const VectorInt& varXvalid, int iech0)
{
  Db* datap = data->clone();
  DbStringFormat* dbfmt =
    DbStringFormat::createFromFlags(false, true, false, false, true);

  // Delete the cross-validation information
  for (int i = 0, nval = (int) varXvalid.size(); i < nval; i++)
    datap->setLocVariable(ELoc::Z, iech0, varXvalid[i], TEST);
  datap->display(dbfmt);

  return datap;
}

/****************************************************************************/
/*!
 ** Testing Bayes option
 **
 *****************************************************************************/
static void _firstTest(Db* data,
                       Db* target,
                       Model* model,
                       ANeigh* neigh,
                       const VectorDouble& means,
                       const VectorDouble& PriorMean,
                       MatrixSquareSymmetric& PriorCov)
{
  // Local parameters
  bool debugPrint = false;
  bool debugSchur = false;
  bool flagBayes    = false;
  int iech0 = 0;
  if (debugPrint) OptDbg::setReference(iech0 + 1);

  // Title
  mestitle(0, "Checking Bayes option");
  message("Compare:\n");
  message("- Kriging with traditional code\n");
  message("- Estimation performed with 'KrigingCalcul'\n");
  message("- Option: Bayesian\n");

  // ---------------------- Using Standard Kriging procedure ---------------
  mestitle(1, "Using Standard Kriging procedure");
  Table table;
  if (!flagBayes)
  {
    kriging(data, target, model, neigh, EKrigOpt::POINT, true, true, true);
    table = target->printOneSample(iech0, {"Kriging*"}, true, true);
    target->deleteColumn("Kriging*");
  }
  else
  {
    kribayes(data, target, model, neigh, PriorMean, PriorCov, true, true);
    table = target->printOneSample(iech0, {"Bayes*"}, true, true);
    target->deleteColumn("Bayes*");
  }
  table.display();

  // ---------------------- Using Schur Class ------------------------------
  mestitle(1, "Using Schur class");

  MatrixSquareSymmetric Sigma00 = model->getSillValues(0);
  MatrixSquareSymmetric Sigma   = model->evalCovMatrixSymmetric(data);
  MatrixRectangular X           = model->evalDriftMatrix(data);
  MatrixRectangular Sigma0      = model->evalCovMatrix(data, target);
  MatrixRectangular X0          = model->evalDriftMatrix(target);
  VectorDouble Z = data->getMultipleValuesActive(VectorInt(), VectorInt(), means);
  KrigingCalcul Kcalc(&Z, &Sigma, &X, &Sigma00, &means);
  Kcalc.setTarget(&Sigma0, &X0);
  if (flagBayes) Kcalc.setBayes(&PriorMean, &PriorCov);

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
  if (debugSchur) Kcalc.printStatus();
}

/****************************************************************************/
/*!
 ** Testing Collocated option
 **
 *****************************************************************************/
static void _secondTest(Db* data, Db* target, Model* model, ANeigh* neigh, const VectorDouble& means)
{
  // Local parameters
  int nvar = model->getVariableNumber();
  if (nvar <= 1) return;
  // Set of ranks of collocated variables
  VectorInt varColCok = {1};
  int iech0           = 0;
  AStringFormat format;
  bool debugSchur = false;

  // Title
  mestitle(0,"Compare the Collocated Option (in Unique Neighborhood):");
  message("- using Standard Kriging on the Complemented Data Set\n");
  message("- using 'KrigingCalcul' on the Complemented Data Set\n");
  message("- using 'KrigingCalcul' on Standard Data Set adding Collocated Option\n");

  // Creating the Complemented Data Set
  VectorDouble valuesTarget(nvar, TEST);
  for (int ivar = 0; ivar < (int)varColCok.size(); ivar++)
    valuesTarget[varColCok[ivar]] = law_gaussian();
  Db* dataP = _dataComplement(data, target, valuesTarget);

  // ---------------------- Using Standard Kriging procedure ---------------
  mestitle(1, "Using Standard Kriging procedure");

  Table table;
  kriging(dataP, target, model, neigh, EKrigOpt::POINT, true, true, true);
  table = target->printOneSample(iech0, {"Kriging*"}, true, true);
  table.display();
  target->deleteColumn("Kriging*");

  // ---------------------- With complemented Data Base ---------------------
  mestitle(1, "With complemented input Data Base");

  MatrixSquareSymmetric Sigma00P = model->getSillValues(0);
  MatrixSquareSymmetric SigmaP   = model->evalCovMatrixSymmetric(dataP);
  MatrixRectangular XP           = model->evalDriftMatrix(dataP);
  MatrixRectangular Sigma0P      = model->evalCovMatrix(dataP, target);
  MatrixRectangular X0P          = model->evalDriftMatrix(target);
  VectorDouble ZP = dataP->getMultipleValuesActive(VectorInt(), VectorInt(), means);

  KrigingCalcul KcalcP(&ZP, &SigmaP, &XP, &Sigma00P, &means);
  KcalcP.setTarget(&Sigma0P, &X0P);

  VH::display("Kriging Value(s)", KcalcP.getEstimation());
  VH::display("Standard Deviation of Estimation Error", KcalcP.getStdv());
  VH::display("Variance of Estimator", KcalcP.getVarianceZstar());
  if (debugSchur) KcalcP.printStatus();

  // ---------------------- With Collocated Option -------------------------
  mestitle(1, "With Collocated Option");

  MatrixSquareSymmetric Sigma00 = model->getSillValues(0);
  MatrixSquareSymmetric Sigma   = model->evalCovMatrixSymmetric(data);
  MatrixRectangular X           = model->evalDriftMatrix(data);
  MatrixRectangular Sigma0      = model->evalCovMatrix(data, target);
  MatrixRectangular X0          = model->evalDriftMatrix(target);
  VectorDouble Z =
    data->getMultipleValuesActive(VectorInt(), VectorInt(), means);
  KrigingCalcul Kcalc(&Z, &Sigma, &X, &Sigma00, &means);
  Kcalc.setTarget(&Sigma0, &X0);
  Kcalc.setColCokUnique(&valuesTarget, &varColCok);

  VH::display("Kriging Value(s)", Kcalc.getEstimation());
  VH::display("Standard Deviation of Estimation Error", Kcalc.getStdv());
  VH::display("Variance of Estimator", Kcalc.getVarianceZstar());
  if (debugSchur) Kcalc.printStatus();

  delete dataP;
}

/****************************************************************************/
/*!
 ** Testing Cross-validation option
 **
 *****************************************************************************/
static void _thirdTest(Db* data, Model* model, const VectorDouble& means)
{
  // Set of ranks of cross-validated information
  VectorInt varXvalid = {0,1};
  int iech0           = 1;
  AStringFormat format;
  bool debugSchur = false;

  // Title
  mestitle(0, "Cross-Validation in Unique Neighborhood");
  message("Compare the Cross-validation Option (in Unique Neighborhood):\n");
  message("- using Standard Kriging on the Deplemented Data Set\n");
  message("- using 'KrigingCalcul' on Initial Set with Cross-validation option\n");

  const VectorVectorInt index = data->getMultipleRanksActive();
  VectorInt rankXvalid        = Db::getMultipleRanks(index, varXvalid, {iech0});

  // Creating the Complemented Data Set
  Db* targetP = data->clone();
  Db* dataP = _dataTargetDeplement(data, varXvalid, iech0);

  // ---------------------- With deplemented Data Base ---------------------
  mestitle(1, "With deplemented input Data Base");

  MatrixSquareSymmetric Sigma00P = model->getSillValues(0);
  MatrixSquareSymmetric SigmaP   = model->evalCovMatrixSymmetric(dataP);
  MatrixRectangular XP           = model->evalDriftMatrix(dataP);
  MatrixRectangular Sigma0P      = model->evalCovMatrix(dataP, targetP, -1, -1, VectorInt(), {iech0});
  MatrixRectangular X0P          = model->evalDriftMatrix(targetP, -1, {iech0});
  VectorDouble ZP =
    dataP->getMultipleValuesActive(VectorInt(), VectorInt(), means);

  KrigingCalcul KcalcP(&ZP, &SigmaP, &XP, &Sigma00P, &means);
  KcalcP.setTarget(&Sigma0P, &X0P);

  VH::display("Kriging Value(s)", KcalcP.getEstimation());
  VH::display("Standard Deviation of Estimation Error", KcalcP.getStdv());
  VH::display("Variance of Estimator", KcalcP.getVarianceZstar());

  if (debugSchur) KcalcP.printStatus();

  // ---------------------- With Cross-validation Option -------------------------
  mestitle(1, "With Cross-Validation Option");

  MatrixSquareSymmetric Sigma00 = model->getSillValues(0);
  MatrixSquareSymmetric Sigma   = model->evalCovMatrixSymmetric(data);
  MatrixRectangular X           = model->evalDriftMatrix(data);
  MatrixRectangular Sigma0      = model->evalCovMatrix(data, targetP);
  MatrixRectangular X0          = model->evalDriftMatrix(targetP);
  VectorDouble Z =
    data->getMultipleValuesActive(VectorInt(), VectorInt(), means);
  
  KrigingCalcul Kcalc(&Z, &Sigma, &X, &Sigma00, &means);
  Kcalc.setXvalidUnique(&rankXvalid);

  VH::display("Kriging Value(s)", Kcalc.getEstimation());
  VH::display("Standard Deviation of Estimation Error", Kcalc.getStdv());
  VH::display("Variance of Estimator", Kcalc.getVarianceZstar());

  if (debugSchur) Kcalc.printStatus();

  delete dataP;
  delete targetP;
}

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
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Global parameters
  int ndim = 2;
  law_set_random_seed(32131);
  AStringFormat format;

  defineDefaultSpace(ESpaceType::RN, ndim);
  OptCst::define(ECst::NTCOL, -1);
  OptCst::define(ECst::NTROW, -1);
  OptDbg::define(EDbg::BAYES);

  // Parameters
  bool debugPrint   = false;
  int nech          = 3;
  int nvar          = 3;
  int nfex          = 0;
  int nbfl          = (nfex + 1) * nvar;
  bool flagSK       = false;
  bool flagBayes    = false;
  if (flagBayes) flagSK = false;
  int mode = 3;

  // Generate the data base
  Db* data = Db::createFillRandom(nech, ndim, nvar, nfex);

  // Generate the target file
  Db* target = Db::createFillRandom(1, ndim, 0, nfex);

  // Create the constant Mean vector
  VectorDouble means(nvar, 0.);
  if (flagSK) means = VH::simulateGaussian(nvar);

  // Create the Model
  Model* model;
  double range = 0.7;
  MatrixSquareSymmetric* sills =
    MatrixSquareSymmetric::createRandomDefinitePositive(nvar);
  model = Model::createFromParam(ECov::SPHERICAL, range, 0., 0., VectorDouble(),
                                 sills->getValues());
  model->setMeans(means);
  if (!flagSK) model->setDriftIRF(0, nfex);

  // Create the Bayesian Priors for Drift coefficients
  VectorDouble PriorMean = VH::simulateGaussian(nbfl);
  MatrixSquareSymmetric PriorCov(nbfl);
  PriorCov.setDiagonal(VH::simulateUniform(nbfl, 0.1, 0.5));

  // Unique Neighborhood
  NeighUnique* neigh = NeighUnique::create();

  // Define the verbose option
  int iech0 = 0;
  if (debugPrint) OptDbg::setReference(iech0 + 1);

  // First Test
  if (mode == 0 || mode == 1)
  {
    Db* dataLocal = data->clone();
    _firstTest(dataLocal, target, model, neigh, means, PriorMean, PriorCov);
    delete dataLocal;
  }

  // Second Test
  if (mode == 0 || mode == 2)
  {
    Db* dataLocal = data->clone();
    _secondTest(dataLocal, target, model, neigh, means);
    delete dataLocal;
  }

  // Third Test
  if (mode == 0 || mode == 3)
  {
    Db* dataLocal = data->clone();
    _thirdTest(dataLocal, model, means);
    delete dataLocal;
  }

  // Free pointers

  delete sills;
  delete neigh;
  delete data;
  delete target;
  delete model;

  return (0);
}
