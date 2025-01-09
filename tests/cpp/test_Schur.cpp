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
  DbStringFormat* dbfmt =
    DbStringFormat::createFromFlags(false, false, false, false, true);
  
  int iech = datap->addSamples(1);
  datap->setSampleCoordinates(iech, target->getSampleCoordinates(0));
  datap->setLocVariables(ELoc::Z, iech, valuesTarget);
  datap->display(dbfmt);
  return datap;
}

static Db* _dataTargetDeplement(Db* data, const VectorInt& varXvalid, int iech0)
{
  Db* datap = data->clone();
  DbStringFormat* dbfmt =
    DbStringFormat::createFromFlags(false, false, false, false, true);

  // Delete the cross-validation information
  for (int i = 0, nval = (int) varXvalid.size(); i < nval; i++)
    datap->setLocVariable(ELoc::Z, iech0, varXvalid[i], TEST);
  datap->display(dbfmt);

  return datap;
}

static Db* _dataAsIs(Db* data)
{
  Db* datap = data->clone();
  DbStringFormat* dbfmt =
    DbStringFormat::createFromFlags(false, false, false, false, true);
  datap->display(dbfmt);
  delete dbfmt;
  
  return datap;
}

/****************************************************************************/
/*!
 ** Testing the Bayesian option
 **
 *****************************************************************************/
static void _firstTest(Db* data,
                       Db* target,
                       ModelGeneric* model,
                       ANeigh* neigh,
                       const VectorDouble& means,
                       const VectorDouble& PriorMean,
                       MatrixSquareSymmetric& PriorCov)
{
  Model* modelc = dynamic_cast<Model*>(model);
  if (!modelc->hasDrift())
  {
    messerr("The 'Model' must have drift defined to check 'Bayesian' option");
    return;
  }
  
  
  // Local parameters
  bool debugPrint = false;
  bool debugSchur = false;
  int iech0 = 0;
  if (debugPrint) OptDbg::setReference(iech0 + 1);

  // Title
  mestitle(0, "Bayes option");
  message("Compare:\n");
  message("- Kriging with traditional code\n");
  message("- Estimation performed with 'KrigingCalcul'\n");
  message("Option: Bayesian\n");

  // Creating the data file
  Db* dataP = _dataAsIs(data);
   
  // ---------------------- Using Standard Kriging procedure ---------------
  mestitle(1, "Using Standard Kriging procedure");
  Table table;
  kribayes(dataP, target, modelc, neigh, PriorMean, PriorCov, true, true);
  table = target->printOneSample(iech0, {"Bayes*"}, true, true);
  target->deleteColumn("Bayes*");
  table.display();

  // ---------------------- Using Schur Class ------------------------------
  mestitle(1, "Using Schur class");

  MatrixSquareSymmetric Sigma00 = modelc->getSillValues(0);
  MatrixSquareSymmetric Sigma   = modelc->evalCovMatrixSymmetric(data);
  MatrixRectangular X           = model->evalDriftMatrix(data);
  MatrixRectangular Sigma0      = modelc->evalCovMatrix(data, target);
  MatrixRectangular X0          = model->evalDriftMatrix(target);
  VectorDouble Z = data->getMultipleValuesActive(VectorInt(), VectorInt(), means);
  KrigingCalcul Kcalc;
  Kcalc.setData(&Z, &means);
  Kcalc.setLHS(&Sigma, &X);
  Kcalc.setVar(&Sigma00);
  Kcalc.setRHS(&Sigma0, &X0);
  Kcalc.setBayes(&PriorMean, &PriorCov);

  VH::display("Prior Mean", PriorMean);
  message("Prior Variance-Covariance Matrix\n");
  PriorCov.display();
  VectorDouble beta = Kcalc.getPostMean();
  if (!beta.empty()) VH::display("Posterior Mean", beta);
  message("Posterior Variance-Covariance Matrix\n");
  Kcalc.getPostCov()->display();

  VH::display("Kriging Value(s)", Kcalc.getEstimation());
  VH::display("Standard Deviation of Estimation Error", Kcalc.getStdv());
  VH::display("Variance of Estimator", Kcalc.getVarianceZstar());
  if (debugSchur) Kcalc.printStatus();

  delete dataP;
}

/****************************************************************************/
/*!
 ** Testing Collocated option
 **
 *****************************************************************************/
static void _secondTest(Db* data, Db* target, ModelGeneric* model, const VectorDouble& means)
{
  Model* modelc = dynamic_cast<Model*>(model);
  // Local parameters
  int nvar = modelc->getVariableNumber();
  VectorInt varColCok = {0, 2}; // Ranks of collcated variables
  bool debugSchur     = false;
  if (nvar <= 1)
  {
    messerr("The collocated Test only makes sense for more than 1 variable");
    return;
  }

  // Title
  mestitle(0,"Collocated Option (in Unique Neighborhood):");
  message("- using 'KrigingCalcul' on the Complemented Data Set\n");
  message(
    "- using 'KrigingCalcul' on Standard Data Set adding Collocated Option\n");
  VH::display("- Collocated Variable ranks", varColCok, false);

  // Creating the Complemented Data Set
  VectorDouble valuesTarget(nvar, TEST);
  for (int ivar = 0; ivar < (int)varColCok.size(); ivar++)
    valuesTarget[varColCok[ivar]] = law_gaussian();
  Db* dataP = _dataComplement(data, target, valuesTarget);

  // ---------------------- With complemented Data Base ---------------------
  mestitle(1, "With Complemented input Data Base");

  MatrixSquareSymmetric Sigma00P = modelc->getSillValues(0);
  MatrixSquareSymmetric SigmaP   = model->evalCovMatrixSymmetric(dataP);
  MatrixRectangular XP           = model->evalDriftMatrix(dataP);
  MatrixRectangular Sigma0P      = model->evalCovMatrix(dataP, target);
  MatrixRectangular X0P          = model->evalDriftMatrix(target);
  VectorDouble ZP = dataP->getMultipleValuesActive(VectorInt(), VectorInt(), means);

  KrigingCalcul KcalcP;
  KcalcP.setData(&ZP, &means);
  KcalcP.setLHS(&SigmaP, &XP);
  KcalcP.setRHS(&Sigma0P, &X0P);
  KcalcP.setVar(&Sigma00P);

  VH::display("Kriging Value(s)", KcalcP.getEstimation());
  VH::display("Standard Deviation of Estimation Error", KcalcP.getStdv());
  VH::display("Variance of Estimator", KcalcP.getVarianceZstar());

  if (debugSchur) KcalcP.printStatus();

  // ---------------------- With Collocated Option -------------------------
  mestitle(1, "With Collocated Option");

  MatrixSquareSymmetric Sigma00 = modelc->getSillValues(0);
  MatrixSquareSymmetric Sigma   = model->evalCovMatrixSymmetricOptim(data);
  MatrixRectangular X           = model->evalDriftMatrix(data);
  MatrixRectangular Sigma0      = model->evalCovMatrixOptim(data, target);
  MatrixRectangular X0          = model->evalDriftMatrix(target);
  VectorDouble Z = data->getMultipleValuesActive(VectorInt(), VectorInt(), means);

  KrigingCalcul Kcalc;
  Kcalc.setData(&Z, &means);
  Kcalc.setLHS(&Sigma, &X);
  Kcalc.setRHS(&Sigma0, &X0);
  Kcalc.setVar(&Sigma00);
  // Subtract the mean (non zero for SK only) from the Collocated values
  VH::subtractInPlace(valuesTarget, means);
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
static void _thirdTest(Db* data, ModelGeneric* model, const VectorDouble& means)
{
  Model* modelc = dynamic_cast<Model*>(model);
  // Set of ranks of cross-validated information
  VectorInt varXvalid = {1,2};
  int iech0           = 1;
  AStringFormat format;
  bool debugSchur = false;

  const VectorVectorInt index = data->getMultipleRanksActive();
  VectorInt rankXvalidEqs = Db::getMultipleSelectedIndices(index, varXvalid, {iech0});
  VectorInt rankXvalidVars = Db::getMultipleSelectedVariables(index, varXvalid, {iech0});

  // Title
  mestitle(0, "Cross-Validation (in Unique Neighborhood)");
  message("Compare the Cross-validation Option (in Unique Neighborhood):\n");
  message("- using Standard Kriging on the Deplemented Data Set\n");
  message("- using 'KrigingCalcul' on Initial Set with Cross-validation option\n");

  // Creating the Complemented Data Set
  Db* targetP = data->clone();
  Db* dataP = _dataTargetDeplement(data, varXvalid, iech0);

  // ----------------------With Deplemented Data Base ---------------------
  mestitle(1, "With Deplemented input Data Base");

  MatrixSquareSymmetric Sigma00P = modelc->getSillValues(0);
  MatrixSquareSymmetric SigmaP   = model->evalCovMatrixSymmetricOptim(dataP);
  MatrixRectangular XP           = model->evalDriftMatrix(dataP);
  MatrixRectangular Sigma0P      = model->evalCovMatrixOptim(dataP, targetP, -1, -1, VectorInt(), {iech0});
  MatrixRectangular X0P          = model->evalDriftMatrix(targetP, -1, {iech0});
  VectorDouble ZP =
    dataP->getMultipleValuesActive(VectorInt(), VectorInt(), means);

  KrigingCalcul KcalcP;
  KcalcP.setData(&ZP, &means);
  KcalcP.setLHS(&SigmaP, &XP);
  KcalcP.setRHS(&Sigma0P, &X0P);
  KcalcP.setVar(&Sigma00P);

  VH::display("Kriging Value(s)", KcalcP.getEstimation());
  VH::display("Standard Deviation of Estimation Error", KcalcP.getStdv());
  VH::display("Variance of Estimator", KcalcP.getVarianceZstar());

  if (debugSchur) KcalcP.printStatus();

  // ---------------------- With Cross-validation Option -------------------------
  mestitle(1, "With Cross-Validation Option");

  MatrixSquareSymmetric Sigma00 = modelc->getSillValues(0);
  MatrixSquareSymmetric Sigma   = model->evalCovMatrixSymmetricOptim(data);
  MatrixRectangular X           = model->evalDriftMatrix(data);
  MatrixRectangular Sigma0      = model->evalCovMatrixOptim(data, targetP);
  MatrixRectangular X0          = model->evalDriftMatrix(targetP);
  VectorDouble Z =
    data->getMultipleValuesActive(VectorInt(), VectorInt(), means);

  KrigingCalcul Kcalc;
  Kcalc.setData(&Z, &means);
  Kcalc.setLHS(&Sigma, &X);
  Kcalc.setVar(&Sigma00);
  Kcalc.setXvalidUnique(&rankXvalidEqs, &rankXvalidVars);

  VH::display("Kriging Value(s)", Kcalc.getEstimation());
  VH::display("Standard Deviation of Estimation Error", Kcalc.getStdv());
  VH::display("Variance of Estimator", Kcalc.getVarianceZstar());

  if (debugSchur) Kcalc.printStatus();

  delete dataP;
  delete targetP;
}

/****************************************************************************/
/*!
 ** Testing Dual option
 **
 ** Note: Means are set to 0 to check SK option
 **
 *****************************************************************************/
static void _fourthTest(Db* data, Db* target, ModelGeneric* model, const VectorDouble& means)
{
  // Title
  mestitle(0, "Estimation using Dual option or not (in Unique Neighborhood):");

  Model* modelc = dynamic_cast<Model*>(model);
  // ---------------------- Without Dual option ---------------------
  mestitle(1, "Without Dual option");

  MatrixSquareSymmetric Sigma00 = modelc->getSillValues(0);
  MatrixSquareSymmetric Sigma   = model->evalCovMatrixSymmetricOptim(data);
  MatrixRectangular X           = model->evalDriftMatrix(data);
  MatrixRectangular Sigma0      = model->evalCovMatrixOptim(data, target);
  MatrixRectangular X0          = model->evalDriftMatrix(target);
  VectorDouble Z = data->getMultipleValuesActive(VectorInt(), VectorInt(), means);

  KrigingCalcul Kcalc1(false);
  Kcalc1.setData(&Z, &means);
  Kcalc1.setLHS(&Sigma, &X);
  Kcalc1.setRHS(&Sigma0, &X0);
  Kcalc1.setVar(&Sigma00);

  VH::display("Kriging Value(s)", Kcalc1.getEstimation());
  VH::display("Standard Deviation of Estimation Error", Kcalc1.getStdv());
  VH::display("Variance of Estimator", Kcalc1.getVarianceZstar());

  // ---------------------- With Dual Option -------------------------
  mestitle(1, "With Dual Option (only Estimation is available)");

  KrigingCalcul Kcalc2(true);
  Kcalc2.setData(&Z, &means);
  Kcalc2.setLHS(&Sigma, &X);
  Kcalc2.setRHS(&Sigma0, &X0);

  VH::display("Kriging Value(s)", Kcalc2.getEstimation());
}

/****************************************************************************/
/*!
 ** Main Program
 **
 ** This test is composed of several parts, comparing the results of
 ** traditional Kriging with the results of the Algebraic calculations
 ** provided within 'KrigingCalcul'.
 ** Different scenarios are elaborated:
 ** 1) Bayesian case
 ** 2) Test on Collocated CoKriging in Unique Neighborhood
 ** 3) Test on Cross-Validation in Unique Neighborhood
 ** 4) Test on Estimation using the Dual option or not
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
  int mode = 0;

  // Generate the data base
  Db* data = Db::createFillRandom(nech, ndim, nvar, nfex);

  // Generate the target file
  Db* target = Db::createFillRandom(1, ndim, 0, nfex);

  // Create the constant Mean vector
  VectorDouble means(nvar, 0.);
  if (flagSK) means = VH::simulateGaussian(nvar);

  // Create the Model
  ModelGeneric* model;
  
  double scale = 0.7;
  MatrixSquareSymmetric* sills =
    MatrixSquareSymmetric::createRandomDefinitePositive(nvar);
  model = Model::createFromParam(ECov::EXPONENTIAL, scale, 0., 0., VectorDouble(),
                                 *sills, VectorDouble(), nullptr, false);
  Model* modelc = dynamic_cast<Model*>(model);
  modelc->setMeans(means);
  if (!flagSK) modelc->setDriftIRF(0, nfex);

  // Create the Bayesian Priors for Drift coefficients
  VectorDouble PriorMean = VH::simulateGaussian(nbfl);
  MatrixSquareSymmetric PriorCov(nbfl);
  PriorCov.setDiagonal(VH::simulateUniform(nbfl, 0.1, 0.5));

  // Unique Neighborhood
  NeighUnique* neigh = NeighUnique::create();

  // Define the verbose option
  int iech0 = 0;
  if (debugPrint) OptDbg::setReference(iech0 + 1);

  // Test on Bayesian
  if (mode == 0 || mode == 1)
  {
    Db* dataLocal = data->clone();
    _firstTest(dataLocal, target, model, neigh, means, PriorMean, PriorCov);
    delete dataLocal;
  }

  // Test on Collocated CoKriging in Unique Neighborhood
  if (mode == 0 || mode == 2)
  {
    Db* dataLocal = data->clone();
    _secondTest(dataLocal, target, model, means);
    delete dataLocal;
  }

  // Test on Cross-Validation in Unique Neighborhood
  if (mode == 0 || mode == 3)
  {
    Db* dataLocal = data->clone();
    _thirdTest(dataLocal, model, means);
    delete dataLocal;
  }

  // Test on Estimation using the Dual option or not
  if (mode == 0 || mode == 4)
  {
    Db* dataLocal = data->clone();
    _fourthTest(dataLocal, target, model, means);
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
