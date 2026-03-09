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
#include "Estimation/Estimations.hpp"

#include "Basic/OptCustom.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Enum/EMorpho.hpp"
#include "Estimation/CalcGlobal.hpp"
#include "Estimation/CalcImage.hpp"
#include "Estimation/CalcKrigingFactors.hpp"
#include "Estimation/CalcKrigingGradient.hpp"
#include "Estimation/CalcKrigingSimpleCase.hpp"
#include "Estimation/CalcSimpleInterpolation.hpp"
#include "Model/Model.hpp"
#include "Model/ModelGeneric.hpp"
#include "Neigh/ANeigh.hpp"
#include "Neigh/NeighBench.hpp"
#include "Neigh/NeighMoving.hpp"
#include "Neigh/NeighUnique.hpp"

namespace gstlrn
{
ANeigh* _defaultNeighborhood(ANeigh* neigh, Db* dbin, Id maxNumber = 500)
{
  // If a neighborhood is already defined, this is the correct solution
  if (neigh != nullptr) return neigh;

  // We are about to define a Unique Neighborhood by default
  // Let us first check that the number of active samples is not too large
  if (dbin != nullptr)
  {
    Id nech = dbin->getNSample(true);
    if (nech > maxNumber)
    {
      messerr("No neighborhood has been defined");
      messerr("The number of active samples (%d) is too large (>%d)", nech, maxNumber);
      messerr("to allow the definition of a Unique Neighborhood by default");
      return nullptr;
    }
  }

  // Create a default unique neighborhood
  auto* neighUnique = new NeighUnique();
  return neighUnique;
}

Global_Result global_arithmetic(Db* dbin,
                                DbGrid* dbgrid,
                                ModelGeneric* model,
                                Id ivar0,
                                bool verbose)
{
  Global_Result gres;
  CalcGlobal global(ivar0, verbose);
  global.setDbin(dbin);
  global.setDbout(dbgrid);
  global.setModelGeneric(model);
  global.setFlagArithmetic(true);

  if (global.run())
    gres = global.getGRes();
  return gres;
}

Global_Result global_kriging(Db* dbin,
                             Db* dbout,
                             ModelGeneric* model,
                             Id ivar0,
                             bool verbose)
{
  Global_Result gres;
  CalcGlobal global(ivar0, verbose);
  global.setDbin(dbin);
  global.setDbout(dbout);
  global.setModelGeneric(model);
  global.setFlagKriging(true);

  if (global.run())
    gres = global.getGRes();
  return gres;
}

/****************************************************************************/
/*!
 **  Kriging (Factorial) a regular grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid     input and output Db grid structure
 ** \param[in]  model      Model structure
 ** \param[in]  neigh      ANeigh structure
 ** \param[in]  flagFFT    True if the FFT version is to be used
 ** \param[in]  verbose    Verbose flag
 ** \param[in]  seed       Seed used for random number generation
 ** \param[in]  namconv    Naming Convention
 **
 *****************************************************************************/
Id krimage(DbGrid* dbgrid,
           Model* model,
           ANeigh* neigh,
           bool flagFFT,
           bool verbose,
           Id seed,
           const NamingConvention& namconv)
{
  CalcImage image;

  image.setDbin(dbgrid);
  image.setDbout(dbgrid);
  image.setModelGeneric(model);
  image.setNeigh(neigh);
  image.setFlagFFT(flagFFT);
  image.setSeed(seed);
  image.setVerbose(verbose);
  image.setNamingConvention(namconv);

  image.setFlagFilter(true);

  Id error = (image.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Smooth a regular grid
 **
 ** \return  Error return code
 **
 ** \param[in]  dbgrid     input and output Db grid structure
 ** \param[in]  neigh      ANeigh structure
 ** \param[in]  type       1 for Uniform; 2 for Gaussian
 ** \param[in]  range      Range (used for Gaussian only)
 ** \param[in]  namconv    Naming Convention
 **
 *****************************************************************************/
Id dbSmoother(DbGrid* dbgrid,
              ANeigh* neigh,
              Id type,
              double range,
              const NamingConvention& namconv)
{
  CalcImage image;

  image.setDbin(dbgrid);
  image.setDbout(dbgrid);
  image.setNeigh(neigh);
  image.setNamingConvention(namconv);

  image.setFlagSmooth(true);
  image.setSmoothType(type);
  image.setSmoothRange(range);

  Id error = (image.run()) ? 0 : 1;
  return error;
}

/**
 * Perform a Morphological operation on an image stored in Db
 * @param dbgrid  Target IN/OUT Db (must be a Grid)
 * @param oper    Type of morphological operation
 * @param vmin    Minimum threshold value
 * @param vmax    Maximum threshold value
 * @param option  Option
 * @param radius  Radius
 * @param verbose Verbose option
 * @param flagDistErode True: Inflate the grain; False: Reduce the grain
 * @param namconv Naming convention
 * @return
 */
GSTLEARN_EXPORT Id dbMorpho(DbGrid* dbgrid,
                            const EMorpho& oper,
                            double vmin,
                            double vmax,
                            Id option,
                            const VectorInt& radius,
                            bool flagDistErode,
                            bool verbose,
                            const NamingConvention& namconv)
{
  CalcImage image;

  image.setDbin(dbgrid);
  image.setDbout(dbgrid);
  image.setNamingConvention(namconv);

  image.setFlagMorpho(true);
  image.setOper(oper);
  image.setVmin(vmin);
  image.setVmax(vmax);
  image.setOption(option);
  image.setRadius(radius);
  image.setDistErode(flagDistErode);
  image.setVerbose(verbose);

  // Particular case of the number of output variables
  Id nvar = 1;
  if (oper == EMorpho::GRADIENT) nvar = dbgrid->getNDim();
  image.setNvarMorpho(nvar);

  Id error = (image.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Standard Kriging
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  model       ModelGeneric structure
 ** \param[in]  neigh       ANeigh structure
 ** \param[in]  flag_est    Option for storing the estimation
 ** \param[in]  flag_std    Option for storing the standard deviation
 ** \param[in]  flag_varz   Option for storing the variance of the estimator
 **                         (only available for stationary model)
 ** \param[in]  krigopt     KrigOpt structure
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
Id kriging(Db* dbin,
           Db* dbout,
           ModelGeneric* model,
           ANeigh* neigh,
           bool flag_est,
           bool flag_std,
           bool flag_varz,
           const KrigOpt& krigopt,
           const NamingConvention& namconv)
{
  auto* neighLocal = _defaultNeighborhood(neigh, dbin);
  auto* neighBench = dynamic_cast<NeighBench*>(neighLocal);
  if (krigopt.getCalcul() == EKrigOpt::POINT &&
      !krigopt.hasColcok() &&
      !krigopt.hasMatLC() &&
      neighBench == nullptr &&
      model->getNVar() == 1 &&
      OptCustom::query("NotOptimSimpleCase", 0) == 0 &&
      !dbout->hasLocator(ELoc::DOM))
  {
    OptCustom::define("Optim", 1);
    CalcKrigingSimpleCase krige(flag_est, flag_std, flag_varz);
    krige.setDbin(dbin);
    krige.setDbout(dbout);
    krige.setModelGeneric(model);
    krige.setNeigh(neighLocal);
    krige.setNamingConvention(namconv);
    Id result = krige.run();
    OptCustom::undefine("Optim");
    return 1 - result;
  }

  CalcKriging krige(flag_est, flag_std, flag_varz);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModelGeneric(model);
  krige.setNeigh(neighLocal);
  krige.setKrigopt(krigopt);
  krige.setNamingConvention(namconv);

  Id error = (krige.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Standard Block Kriging with variable cell dimension
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  model       ModelGeneric structure
 ** \param[in]  neigh       ANeigh structure
 ** \param[in]  flag_est    Option for the storing the estimation
 ** \param[in]  flag_std    Option for the storing the standard deviation
 ** \param[in]  krigopt     KrigOpt structure
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
Id krigcell(Db* dbin,
            Db* dbout,
            ModelGeneric* model,
            ANeigh* neigh,
            bool flag_est,
            bool flag_std,
            const KrigOpt& krigopt,
            const NamingConvention& namconv)
{
  auto* neighLocal = _defaultNeighborhood(neigh, dbin);
  CalcKriging krige(flag_est, flag_std, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModelGeneric(model);
  krige.setNeigh(neighLocal);
  krige.setKrigopt(krigopt);
  krige.setNamingConvention(namconv);

  Id error = (krige.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Estimation with Bayesian Drift
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure
 ** \param[in]  dbout      output Db structure
 ** \param[in]  model      ModelGeneric structure
 ** \param[in]  neigh      ANeigh structure
 ** \param[in]  flag_est   Pointer for the storing the estimation
 ** \param[in]  flag_std   Pointer for the storing the standard deviation
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
Id kribayes(Db* dbin,
            Db* dbout,
            ModelGeneric* model,
            ANeigh* neigh,
            bool flag_est,
            bool flag_std,
            const NamingConvention& namconv)
{
  auto* neighLocal = _defaultNeighborhood(neigh, dbin);
  CalcKriging krige(flag_est, flag_std, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModelGeneric(model);
  krige.setNeigh(neighLocal);
  krige.setNamingConvention(namconv);

  krige.setFlagBayes(true);

  Id error = (krige.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Perform kriging and return the calculation elements
 **
 ** \return  A Krigtest_Res structure
 **
 ** \param[in]  dbin        input Db structure
 ** \param[in]  dbout       output Db structure
 ** \param[in]  model       ModelGeneric structure
 ** \param[in]  neigh       ANeigh structure
 ** \param[in]  iech0       Rank of the target sample
 ** \param[in]  krigopt     KrigOpt structure
 ** \param[in]  verbose     When TRUE, the full debugging flag is switched ON
 **                         (the current status is reset after the run)
 **
 *****************************************************************************/
Krigtest_Res krigtest(Db* dbin,
                      Db* dbout,
                      ModelGeneric* model,
                      ANeigh* neigh,
                      Id iech0,
                      const KrigOpt& krigopt,
                      bool verbose)
{
  auto* neighLocal = _defaultNeighborhood(neigh, dbin);
  CalcKriging krige(true, true, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModelGeneric(model);
  krige.setNeigh(neighLocal);
  krige.setKrigopt(krigopt);
  krige.setIechSingleTarget(iech0);
  krige.setVerboseSingleTarget(verbose);

  (void)krige.run();

  return krige.getKtest();
}

/****************************************************************************/
/*!
 **  Punctual Kriging in the Anamorphosed Gaussian Model
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure
 ** \param[in]  dbout      output Db structure
 ** \param[in]  model      ModelGeneric structure
 ** \param[in]  neigh      ANeigh structure
 ** \param[in]  anam       AAnam structure
 ** \param[in]  namconv    Naming convention
 **
 *****************************************************************************/
Id kriggam(Db* dbin,
           Db* dbout,
           ModelGeneric* model,
           ANeigh* neigh,
           AAnam* anam,
           const NamingConvention& namconv)
{
  CalcKriging krige(true, true, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModelGeneric(model);
  krige.setNeigh(neigh);
  krige.setNamingConvention(namconv);

  krige.setFlagGam(true);
  krige.setAnam(anam);

  Id error = (krige.run()) ? 0 : 1;
  return error;
}

/**
 * Standard Cross-Validation
 *
 * @param db Db structure
 * @param model ModelGeneric structure
 * @param neigh ANeigh structure
 * @param flag_kfold True if a code (K-FOLD) is used
 * @param flag_xvalid_est Option for storing the estimation: 1 for Z*-Z; -1 for
 * Z*; 0 not stored
 * @param flag_xvalid_std Option for storing the standard deviation: 1:for
 * (Z*-Z)/S; -1 for S; 0 not stored
 * @param flag_xvalid_varz Option for storing the variance of the estimator: 1
 * to store and 0 not stored
 * @param krigopt KrigOpt structure
 * @param namconv Naming Convention
 * @return Error return code
 */
Id xvalid(Db* db,
          ModelGeneric* model,
          ANeigh* neigh,
          bool flag_kfold,
          Id flag_xvalid_est,
          Id flag_xvalid_std,
          Id flag_xvalid_varz,
          const KrigOpt& krigopt,
          const NamingConvention& namconv)
{
  auto* neighLocal = _defaultNeighborhood(neigh, db);
  CalcKriging krige(flag_xvalid_est != 0, flag_xvalid_std != 0,
                    flag_xvalid_varz != 0);
  krige.setDbin(db);
  krige.setDbout(db);
  krige.setModelGeneric(model);
  krige.setNeigh(neighLocal);
  krige.setNamingConvention(namconv);

  krige.setFlagXvalid(true);
  krige.setFlagXvalidEst(flag_xvalid_est);
  krige.setFlagXvalidStd(flag_xvalid_std);
  krige.setFlagXvalidVarZ(flag_xvalid_varz);
  krige.setFlagKfold(flag_kfold);
  krige.setKrigopt(krigopt);

  Id error = (krige.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Check the Neighborhood
 **
 ** \return  Error return code (0: success, 1: error)
 **
 ** \param[in]  dbin       input Db structure
 ** \param[in]  dbout      output Db structure
 ** \param[in]  model      ModelGeneric structure (optional)
 ** \param[in]  neigh      ANeigh structure
 ** \param[in]  namconv    Naming Convention
 **
 ** \remark This procedure creates the following arrays:
 ** \remark 1 - The number of selected samples
 ** \remark 2 - The maximum neighborhood distance
 ** \remark 3 - The minimum neighborhood distance
 ** \remark 4 - The number of non-empty sectors
 ** \remark 5 - The number of consecutive empty sectors
 **
 *****************************************************************************/
Id test_neigh(Db* dbin,
              Db* dbout,
              ModelGeneric* model,
              ANeigh* neigh,
              const NamingConvention& namconv)
{
  auto* neighLocal = _defaultNeighborhood(neigh, dbin);
  CalcKriging krige(false, false, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModelGeneric(model);
  krige.setNeigh(neighLocal);
  krige.setNamingConvention(namconv);

  krige.setFlagNeighOnly(true);

  Id error = (krige.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Disjunctive Kriging
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure (containing the factors)
 ** \param[in]  dbout      output Grid Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neigh      ANeigh structure

 ** \param[in]  flag_est   Option for the storing the estimation
 ** \param[in]  flag_std   Option for the storing the standard deviation
 ** \param[in]  krigopt    Krigopt structure
 ** \param[in]  namconv    Naming convention
 **
 ** \remark When the change of support is defined through the Anamorphosis
 ** \remark the 'calcul' option must be set to POINT and 'ndisc' does not
 ** \remark have to be defined
 **
 *****************************************************************************/
Id krigingFactors(Db* dbin,
                  Db* dbout,
                  Model* model,
                  ANeigh* neigh,
                  bool flag_est,
                  bool flag_std,
                  const KrigOpt& krigopt,
                  const NamingConvention& namconv)
{
  CalcKrigingFactors krige(flag_est, flag_std);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModelGeneric(model);
  krige.setNeigh(neigh);
  krige.setKrigopt(krigopt);
  krige.setNamingConvention(namconv);

  krige.setIuidFactors(dbin->getUIDsByLocator(ELoc::Z));

  Id error = (krige.run()) ? 0 : 1;
  return error;
}

Id krigingGradient(Db* dbin,
                   Db* dbout,
                   ModelGeneric* model,
                   ANeigh* neigh,
                   bool flag_est,
                   bool flag_std,
                   double ball_radius,
                   bool flagForceNumeric,
                   const NamingConvention& namconv)
{
  CalcKrigingGradient krigeGradient(flag_est, flag_std, ball_radius);
  krigeGradient.setDbin(dbin);
  krigeGradient.setDbout(dbout);
  krigeGradient.setModelGeneric(model);
  krigeGradient.setNeigh(neigh);
  krigeGradient.setFlagForceNumeric(flagForceNumeric);
  krigeGradient.setNamingConvention(namconv);

  Id error = (krigeGradient.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Inverse distance estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  exponent    exponent of the inverse distance
 ** \param[in]  flag_expand True for expansion option (if dbin is Grid)
 ** \param[in]  dmax        Maximum search radius (if dbin is Points)
 ** \param[in]  flag_est    True if the estimation must be calculated
 ** \param[in]  flag_std    True if the St. Dev. must be calculated
 ** \param[in]  model       Model structure (used for St. Dev.)
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
Id inverseDistance(Db* dbin,
                   Db* dbout,
                   double exponent,
                   bool flag_expand,
                   double dmax,
                   bool flag_est,
                   bool flag_std,
                   Model* model,
                   const NamingConvention& namconv)
{
  CalcSimpleInterpolation interpol;
  interpol.setDbin(dbin);
  interpol.setDbout(dbout);
  interpol.setFlagEst(flag_est);
  interpol.setFlagStd(flag_std);
  interpol.setModelGeneric(model);
  interpol.setNamingConvention(namconv);

  interpol.setFlagInvDist(true);
  interpol.setExponent(exponent);
  interpol.setFlagExpand(flag_expand);
  interpol.setDmax(dmax);

  Id error = (interpol.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Moving Average estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  neigh       ANeigh structure
 ** \param[in]  flag_est    True if the estimation must be calculated
 ** \param[in]  flag_std    True if the St. Dev. must be calculated
 ** \param[in]  model       Model structure (used for St. Dev.)
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
GSTLEARN_EXPORT Id movingAverage(Db* dbin,
                                 Db* dbout,
                                 ANeigh* neigh,
                                 bool flag_est,
                                 bool flag_std,
                                 Model* model,
                                 const NamingConvention& namconv)
{
  CalcSimpleInterpolation interpol;
  interpol.setDbin(dbin);
  interpol.setDbout(dbout);
  interpol.setNeigh(neigh);
  interpol.setFlagEst(flag_est);
  interpol.setFlagStd(flag_std);
  interpol.setModelGeneric(model);
  interpol.setNamingConvention(namconv);

  interpol.setFlagMovAve(true);

  Id error = (interpol.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Moving Median estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  neigh       ANeigh structure
 ** \param[in]  flag_est    True if the estimation must be calculated
 ** \param[in]  flag_std    True if the St. Dev. must be calculated
 ** \param[in]  model       Model structure (used for St. Dev.)
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
GSTLEARN_EXPORT Id movingMedian(Db* dbin,
                                Db* dbout,
                                ANeigh* neigh,
                                bool flag_est,
                                bool flag_std,
                                Model* model,
                                const NamingConvention& namconv)
{
  CalcSimpleInterpolation interpol;
  interpol.setDbin(dbin);
  interpol.setDbout(dbout);
  interpol.setNeigh(neigh);
  interpol.setFlagEst(flag_est);
  interpol.setFlagStd(flag_std);
  interpol.setModelGeneric(model);
  interpol.setNamingConvention(namconv);

  interpol.setFlagMovMed(true);

  Id error = (interpol.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Nearest Neighbor estimation
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  flag_est    True if the estimation must be calculated
 ** \param[in]  flag_std    True if the St. Dev. must be calculated
 ** \param[in]  model       Model structure (used for St. Dev.)
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
GSTLEARN_EXPORT Id nearestNeighbor(Db* dbin,
                                   Db* dbout,
                                   bool flag_est,
                                   bool flag_std,
                                   Model* model,
                                   const NamingConvention& namconv)
{
  CalcSimpleInterpolation interpol;
  interpol.setDbin(dbin);
  interpol.setDbout(dbout);
  interpol.setFlagEst(flag_est);
  interpol.setFlagStd(flag_std);
  interpol.setModelGeneric(model);

  NeighMoving neighM(false, 1, 1.e6);
  interpol.setNeigh(&neighM);
  interpol.setNamingConvention(namconv);

  interpol.setFlagNearest(true);

  Id error = (interpol.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Polynomial estimation using Least Squares
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output Db structure
 ** \param[in]  neigh       ANeigh structure
 ** \param[in]  order       Order of the polynomial
 ** \param[in]  namconv     Naming Convention
 **
 *****************************************************************************/
GSTLEARN_EXPORT Id leastSquares(Db* dbin,
                                Db* dbout,
                                ANeigh* neigh,
                                Id order,
                                const NamingConvention& namconv)
{
  CalcSimpleInterpolation interpol;
  interpol.setDbin(dbin);
  interpol.setDbout(dbout);
  interpol.setNeigh(neigh);
  interpol.setNamingConvention(namconv);

  interpol.setFlagLstSqr(true);
  interpol.setOrder(order);

  Id error = (interpol.run()) ? 0 : 1;
  return error;
}
} // namespace gstlrn