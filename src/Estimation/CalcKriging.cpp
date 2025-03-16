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
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Enum/EKrigOpt.hpp"
#include "Estimation/CalcKrigingSimpleCase.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Basic/OptDbg.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"

#include <math.h>

CalcKriging::CalcKriging(bool flag_est, bool flag_std, bool flag_varZ)
    : ACalcInterpolator(),
    _flagEst(flag_est),
    _flagStd(flag_std),
    _flagVarZ(flag_varZ),
    _calcul(EKrigOpt::POINT),
    _ndiscs(),
    _rankColCok(),
    _matLC(nullptr),
    _flagDGM(false),
    _nameCoord(),
    _flagBayes(false),
    _priorMean(),
    _priorCov(),
    _iechSingleTarget(-1),
    _verboseSingleTarget(false),
    _flagPerCell(false),
    _flagGam(false),
    _anam(nullptr),
    _flagXvalid(false),
    _flagKfold(false),
    _flagXvalidEst(0),
    _flagXvalidStd(0),
    _flagXvalidVarZ(0),
    _flagNeighOnly(false),
    _nbNeigh(5),
    _iptrEst(-1),
    _iptrStd(-1),
    _iptrVarZ(-1),
    _iptrNeigh(-1)
{
}

CalcKriging::~CalcKriging()
{
}

void CalcKriging::setCalcul(const EKrigOpt &calcul)
{
  _calcul = calcul;

  // Temporary code for incorporating DGM as a EKrigOpt option

  if (_calcul == EKrigOpt::DGM)
    setFlagDgm(true);
  else
    setFlagDgm(false);
}

bool CalcKriging::_check()
{
  if (! ACalcInterpolator::_check()) return false;

  if (! hasDbin()) return false;
  if (! hasDbout()) return false;
  if (! hasModel()) return false;
  if (! hasNeigh()) return false;
  if (getNeigh()->getType() == ENeigh::IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    return 1;
  }

  if (_flagVarZ)
  {
    if (getModel()->isNoStat())
    {
      messerr("Variance of Estimator is limited to Stationary Covariance");
      return false;
    }
  }
  if (_flagDGM)
  {
    if (! getDbout()->isGrid())
    {
      messerr("For DGM option, the argument 'dbout'  should be a Grid");
      return false;
    }
    const Model* model = dynamic_cast<const Model*>(getModel());
    if (model == nullptr)
    {
      messerr("The 'model' must be a Model (not a ModelGeneric)");
      return false;
    }
    if (! model->hasAnam())
    {
      messerr("For DGM option, the Model must have an Anamorphosis attached");
      return false;
    }
    if (! model->isChangeSupportDefined())
    {
      messerr("DGM option requires a Change of Support to be defined");
      return false;
    }
  }
  return true;
}

bool CalcKriging::_preprocess()
{
  if (!ACalcInterpolator::_preprocess()) return false;

  if (_matLC != nullptr) _setNvar(_matLC->getNRows(), true);

  int status = 1;
  if (_iechSingleTarget >= 0) status = 2;

  if (_flagEst)
  {
    _iptrEst = _addVariableDb(2, status, ELoc::UNKNOWN, 0, _getNVar(), TEST);
    if (_iptrEst < 0) return false;
  }
  if (_flagStd)
  {
    _iptrStd = _addVariableDb(2, status, ELoc::UNKNOWN, 0, _getNVar(), TEST);
    if (_iptrStd < 0) return false;
  }
  if (_flagVarZ)
  {
    _iptrVarZ = _addVariableDb(2, status, ELoc::UNKNOWN, 0, _getNVar(), TEST);
    if (_iptrVarZ < 0) return false;
  }
  if (_flagNeighOnly)
  {
    _iptrNeigh = _addVariableDb(2, status, ELoc::UNKNOWN, 0, _nbNeigh, TEST);
    if (_iptrNeigh < 0) return false;
  }

  // Centering the Data (for DGM)

  if (_flagDGM)
  {
    // Centering (only if the output file is a Grid)
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
    if (dbgrid != nullptr)
    {
      // Duplicating the coordinate variable names before centering
      _nameCoord = getDbin()->getNamesByLocator(ELoc::X);
      if (_centerDataToGrid(dbgrid)) return false;
    }
  }

  return true;
}

bool CalcKriging::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  int nvar = _getNVar();
  if (_flagXvalid)
  {
    if (_flagXvalidStd > 0)
      _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrStd, "stderr", 1,
                      false);
    else if (_flagXvalidStd < 0)
      _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrStd, "stdev", 1,
                      false);

    if (_flagXvalidEst > 0)
      _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrEst, "esterr", 1);
    else if (_flagXvalidEst < 0)
      _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrEst, "estim", 1);

    if (_flagXvalidVarZ != 0)
      _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrVarZ, "varz", 1);
  }
  else if (_flagNeighOnly)
  {
    _renameVariable(2, VectorString(), ELoc::Z, 1, _iptrNeigh, "Number", 1);
    _renameVariable(2, VectorString(), ELoc::Z, 1, _iptrNeigh + 1, "MaxDist",
                    1);
    _renameVariable(2, VectorString(), ELoc::Z, 1, _iptrNeigh + 2, "MinDist",
                    1);
    _renameVariable(2, VectorString(), ELoc::Z, 1, _iptrNeigh + 3, "NbNESect",
                    1);
    _renameVariable(2, VectorString(), ELoc::Z, 1, _iptrNeigh + 4, "NbCESect",
                    1);
  }
  else if (_flagDGM)
  {
    if (!_nameCoord.empty()) getDbin()->setLocators(_nameCoord, ELoc::X, 0);

    _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrVarZ, "varz", 1);
    _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrStd, "stdev", 1);
    _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrEst, "estim", 1);
  }
  else
  {
    if (_matLC == nullptr)
    {
      _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrVarZ, "varz", 1);
      _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrStd, "stdev", 1);
      _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrEst, "estim", 1);
    }
    else
    {
      _renameVariable(2, {"LC"}, ELoc::UNKNOWN, nvar, _iptrVarZ, "varz", 1);
      _renameVariable(2, {"LC"}, ELoc::UNKNOWN, nvar, _iptrStd, "stdev", 1);
      _renameVariable(2, {"LC"}, ELoc::UNKNOWN, nvar, _iptrEst, "estim", 1);
    }
  }

  return true;
}

void CalcKriging::_rollback()
{
  _cleanVariableDb(1);
}

void CalcKriging::_storeResultsForExport(const KrigingSystem& ksys)
{
  _ktest.ndim  = ksys.getNDim();
  _ktest.nvar  = ksys.getNVar();
  _ktest.nech  = ksys.getNech();
  _ktest.CSize = ksys.getCovSize();
  _ktest.DSize = ksys.getDriftSize();
  _ktest.nrhs  = ksys.getNrhs();
  _ktest.nbgh  = ksys.getSampleNbgh();
  _ktest.xyz   = ksys.getSampleCoordinates();
  _ktest.data  = ksys.getSampleData();
  _ktest.lhs   = ksys.getLHS();
  _ktest.lhsF  = ksys.getLHSF();
  _ktest.rhs   = ksys.getRHS();
  _ktest.rhsF  = ksys.getRHSF();
  _ktest.wgt   = ksys.getWeights();
  _ktest.mu    = ksys.getMu();
  _ktest.var   = ksys.getVariance();
}

/****************************************************************************/
/*!
 **  Standard Kriging
 **
 ** \return  Error return code
 **
 *****************************************************************************/
bool CalcKriging::_run()
{
  /* Setting options */

  KrigingSystem ksys(getDbin(), getDbout(), getModel(), getNeigh());
  if (ksys.updKrigOptEstim(_iptrEst, _iptrStd, _iptrVarZ)) return false;
  if (ksys.setKrigOptCalcul(_calcul, _ndiscs, _flagPerCell)) return false;
  if (ksys.setKrigOptColCok(_rankColCok)) return false;
  if (ksys.setKrigOptMatLC(_matLC)) return false;
  if (_flagDGM)
  {
    if (ksys.setKrigOptDGM(true)) return false;
  }
  if (_flagBayes)
  {
    ksys.setKrigOptBayes(true, _priorMean, _priorCov);
  }
  if (_flagGam)
  {
    if (ksys.setKrigOptAnamophosis(_anam)) return false;
  }
  if (_flagXvalid)
  {
    if (ksys.setKrigOptXValid(true, _flagKfold, _flagXvalidEst > 0,
                              _flagXvalidStd > 0, _flagXvalidVarZ != 0))
      return false;
  }
  if (_flagNeighOnly)
  {
    if (ksys.updKrigOptNeighOnly(_iptrNeigh)) return false;
  }
  if (!ksys.isReady()) return false;

  /***************************************/
  /* Loop on the targets to be processed */
  /***************************************/

  for (int iech_out = 0, nech_out = getDbout()->getNSample(); iech_out < nech_out; iech_out++)
  {
    if (_iechSingleTarget > 0)
    {
      if (iech_out != _iechSingleTarget) continue;
      if (_verboseSingleTarget) OptDbg::defineAll();
    }
    else
    {
      mes_process("Kriging sample", getDbout()->getNSample(), iech_out);
    }

    bool error = ksys.estimate(iech_out);

    if (_iechSingleTarget > 0)
    {
      if (_verboseSingleTarget) OptDbg::undefineAll();
    }
    if (error) return false;
  }

  // Store the results in an API structure (only if flagSingleTarget)

  if (_iechSingleTarget >= 0) _storeResultsForExport(ksys);

  ksys.conclusion();

  return true;
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
 ** \param[in]  calcul      Kriging calculation option (EKrigOpt)
 ** \param[in]  ndiscs      Array giving the discretization counts
 ** \param[in]  flag_est    Option for storing the estimation
 ** \param[in]  flag_std    Option for storing the standard deviation
 ** \param[in]  flag_varz   Option for storing the variance of the estimator
 **                         (only available for stationary model)
 ** \param[in]  rank_colcok Option for running Collocated Cokriging
 ** \param[in]  matLC       Matrix of linear combination (or NULL)
 **                         (Dimension: nvarLC * model->getNVar())
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int kriging(Db* dbin,
            Db* dbout,
            ModelGeneric* model,
            ANeigh* neigh,
            const EKrigOpt& calcul,
            bool flag_est,
            bool flag_std,
            bool flag_varz,
            const VectorInt& ndiscs,
            const VectorInt& rank_colcok,
            const MatrixRectangular* matLC,
            const NamingConvention& namconv)
{
  NeighUnique* neighUnique = dynamic_cast<NeighUnique*>(neigh);
  if (calcul == EKrigOpt::POINT && rank_colcok.empty() && 
      matLC == nullptr && neighUnique != nullptr &&
      model->getNVar() == 1)
  { 
  CalcKrigingSimpleCase krige(flag_est, flag_std, flag_varz);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeigh(neigh);
  krige.setCalcul(calcul);
  krige.setNamingConvention(namconv);
  return 1 - krige.run();
  }

  CalcKriging krige(flag_est, flag_std, flag_varz);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeigh(neigh);
  krige.setNamingConvention(namconv);

  krige.setCalcul(calcul);
  krige.setNdisc(ndiscs);
  krige.setRankColCok(rank_colcok);
  krige.setMatLC(matLC);

  // Run the calculator
  int error = 1 - krige.run();

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
 ** \param[in]  ndiscs      Array giving the discretization counts
 ** \param[in]  flag_est    Option for the storing the estimation
 ** \param[in]  flag_std    Option for the storing the standard deviation
 ** \param[in]  rank_colcok Option for running Collocated Cokriging
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int krigcell(Db* dbin,
             Db* dbout,
             ModelGeneric* model,
             ANeigh* neigh,
             bool flag_est,
             bool flag_std,
             const VectorInt& ndiscs,
             const VectorInt& rank_colcok,
             const NamingConvention& namconv)
{
  CalcKriging krige(flag_est, flag_std, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeigh(neigh);
  krige.setNamingConvention(namconv);

  krige.setCalcul(EKrigOpt::BLOCK);
  krige.setNdisc(ndiscs);
  krige.setRankColCok(rank_colcok);
  krige.setFlagPerCell(true);

  // Run the calculator
  int error = (krige.run()) ? 0 : 1;
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
 ** \param[in]  prior_mean Array giving the prior means for the drift terms
 ** \param[in]  prior_cov  Array containing the prior covariance matrix
 **                        for the drift terms
 ** \param[in]  flag_est   Pointer for the storing the estimation
 ** \param[in]  flag_std   Pointer for the storing the standard deviation
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int kribayes(Db* dbin,
             Db* dbout,
             ModelGeneric* model,
             ANeigh* neigh,
             const VectorDouble& prior_mean,
             const MatrixSquareSymmetric& prior_cov,
             bool flag_est,
             bool flag_std,
             const NamingConvention& namconv)
{
  CalcKriging krige(flag_est, flag_std, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeigh(neigh);
  krige.setNamingConvention(namconv);

  krige.setFlagBayes(true);
  krige.setPriorMean(prior_mean);
  krige.setPriorCov(prior_cov);

  // Run the calculator
  int error = (krige.run()) ? 0 : 1;
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
 ** \param[in]  calcul      Kriging calculation option (EKrigOpt)
 ** \param[in]  ndiscs      Array giving the discretization counts
 ** \param[in]  flagPerCell Use local block extensions (when defined)
 ** \param[in]  verbose     When TRUE, the full debugging flag is switched ON
 **                         (the current status is reset after the run)
 **
 *****************************************************************************/
Krigtest_Res krigtest(Db* dbin,
                      Db* dbout,
                      ModelGeneric* model,
                      ANeigh* neigh,
                      int iech0,
                      const EKrigOpt& calcul,
                      const VectorInt& ndiscs,
                      bool flagPerCell,
                      bool verbose)
{
  CalcKriging krige(true, true, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeigh(neigh);

  krige.setCalcul(calcul);
  krige.setNdisc(ndiscs);
  krige.setIechSingleTarget(iech0);
  krige.setVerboseSingleTarget(verbose);
  krige.setFlagPerCell(flagPerCell);

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
int kriggam(Db* dbin,
            Db* dbout,
            ModelGeneric* model,
            ANeigh* neigh,
            AAnam* anam,
            const NamingConvention& namconv)
{
  CalcKriging krige(true, true, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeigh(neigh);
  krige.setNamingConvention(namconv);

  krige.setFlagGam(true);
  krige.setAnam(anam);

  // Run the calculator
  int error = (krige.run()) ? 0 : 1;
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
 * @param rank_colcok Option for running Collocated Cokriging
 * @param namconv Naming Convention
 * @return Error return code
 */
int xvalid(Db* db,
           ModelGeneric* model,
           ANeigh* neigh,
           bool flag_kfold,
           int flag_xvalid_est,
           int flag_xvalid_std,
           int flag_xvalid_varz,
           const VectorInt& rank_colcok,
           const NamingConvention& namconv)
{
  CalcKriging krige(flag_xvalid_est != 0, flag_xvalid_std != 0,
                    flag_xvalid_varz != 0);
  krige.setDbin(db);
  krige.setDbout(db);
  krige.setModel(model);
  krige.setNeigh(neigh);
  krige.setNamingConvention(namconv);

  krige.setFlagXvalid(true);
  krige.setFlagXvalidEst(flag_xvalid_est);
  krige.setFlagXvalidStd(flag_xvalid_std);
  krige.setFlagXvalidVarZ(flag_xvalid_varz);
  krige.setFlagKfold(flag_kfold);
  krige.setRankColCok(rank_colcok);

  // Run the calculator
  int error = (krige.run()) ? 0 : 1;
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
int test_neigh(Db* dbin,
               Db* dbout,
               ModelGeneric* model,
               ANeigh* neigh,
               const NamingConvention& namconv)
{
  CalcKriging krige(false, false, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeigh(neigh);
  krige.setNamingConvention(namconv);

  krige.setFlagNeighOnly(true);

  // Run the calculator
  int error = (krige.run()) ? 0 : 1;
  return error;
}
