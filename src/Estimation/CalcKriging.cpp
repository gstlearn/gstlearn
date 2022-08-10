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
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Estimation/KrigingSystem.hpp"

#include <math.h>

CalcKriging::CalcKriging(bool flag_est, bool flag_std, bool flag_varZ)
    : ACalcInterpolator(),
    _flagEst(flag_est),
    _flagStd(flag_std),
    _flagVarZ(flag_varZ),
    _calcul(EKrigOpt::PONCTUAL),
    _ndisc(),
    _rankColCok(),
    _matCL(),
    _flagDGM(false),
    _rCoeff(1.),
    _flagBayes(false),
    _priorMean(),
    _priorCov(),
    _flagProf(false),
    _flagSingleTarget(false),
    _iechSingleTarget(-1),
    _flagPerCell(false),
    _iptrEst(-1),
    _iptrStd(-1),
    _iptrVarZ(-1)
{
}

CalcKriging::~CalcKriging()
{
}

bool CalcKriging::_check()
{
  if (! ACalcInterpolator::_check()) return false;

  if (! hasDbin()) return false;
  if (! hasDbout()) return false;
  if (! hasModel()) return false;
  if (! hasNeighParam()) return false;
  if (getNeighparam()->getType() == ENeigh::IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    return 1;
  }
  if (_flagDGM && ! getDbout()->isGrid())
  {
    messerr("For DGM option, the argument 'dbout'  should be a Grid");
    return false;
  }
  return true;
}

bool CalcKriging::_preprocess()
{
  int status = 1;
  if (_flagSingleTarget) status = 2;

  if (_flagEst)
  {
    _iptrEst = _addVariableDb(2, status, ELoc::UNKNOWN, _getNVar(), 0.);
    if (_iptrEst < 0) return false;
  }
  if (_flagStd)
  {
    _iptrStd = _addVariableDb(2, status, ELoc::UNKNOWN, _getNVar(), 0.);
    if (_iptrStd < 0) return false;
  }
  if (_flagVarZ)
  {
    _iptrVarZ = _addVariableDb(2, status, ELoc::UNKNOWN, _getNVar(), 0.);
    if (_iptrVarZ < 0) return false;
  }

  // Centering the Data (for DGM)

  if (_flagDGM)
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
    if (db_center_point_to_grid(getDbin(), dbgrid)) return false;
  }

  return true;
}

bool CalcKriging::_postprocess()
{
  int nvar = _getNVar();
  _renameVariable(nvar, _iptrVarZ, "varz", 1);
  _renameVariable(nvar, _iptrStd, "stdev", 1);
  _renameVariable(nvar, _iptrEst, "estim", 1);
  return true;
}

void CalcKriging::_rollback()
{
  _cleanVariableDb(1);
}

int CalcKriging::_getNVar() const
{
  int nvar = _matCL.empty() ? getModel()->getVariableNumber() : _matCL.size();
  return nvar;
}

void CalcKriging::_storeResultsForExport(const KrigingSystem& ksys)
{
  /* Extract relevant information */

  _ktest.ndim = ksys.getNDim();
  _ktest.nech = ksys.getNRed();
  _ktest.nrhs = 1;
  _ktest.neq  = ksys.getNeq();
  _ktest.nbgh = ksys.getSampleIndices();
  _ktest.xyz  = ksys.getSampleCoordinates();
  _ktest.data = ksys.getSampleData();
  _ktest.zam  = ksys.getZam();
  _ktest.lhs  = ksys.getLHS();
  _ktest.rhs  = ksys.getRHSC();
  _ktest.wgt  = ksys.getWeights();
  _ktest.var  = ksys.getVariance();
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

  KrigingSystem ksys(getDbin(), getDbout(), getModel(), getNeighparam());
  if (ksys.updKrigOptEstim(_iptrEst, _iptrStd, _iptrVarZ)) return false;
  if (ksys.setKrigOptCalcul(_calcul, _ndisc, _flagPerCell)) return false;
  if (ksys.setKrigOptColCok(_rankColCok)) return false;
  if (ksys.setKrigOptMatCL(_matCL)) return false;
  if (_flagDGM)
  {
    if (ksys.setKrigOptDGM(true, _rCoeff)) return false;
  }
  if (_flagBayes)
  {
    ksys.setKrigOptBayes(true, _priorMean, _priorCov);
  }
  if (_flagProf)
  {
    if (ksys.setKrigoptCode(true)) return false;
  }
  if (! ksys.isReady()) return false;

  /* Loop on the targets to be processed */

  for (int iech_out = 0; iech_out < getDbout()->getSampleNumber(); iech_out++)
  {
    if (_flagSingleTarget)
    {
      if (iech_out != _iechSingleTarget) continue;
    }
    else
    {
      mes_process("Kriging sample", getDbout()->getSampleNumber(), iech_out);
    }
    if (ksys.estimate(iech_out)) return false;
  }

  // Store the results in an API structure (only if flagSingleTarget)

  if (_flagSingleTarget)
    _storeResultsForExport(ksys);

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
 ** \param[in]  model       Model structure
 ** \param[in]  neighparam  ANeighParam structure
 ** \param[in]  calcul      Kriging calculation option (EKrigOpt)
 ** \param[in]  ndisc       Array giving the discretization counts
 ** \param[in]  flag_est    Option for storing the estimation
 ** \param[in]  flag_std    Option for storing the standard deviation
 ** \param[in]  flag_varz   Option for storing the variance of the estimator
 ** \param[in]  rank_colcok Option for running Collocated Cokriging
 ** \param[in]  matCL       Matrix of linear combination (or NULL)
 **                         (Dimension: nvarCL * model->getNVar())
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int kriging(Db *dbin,
            Db *dbout,
            Model *model,
            ANeighParam *neighparam,
            const EKrigOpt &calcul,
            bool flag_est,
            bool flag_std,
            bool flag_varz,
            VectorInt ndisc,
            VectorInt rank_colcok,
            VectorVectorDouble matCL,
            const NamingConvention& namconv)
{
  CalcKriging krige(flag_est, flag_std, flag_varz);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeighparam(neighparam);
  krige.setNamingConvention(namconv);

  krige.setCalcul(calcul);
  krige.setNdisc(ndisc);
  krige.setRankColCok(rank_colcok);
  krige.setMatCl(matCL);

  // Run the calculator
  int error = (krige.run()) ? 0 : 1;
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
 ** \param[in]  model       Model structure
 ** \param[in]  neighparam  ANeighParam structure
 ** \param[in]  ndisc       Array giving the discretization counts
 ** \param[in]  flag_est    Option for the storing the estimation
 ** \param[in]  flag_std    Option for the storing the standard deviation
 ** \param[in]  rank_colcok Option for running Collocated Cokriging
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int krigcell(Db *dbin,
             Db *dbout,
             Model *model,
             ANeighParam *neighparam,
             bool flag_est,
             bool flag_std,
             VectorInt ndisc,
             VectorInt rank_colcok,
             const NamingConvention& namconv)
{
  CalcKriging krige(flag_est, flag_std, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeighparam(neighparam);
  krige.setNamingConvention(namconv);

  krige.setCalcul(EKrigOpt::BLOCK);
  krige.setNdisc(ndisc);
  krige.setRankColCok(rank_colcok);
  krige.setFlagPerCell(true);

  // Run the calculator
  int error = (krige.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Kriging in the Gaussian Discrete Model
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin        Input Db structure
 ** \param[in]  dbout       Output DbGrid structure
 ** \param[in]  model       Model structure
 ** \param[in]  neighparam  ANeighParam structure
 ** \param[in]  flag_est    Option for storing the estimation
 ** \param[in]  flag_std    Option for storing the standard deviation
 ** \param[in]  flag_varz   Option for storing the variance of the estimator
 ** \param[in]  rval        Change of support coefficient
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int krigdgm(Db *dbin,
            DbGrid *dbout,
            Model *model,
            ANeighParam *neighparam,
            bool flag_est,
            bool flag_std,
            bool flag_varz,
            double rval,
            const NamingConvention& namconv)
 {
  CalcKriging krige(flag_est, flag_std, flag_varz);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeighparam(neighparam);
  krige.setNamingConvention(namconv);

  krige.setFlagDgm(true);
  krige.setRCoeff(rval);

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
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  prior_mean Array giving the prior means for the drift terms
 ** \param[in]  prior_cov  Array containing the prior covariance matrix
 **                        for the drift terms
 ** \param[in]  flag_est   Pointer for the storing the estimation
 ** \param[in]  flag_std   Pointer for the storing the standard deviation
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int kribayes(Db *dbin,
             Db *dbout,
             Model *model,
             ANeighParam *neighparam,
             const VectorDouble& prior_mean,
             const VectorDouble& prior_cov,
             bool flag_est,
             bool flag_std,
             const NamingConvention& namconv)
{
  CalcKriging krige(flag_est, flag_std, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeighparam(neighparam);
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
 **  Punctual Kriging based on profiles
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure
 ** \param[in]  dbout      output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  flag_est   Option for the storing the estimation
 ** \param[in]  flag_std   Option for the storing the standard deviation
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int krigprof(Db *dbin,
             Db *dbout,
             Model *model,
             ANeighParam *neighparam,
             bool flag_est,
             bool flag_std,
             const NamingConvention& namconv)
{
  CalcKriging krige(flag_est, flag_std, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeighparam(neighparam);
  krige.setNamingConvention(namconv);

  krige.setFlagProf(true);

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
 ** \param[in]  dbin       input Db structure
 ** \param[in]  dbout      output Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  iech0      Rank of the target sample
 ** \param[in]  calcul     Kriging calculation option (EKrigOpt)
 ** \param[in]  ndisc      Array giving the discretization counts
 **
 *****************************************************************************/
Krigtest_Res krigtest(Db *dbin,
                      Db *dbout,
                      Model *model,
                      ANeighParam *neighparam,
                      int iech0,
                      const EKrigOpt &calcul,
                      VectorInt ndisc)
{
  CalcKriging krige(true, true, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeighparam(neighparam);

  krige.setCalcul(calcul);
  krige.setNdisc(ndisc);
  krige.setFlagSingleTarget(true);
  krige.setIechSingleTarget(iech0);

  (void) krige.run();

  return krige.getKtest();
}

/**
 * Dummy function inserted here to check VectorDouble argument.
 * @param toto
 * @return
 */
int krigingBidon(const VectorDouble& toto)
{
  if (! toto.empty()) ut_vector_display("mon titre", toto);
  int error = 0;
  return error;
}
