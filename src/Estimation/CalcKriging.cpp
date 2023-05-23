/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clauses                                                     */
/*                                                                            */
/******************************************************************************/
#include "geoslib_f.h"

#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Estimation/CalcKriging.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Basic/OptDbg.hpp"

#include "Matrix/MatrixEigen.hpp"
#include <Eigen/Dense>
#include <Eigen/Cholesky>

#include <math.h>

CalcKriging::CalcKriging(bool flag_est, bool flag_std, bool flag_varZ)
    : ACalcInterpolator(),
    _flagEst(flag_est),
    _flagStd(flag_std),
    _flagVarZ(flag_varZ),
    _calcul(EKrigOpt::POINT),
    _ndisc(),
    _rankColCok(),
    _matCL(),
    _flagDGM(false),
    _nameCoord(),
    _flagBayes(false),
    _priorMean(),
    _priorCov(),
    _flagProf(false),
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
  if (! hasNeighParam()) return false;
  if (getNeighparam()->getType() == ENeigh::IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    return 1;
  }
  if (_flagVarZ)
  {
    if (! getModel()->isStationary())
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
    if (! getModel()->hasAnam())
    {
      messerr("For DGM option, the Model must have an Anamorphosis attached");
      return false;
    }
    if (! getModel()->isChangeSupportDefined())
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
    DbGrid *dbgrid = dynamic_cast<DbGrid*>(getDbout());
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
      _renameVariable(2, nvar, _iptrStd, "stderr", 1, false);
    else if (_flagXvalidStd < 0)
      _renameVariable(2, nvar, _iptrStd, "stdev", 1, false);

    if (_flagXvalidEst > 0)
      _renameVariable(2, nvar, _iptrEst, "esterr", 1);
    else if (_flagXvalidEst < 0)
      _renameVariable(2, nvar, _iptrEst, "estim", 1);

    if (_flagXvalidVarZ != 0)
      _renameVariable(2, nvar, _iptrVarZ, "varz", 1);
  }
  else if (_flagNeighOnly)
  {
    _renameVariable(2, 1, _iptrNeigh  , "Number", 1);
    _renameVariable(2, 1, _iptrNeigh+1, "MaxDist", 1);
    _renameVariable(2, 1, _iptrNeigh+2, "MinDist", 1);
    _renameVariable(2, 1, _iptrNeigh+3, "NbNESect", 1);
    _renameVariable(2, 1, _iptrNeigh+4, "NbCESect", 1);
  }
  else if (_flagDGM)
  {
    if (!_nameCoord.empty())
      getDbin()->setLocators(_nameCoord, ELoc::X);

    _renameVariable(2, nvar, _iptrVarZ, "varz", 1);
    _renameVariable(2, nvar, _iptrStd, "stdev", 1);
    _renameVariable(2, nvar, _iptrEst, "estim", 1);
  }
  else
  {
    _renameVariable(2, nvar, _iptrVarZ, "varz", 1);
    _renameVariable(2, nvar, _iptrStd, "stdev", 1);
    _renameVariable(2, nvar, _iptrEst, "estim", 1);
  }

  return true;
}

void CalcKriging::_rollback()
{
  _cleanVariableDb(1);
}

int CalcKriging::_getNVar() const
{
  int nvar = (_matCL.empty() || _matCL[0].empty()) ? getModel()->getVariableNumber() : (int) _matCL.size();
  return nvar;
}

void CalcKriging::_storeResultsForExport(const KrigingSystem& ksys)
{
  int ndim = ksys.getNDim();

  /* Extract relevant information */

  _ktest.ndim = ndim;
  _ktest.nvar = 1;
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
    if (ksys.setKrigOptDGM(true)) return false;
  }
  if (_flagBayes)
  {
    ksys.setKrigOptBayes(true, _priorMean, _priorCov);
  }
  if (_flagProf)
  {
    if (ksys.setKrigoptCode(true)) return false;
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
  if (! ksys.isReady()) return false;

  /* Loop on the targets to be processed */

  for (int iech_out = 0; iech_out < getDbout()->getSampleNumber(); iech_out++)
  {
    if (_iechSingleTarget > 0)
    {
      if (iech_out != _iechSingleTarget) continue;
      if (_verboseSingleTarget) OptDbg::defineAll();
    }
    else
    {
      mes_process("Kriging sample", getDbout()->getSampleNumber(), iech_out);
    }

    bool error = ksys.estimate(iech_out);

    if (_iechSingleTarget > 0)
    {
      if (_verboseSingleTarget) OptDbg::undefineAll();
    }
    if (error) return false;
  }

  // Store the results in an API structure (only if flagSingleTarget)

  if (_iechSingleTarget >= 0)
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
 **                         (only available for stationary model)
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
 ** \param[in]  dbin        input Db structure
 ** \param[in]  dbout       output Db structure
 ** \param[in]  model       Model structure
 ** \param[in]  neighparam  ANeighParam structure
 ** \param[in]  iech0       Rank of the target sample
 ** \param[in]  calcul      Kriging calculation option (EKrigOpt)
 ** \param[in]  ndisc       Array giving the discretization counts
 ** \param[in]  flagPerCell Use local block extensions (when defined)
 ** \param[in]  verbose     When TRUE, the full debugging flag is switched ON
 **                         (the current status is reset after the run)
 **
 *****************************************************************************/
Krigtest_Res krigtest(Db *dbin,
                      Db *dbout,
                      Model *model,
                      ANeighParam *neighparam,
                      int iech0,
                      const EKrigOpt &calcul,
                      VectorInt ndisc,
                      bool flagPerCell,
                      bool verbose)
{
  CalcKriging krige(true, true, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeighparam(neighparam);

  krige.setCalcul(calcul);
  krige.setNdisc(ndisc);
  krige.setIechSingleTarget(iech0);
  krige.setVerboseSingleTarget(verbose);
  krige.setFlagPerCell(flagPerCell);

  (void) krige.run();

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
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  anam       AAnam structure
 ** \param[in]  namconv    Naming convention
 **
 *****************************************************************************/
int kriggam(Db *dbin,
            Db *dbout,
            Model *model,
            ANeighParam *neighparam,
            AAnam *anam,
            const NamingConvention& namconv)
{
  CalcKriging krige(true, true, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeighparam(neighparam);
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
 * @param model Model structure
 * @param neighparam ANeighParam structure
 * @param flag_kfold True if a code (K-FOLD) is used
 * @param flag_xvalid_est Option for storing the estimation: 1 for Z*-Z; -1 for Z*
 * @param flag_xvalid_std Option for storing the standard deviation: 1:for (Z*-Z)/S; -1 for S
 * @param flag_xvalid_varz Option for storing the variance of the estimator
 * @param rank_colcok Option for running Collocated Cokriging
 * @param namconv Naming Convention
 * @return Error return code
 */
int xvalid(Db *db,
           Model *model,
           ANeighParam *neighparam,
           bool flag_kfold,
           int flag_xvalid_est,
           int flag_xvalid_std,
           int flag_xvalid_varz,
           VectorInt rank_colcok,
           const NamingConvention& namconv)
{
  CalcKriging krige(flag_xvalid_est != 0,
                    flag_xvalid_std != 0,
                    flag_xvalid_varz != 0);
  krige.setDbin(db);
  krige.setDbout(db);
  krige.setModel(model);
  krige.setNeighparam(neighparam);
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
 ** \param[in]  model      Model structure (optional)
 ** \param[in]  neighparam ANeighParam structure
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
int test_neigh(Db *dbin,
               Db *dbout,
               Model *model,
               ANeighParam *neighparam,
               const NamingConvention &namconv)
{
  CalcKriging krige(false, false, false);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeighparam(neighparam);
  krige.setNamingConvention(namconv);
  krige.setFlagNeighOnly(true);

  // Run the calculator
  int error = (krige.run()) ? 0 : 1;
  return error;
}



void krigingExperimental(const Db *dbin,
                        Db *dbout,
                        Model *model,
						bool flag_est,
                        bool flag_std,
						bool flag_varz,
                        const NamingConvention& namconv)
{

	int nech = dbin->getSampleNumber(true);
	int nechout = dbout->getSampleNumber(true);
	int nloc = 0;

	VectorDouble res = VectorDouble(nechout);
	VectorDouble z = dbin->getColumnByLocator(ELoc::Z, 0, true, false);

	double mean = model->getMean(0);
	MatrixEigen driftsE;
	VectorDouble coeffs;

	if ( model->getDriftNumber() == 0)
	{
		VH::addConstant(z, -mean);
	}
	else
	{
		VectorVectorDouble drifts = model->getDrifts(dbin, true);
		coeffs = VectorDouble(model->getDriftNumber());
		driftsE = MatrixEigen(drifts,true);
	}

	auto C  = model->evalCovMatrixEigen(dbin);
	auto C0 = model->evalCovMatrixEigen(dbin,dbout);


	VectorDouble drift;
	VectorDouble vars;
	VectorDouble varest;

	if ( model->getDriftNumber() > 0)
	{
		VectorDouble tempv = VectorDouble(model->getDriftNumber());
		auto tempDriftsE = C.solve(driftsE); // C^{-1} F
		auto temp = MatrixEigen::prodT1(driftsE, tempDriftsE); // F'C^{-1}F   nd x nd
		tempDriftsE.prodTMatVecInPlace(z,tempv);  //                F'C^{-1}z nd x 1
		temp.solve(tempv,coeffs);	// beta = (F'C^{-1}F)^{-1}F'C^{-1}z				  nd x 1
		drift = model->evalDrifts(dbin, coeffs,0, true);		// F beta = n x 1
		VH::subtractInPlace(z, drift);							// (z - F beta) = (I - F (F'C^{-1}F)^{-1}F'C^{-1})z
																// (I - F (F'C^{-1}F)^{-1}F)
	}

	if (!flag_std && !flag_varz)
	{
		VectorDouble dual = VectorDouble(nech);
		C.solve(z,dual);
		C0.prodTMatVecInPlace(dual,res);
	}
	else
	{
		varest.resize(dbout->getSampleNumber(true));
		auto weights = C.solve(C0);
		weights.prodTMatVecInPlace(z,res);
		auto varmat = MatrixEigen::productPointwise(weights, C0);
		varmat.sumColsInPlace(varest);


	}

	if (model->getDriftNumber() == 0)
	{
		VH::addConstant(res, mean);
	}
	else
	{
		drift = model->evalDrifts(dbout, coeffs,0, true);
		VH::addInPlace(res, drift);
	}
	if (flag_est)
	{
		dbout->addColumns(res, "r_estim", ELoc::Z, nloc, true, 0.,0);
		nloc++;
	}

	if (flag_varz)
	{
		dbout->addColumns(varest,  "r_varz", ELoc::Z, nloc, true, 0.,0);
		nloc++;
	}

}



void krigingExperimentalEigen(const Db *dbin,
                        Db *dbout,
                        Model *model,
						bool flag_est,
                        bool flag_std,
						bool flag_varz,
                        const NamingConvention& namconv)
{

	int nechout = dbout->getSampleNumber(true);
	int nloc = 0;

	VectorDouble res = VectorDouble(nechout);
	Eigen::Map<Eigen::VectorXd> resE(res.data(),res.size());

	VectorDouble z = dbin->getColumnByLocator(ELoc::Z, 0, true, false);
	Eigen::Map<Eigen::VectorXd> zE(z.data(),z.size());

	double mean = model->getMean(0);

	Eigen::VectorXd coeffs;
	VectorVectorDouble drifts;
	Eigen::MatrixXd driftsE;

	if ( model->getDriftNumber() == 0)
	{
		zE.array()-=mean;
	}
	else
	{
		//// TODO : écrire Model::getDriftsEigen
		drifts = model->getDrifts(dbin, true);
		driftsE = Eigen::MatrixXd((int)drifts[0].size(),(int)drifts.size());
			for(int i = 0; i < (int)drifts.size();i++)
			{
				for (int j = 0; j < (int) drifts[0].size(); j++)
				{
					driftsE(j,i) = drifts[i][j];
				}
			}

		coeffs = Eigen::VectorXd(model->getDriftNumber());
	}

	auto Cm = model->evalCovMatrixEigen(dbin);
	auto C0m =  model->evalCovMatrixEigen(dbin,dbout);
	const auto C  = Cm.getMatrix();
	auto C0 =C0m.getMatrixM();

	Eigen::MatrixXd invC;
	Eigen::LLT<Eigen::MatrixXd> factor;

	Eigen::VectorXd driftE;
	Eigen::MatrixXd tempDriftsE;

	VectorDouble vars;
	VectorDouble varest;

	bool computeVars = flag_std || flag_varz;
	if (computeVars)
	{
		invC = C->inverse();
	}
	else
	{
		factor = C->llt();
	}

	Eigen::MatrixXd invFinvCF;
	if ( model->getDriftNumber() > 0)
	{

		if (computeVars)
			tempDriftsE = invC * driftsE; // C^{-1} F
		else
			tempDriftsE = factor.solve(driftsE);
		auto temp = driftsE.transpose() * tempDriftsE; // F'C^{-1}F   nd x nd
		auto tempv = tempDriftsE.transpose() * zE;  //                F'C^{-1}z nd x 1
		invFinvCF = temp.inverse(); //(F'C^{-1}F)^{-1}
		coeffs = invFinvCF * tempv;	// beta = (F'C^{-1}F)^{-1}F'C^{-1}z				  nd x 1
		driftE = driftsE * coeffs;		// F beta = n x 1
		zE.array()-=driftE.array();		// (z - F beta)

	}

	if (flag_varz)
	{
		varest.resize(dbout->getSampleNumber(true));
	}


	Eigen::Map<Eigen::VectorXd> varestE(varest.data(),varest.size());

	if (!computeVars)
	{
		auto dual = factor.solve(zE);
		resE.noalias() += C0->transpose() * dual;
	}
	else
	{

		vars.resize(dbout->getSampleNumber(true));
		Eigen::Map<Eigen::VectorXd> varsE(vars.data(),vars.size());
		auto weights = invC * *C0;
		resE = weights.transpose() * zE;
		if(model->getDriftNumber() > 0)
		{
			C0->array() -= (((driftsE * invFinvCF) * driftsE.transpose()) * weights).array();
		}

		auto varmat = weights.array() * C0->array();
		varestE.array() = varmat.colwise().sum().array();
	}

	if (model->getDriftNumber() == 0)
	{
		resE.array()+=mean;
	}
	else
	{
		drifts = model->getDrifts(dbout, true);
		driftsE.resize((int)drifts[0].size(),(int)drifts.size());
		for(int i = 0; i < (int)drifts.size();i++)
		{
			for (int j = 0; j < (int) drifts[0].size(); j++)
			{
				driftsE(j,i) = drifts[i][j];
			}
		}

		driftE = driftsE * coeffs;
		resE.array() += driftE.array();

		if (model->getDriftNumber() > 0)
		{
			auto tempA = driftsE.array() * ( driftsE * invFinvCF).array();
			varestE.array()+= tempA.rowwise().sum();
		}
	}
	if (flag_est)
	{
		dbout->addColumns(res, "r_estim", ELoc::Z, nloc, true, 0.,0);
		nloc++;
	}

	if (flag_varz)
	{
		dbout->addColumns(varest,  "r_varz", ELoc::Z, nloc, true, 0.,0);
		nloc++;
	}

}


