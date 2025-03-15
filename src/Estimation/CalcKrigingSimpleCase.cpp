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
#include "Estimation/CalcKrigingSimpleCase.hpp"
#include "Enum/EKrigOpt.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Basic/OptDbg.hpp"
#include "Model/Model.hpp"
#include "Neigh/NeighUnique.hpp"

#include <math.h>

CalcKrigingSimpleCase::CalcKrigingSimpleCase(bool flag_est, bool flag_std, bool flag_varZ)
    : ACalcInterpolator(),
    _flagEst(flag_est),
    _flagStd(flag_std),
    _flagVarZ(flag_varZ),
    _calcul(EKrigOpt::POINT),
    _nameCoord(),
    _iechSingleTarget(-1),
    _flagNeighOnly(false),
    _nbNeigh(5),
    _iptrEst(-1),
    _iptrStd(-1),
    _iptrVarZ(-1),
    _iptrNeigh(-1)
{
}

CalcKrigingSimpleCase::~CalcKrigingSimpleCase()
{
}

void CalcKrigingSimpleCase::setCalcul(const EKrigOpt &calcul)
{
  _calcul = calcul;
}

bool CalcKrigingSimpleCase::_check()
{
  if (! ACalcInterpolator::_check()) return false;

  if (! hasDbin()) return false;
  if (! hasDbout()) return false;
  if (! hasModel()) return false;
  if (! hasNeigh()) return false;


  if (_flagVarZ)
  {
    if (getModel()->isNoStat())
    {
      messerr("Variance of Estimator is limited to Stationary Covariance");
      return false;
    }
  }
 
  return true;
}

bool CalcKrigingSimpleCase::_preprocess()
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

  return true;
}

bool CalcKrigingSimpleCase::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  int nvar = _getNVar();
  
  if (_flagNeighOnly)
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
  else
  {
      _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrVarZ, "varz", 1);
      _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrStd, "stdev", 1);
      _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrEst, "estim", 1);
  }
  return true;
}

void CalcKrigingSimpleCase::_rollback()
{
  _cleanVariableDb(1);
}

void CalcKrigingSimpleCase::_storeResultsForExport(const KrigingSystem& ksys)
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
bool CalcKrigingSimpleCase::_run()
{
  /* Setting options */

  KrigingSystem ksys(getDbin(), getDbout(), getModel(), getNeigh());
  if (ksys.updKrigOptEstim(_iptrEst, _iptrStd, _iptrVarZ)) return false;
  if (ksys.setKrigOptCalcul(_calcul, 0, false)) return false;

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
    }
    else
    {
      mes_process("Kriging sample", getDbout()->getNSample(), iech_out);
    }

    bool error = ksys.estimate(iech_out);

    if (error) return false;
  }

  // Store the results in an API structure (only if flagSingleTarget)

  if (_iechSingleTarget >= 0) _storeResultsForExport(ksys);

  ksys.conclusion();

  return true;
}
