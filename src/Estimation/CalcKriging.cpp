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
#include "Estimation/CalcKriging.hpp"
#include "Basic/OptDbg.hpp"
#include "Db/Db.hpp"
#include "Db/DbGrid.hpp"
#include "Estimation/CalcKrigingSimpleCase.hpp"
#include "Estimation/KrigOpt.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeigh.hpp"
#include "Neigh/NeighBench.hpp"
#include "Neigh/NeighUnique.hpp"

#include <cmath>

namespace gstlrn
{

  CalcKriging::CalcKriging(bool flag_est, bool flag_std, bool flag_varZ)
    : ACalcInterpolator()
    , _flagEst(flag_est)
    , _flagStd(flag_std)
    , _flagVarZ(flag_varZ)
    , _nameCoord()
    , _flagBayes(false)
    , _iechSingleTarget(-1)
    , _verboseSingleTarget(false)
    , _flagGam(false)
    , _anam(nullptr)
    , _flagXvalid(false)
    , _flagKfold(false)
    , _flagXvalidEst(0)
    , _flagXvalidStd(0)
    , _flagXvalidVarZ(0)
    , _flagNeighOnly(false)
    , _nbNeigh(5)
    , _iptrEst(-1)
    , _iptrStd(-1)
    , _iptrVarZ(-1)
    , _iptrNeigh(-1)
  {
  }

  CalcKriging::~CalcKriging() {}

  bool CalcKriging::_check()
  {
    if (!ACalcInterpolator::_check()) return false;

    if (!hasDbin()) return false;
    if (!hasDbout()) return false;
    if (!hasModelGeneric()) return false;
    if (!hasNeigh()) return false;
    if (getNeigh()->getType() == ENeigh::IMAGE)
    {
      messerr("This tool cannot function with an IMAGE neighborhood");
      return 1;
    }

    if (_flagVarZ)
    {
      if (getModelGeneric()->isNoStat())
      {
        messerr("Variance of Estimator is limited to Stationary Covariance");
        return false;
      }
    }
    return true;
  }

  bool CalcKriging::_preprocess()
  {
    if (!ACalcInterpolator::_preprocess()) return false;

    if (getKrigopt().hasMatLC()) _setNvar(getKrigopt().getMatLCNRows(), true);

    Id status = 1;
    if (_iechSingleTarget >= 0) status = 2;

    if (_flagEst)
    {
      _iptrEst =
        _addVariableDb(2, status, ELoc::UNDEFINED, 0, _getNVar(), TEST);
      if (_iptrEst < 0) return false;
    }
    if (_flagStd)
    {
      _iptrStd =
        _addVariableDb(2, status, ELoc::UNDEFINED, 0, _getNVar(), TEST);
      if (_iptrStd < 0) return false;
    }
    if (_flagVarZ)
    {
      _iptrVarZ =
        _addVariableDb(2, status, ELoc::UNDEFINED, 0, _getNVar(), TEST);
      if (_iptrVarZ < 0) return false;
    }
    if (_flagNeighOnly)
    {
      _iptrNeigh =
        _addVariableDb(2, status, ELoc::UNDEFINED, 0, _nbNeigh, TEST);
      if (_iptrNeigh < 0) return false;
    }

    // Centering the Data (for DGM only)

    if (getKrigopt().hasFlagDGM())
    {
      // Centering (only if the output file is a Grid)
      auto* dbgrid = dynamic_cast<DbGrid*>(getDbout());
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

    auto nvar = _getNVar();
    if (_flagXvalid)
    {
      if (_flagXvalidStd > 0)
        _renameVariable(
          2, VectorString(), ELoc::Z, nvar, _iptrStd, "stderr", 1, false);
      else if (_flagXvalidStd < 0)
        _renameVariable(
          2, VectorString(), ELoc::Z, nvar, _iptrStd, "stdev", 1, false);

      if (_flagXvalidEst > 0)
        _renameVariable(
          2, VectorString(), ELoc::Z, nvar, _iptrEst, "esterr", 1);
      else if (_flagXvalidEst < 0)
        _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrEst, "estim", 1);

      if (_flagXvalidVarZ != 0)
        _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrVarZ, "varz", 1);
    }
    else if (_flagNeighOnly)
    {
      _renameVariable(2, VectorString(), ELoc::Z, 1, _iptrNeigh, "Number", 1);
      _renameVariable(
        2, VectorString(), ELoc::Z, 1, _iptrNeigh + 1, "MaxDist", 1);
      _renameVariable(
        2, VectorString(), ELoc::Z, 1, _iptrNeigh + 2, "MinDist", 1);
      _renameVariable(
        2, VectorString(), ELoc::Z, 1, _iptrNeigh + 3, "NbNESect", 1);
      _renameVariable(
        2, VectorString(), ELoc::Z, 1, _iptrNeigh + 4, "NbCESect", 1);
    }
    else if (getKrigopt().hasFlagDGM())
    {
      if (!_nameCoord.empty()) getDbin()->setLocators(_nameCoord, ELoc::X, 0);

      _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrVarZ, "varz", 1);
      _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrStd, "stdev", 1);
      _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrEst, "estim", 1);
    }
    else
    {
      if (!getKrigopt().hasMatLC())
      {
        _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrVarZ, "varz", 1);
        _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrStd, "stdev", 1);
        _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrEst, "estim", 1);
      }
      else
      {
        _renameVariable(2, {"LC"}, ELoc::UNDEFINED, nvar, _iptrVarZ, "varz", 1);
        _renameVariable(2, {"LC"}, ELoc::UNDEFINED, nvar, _iptrStd, "stdev", 1);
        _renameVariable(2, {"LC"}, ELoc::UNDEFINED, nvar, _iptrEst, "estim", 1);
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
    _ktest.ndim = ksys.getNDim();
    _ktest.nvar = ksys.getNVar();
    _ktest.nech = ksys.getNech();
    _ktest.CSize = ksys.getCovSize();
    _ktest.DSize = ksys.getDriftSize();
    _ktest.nrhs = ksys.getNrhs();
    _ktest.nbgh = ksys.getSampleNbgh();
    _ktest.xyz = ksys.getSampleCoordinates();
    _ktest.data = ksys.getSampleData();
    _ktest.lhs = ksys.getLHS();
    _ktest.lhsF = ksys.getLHSF();
    _ktest.rhs = ksys.getRHS();
    _ktest.rhsF = ksys.getRHSF();
    _ktest.wgt = ksys.getWeights();
    _ktest.mu = ksys.getMu();
    _ktest.var = ksys.getVariance();
  }

  /**
   * @brief The Kriging weights are added to the input Db
   * They are multiplied by 100 for better legibility
   *
   * @param namconv
   */
  void CalcKriging::addWeights(const NamingConvention& namconv)
  {
    Id nvar = _getNVar();
    Id nech = _ktest.nech;
    auto names = getDbin()->getNamesByLocator(ELoc::Z);
    auto iuid = getDbin()->addColumnsByConstant(nvar * nvar, TEST, "Weight");

    // Loop on the variables
    for (Id ivar = 0; ivar < nvar; ++ivar)
    {
      // Loop on the data variables
      Id lec = 0;
      for (Id jvar = 0; jvar < nvar; ++jvar)

        // Loop on the samples involved in neighborhood
        for (Id iech = 0; iech < nech; ++iech)
        {
          Id jech = _ktest.nbgh[iech];

          auto value = getDbin()->getZVariable(jech, jvar);
          if (FFFF(value)) continue;

          getDbin()->setValueByUID(
            jech, iuid + jvar + nvar * ivar,
            100. * _ktest.wgt.getValue(lec++, ivar));
        }
    }
    namconv.setNamesAndLocators(
      getDbin(), names, ELoc::Z, nvar, getDbin(), iuid, String(), nvar, false);
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
    KrigingSystem ksys(
      getDbin(), getDbout(), getModelGeneric(), getNeigh(), getKrigopt());
    if (ksys.updKrigOptEstim(_iptrEst, _iptrStd, _iptrVarZ)) return false;
    if (_flagBayes)
    {
      ksys.setKrigOptBayes(true);
    }
    if (_flagGam)
    {
      if (ksys.setKrigOptAnamophosis(_anam)) return false;
    }
    if (_flagXvalid)
    {
      if (ksys.setKrigOptXValid(
            true, _flagKfold, _flagXvalidEst > 0, _flagXvalidStd > 0,
            _flagXvalidVarZ != 0))
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

    for (Id iech_out = 0, nech_out = getDbout()->getNSample();
         iech_out < nech_out; iech_out++)
    {
      if (_iechSingleTarget >= 0)
      {
        if (iech_out != _iechSingleTarget) continue;
        if (_verboseSingleTarget) OptDbg::defineAll();
      }
      else
      {
        mes_process("Kriging sample", getDbout()->getNSample(), iech_out);
      }

      bool error = ksys.estimate(iech_out);

      if (_iechSingleTarget >= 0)
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

} // namespace gstlrn
