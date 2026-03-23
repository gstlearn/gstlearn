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
#include "Estimation/CalcKrigingSimpleCase.hpp"
#include "Basic/OptCustom.hpp"
#include "Basic/VectorNumT.hpp"
#include "Db/Db.hpp"
#include "Estimation/KrigingAlgebraSimpleCase.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Estimation/KrigingSystemSimpleCase.hpp"
#include "Model/Model.hpp"
#include "Neigh/ANeigh.hpp"
#include "Neigh/NeighUnique.hpp"

#include <cmath>
#include <omp.h>

namespace gstlrn
{
  CalcKrigingSimpleCase::CalcKrigingSimpleCase(bool flag_est,
                                               bool flag_std,
                                               bool flag_varZ)
    : ACalcInterpolator()
    , _flagEst(flag_est)
    , _flagStd(flag_std)
    , _flagVarZ(flag_varZ)
    , _nameCoord()
    , _iechSingleTarget(-1)
    , _iptrEst(-1)
    , _iptrStd(-1)
    , _iptrVarZ(-1)
  {
  }

  CalcKrigingSimpleCase::~CalcKrigingSimpleCase() {}

  bool CalcKrigingSimpleCase::_check()
  {
    if (!ACalcInterpolator::_check()) return false;

    if (!hasDbin()) return false;
    if (!hasDbout()) return false;
    if (!hasModelGeneric()) return false;
    if (!hasNeigh()) return false;

    if (_flagVarZ)
    {
      if (getModelGeneric()->isNoStat())
      {
        messerr(
          "Variance of Estimator is limited to Stationary Covariance"); // Why?
        messerr(
          "Variance of Estimator is limited to Stationary Covariance"); // Why?
        return false;
      }
    }

    return true;
  }

  bool CalcKrigingSimpleCase::_preprocess()
  {
    if (!ACalcInterpolator::_preprocess()) return false;

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

    return true;
  }

  bool CalcKrigingSimpleCase::_postprocess()
  {
    /* Free the temporary variables */
    _cleanVariableDb(2);

    auto nvar = _getNVar();

    _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrVarZ, "varz", 1);
    _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrStd, "stdev", 1);
    _renameVariable(2, VectorString(), ELoc::Z, nvar, _iptrEst, "estim", 1);

    return true;
  }

  void CalcKrigingSimpleCase::_rollback()
  {
    _cleanVariableDb(1);
  }

  void CalcKrigingSimpleCase::_storeResultsForExport(
    const KrigingSystemSimpleCase& ksys,
    KrigingAlgebraSimpleCase& algebra,
    Id iechout)
  {
    _ktest.ndim = ksys.getNDim();
    _ktest.nvar = ksys.getNVar();
    _ktest.xyz = ksys.getSampleCoordinates(algebra, iechout);
    //_ktest.lhs   = ksys.getLHS();
    _ktest.wgt = gstlrn::KrigingSystemSimpleCase::getWeights(algebra);
    _ktest.mu = gstlrn::KrigingSystemSimpleCase::getMu(algebra);
    // _ktest.var   = ksys.getVariance();
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

    KrigingSystemSimpleCase ksys(getDbin(),
                                 getDbout(),
                                 getModelGeneric(),
                                 getNeigh());
    if (ksys.updKrigOptEstim(_iptrEst, _iptrStd, _iptrVarZ)) return false;
    VectorDouble tabwork(getDbin()->getNSample());
    if (!ksys.isReady(tabwork)) return false;

    /***************************************/
    /* Loop on the targets to be processed */
    /***************************************/

    KrigingAlgebraSimpleCase algebra(ksys.getAlgebra());
    bool use_parallel = !getModelGeneric()->isNoStat();
    Id nech_out = getDbout()->getNSample();
    auto nbthread = static_cast<I32>(
      OptCustom::query("ompthreads",
                       1)); // TODO : would like to use more threads

    SpacePoint pin(getModelGeneric()->getSpace());
    SpacePoint pout(getModelGeneric()->getSpace());
    ModelGeneric model(*ksys.getModel());
    auto ndim = getModelGeneric()->getSpace()->getNDim();
    const VectorVectorDouble coords = getDbout()->getAllCoordinates();
    static ANeigh* neigh = nullptr;

    const double* selcol = getDbout()->getColumnPtr(ELoc::SEL, 0);
    bool hassel = getDbout()->hasLocator(ELoc::SEL);
#pragma omp threadprivate(neigh)
#pragma omp parallel for firstprivate(pin, pout, tabwork, algebra, model)      \
  schedule(guided) num_threads(nbthread) if (use_parallel)
    for (Id iech_out = 0; iech_out < nech_out; iech_out++)
    {
      if (hassel && !selcol[iech_out]) continue; // Skip non-selected targets
      if (neigh == nullptr)
      {
        neigh = static_cast<ANeigh*>(getNeigh()->clone());
        getDbout()->initThread();
      }
      else
      {
        neigh->reset();
      }
      // TODO : encapsulate in Db (threadsafe)
      for (Id idim = 0; idim < static_cast<Id>(ndim); idim++)
      {
        pin.setCoord(idim, coords[idim][iech_out]);
      }
      ksys.estimate(iech_out, pin, pout, tabwork, algebra, model, neigh);

      // Store the results in an API structure (only if flagSingleTarget)
      if (_iechSingleTarget >= 0)
        _storeResultsForExport(ksys, algebra, iech_out);
    }

#pragma omp parallel num_threads(nbthread) if (use_parallel)
    {
      delete neigh;
      neigh = nullptr;
    }
    ksys.conclusion();

    return true;
  }
} // namespace gstlrn
