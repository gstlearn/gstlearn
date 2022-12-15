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
#include "geoslib_old_f.h"

#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Estimation/CalcFactorKriging.hpp"
#include "Anamorphosis/AAnam.hpp"

CalcFactorKriging::CalcFactorKriging(bool flag_est, bool flag_std)
    : ACalcInterpolator(),
      _flagEst(flag_est),
      _flagStd(flag_std),
      _calcul(EKrigOpt::PONCTUAL),
      _ndisc(),
      _iptrEst(-1),
      _iptrStd(-1),
      _iuidFactors()
{
}

CalcFactorKriging::~CalcFactorKriging()
{
}

bool CalcFactorKriging::_check()
{
  // Turn the problem to Monovariate before checking consistency with 'Model'
  getDbin()->clearLocators(ELoc::Z);
  getDbin()->setLocatorByUID(_iuidFactors[0], ELoc::Z);

  if (! ACalcInterpolator::_check()) return false;

  if (! hasDbin()) return false;
  if (! hasDbout()) return false;
  if (! hasModel()) return false;
  if (! hasNeighParam()) return false;

  if (! getDbout()->isGrid())
  {
    messerr("This application is limited to output Grid file");
    return false;
  }
  if (getNeighparam()->getType() == ENeigh::IMAGE)
  {
    messerr("This tool cannot function with an IMAGE neighborhood");
    return false;
  }
  if (getModel()->getVariableNumber() != 1)
  {
    messerr("This application is limited to the monovariate Model case");
    return false;
  }
  if (! getModel()->hasAnam())
  {
    messerr("Argument 'model' should has an Anamorphosis attached");
    return false;
  }
  // If change of support is defined through the anamorphosis,
  // the calculation option (EKrigOpt) should be set to PONCTUAL
  // in order to avoid additional block randomization
  if (_calcul == EKrigOpt::BLOCK && _ndisc.empty())
  {
    messerr("For Block estimate, you must specify the discretization");
    return false;
  }
  return true;
}

bool CalcFactorKriging::_preprocess()
{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  const AAnam* anam = getModel()->getAnam();

  // Check if the change of support is defined in the Anamorphosis
  bool flag_change_support = anam->isChangeSupportDefined();

  // Centering the information (onyl when a change of support is defined)
  if (flag_change_support)
  {
    if (_ndisc.empty())
    {
      // Center the information in the blocks of the output grid
      if (db_center_point_to_grid(getDbin(), dbgrid)) return false;
    }
    if (! _ndisc.empty())
    {
      // Center the information in sub-blocks when the output grid defines panels
      DbGrid* dbsmu = db_create_grid_divider(dbgrid, _ndisc, 1);
      int error = db_center_point_to_grid(getDbin(), dbsmu);
      delete dbsmu;
      if (error) return false;
    }
  }

  if (_flagEst)
  {
    _iptrEst = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, _getNFactors(), 0.);
    if (_iptrEst < 0) return false;
  }
  if (_flagStd)
  {
    _iptrStd = _addVariableDb(2, 1, ELoc::UNKNOWN, 0, _getNFactors(), 0.);
    if (_iptrStd < 0) return false;
  }
  return true;
}

bool CalcFactorKriging::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  getDbin()->setLocatorsByUID(_iuidFactors, ELoc::Z);

  int nfactor = _getNFactors();
  _renameVariable(2, nfactor, _iptrStd, "stdev", 1);
  _renameVariable(2, nfactor, _iptrEst, "estim", 1);

  return true;
}

void CalcFactorKriging::_rollback()
{
  _cleanVariableDb(1);
}

int CalcFactorKriging::_getNFactors() const
{
  return (int) _iuidFactors.size();
}

/****************************************************************************/
/*!
 **  Standard Kriging
 **
 ** \return  Error return code
 **
 *****************************************************************************/
bool CalcFactorKriging::_run()
{
  KrigingSystem ksys(getDbin(), getDbout(), getModel(), getNeighparam());
  if (ksys.updKrigOptEstim(_iptrEst, _iptrStd, -1)) return 1;
  if (ksys.setKrigOptCalcul(_calcul, _ndisc)) return 1;
  if (ksys.setKrigOptFactorKriging(true)) return 1;
  if (! ksys.isReady()) return 1;

  /* Loop on the targets to be processed */

  for (int iech_out = 0; iech_out < getDbout()->getSampleNumber(); iech_out++)
  {
    mes_process("Disjunctive Kriging for cell",
                getDbout()->getSampleNumber(),iech_out);

    for (int iclass = 1; iclass <= _getNFactors(); iclass++)
    {
      int jptr_est = (_flagEst) ? _iptrEst + iclass - 1 : -1;
      int jptr_std = (_flagStd) ? _iptrStd + iclass - 1 : -1;
      getDbin()->clearLocators(ELoc::Z);
      getDbin()->setLocatorByUID(_iuidFactors[iclass - 1], ELoc::Z);
      if (ksys.updKrigOptEstim(jptr_est, jptr_std, -1)) return 1;
      if (ksys.updKrigOptIclass(iclass, _getNFactors())) return 1;
      if (ksys.estimate(iech_out)) return 1;
    }
  }
  return true;
}

/****************************************************************************/
/*!
 **  Disjunctive Kriging
 **
 ** \return  Error return code
 **
 ** \param[in]  dbin       input Db structure (containing the factors)
 ** \param[in]  dbgrid     output Grid Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neighparam ANeighParam structure
 ** \param[in]  calcul     Type of estimate (from EKrigopt)
 ** \param[in]  ndisc      Discretization parameters (or empty)
 ** \param[in]  flag_est   Option for the storing the estimation
 ** \param[in]  flag_std   Option for the storing the standard deviation
 ** \param[in]  namconv    Naming convention
 **
 ** \remark When the change of support is defined through the Anamorphosis
 ** \remark the 'calcul' option must be set to PONCTUAL and 'ndisc' does not
 ** \remark have to be defined
 **
 *****************************************************************************/
int KrigingFactors(Db *dbin,
                   DbGrid *dbgrid,
                   Model *model,
                   ANeighParam *neighparam,
                   const EKrigOpt &calcul,
                   const VectorInt &ndisc,
                   bool flag_est,
                   bool flag_std,
                   const NamingConvention &namconv)
{
  CalcFactorKriging krige(flag_est, flag_std);
  krige.setDbin(dbin);
  krige.setDbout(dbgrid);
  krige.setModel(model);
  krige.setNeighparam(neighparam);
  krige.setNamingConvention(namconv);

  krige.setCalcul(calcul);
  krige.setNdisc(ndisc);
  krige.setIuidFactors(dbin->getUIDsByLocator(ELoc::Z));

  // Run the calculator
  int error = (krige.run()) ? 0 : 1;
  return error;
}
