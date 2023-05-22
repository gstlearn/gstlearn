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
#include "geoslib_old_f.h"
#include "geoslib_f.h"

#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Estimation/CalcKrigingFactors.hpp"
#include "Anamorphosis/AAnam.hpp"

CalcKrigingFactors::CalcKrigingFactors(bool flag_est, bool flag_std)
    : ACalcInterpolator(),
      _flagEst(flag_est),
      _flagStd(flag_std),
      _calcul(EKrigOpt::POINT),
      _ndisc(),
      _nameCoord(),
      _iptrEst(-1),
      _iptrStd(-1),
      _iuidFactors()
{
}

CalcKrigingFactors::~CalcKrigingFactors()
{
}

bool CalcKrigingFactors::_check()
{
  // Turn the problem to Monovariate before checking consistency with 'Model'
  getDbin()->clearLocators(ELoc::Z);
  getDbin()->setLocatorByUID(_iuidFactors[0], ELoc::Z);

  if (! ACalcInterpolator::_check()) return false;

  if (! hasDbin()) return false;
  if (! hasDbout()) return false;
  if (! hasModel()) return false;
  if (! hasNeigh()) return false;

  if (getNeigh()->getType() == ENeigh::IMAGE)
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
  // the calculation option (EKrigOpt) should be set to POINT
  // in order to avoid additional block randomization
  if (_calcul == EKrigOpt::BLOCK && _ndisc.empty())
  {
    messerr("For Block estimate, you must specify the discretization");
    return false;
  }
  return true;
}

bool CalcKrigingFactors::_hasChangeSupport() const
{
  const AAnam* anam = getModel()->getAnam();
  if (anam == nullptr) return false;

  // Check if the change of support is defined in the Anamorphosis
  return anam->isChangeSupportDefined();
}

bool CalcKrigingFactors::_preprocess()
{
  // Centering the information (only when a change of support is defined)
  if (_hasChangeSupport())
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
    if (dbgrid == nullptr)
    {
      messerr("Due to change of support, 'dbout' should be a Grid");
      return false;
    }
    if (_ndisc.empty())
    {
      // Center the information in the blocks of the output grid
      // Duplicating the coordinate variable names before centering
      _nameCoord = getDbin()->getNamesByLocator(ELoc::X);
      int error = _centerDataToGrid(dbgrid);
      if (error) return false;
    }
    if (! _ndisc.empty())
    {
      // Center the information in sub-blocks when the output grid defines panels
      DbGrid* dbsmu = db_create_grid_divider(dbgrid, _ndisc, 1);
      _nameCoord = getDbin()->getNamesByLocator(ELoc::X);
      int error = _centerDataToGrid(dbsmu);
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

bool CalcKrigingFactors::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  getDbin()->setLocatorsByUID(_iuidFactors, ELoc::Z);

  int nfactor = _getNFactors();
  _renameVariable(2, nfactor, _iptrStd, "stdev", 1);
  _renameVariable(2, nfactor, _iptrEst, "estim", 1);

  // Centering the information (only when a change of support is defined)
  if (_hasChangeSupport() && ! _nameCoord.empty())
    getDbin()->setLocators(_nameCoord, ELoc::X);

  return true;
}

void CalcKrigingFactors::_rollback()
{
  _cleanVariableDb(1);
}

int CalcKrigingFactors::_getNFactors() const
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
bool CalcKrigingFactors::_run()
{
  KrigingSystem ksys(getDbin(), getDbout(), getModel(), getNeigh());
  if (ksys.updKrigOptEstim(_iptrEst, _iptrStd, -1)) return 1;
  if (ksys.setKrigOptCalcul(_calcul, _ndisc)) return 1;
  if (ksys.setKrigOptFactorKriging(true)) return 1;
  if (! ksys.isReady()) return 1;

  // Loop on the targets to be processed

  int ntotal = getDbout()->getSampleNumber() * _getNFactors();
  int iproc = 0;
  for (int iclass = 1; iclass <= _getNFactors(); iclass++)
  {
    int jptr_est = (_flagEst) ? _iptrEst + iclass - 1 : -1;
    int jptr_std = (_flagStd) ? _iptrStd + iclass - 1 : -1;
    getDbin()->clearLocators(ELoc::Z);
    getDbin()->setLocatorByUID(_iuidFactors[iclass - 1], ELoc::Z);
    if (ksys.updKrigOptEstim(jptr_est, jptr_std, -1)) return 1;
    if (ksys.updKrigOptIclass(iclass, _getNFactors())) return 1;

    for (int iech_out = 0; iech_out < getDbout()->getSampleNumber(); iech_out++)
    {
      mes_process("Disjunctive Kriging for cell", ntotal, iproc);
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
 ** \param[in]  dbout      output Grid Db structure
 ** \param[in]  model      Model structure
 ** \param[in]  neigh      ANeigh structure
 ** \param[in]  calcul     Type of estimate (from EKrigopt)
 ** \param[in]  ndisc      Discretization parameters (or empty)
 ** \param[in]  flag_est   Option for the storing the estimation
 ** \param[in]  flag_std   Option for the storing the standard deviation
 ** \param[in]  namconv    Naming convention
 **
 ** \remark When the change of support is defined through the Anamorphosis
 ** \remark the 'calcul' option must be set to POINT and 'ndisc' does not
 ** \remark have to be defined
 **
 *****************************************************************************/
int krigingFactors(Db *dbin,
                   Db *dbout,
                   Model *model,
                   ANeigh *neigh,
                   const EKrigOpt &calcul,
                   const VectorInt &ndisc,
                   bool flag_est,
                   bool flag_std,
                   const NamingConvention &namconv)
{
  CalcKrigingFactors krige(flag_est, flag_std);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeigh(neigh);
  krige.setNamingConvention(namconv);

  krige.setCalcul(calcul);
  krige.setNdisc(ndisc);
  krige.setIuidFactors(dbin->getUIDsByLocator(ELoc::Z));

  // Run the calculator
  int error = (krige.run()) ? 0 : 1;
  return error;
}
