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
#include "Estimation/KrigingSystem.hpp"
#include "Estimation/CalcKrigingFactors.hpp"
#include "Anamorphosis/AAnam.hpp"
#include "Model/Model.hpp"

CalcKrigingFactors::CalcKrigingFactors(bool flag_est, bool flag_std)
  : ACalcInterpolator()
  , _flagEst(flag_est)
  , _flagStd(flag_std)
  , _nameCoord()
  , _iptrEst(-1)
  , _iptrStd(-1)
  , _iuidFactors()
  , _modelLocal(nullptr)
{
}

CalcKrigingFactors::~CalcKrigingFactors()
{
}

bool CalcKrigingFactors::_check()
{
  // Turn the problem to Monovariate before checking consistency with 'Model'
  getDbin()->clearLocators(ELoc::Z);
  getDbin()->setLocatorByUID(_iuidFactors[0], ELoc::Z, 0);

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

  _modelLocal = dynamic_cast<Model*>(getModel());
  if (_modelLocal == nullptr)
  {
    messerr("The model must be of type Model (not ModelGeneric)");
    return false;
  }
  if (_modelLocal->getNVar() != 1)
  {
    messerr("This application is limited to the monovariate Model case");
    return false;
  }
  if (! _modelLocal->hasAnam())
  {
    messerr("Argument 'model' should has an Anamorphosis attached");
    return false;
  }
  // If change of support is defined through the anamorphosis,
  // the calculation option (EKrigOpt) should be set to POINT
  // in order to avoid additional block randomization
  if (getKrigopt().getCalcul() == EKrigOpt::BLOCK &&
      ! getKrigopt().hasDiscs())
  {
    messerr("For Block estimate, you must specify the discretization");
    return false;
  }
  return true;
}

bool CalcKrigingFactors::_hasChangeSupport() const
{
  const AAnam* anam = _modelLocal->getAnam();
  if (anam == nullptr) return false;

  // Check if the change of support is defined in the Anamorphosis
  return anam->isChangeSupportDefined();
}

bool CalcKrigingFactors::_preprocess()
{
  if (!ACalcInterpolator::_preprocess()) return false;
  
  // Centering the information (only when a change of support is defined)
  if (_hasChangeSupport())
  {
    DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
    if (dbgrid == nullptr)
    {
      messerr("Due to change of support, 'dbout' should be a Grid");
      return false;
    }
    if (! getKrigopt().hasDiscs())
    {
      // Center the information in the blocks of the output grid
      // Duplicating the coordinate variable names before centering
      _nameCoord = getDbin()->getNamesByLocator(ELoc::X);
      int error = _centerDataToGrid(dbgrid);
      if (error) return false;
    }
    if (getKrigopt().hasDiscs())
    {
      // Center the information in sub-blocks when the output grid defines panels
      VectorInt ndiscs = getKrigopt().getDiscs();
      DbGrid* dbsmu = DbGrid::createDivider(dbgrid, ndiscs, 1);
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

  getDbin()->setLocatorsByUID(_iuidFactors, ELoc::Z, 0);

  int nfactor = _getNFactors();
  _renameVariable(2, VectorString(), ELoc::Z, nfactor, _iptrStd, "stdev", 1);
  _renameVariable(2, VectorString(), ELoc::Z, nfactor, _iptrEst, "estim", 1);

  // Centering the information (only when a change of support is defined)
  if (_hasChangeSupport() && ! _nameCoord.empty())
    getDbin()->setLocators(_nameCoord, ELoc::X, 0);

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
  KrigingSystem ksys(getDbin(), getDbout(), getModel(), getNeigh(), getKrigopt());
  if (ksys.updKrigOptEstim(_iptrEst, _iptrStd, -1)) return 1;
  if (ksys.setKrigOptFactorKriging(true)) return 1;
  if (! ksys.isReady()) return 1;

  // Loop on the targets to be processed

  int ntotal = getDbout()->getNSample() * _getNFactors();
  int iproc = 0;
  for (int iclass = 1; iclass <= _getNFactors(); iclass++)
  {
    int jptr_est = (_flagEst) ? _iptrEst + iclass - 1 : -1;
    int jptr_std = (_flagStd) ? _iptrStd + iclass - 1 : -1;
    getDbin()->clearLocators(ELoc::Z);
    getDbin()->setLocatorByUID(_iuidFactors[iclass - 1], ELoc::Z, 0);
    if (ksys.updKrigOptEstim(jptr_est, jptr_std, -1)) return 1;
    if (ksys.updKrigOptIclass(iclass, _getNFactors())) return 1;

    for (int iech_out = 0; iech_out < getDbout()->getNSample(); iech_out++)
    {
      mes_process("Disjunctive Kriging for cell", ntotal, iproc);
      if (ksys.estimate(iech_out)) return 1;
    }
  }
  ksys.conclusion();

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
int krigingFactors(Db *dbin,
                   Db *dbout,
                   Model *model,
                   ANeigh *neigh,
                   bool flag_est,
                   bool flag_std,
                   const KrigOpt& krigopt,
                   const NamingConvention &namconv)
{
  CalcKrigingFactors krige(flag_est, flag_std);
  krige.setDbin(dbin);
  krige.setDbout(dbout);
  krige.setModel(model);
  krige.setNeigh(neigh);
  krige.setKrigopt(krigopt);
  krige.setNamingConvention(namconv);

  krige.setIuidFactors(dbin->getUIDsByLocator(ELoc::Z));

  // Run the calculator
  int error = (krige.run()) ? 0 : 1;
  return error;
}
