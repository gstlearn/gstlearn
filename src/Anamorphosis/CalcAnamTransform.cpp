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
#include "geoslib_f_private.h"
#include "geoslib_old_f.h"

#include "Db/Db.hpp"
#include "Basic/VectorHelper.hpp"
#include "Anamorphosis/CalcAnamTransform.hpp"
#include "Anamorphosis/AnamContinuous.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"

#include <math.h>

CalcAnamTransform::CalcAnamTransform(AAnam* anam)
    : ACalcDbVarCreator(),
      _iattVar(-1),
      _iattFac(-1),
      _iattSel(-1),
      _flagVars(false),
      _flagToFactors(false),
      _flagDisjKrig(false),
      _flagCondExp(false),
      _flagUniCond(false),
      _flagZToY(true),
      _flagNormalScore(false),
      _ifacs(),
      _iptrEst(),
      _iptrStd(),
      _nbsimu(0),
      _flagOK(false),
      _proba(TEST),
      _anam(anam),
      _selectivity(nullptr)
{
}

CalcAnamTransform::~CalcAnamTransform()
{
}

bool CalcAnamTransform::_hasAnam(const EAnam& anamType) const
{
  if (_anam == nullptr)
  {
    messerr("The argument 'anam' must be defined");
    return false;
  }
  if (anamType == EAnam::UNKNOWN) return true;
  if (anamType != _anam->getType())
  {
    messerr("The argument 'anam'  should be of type");
    return false;
  }
  return true;
}

bool CalcAnamTransform::_hasInputVarDefined(int mode) const
{
  // Check that the vector is not empty
  if (_iptrEst.empty())
  {
    messerr("'db' should contain an Estimate variable");
    return false;
  }

  // Check that the UID are valid
  for (int i = 0; i < (int) _iptrEst.size(); i++)
    if (_iptrEst[i] < 0)
    {
      messerr("An estimation variable is not correctly defined");
      return false;
    }

  // Check that the vector is not empty
  if (_iptrStd.empty())
  {
    if (mode == 0)
      messerr("'db' should contain an St.Dev of Estimation Error variable");
    else
      messerr("'db' should contain an Variance of Estimation Error variable");
    return false;
  }

  // Check that the UID are valid
  for (int i = 0; i < (int) _iptrStd.size(); i++)
    if (_iptrStd[i] < 0)
    {
      if (mode == 0)
        messerr("A St. Dev. variable is not correctly defined");
      else
        messerr("A Variance variable is not correctly defined");
      return false;
    }

  return true;
}

bool CalcAnamTransform::_hasSelectivity() const
{
  int ncuts = _selectivity->getNCuts();
  if (ncuts <= 0 && ! _selectivity->isOnlyZDefined())
  {
    messerr("You must define some cutoff values");
    return false;
  }
  int nvarout = _selectivity->getVariableNumber();
  if (nvarout <= 0)
  {
    messerr("No recovery function is defined");
    return false;
  }
  return true;
}

bool CalcAnamTransform::_hasVariableNumber(bool equal1) const
{
  int number = getDb()->getLocatorNumber(ELoc::Z);
  if (! equal1)
  {
    if (number <= 0)
    {
      messerr("The argument 'db'  must have some variable(s) defined");
      return false;
    }
  }
  else
  {
    if (number != 1)
    {
      messerr("The argument 'db'  must have a single variable defined");
      return false;
    }
  }
  return true;
}

bool CalcAnamTransform::_check()
{
  if (!  hasDb()) return false;
  if (! _hasAnam()) return false;
  if (! _hasVariableNumber()) return false;

  if (_flagVars)
  {
    AnamContinuous* anamC = dynamic_cast<AnamContinuous*>(_anam);
    if (anamC == nullptr)
    {
      messerr("The argument 'anam'  must be of type AnamContinuous");
      return false;
    }
    return true;
  }

  if (_flagToFactors)
  {
    if (! _hasVariableNumber(true)) return false;
    int nmax = _anam->getNFactor();
    int nfact = _getNfact();
    for (int ifac = 0; ifac < nfact; ifac++)
      if (_ifacs[ifac] < 1 || _ifacs[ifac] > nmax)
      {
        messerr("Error in the rank of the factor(%d): it should lie in [1,%d]",
                _ifacs[ifac], nmax);
        return false;
      }
    return true;
  }

  if (_flagDisjKrig)
  {
    if (! _hasAnam(EAnam::HERMITIAN)) return false;
    if (! _hasInputVarDefined()) return false;
    if (! _hasSelectivity()) return false;
    return true;
  }

  if (_flagCondExp)
  {
    if (! _hasAnam(EAnam::HERMITIAN)) return false;
    if (! _hasInputVarDefined()) return false;
    if (! _hasSelectivity()) return false;
    return true;
  }

  if (_flagUniCond)
  {
    if (! _hasAnam(EAnam::HERMITIAN)) return false;
    if (! _hasInputVarDefined(1)) return false;
    if (! _hasSelectivity()) return false;
    if (_selectivity->isUsed(ESelectivity::Z))
    {
      messerr("The recovery option 'Z' is not available in this function");
      return false;
    }
    return true;
  }

  messerr("No Transformation option has been defined");
  return false;
}

bool CalcAnamTransform::_preprocess()
{
  if (_flagVars)
  {
    int nvar = _getNVar();
    _iattVar = getDb()->addColumnsByConstant(nvar);
    return (_iattVar >= 0);
  }

  if (_flagToFactors)
  {
    int nfact = _getNfact();
    _iattFac = getDb()->addColumnsByConstant(nfact);
    return (_iattFac >= 0);
  }

  if (_flagDisjKrig)
  {
    int nvarout = _getNSel();
    _iattSel = getDb()->addColumnsByConstant(nvarout, TEST);
    if (_iattSel < 0) return 1;
    return true;
  }

  if (_flagCondExp)
  {
    int nvarout = _getNSel();
    _iattSel = getDb()->addColumnsByConstant(nvarout, TEST);
    if (_iattSel < 0) return 1;
    return true;
  }

  if (_flagUniCond)
  {
    int nvarout = _getNSel();
    _iattSel = getDb()->addColumnsByConstant(nvarout, TEST);
    if (_iattSel < 0) return 1;
    return true;
  }

  return false;
}

bool CalcAnamTransform::_postprocess()
{
  /* Free the temporary variables */
  _cleanVariableDb(2);

  if (_flagVars)
  {
    int nvar = _getNVar();
    _renameVariable(nvar, _iattVar, ELoc::Z, String(), 1);
    return true;
  }

  if (_flagToFactors)
  {
    int nfact = _getNfact();
    _renameVariable(1, _iattFac, ELoc::Z, String(), nfact);
    return true;
  }

  if (_flagDisjKrig)
  {
    int nsel = _getNSel();
    for (int i = 0; i < nsel; i++)
      _renameVariable(1, _iattSel+i, ELoc::UNKNOWN, _selectivity->getVariableName(i), 1);
    return true;
  }

  if (_flagCondExp)
  {
    int nsel = _getNSel();
    for (int i = 0; i < nsel; i++)
      _renameVariable(1, _iattSel+i, ELoc::Z, _selectivity->getVariableName(i), 1);
    return true;
  }

  if (_flagUniCond)
  {
    int nsel = _getNSel();
    for (int i = 0; i < nsel; i++)
      _renameVariable(1, _iattSel+i, ELoc::Z, _selectivity->getVariableName(i), 1);
    return true;
  }

  return true;
}

void CalcAnamTransform::_rollback()
{
  _cleanVariableDb(1);
}

bool CalcAnamTransform::_run()
{
  if (_flagVars)
  {
    if (_flagZToY)
    {

      // Transform from Raw to Gaussian
      if (_flagNormalScore)
      {
        // Transform using Normal Score
        if (_ZToYByNormalScore()) return true;
      }
      else
      {

        // Transform using Hermite expansion
        if (_ZToYByHermite()) return true;
      }
    }
    else
    {

      // Transform from Gaussian to Raw
      if (_YToZByHermite()) return true;
    }
    return false;
  }

  if (_flagToFactors)
  {
    return _ZToFactors();
  }

  if (_flagDisjKrig)
  {
    return _FactorsToSelectivity();
  }

  if (_flagCondExp)
  {
    return (!_conditionalExpectation(getDb(), _anam, _selectivity, _iattSel,
                                 _iptrEst[0], _iptrStd[0], _flagOK, _proba,
                                 _nbsimu));
  }

  if (_flagUniCond)
  {
    AnamHermite* anam_hermite = dynamic_cast<AnamHermite*>(_anam);
    return (!_uniformConditioning(getDb(), anam_hermite, _selectivity, _iattSel,
                              _iptrEst[0], _iptrStd[0]));
  }

  return false;
}

bool CalcAnamTransform::_ZToYByNormalScore()
{
  int nech = getDb()->getSampleNumber();
  int nvar = _getNVar();

  // Reading the optional weight

  VectorDouble wt;
  if (getDb()->hasLocVariable(ELoc::W))
    wt = getDb()->getColumnByLocator(ELoc::W);
  else
    wt.resize(nech, 1.);

  // Loop on the variables

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble data = getDb()->getColumnByLocator(ELoc::Z, ivar);
    VectorDouble vec = VH::normalScore(data, wt);
    if (! vec.empty())
      getDb()->setColumnByUID(vec, _iattVar + ivar);
  }
  return true;
}

bool CalcAnamTransform::_ZToYByHermite()
{
  int nvar = _getNVar();
  const AnamContinuous* anamC = dynamic_cast<const AnamContinuous*>(getAnam());

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble z = getDb()->getColumnByLocator(ELoc::Z, ivar, true);
    if (z.size() <= 0) continue;
    VectorDouble y = anamC->rawToGaussianVector(z);
    getDb()->setColumnByUID(y, _iattVar + ivar, true);
  }
  return true;
}

bool CalcAnamTransform::_YToZByHermite()
{
  int nvar = _getNVar();
  const AnamContinuous* anamC = dynamic_cast<const AnamContinuous*>(getAnam());

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble y = getDb()->getColumnByLocator(ELoc::Z, ivar, true);
    if (y.size() <= 0) continue;
    VectorDouble z = anamC->gaussianToRawVector(y);
    getDb()->setColumnByUID(z, _iattVar + ivar, true);
  }
  return true;
}

bool CalcAnamTransform::_ZToFactors()
{
  int nfact = _getNfact();

  for (int iech = 0; iech < getDb()->getSampleNumber(); iech++)
  {
    if (! getDb()->isActive(iech)) continue;
    double zval = getDb()->getZVariable(iech, 0);
    if (FFFF(zval)) continue;
    VectorDouble factors = _anam->z2factor(zval, _ifacs);
    if (factors.empty()) continue;
    for (int ifac = 0; ifac < nfact; ifac++)
      getDb()->setArray(iech, _iattFac + ifac, factors[ifac]);
  }
  return true;
}

bool CalcAnamTransform::_FactorsToSelectivity()
{
  AnamHermite    *anam_hermite     = dynamic_cast<AnamHermite*>(_anam);
  AnamDiscreteDD *anam_discrete_DD = dynamic_cast<AnamDiscreteDD*>(_anam);
  AnamDiscreteIR *anam_discrete_IR = dynamic_cast<AnamDiscreteIR*>(_anam);

  /* Dispatch according to the type of Anamorphosis */

  switch (_anam->getType().toEnum())
  {
    case EAnam::E_HERMITIAN:
      anam_hermite->factor2Selectivity(getDb(), _selectivity, _iptrEst,
                                       _iptrStd, _iattSel);
      return true;

    case EAnam::E_DISCRETE_DD:
      anam_discrete_DD->factor2Selectivity(getDb(), _selectivity, _iptrEst,
                                           _iptrStd, _iattSel);
      return true;

    case EAnam::E_DISCRETE_IR:
      anam_discrete_IR->factor2Selectivity(getDb(), _selectivity, _iptrEst,
                                           _iptrStd, _iattSel);
      return true;

    default:
      messerr("This method is not programmed yet for this anamorphosis");
      return false;
  }
}


/*****************************************************************************/
/*!
 **  Calculate the recoveries (z,T,Q,m,B) starting from the factors
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Point anamorphosis
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  name_est     Array of variable names for factor estimation
 ** \param[in]  name_std     Array of variable names for factor St. Dev.
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
int DisjunctiveKriging(Db *db,
                        AAnam *anam,
                        Selectivity *selectivity,
                        const VectorString &name_est,
                        const VectorString &name_std,
                        const NamingConvention &namconv)
{
  CalcAnamTransform transfo(anam);
  transfo.setDb(db);
  transfo.setSelectivity(selectivity);
  transfo.setIptrEst(db->getUIDs(name_est));
  transfo.setIptrStd(db->getUIDs(name_std));
  transfo.setFlagDisjKrig(true);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
}

/*****************************************************************************/
/*!
 **  Calculate the Conditional Expectation
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Point anamorphosis
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  name_est     Name of the Kriging estimate
 ** \param[in]  name_std     Name of the Kriging St. deviation
 ** \param[in]  flag_OK      1 if kriging has ben performed in Ordinary Kriging
 ** \param[in]  proba        Probability
 ** \param[in]  nbsimu       Number of Simulation outcomes
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
int ConditionalExpectation(Db *db,
                           AAnam *anam,
                           Selectivity *selectivity,
                           const String &name_est,
                           const String &name_std,
                           bool flag_OK,
                           double proba,
                           int nbsimu,
                           const NamingConvention &namconv)
{
  CalcAnamTransform transfo(anam);
  transfo.setDb(db);
  transfo.setSelectivity(selectivity);
  transfo.setIptrEst({db->getUID(name_est)});
  transfo.setIptrStd({db->getUID(name_std)});
  transfo.setFlagCondExp(true);
  transfo.setFlagOk(flag_OK);
  transfo.setNbsimu(nbsimu);
  transfo.setProba(proba);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
}

/*****************************************************************************/
/*!
 **  Calculate the Uniform Conditioning
 **
 ** \return  Error return code
 **
 ** \param[in]  db           Db structure containing the factors (Z-locators)
 ** \param[in]  anam         Point anamorphosis
 ** \param[in]  selectivity  Selectivity structure
 ** \param[in]  name_est     Name of the Kriging estimate
 ** \param[in]  name_varz    Name of the Variance of Kriging estimate
 ** \param[in]  namconv      Naming Convention
 **
 ** \remark We need the variance of Estimation Error... even if it will be
 ** \remark temporarily stored in a member names iptrStd.
 **
 *****************************************************************************/
int UniformConditioning(Db *db,
                        AAnam *anam,
                        Selectivity *selectivity,
                        const String &name_est,
                        const String &name_varz,
                        const NamingConvention &namconv)
{
  CalcAnamTransform transfo(anam);
  transfo.setDb(db);
  transfo.setSelectivity(selectivity);
  transfo.setIptrEst({db->getUID(name_est)});
  transfo.setIptrStd({db->getUID(name_varz)});
  transfo.setFlagUniCond(true);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
}
