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
#include "Basic/Law.hpp"
#include "Basic/Table.hpp"
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
      _flagFromFactors(false),
      _flagZToY(true),
      _flagNormalScore(false),
      _ifacs(),
      _iptrEst(),
      _iptrStd(),
      _anam(anam),
      _selectivity(nullptr)
{
}

CalcAnamTransform::~CalcAnamTransform()
{
}

bool CalcAnamTransform::_check()
{
  if (! hasDb()) return false;

  if (_anam == nullptr)
  {
    messerr("The argument 'anam'  must be defined");
    return false;
  }
  int number = getDb()->getLocatorNumber(ELoc::Z);
  if (number <= 0)
  {
    messerr("The argument 'db'  must have some variable(s) defined");
    return false;
  }

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
    if (number != 1)
    {
      messerr("The argument 'db' must contain a single variable");
      return false;
    }
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

  if (_flagFromFactors)
  {
    if (_iptrEst.empty() && _iptrStd.empty())
    {
      messerr("No factor is defined");
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

  messerr("No Transformation option has been defined");
  return false;
}

bool CalcAnamTransform::_preprocess()
{
  if (_flagVars)
  {
    int nvar = _getNVar();
    _iattVar = getDb()->addColumnsByConstant(nvar);
    if (_iattVar < 0) return false;
    return true;
  }

  if (_flagToFactors)
  {
    int nfact = _getNfact();
    _iattFac = getDb()->addColumnsByConstant(nfact);
    if (_iattFac < 0) return false;
    return true;
  }

  if (_flagFromFactors)
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
  if (_flagVars)
  {
    int nvar = _getNVar();
    _renameVariable(ELoc::Z, nvar, _iattVar, String(), 1);
    return true;
  }

  if (_flagToFactors)
  {
    int nfact = _getNfact();
    _renameVariable(ELoc::Z, 1, _iattFac, String(), nfact);
    return true;
  }

  if (_flagFromFactors)
  {
    int nsel = _getNSel();
    _renameVariable(ELoc::Z, 1, _iattSel, String(), nsel);
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
    if (_ZToFactors()) return true;
    return false;
  }

  if (_flagFromFactors)
  {
    if (_FactorsToSelectivity()) return true;
    return false;
  }
  return false;
}

bool CalcAnamTransform::_ZToYByNormalScore()
{
  int nech = getDb()->getSampleNumber();
  int nvar = _getNVar();

  // Reading the optional weight

  VectorDouble wt;
  if (getDb()->hasWeight())
    wt = getDb()->getColumnByLocator(ELoc::W);
  else
    wt.resize(nech, 1.);

  // Loop on the variables

  for (int ivar = 0; ivar < nvar; ivar++)
  {
    VectorDouble data = getDb()->getColumnByLocator(ELoc::Z, ivar);

    // Check that all weights are positive
    double wtotal = 0.;
    for (int iech = 0; iech < nech; iech++)
    {
      bool defined = false;
      if (getDb()->isActive(iech) && !FFFF(data[iech]))
      {
        defined = true;
        if (wt[iech] < 0.)
        {
          messerr("The weight of sample (%d) is negative (%lf)", iech + 1,
                  wt[iech]);
          return false;
        }
        wtotal += wt[iech];
      }
      if (!defined) data[iech] = TEST;
    }
    if (wtotal <= 0.)
    {
      messerr("The sum of the weights is not positive");
      return false;
    }

    // Get the list of indices sorted by increasing values of data
    VectorInt idx = ut_vector_sort_indices(data);

    // Loop on the samples
    double wlocal = 0.;
    for (int iech = 0; iech < nech; iech++)
    {
      int jech = idx[iech];
      wlocal += wt[jech];
      double z = wlocal / ((double) nech + 1.);
      double gaus = law_invcdf_gaussian(z);
      getDb()->setArray(iech, _iattVar + ivar, gaus);
    }
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
    VectorDouble y = anamC->RawToGaussianVector(z);
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
    VectorDouble y = getDb()->getColumnByLocator(ELoc::Z, ivar);
    if (y.size() <= 0) continue;
    VectorDouble z = anamC->GaussianToRawVector(y);
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
    double zval = getDb()->getVariable(iech, 0);
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
  AnamHermite *anam_hermite = dynamic_cast<AnamHermite*>(_anam);
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
  return false;
}

int RawToGaussian(Db *db,
                  AAnam* anam,
                  const ELoc &locatorType,
                  const NamingConvention &namconv)
{
  CalcAnamTransform transfo(anam);
  transfo.setFlagVars(true);
  transfo.setFlagZToY(true);
  transfo.setFlagNormalScore(false);
  transfo.setDb(db);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
}

int RawToGaussian(Db *db,
                  AAnam* anam,
                  const String &name,
                  const NamingConvention &namconv)
{
  if (db == nullptr) return 1;
  db->setLocator(name, ELoc::Z);

  CalcAnamTransform transfo(anam);
  transfo.setFlagVars(true);
  transfo.setFlagZToY(true);
  transfo.setFlagNormalScore(false);
  transfo.setDb(db);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
}

int GaussianToRaw(Db *db,
                  AAnam *anam,
                  const ELoc &locatorType,
                  const NamingConvention &namconv)
{
  CalcAnamTransform transfo(anam);
  transfo.setFlagZToY(true);
  transfo.setFlagNormalScore(false);
  transfo.setDb(db);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
}

int GaussianToRaw(Db *db,
                  AAnam *anam,
                  const String& name,
                  const NamingConvention &namconv)
{
  if (db == nullptr) return 1;
  db->setLocator(name, ELoc::Z);

  CalcAnamTransform transfo(anam);
  transfo.setFlagVars(true);
  transfo.setFlagZToY(true);
  transfo.setFlagNormalScore(false);
  transfo.setDb(db);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
}

/****************************************************************************/
/*!
 **  Transform the target variable inti Gaussian by Normal Score
 **
 ** \return  Error return code
 **
 ** \param[in]  db         Db Structure
 ** \param[in]  namconv    Naming convention
 **
 *****************************************************************************/
int NormalScore(Db *db, const NamingConvention &namconv)
{
  AnamHermite anam = AnamHermite(1);
  CalcAnamTransform transfo(&anam);
  transfo.setDb(db);
  transfo.setFlagVars(true);
  transfo.setFlagZToY(true);
  transfo.setFlagNormalScore(true);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
}

/*****************************************************************************/
/*!
 **  Calculate the factors corresponding to an input data vector
 **
 ** \return  Error return code
 **
 ** \param[in]  db          Db structure
 ** \param[in]  anam        Anamorphosis structure
 ** \param[in]  ifacs       Array of factor ranks
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int RawToFactor(Db *db,
                AAnam *anam,
                const VectorInt &ifacs,
                const NamingConvention &namconv)
{
  CalcAnamTransform transfo(anam);
  transfo.setDb(db);
  transfo.setFlagToFactors(true);
  transfo.setIfacs(ifacs);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
}

/*****************************************************************************/
/*!
 **  Calculate the factors corresponding to an input data vector
 **
 ** \return  Error return code
 **
 ** \param[in]  db          Db structure
 ** \param[in]  anam        Anamorphosis structure
 ** \param[in]  nfactor     Number of first factors
 ** \param[in]  namconv     Naming convention
 **
 *****************************************************************************/
int RawToFactor(Db *db,
                AAnam *anam,
                int nfactor,
                const NamingConvention &namconv)
{
  CalcAnamTransform transfo(anam);
  transfo.setDb(db);
  transfo.setFlagToFactors(true);
  VectorInt ifacs = ut_ivector_sequence(nfactor, 1);
  transfo.setIfacs(ifacs);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
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
 ** \param[in]  names_est    Array of variable names for factor estimation
 ** \param[in]  names_std    Array of variable names for factor St. Dev.
 ** \param[in]  namconv      Naming convention
 **
 *****************************************************************************/
int FactorToSelectivity(Db *db,
                        AAnam *anam,
                        Selectivity *selectivity,
                        const VectorString &names_est,
                        const VectorString &names_std,
                        const NamingConvention &namconv)
{
  CalcAnamTransform transfo(anam);
  transfo.setDb(db);
  transfo.setSelectivity(selectivity);
  transfo.setIptrEst(db->getUIDs(names_est));
  transfo.setIptrStd(db->getUIDs(names_std));
  transfo.setFlagFromFactors(true);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
}

