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
#include "Anamorphosis/CalcAnamTransform.hpp"
#include "Anamorphosis/AnamContinuous.hpp"
#include "Anamorphosis/AnamHermite.hpp"

#include <math.h>

CalcAnamTransform::CalcAnamTransform(AAnam* anam)
    : ACalcDbVarCreator(),
      _iatt(-1),
      _flagZToY(true),
      _flagNormalScore(false),
      _ifacs(),
      _anam(anam)
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

  if (! _toFactors())
  {
    AnamContinuous* anamC = dynamic_cast<AnamContinuous*>(_anam);
    if (anamC == nullptr)
    {
      messerr("The argument 'anam'  must be of type AnamContinuous");
      return false;
    }
  }
  else
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
  }
  return true;
}

bool CalcAnamTransform::_preprocess()
{
  if (! _toFactors())
  {
    int nvar = _getNVar();
    _iatt = getDb()->addColumnsByConstant(nvar);
    if (_iatt < 0) return false;
  }
  else
  {
    int nfact = _getNfact();
    _iatt = getDb()->addColumnsByConstant(nfact);
    if (_iatt < 0) return false;
  }
  return true;
}

bool CalcAnamTransform::_postprocess()
{
  if (! _toFactors())
  {
    int nvar = _getNVar();
  _renameVariable(ELoc::Z, nvar, _iatt, String(), 1);
  }
  else
  {
    int nfact = _getNfact();
  _renameVariable(ELoc::Z, 1, _iatt, String(), nfact);
  }

  return true;
}

void CalcAnamTransform::_rollback()
{
  _cleanVariableDb(1);
}

bool CalcAnamTransform::_run()
{
  if (!_toFactors())
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
  }
  else
  {
    if (_ZToFactors()) return true;
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
      getDb()->setArray(iech, _iatt + ivar, gaus);
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
    getDb()->setColumnByUID(y, _iatt + ivar, true);
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
    getDb()->setColumnByUID(z, _iatt + ivar, true);
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
      getDb()->setArray(iech, _iatt + ifac, factors[ifac]);
  }
  return true;
}

int RawToGaussian(Db *db,
                  AAnam* anam,
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

int RawToGaussian(Db *db,
                  AAnam* anam,
                  const String &name,
                  const NamingConvention &namconv)
{
  if (db == nullptr) return 1;
  db->setLocator(name, ELoc::Z);

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
  VectorInt ifacs = ut_ivector_sequence(nfactor, 1);
  transfo.setIfacs(ifacs);
  transfo.setNamingConvention(namconv);

  // Run the calculator
  int error = (transfo.run()) ? 0 : 1;
  return error;
}

