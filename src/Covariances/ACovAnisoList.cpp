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
#include "Covariances/ACovAnisoList.hpp"

#include "Space/ASpace.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "geoslib_f.h"
#include <math.h>

ACovAnisoList::ACovAnisoList(const ASpace* space)
: ACov(space),
  _covs(),
  _filtered()
{
}

ACovAnisoList::ACovAnisoList(const ACovAnisoList &r)
: ACov(r),
  _covs(),
  _filtered(r._filtered)
{
  for (auto e: r._covs)
  {
    _covs.push_back(dynamic_cast<CovAniso*>(e->clone()));
  }
}

ACovAnisoList& ACovAnisoList::operator=(const ACovAnisoList &r)
{
  if (this != &r)
  {
    ACov::operator=(r);
    for (auto e: r._covs)
    {
      _covs.push_back(dynamic_cast<CovAniso*>(e->clone()));
    }
    _filtered = r._filtered;
  }
  return *this;
}

ACovAnisoList::~ACovAnisoList()
{
  delAllCov();
}

void ACovAnisoList::addCovList(const ACovAnisoList* covs)
{
  int ncov = covs->getCovNumber();
  for (int icov = 0; icov < ncov; icov++)
    addCov(covs->getCova(icov));
}

void ACovAnisoList::addCov(const CovAniso* cov)
{
  if (getCovNumber() > 0)
  {
    // A covariance has already been considered.
    // Check that the current Context is similar to the one of the newly
    // added covariance

    if (! cov->getContext().isEqual(_covs[0]->getContext()))
    {
      messerr("Error: Covariances in the same ACovAnisoList should share the same Context");
      messerr("Operation is cancelled");
      return;
    }
  }
  _covs.push_back(dynamic_cast<CovAniso*>(cov->clone()));
  _filtered.push_back(false);
}

void ACovAnisoList::delCov(unsigned int i)
{
  if (! _isCovarianceIndexValid(i)) return;
  delete _covs[i];
  _covs.erase(_covs.begin() + i);
  _filtered.erase(_filtered.begin() + i);
}

void ACovAnisoList::delAllCov()
{
  for (auto e: _covs)
  {
    delete e;
  }
  _covs.clear();
  _filtered.clear();
}

void ACovAnisoList::setFiltered(unsigned int i, bool filtered)
{
  if (! _isCovarianceIndexValid(i)) return;
  _filtered[i] = filtered;
}

bool ACovAnisoList::isConsistent(const ASpace* /*space*/) const
{
  /// TODO : ACovAnisoList::isConsistent
  return true;
}

int ACovAnisoList::getNVariables() const
{
  if (getCovNumber() > 0)
    return _covs[0]->getNVariables();
  return 0;
}

double ACovAnisoList::eval0(int ivar, int jvar, const CovCalcMode& mode) const
{
  double cov = 0.;
  for (unsigned int i=0, n=getCovNumber(); i<n; i++)
  {
    if (mode.getMember() != ECalcMember::LHS && isFiltered(i))
      continue;
    if (mode.isFilterNugget() && getType(i) == ECov::NUGGET)
      continue;
    if (mode.getKeepOnlyCovIdx() == i)
      return _covs[i]->eval0(ivar, jvar, mode);
    cov += _covs[i]->eval0(ivar, jvar, mode);
  }

  // Normalization
  if (mode.getNormalized())
  {
    cov /= _getNormalizationFactor(ivar, jvar, mode);
  }

  return cov;
}

double ACovAnisoList::eval(int ivar,
                           int jvar,
                           const SpacePoint& p1,
                           const SpacePoint& p2,
                           const CovCalcMode& mode) const
{
  double cov = 0.;
  for (unsigned int i=0, n=getCovNumber(); i<n; i++)
  {
    if (mode.getMember() != ECalcMember::LHS && isFiltered(i))
      continue;
    if (mode.isFilterNugget() && getType(i) == ECov::NUGGET)
      continue;
    if (mode.getKeepOnlyCovIdx() == i)
      return _covs[i]->eval(ivar, jvar, p1, p2, mode);
    cov += _covs[i]->eval(ivar, jvar, p1, p2, mode);
  }

  // Normalization
  if (mode.getNormalized())
  {
    cov /= _getNormalizationFactor(ivar, jvar, mode);
  }

  return cov;
}

String ACovAnisoList::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  for (int icov = 0; icov < getCovNumber(); icov++)
  {
    sstr << getCova(icov)->toString();
    if (isFiltered(icov))
      sstr << "  (This component is Filtered)" << std::endl;
  }

  // Display the Total Sill (optional)
  if (isStationary())
  {
    if (getNVariables() <= 1)
    {
      sstr << "Total Sill     = " << toDouble(getTotalSill(0,0));
    }
    else
    {
      sstr << toMatrix("Total Sill",VectorString(),VectorString(),0,
                       getNVariables(),getNVariables(),
                       getTotalSill().getValues());
    }
  }
  sstr << std::endl;
  return sstr.str();
}

int ACovAnisoList::getCovNumber() const
{
  return static_cast<int> (_covs.size());
}

bool ACovAnisoList::isFiltered(unsigned int i) const
{
  if (! _isCovarianceIndexValid(i)) return false;
  return _filtered[i];
}

bool ACovAnisoList::hasRange() const
{
  for (unsigned int i=0, n=getCovNumber(); i<n; i++)
  {
    if (!getCova(i)->hasRange())
      return false;
  }
  return true;
}

bool ACovAnisoList::isStationary() const
{
  for (unsigned int i=0, n=getCovNumber(); i<n; i++)
  {
    if (getCova(i)->getMinOrder() >= 0)
      return false;
  }
  return true;
}

const CovAniso* ACovAnisoList::getCova(unsigned int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return nullptr;
  return _covs[icov];
}
CovAniso* ACovAnisoList::getCova(unsigned int icov)
{
  if (! _isCovarianceIndexValid(icov)) return nullptr;
  return _covs[icov];
}
const ECov& ACovAnisoList::getType(unsigned int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return ECov::UNKNOWN;
  return _covs[icov]->getType();
}
String ACovAnisoList::getCovName(unsigned int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return String();
  return _covs[icov]->getCovName();
}
double ACovAnisoList::getParam(unsigned int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return 0.;
  return _covs[icov]->getParam();
}
const MatrixSquareSymmetric& ACovAnisoList::getSill(unsigned int icov) const
{
  return _covs[icov]->getSill();
}
double ACovAnisoList::getSill(unsigned int icov, int ivar, int jvar) const
{
  if(! _isCovarianceIndexValid(icov)) return 0.;
  return _covs[icov]->getSill(ivar, jvar);
}
void ACovAnisoList::setSill(unsigned int icov, int ivar, int jvar, double value)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covs[icov]->setSill(ivar, jvar, value);
}

void ACovAnisoList::setType(unsigned int icov, const ECov& type)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covs[icov]->setType(type);
}

int ACovAnisoList::getGradParamNumber(unsigned int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return 0;
  return _covs[icov]->getGradParamNumber();
}

/**
 * Calculate the total sill of the model for given pair of variables
 * Returns 0 as soon as one structure has no range
 * @param ivar Rank of the first variable
 * @param jvar Rank of the second variable
 */
double ACovAnisoList::getTotalSill(int ivar, int jvar) const
{
  return eval0(ivar, jvar);
}

MatrixSquareGeneral ACovAnisoList::getTotalSill() const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      mat.setValue(ivar,jvar,eval0(ivar,jvar));
  return mat;
}

bool ACovAnisoList::_isCovarianceIndexValid(unsigned int i) const
{
  if (i >= (unsigned int) getCovNumber())
  {
    mesArg("Covariance Index",i,getCovNumber());
    return false;
  }
  return true;
}

/**
 * Returns the largest range (in any direction in the heterotopic case)
 * @return
 */
double ACovAnisoList::getMaximumDistance() const

{
  double maxdist = 0.;
  for (int icov = 0; icov < getCovNumber(); icov++)
  {
    const CovAniso* cova = getCova(icov);
    if (! cova->hasRange()) continue;
    double range = cova->getRange();
    if (range > maxdist) maxdist = range;
  }
  return maxdist;
}

double ACovAnisoList::_getNormalizationFactor(int ivar,
                                       int jvar,
                                       const CovCalcMode& mode) const
{
  double c00 = eval0(ivar, ivar, mode);
  if (c00 < 0 || FFFF(c00))
  my_throw("Normalization required but C00 not defined");
  if (ivar != jvar)
  {
    double c11 = eval0(jvar, jvar, mode);
    if (c11 < 0 || FFFF(c11))
    my_throw("Normalization required but C00 is not defined");
    c00 = sqrt(c00 * c11);
  }
  return c00;
}
