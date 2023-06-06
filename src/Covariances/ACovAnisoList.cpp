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
#include "Covariances/ACovAnisoList.hpp"

#include "Space/ASpace.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Covariances/CovLMGradient.hpp"

#include <math.h>
#include <vector>

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
    _covs.push_back(e->clone());
}

ACovAnisoList& ACovAnisoList::operator=(const ACovAnisoList &r)
{
  if (this != &r)
  {
    ACov::operator=(r);
    for (auto e: r._covs)
      _covs.push_back(e->clone());
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
  _covs.push_back(cov->clone());
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

  if (mode.getKeepOnlyCovIdx() >= 0)
    cov = _covs[mode.getKeepOnlyCovIdx()]->eval0(ivar, jvar, mode);
  else
  {
    for (int i=0, n=getCovNumber(); i<n; i++)
    {
      if (mode.getMember() != ECalcMember::LHS && isFiltered(i))
        continue;
      if (mode.getCovFiltered(i))
        continue;
      cov += _covs[i]->eval0(ivar, jvar, mode);
    }
  }
  return cov;
}

void ACovAnisoList::evalOptim(const SpacePoint &p1,
                              VectorDouble &res,
                              VectorDouble &temp,
                              SpacePoint &pttr,
                              int ivar,
                              int jvar,
                              const CovCalcMode &mode) const
{
  for (auto &e : res)
    e = 0;
  for (int i = 0, n = getCovNumber(); i < n; i++)
  {
    _covs[i]->evalOptim(p1, res, temp, pttr, ivar, jvar, mode);
  }
}

double ACovAnisoList::eval(const SpacePoint& p1,
                           const SpacePoint& p2,
                           int ivar,
                           int jvar,
                           const CovCalcMode& mode) const
{
  double cov = 0.;
  if (mode.getKeepOnlyCovIdx() >= 0)
    cov = _covs[mode.getKeepOnlyCovIdx()]->eval(p1, p2, ivar, jvar, mode);
  else
  {
    for (unsigned int i=0, n=getCovNumber(); i<n; i++)
    {
      if (mode.getMember() != ECalcMember::LHS && isFiltered(i))
        continue;
      if (mode.getCovFiltered(i))
        continue;
      cov += _covs[i]->eval(p1, p2, ivar, jvar, mode);
    }
  }
  return cov;
}

double ACovAnisoList::evalBasic(const SpacePoint &p1,
                                const SpacePoint &p2,
                                int ivar,
                                int jvar,
                                const CovCalcMode &mode) const
{
  double cov = 0.;
  for (unsigned int i=0, n=getCovNumber(); i<n; i++)
    cov += _covs[i]->eval(p1, p2, ivar, jvar, mode);
   return cov;
}

String ACovAnisoList::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (getCovNumber() <= 0) return sstr.str();

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

CovAniso ACovAnisoList::extractCova(int icov) const
{
  return *(_covs[icov]);
}

int ACovAnisoList::getMinOrder() const
{
  int nmini = -1;
  for (unsigned i = 0, n = getCovNumber(); i<n; i++)
  {
    int locmini = _covs[i]->getMinOrder();
    if (locmini > nmini) nmini = locmini;
  }
  return nmini;
}

const CovAniso* ACovAnisoList::getCova(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return nullptr;
  return _covs[icov];
}
CovAniso* ACovAnisoList::getCova(int icov)
{
  if (! _isCovarianceIndexValid(icov)) return nullptr;
  return _covs[icov];
}
void ACovAnisoList::setCova(int icov, CovAniso* covs)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covs[icov] = covs;
}
const ECov& ACovAnisoList::getType(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return ECov::UNKNOWN;
  return _covs[icov]->getType();
}
String ACovAnisoList::getCovName(int icov) const
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
 * Returns TEST as soon as one structure has no sill
 * @param ivar Rank of the first variable
 * @param jvar Rank of the second variable
 */
double ACovAnisoList::getTotalSill(int ivar, int jvar) const
{
  double sill_total = 0.;
  for (int icov = 0; icov < getCovNumber(); icov++)
  {
    const CovAniso* cova = getCova(icov);
    if (cova->getMinOrder() >= 0) return TEST;
    sill_total += cova->getSill(ivar, jvar);
  }
  return sill_total;
}

MatrixSquareGeneral ACovAnisoList::getTotalSill() const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      mat.setValue(ivar,jvar,getTotalSill(ivar,jvar));
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

void ACovAnisoList::copyCovContext(const CovContext& ctxt)
{
  int number = (int) _covs.size();
  for (int i = 0; i < number; i++)
    _covs[i]->copyCovContext(ctxt);
}

void ACovAnisoList::setAllFiltered(bool status)
{
  int number = (int) _covs.size();
  for (int i = 0; i < number; i++)
    _filtered[i] = status;
}

void ACovAnisoList::normalize(double sill, int ivar, int jvar)
{
  double covval = 0.;
  for (unsigned int i=0, n=getCovNumber(); i<n; i++)
    covval += _covs[i]->eval0(ivar, jvar);

  if (covval <= 0. || covval == sill) return;
  double ratio = sill / covval;

  for (unsigned int i=0, n=getCovNumber(); i<n; i++)
  {
    CovAniso* cov = _covs[i];
    cov->setSill(cov->getSill(ivar, jvar) * ratio);
  }
}

bool ACovAnisoList::hasNugget() const
{
  for (int is = 0; is < getCovNumber(); is++)
  {
    if (getType(is) == ECov::NUGGET) return true;
  }
  return false;
}

void ACovAnisoList::preProcess(const std::vector<SpacePoint>& vec) const
{
	for (int is = 0; is < getCovNumber(); is++)
		_covs[is]->preProcess(vec);
}

void ACovAnisoList::cleanPreProcessInfo() const
{
	for (int is = 0; is < getCovNumber(); is++)
		_covs[is]->cleanPreProcessInfo();
}

const ACovAnisoList* ACovAnisoList::reduce(const VectorInt &validVars) const
{
  ACovAnisoList* newcovlist = this->clone();

  for (int is = 0; is < getCovNumber(); is++)
  {
    CovAniso* covs = newcovlist->getCova(is);
    newcovlist->setCova(is,covs->reduce(validVars));
  }
  return newcovlist;
}
