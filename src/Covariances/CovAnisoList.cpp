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
#include "Covariances/CovAnisoList.hpp"
#include "Enum/EModelProperty.hpp"

#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/CovList.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Space/ASpace.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "Covariances/CovLMCConvolution.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Db/Db.hpp"
#include "Space/SpacePoint.hpp"
#include "Anamorphosis/AnamHermite.hpp"

#include <math.h>
#include <vector>


CovAnisoList::CovAnisoList(const CovContext& ctxt)
: CovList(ctxt),
  _covAnisos()
{
}

CovAnisoList::CovAnisoList(const CovAnisoList &r)
: CovList(r.getContext()),
  _covAnisos()
{
  for (auto* e: r._covAnisos)
  {
    _pushCov(e->clone());
  }
  _filtered = r._filtered;
  _updateLists();
}

CovAnisoList& CovAnisoList::operator=(const CovAnisoList &r)
{
  if (this != &r)
  {
    _ctxt = r._ctxt;
    for (auto *e: r._covAnisos)
    {
     _pushCov(e->clone());
    }
    _filtered = r._filtered;
    _updateLists();
  }
  return *this;
}

CovAnisoList::~CovAnisoList()
{
  delAllCov();
}

void CovAnisoList::addCovList(const CovAnisoList* covs)
{
  for (int icov = 0, ncov = covs->getNCov(); icov < ncov; icov++)
    addCov(covs->getCova(icov));
}

void CovAnisoList::addCov(const CovBase* cov)
{
  const CovAniso* covaniso = dynamic_cast<const CovAniso*>(cov);
  if (covaniso != nullptr)
    addCovAniso(covaniso);
  else
    messerr("Error: CovAnisoList::addCov: Covariance is not of type CovAniso");
}

void CovAnisoList::addCovAniso(const CovAniso* cov)
{
  if (getNCov() == 0)
  {
    setNVar(cov->getNVar());
  }
  else
  {
    // A covariance has already been considered.
    // Check that the current Context is similar to the one of the newly
    // added covariance

    if (! cov->getContext().isEqual(_covs[0]->getContext()))
    {
      messerr("Error: Covariances in the same CovAnisoList should share the same Context");
      messerr("Operation is cancelled");
      return;
    }
  }
  _pushCov(cov);
  _filtered.push_back(false);
  _updateLists();
}

void CovAnisoList::_pushCov(const CovAniso* cov)
{
   CovAniso* covcopy = cov->clone();
  _covs.push_back(covcopy);
  _covAnisos.push_back(covcopy);
}

void CovAnisoList::_delCov(int icov)
{
   _covAnisos.erase(_covAnisos.begin() + icov);
}

void CovAnisoList::_delAllCov()
{
  _covAnisos.clear();
}

bool CovAnisoList::isNoStat() const
{
  bool nostat = false;
  for (const auto &e :_covs)
  {
    nostat = nostat || e->isNoStat();
  }
  return nostat;
}
void CovAnisoList::setCovFiltered(int icov, bool filtered)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _filtered[icov] = filtered;
  _updateLists();
}

bool CovAnisoList::isConsistent(const ASpace* /*space*/) const
{
  /// TODO : CovAnisoList::isConsistent
  return true;
}

int CovAnisoList::getNVar() const
{
  if (getNCov() > 0)
    return _covs[0]->getNVar();
  return 0;
}


double CovAnisoList::eval0(int ivar, int jvar, const CovCalcMode* mode) const
{
  double cov      = 0.;
  const VectorInt& list = _getListActiveCovariances(mode);
  for (int i = 0, n = (int)list.size(); i < n; i++)
  {
    int j = list[i];
    cov += _covs[j]->eval0(ivar, jvar, mode);
  }
  return cov;
}

/**
 * Calculate the Matrix of covariance for zero distance
 * @param mat   Covariance matrix (Dimension: nvar * nvar)
 * @param mode  Calculation Options
 *
 * @remarks: Matrix 'mat' should be dimensioned and initialized beforehand
 */
void CovAnisoList::addEval0CovMatBiPointInPlace(MatrixSquareGeneral& mat,
                                                const CovCalcMode* mode) const
{
  const VectorInt& list = _getListActiveCovariances(mode);
  for (int i = 0, n = (int)list.size(); i < n; i++)
  {
    int j = list[i];
    _covs[j]->addEval0CovMatBiPointInPlace(mat, mode);
  }
}

void CovAnisoList::optimizationSetTargetByIndex(int iech) const
{
  for (int is = 0, ns = getNCov(); is < ns; is++)
    _covs[is]->optimizationSetTargetByIndex(iech);
}

/**
 * Calculate the Matrix of covariance between two space points
 * @param p1 Reference of the first space point
 * @param p2 Reference of the second space point
 * @param mat   Covariance matrix (Dimension: nvar * nvar)
 * @param mode  Calculation Options
 *
 * @remarks: Matrix 'mat' should be dimensioned and initialized beforehand
 */
void CovAnisoList::_addEvalCovMatBiPointInPlace(MatrixSquareGeneral& mat,
                                                const SpacePoint& p1,
                                                const SpacePoint& p2,
                                                const CovCalcMode* mode) const
{
  const VectorInt& list = _getListActiveCovariances(mode);
  for (int i = 0, n = (int)list.size(); i < n; i++)
  {
    int j = list[i];
    _covs[j]->addEvalCovMatBiPointInPlace(mat, p1, p2, mode);
  }
}

String CovAnisoList::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (getNCov() <= 0) return sstr.str();

  for (int icov = 0, ncov = getNCov(); icov < ncov; icov++)
  {
    sstr << getCova(icov)->toString();
    if (isFiltered(icov)) sstr << "  (This component is Filtered)" << std::endl;
  }

  // Display the Total Sill (optional)
  if (isStationary())
  {
    if (getNVar() <= 1)
    {
      sstr << "Total Sill     = " << toDouble(getTotalSill(0, 0));
    }
    else
    {
      sstr << toMatrix("Total Sill", VectorString(), VectorString(), 0, getNVar(),
                       getNVar(), getTotalSills().getValues());
    }
  }
  sstr << std::endl;

  return sstr.str();
}

int CovAnisoList::getNCov(bool skipNugget) const
{
  int ncov = (int)_covs.size();
  if (!skipNugget) return ncov;

  int nstruc = 0;
  for (int icov = 0; icov < ncov; icov++)
  {
    if (getCova(icov)->getType() != ECov::NUGGET) nstruc++;
  }
  return nstruc;
}

bool CovAnisoList::isFiltered(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return false;
  return _filtered[icov];
}

bool CovAnisoList::hasRange() const
{
  for (int i = 0, n = getNCov(); i < n; i++)
  {
    if (!getCova(i)->hasRange()) return false;
  }
  return true;
}

bool CovAnisoList::isStationary() const
{
  for (int i = 0, n = getNCov(); i < n; i++)
  {
    if (getCova(i)->getMinOrder() >= 0) return false;
  }
  return true;
}

CovAniso CovAnisoList::extractCova(int icov) const
{
  return *(_covAnisos[icov]);
}

/**
 * @return The Minimum IRF-order induced by the covariances
 */
int CovAnisoList::getCovMinIRFOrder() const
{
  int nmini = -1;
  for (unsigned i = 0, n = getNCov(); i < n; i++)
  {
    int locmini = _covAnisos[i]->getMinOrder();
    if (locmini > nmini) nmini = locmini;
  }
  return nmini;
}

const CovAniso* CovAnisoList::getCova(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return nullptr;
  return _covAnisos[icov];
}
CovAniso* CovAnisoList::getCova(int icov)
{
  if (!_isCovarianceIndexValid(icov)) return nullptr;
  return _covAnisos[icov];
}
void CovAnisoList::setCovAniso(int icov, CovAniso* covs)
{
  if (!_isCovarianceIndexValid(icov)) return;
  _covAnisos[icov] = covs;
  _covs[icov]      = covs;
}
const ECov& CovAnisoList::getCovType(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return ECov::UNKNOWN;
  return _covAnisos[icov]->getType();
}
String CovAnisoList::getCovName(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return String();
  return _covAnisos[icov]->getCovName();
}
double CovAnisoList::getParam(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return 0.;
  return _covAnisos[icov]->getParam();
}
double CovAnisoList::getRange(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return 0.;
  return _covAnisos[icov]->getRange();
}
VectorDouble CovAnisoList::getRanges(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return 0.;
  return _covAnisos[icov]->getRanges();
}
VectorDouble CovAnisoList::getAngles(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return 0.;
  return _covAnisos[icov]->getAnisoAngles();
}
void CovAnisoList::setSill(int icov, int ivar, int jvar, double value)
{
  if (!_isCovarianceIndexValid(icov)) return;
  _covs[icov]->setSill(ivar, jvar, value);
}
void CovAnisoList::setRangeIsotropic(int icov, double range)
{
  if (!_isCovarianceIndexValid(icov)) return;
  _covAnisos[icov]->setRangeIsotropic(range);
}
void CovAnisoList::setParam(int icov, double value)
{
  if (!_isCovarianceIndexValid(icov)) return;
  _covAnisos[icov]->setParam(value);
}
void CovAnisoList::setMarkovCoeffs(int icov, const VectorDouble& coeffs)
{
  if (!_isCovarianceIndexValid(icov)) return;
  _covAnisos[icov]->setMarkovCoeffs(coeffs);
}
void CovAnisoList::setType(int icov, const ECov& type)
{
  if (!_isCovarianceIndexValid(icov)) return;
  _covAnisos[icov]->setType(type);
}

int CovAnisoList::getNGradParam(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return 0;
  return _covAnisos[icov]->getNGradParam();
}

/**
 * Calculate the total sill of the model for given pair of variables
 * Returns TEST as soon as one structure has no sill
 * @param ivar Rank of the first variable
 * @param jvar Rank of the second variable
 */
double CovAnisoList::getTotalSill(int ivar, int jvar) const
{
  double sill_total = 0.;
  for (int icov = 0, ncov = getNCov(); icov < ncov; icov++)
  {
    const CovAniso* cova = getCova(icov);
    if (cova->getMinOrder() >= 0) return TEST;
    sill_total += cova->getSill(ivar, jvar);
  }
  return sill_total;
}

bool CovAnisoList::_isCovarianceIndexValid(int icov) const
{
  return checkArg("Covariance Index", icov, getNCov());
}

/**
 * Returns the largest range (in any direction in the heterotopic case)
 * @return
 */
double CovAnisoList::getMaximumDistance() const

{
  double maxdist = 0.;
  for (int icov = 0, ncov = getNCov(); icov < ncov; icov++)
  {
    const CovAniso* cova = getCova(icov);
    if (!cova->hasRange()) continue;
    double range = cova->getRange();
    if (range > maxdist) maxdist = range;
  }
  return maxdist;
}

void CovAnisoList::_setContext(const CovContext& ctxt)
{
  for (auto& e: _covAnisos)
  {
    e->setContext(ctxt);
  }
}

void CovAnisoList::copyCovContext(const CovContext& ctxt)
{
  int number = (int)_covAnisos.size();
  for (int i = 0; i < number; i++) _covAnisos[i]->copyCovContext(ctxt);
}

void CovAnisoList::normalize(double sill, int ivar, int jvar)
{
  double covval = 0.;
  for (int i = 0, n = getNCov(); i < n; i++) covval += _covs[i]->eval0(ivar, jvar);

  if (covval <= 0. || isEqual(covval, sill)) return;
  double ratio = sill / covval;

  for (int i = 0, n = getNCov(); i < n; i++)
  {
    CovAniso* cov = _covAnisos[i];
    cov->setSill(cov->getSill(ivar, jvar) * ratio);
  }
}

bool CovAnisoList::hasNugget() const
{
  for (int is = 0, ns = getNCov(); is < ns; is++)
  {
    if (getCovType(is) == ECov::NUGGET) return true;
  }
  return false;
}

int CovAnisoList::getRankNugget() const
{
  for (int is = 0, ns = getNCov(); is < ns; is++)
  {
    if (getCovType(is) == ECov::NUGGET) return is;
  }
  return -1;
}

const CovAnisoList* CovAnisoList::createReduce(const VectorInt& validVars) const
{
  CovAnisoList* newcovlist = this->clone();

  for (int is = 0, ns = getNCov(); is < ns; is++)
  {
    CovAniso* covs = newcovlist->getCova(is);
    newcovlist->setCovAniso(is, covs->createReduce(validVars));
  }
  return newcovlist;
}

/**
 * Returns the Ball radius (from the first covariance of _covaList)
 * @return Value of the Ball Radius (if defined, i.e. for Numerical Gradient calculation)
 */
double CovAnisoList::getBallRadius() const
{
  // Check is performed on the first covariance
  const CovAniso* cova = getCova(0);
  double ball_radius   = cova->getBallRadius();
  if (!FFFF(ball_radius)) return ball_radius;
  return 0.;
}

int CovAnisoList::hasExternalCov() const
{
  for (int icov = 0; icov < (int)getNCov(); icov++)
  {
    if (getCovType(icov) == ECov::FUNCTION) return 1;
  }
  return 0;
}

bool CovAnisoList::isChangeSupportDefined() const
{
  if (getAnam() == nullptr)
  {
    return false;
  }
  return getAnam()->isChangeSupportDefined();
}

const AnamHermite* CovAnisoList::getAnamHermite() const
{
  const AAnam* anam = getAnam();
  if (anam == nullptr) return nullptr;
  const AnamHermite* anamH = dynamic_cast<const AnamHermite*>(anam);
  return anamH;
}

const EModelProperty& CovAnisoList::getCovMode() const
{
  if (dynamic_cast<const CovLMCTapering*>(this) != nullptr) return EModelProperty::TAPE;

  if (dynamic_cast<const CovLMCConvolution*>(this) != nullptr)
    return EModelProperty::CONV;

  if (dynamic_cast<const CovLMCAnamorphosis*>(this) != nullptr)
    return EModelProperty::ANAM;

  if (dynamic_cast<const CovLMGradient*>(this) != nullptr) return EModelProperty::GRAD;

  return EModelProperty::NONE;
}

void CovAnisoList::setOptimEnabled(bool status)
{
  for (auto& e: _covAnisos)
  {
    e->setOptimEnabled(status);
  }
}