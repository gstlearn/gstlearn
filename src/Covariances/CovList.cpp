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
#include "Covariances/CovList.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CovBase.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovContext.hpp"
#include "Enum/ECalcMember.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Space/ASpace.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "Db/Db.hpp"
#include "Space/SpacePoint.hpp"
#include "Basic/VectorHelper.hpp"
#include "geoslib_define.h"

#include <math.h>
#include <vector>

CovList::CovList(const CovContext& ctxt)
  : ACov(ctxt)
  , _covs()
  , _filtered()
  , _allActiveCov(true)
  , _allActiveCovList()
  , _activeCovList()
{
}

CovList::CovList(const CovList& r)
  : ACov(r)
{
  for (const auto* e: r._covs)
  {
    _covs.push_back((CovBase*)e->clone());
  }
  _filtered         = r._filtered;
  _allActiveCov     = r._allActiveCov;
  _allActiveCovList = r._allActiveCovList;
  _activeCovList    = r._activeCovList;
  _updateLists();
}

CovList& CovList::operator=(const CovList& r)

{
  if (this != &r)
  {
    ACov::operator=(r);
    for (const auto* e: r._covs)
    {
      _covs.push_back((CovBase*)e->clone());
    }
    _filtered         = r._filtered;
    _allActiveCov     = r._allActiveCov;
    _allActiveCovList = r._allActiveCovList;
    _activeCovList    = r._activeCovList;
    _updateLists();
  }
  return *this;
}

CovList::~CovList()
{
  delAllCov();
}

void CovList::addCovList(const CovList* covs)
{
  for (int icov = 0, ncov = covs->getNCov(); icov < ncov; icov++)
    addCov(covs->getCov(icov));
}

void CovList::addCov(const CovBase* cov)
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

    if (!cov->getContext().isEqual(_covs[0]->getContext()))
    {
      messerr("Error: Covariances in the same CovList should share the same Context");
      messerr("Operation is cancelled");
      return;
    }
  }
  _covs.push_back((CovBase*)cov->clone());
  _filtered.push_back(false);
  _updateLists();
}

void CovList::_updateLists()
{
  int ncov          = getNCov();
  _allActiveCovList = VH::sequence(ncov);

  _activeCovList.clear();
  for (int icov = 0; icov < ncov; icov++)
    if (!_filtered[icov]) _activeCovList.push_back(icov);

  _allActiveCov = _activeCovList.size() == _allActiveCovList.size();
}

void CovList::delCov(int icov)
{
  if (!_isCovarianceIndexValid(icov)) return;
  delete _covs[icov];
  _covs.erase(_covs.begin() + icov);
  _filtered.erase(_filtered.begin() + icov);
  _delCov(icov);
  _updateLists();
}

void CovList::delAllCov()
{
  for (auto& e: _covs)
  {
    delete e;
  }
  _covs.clear();
  _filtered.clear();
  _allActiveCovList.clear();
  _activeCovList.clear();
  _delAllCov();
}

bool CovList::isNoStat() const
{
  bool nostat = false;
  for (const auto& e: _covs)
  {
    nostat = nostat || e->isNoStat();
  }
  return nostat;
}

bool CovList::isConsistent(const ASpace* /*space*/) const
{
  /// TODO : CovList::isConsistent
  return true;
}

int CovList::getNVar() const
{
  if (getNCov() > 0)
    return _covs[0]->getNVar();
  return 0;
}

/**
 * @brief Check if all covariances should be processed
 *
 * @param mode  CovCalcMode structure
 * @return True  all covariances should be treated;
 * @return False only the non-filtered ones are returned
 */
const VectorInt& CovList::_getListActiveCovariances(const CovCalcMode* mode) const
{
  if (mode == nullptr) return _allActiveCovList;
  if (_allActiveCov) return _allActiveCovList;
  if (mode->getMember() != ECalcMember::LHS) return _activeCovList;
  return _allActiveCovList;
}

int CovList::addEvalCovVecRHSInPlace(vect vect,
                                     const VectorInt& index1,
                                     SpacePoint& pin,
                                     SpacePoint& pout,
                                     const int iech2) const
{
  CovCalcMode mode(ECalcMember::RHS);
  const VectorInt& list = _getListActiveCovariances(&mode);
  for (const auto& j: list.getVector())
    _covs[j]->addEvalCovVecRHSInPlace(vect, index1, pin, pout, iech2);
  return 0;
}

double CovList::eval0(int ivar, int jvar, const CovCalcMode* mode) const
{
  double cov            = 0.;
  const VectorInt& list = _getListActiveCovariances(mode);
  for (const auto& j: list.getVector())
    cov += _covs[j]->eval0(ivar, jvar, mode);
  return cov;
}

void CovList::_optimizationSetTarget(SpacePoint& pt) const
{
  for (const auto& e: _covs)
    e->optimizationSetTarget(pt);
}

void CovList::setOptimEnabled(bool flag) const
{
  ACov::setOptimEnabled(flag);
  for (const auto& e: _covs)
    e->setOptimEnabled(flag);
}

double CovList::_eval(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ivar,
                      int jvar,
                      const CovCalcMode* mode) const
{
  double cov            = 0.;
  const VectorInt& list = _getListActiveCovariances(mode);
  for (const auto& j: list.getVector())
    cov += _covs[j]->evalCov(p1, p2, ivar, jvar, mode);
  return cov;
}

void CovList::_load(const SpacePoint& p, bool case1) const
{
  const VectorInt& list = _getListActiveCovariances(nullptr);
  for (const auto& j: list.getVector())
    _covs[j]->load(p, case1);
}

String CovList::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (getNCov() <= 0) return sstr.str();

  for (int icov = 0, ncov = getNCov(); icov < ncov; icov++)
  {
    sstr << getCov(icov)->toString();
    if (isFiltered(icov))
      sstr << "  (This component is Filtered)" << std::endl;
  }
  sstr << std::endl;

  return sstr.str();
}

int CovList::getNCov() const
{
  int ncov = (int)_covs.size();
  return ncov;
}

bool CovList::isFiltered(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return false;
  return _filtered[icov];
}

bool CovList::isAllActiveCovList() const
{
  for (int i = 0, n = getNCov(); i < n; i++)
  {
    if (_filtered[i]) return false;
  }
  return true;
}

const CovBase* CovList::getCov(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return nullptr;
  return _covs[icov];
}

void CovList::setCov(int icov, const CovBase* covs)
{
  if (!_isCovarianceIndexValid(icov)) return;
  delete _covs[icov];
  _covs[icov] = (CovBase*)covs->clone();
}

const ECov& CovList::getCovType(int icov) const
{
  DECLARE_UNUSED(icov)
  return ECov::UNKNOWN;
}

String CovList::getCovName(int icov) const
{
  DECLARE_UNUSED(icov)
  ECov unknown = ECov::UNKNOWN;
  return std::string(unknown.getKey());
}

const MatrixSquareSymmetric& CovList::getSills(int icov) const
{
  return _covs[icov]->getSill();
}
double CovList::getSill(int icov, int ivar, int jvar) const
{
  if (!_isCovarianceIndexValid(icov)) return 0.;
  return _covs[icov]->getSill(ivar, jvar);
}
void CovList::setSill(int icov, int ivar, int jvar, double value)
{
  if (!_isCovarianceIndexValid(icov)) return;
  _covs[icov]->setSill(ivar, jvar, value);
}
void CovList::setSills(int icov, const MatrixSquareSymmetric& sills)
{
  if (!_isCovarianceIndexValid(icov)) return;
  _covs[icov]->setSill(sills);
}
/**
 * Calculate the total sill of the model for given pair of variables
 * @param ivar Rank of the first variable
 * @param jvar Rank of the second variable
 */
double CovList::getTotalSill(int ivar, int jvar) const
{
  double sill_total = 0.;
  for (int icov = 0, ncov = getNCov(); icov < ncov; icov++)
  {
    const CovBase* cova = getCov(icov);
    sill_total += cova->getSill(ivar, jvar);
  }
  return sill_total;
}

MatrixSquareSymmetric CovList::getTotalSills() const
{
  int nvar = getNVar();
  MatrixSquareSymmetric mat(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
      mat.setValue(ivar, jvar, getTotalSill(ivar, jvar));
  return mat;
}

bool CovList::_isCovarianceIndexValid(int icov) const
{
  return checkArg("Covariance Index", icov, getNCov());
}

void CovList::_optimizationPreProcess(int mode, const std::vector<SpacePoint>& ps) const
{
  for (const auto& e: _covs)
    e->optimizationPreProcess(mode, ps);
}

SpacePoint& CovList::_optimizationLoadInPlace(int iech, int mode, int rank) const
{
  for (int is = 1, ns = getNCov(); is < ns; is++)
    (void)_covs[is]->optimizationLoadInPlace(iech, mode, rank);
  return _covs[0]->optimizationLoadInPlace(iech, mode, rank);
}

void CovList::_optimizationPostProcess() const
{
  for (const auto& e: _covs)
    e->optimizationPostProcess();
}

void CovList::_manage(const Db* db1, const Db* db2) const
{
  for (const auto& e: _covs)
    e->manage(db1, db2);
}

/**
 * Update the Model according to the Non-stationary parameters
 * @param icas1 Type of first Db: 1 for Input; 2 for Output
 * @param iech1 Rank of the target within Db1 (or -1)
 * @param icas2 Type of first Db: 1 for Input; 2 for Output
 * @param iech2 Rank of the target within Dbout (or -2)
 */

void CovList::updateCovByPoints(int icas1, int iech1, int icas2, int iech2) const
{
  for (const auto& e: _covs)
    e->updateCovByPoints(icas1, iech1, icas2, iech2);
}

void CovList::setActiveCovListFromOne(int keepOnlyCovIdx) const
{
  _allActiveCov = true;
  _activeCovList.clear();
  if (keepOnlyCovIdx >= 0)
  {
    _activeCovList.push_back(keepOnlyCovIdx);
    _allActiveCov = false;
  }
}

/**
 * Set the list of active covariances from an interval
 * @param inddeb Lower bound of the interval (included)
 * @param indto  Upper bound of the interval (excluded)
 */
void CovList::setActiveCovListFromInterval(int inddeb, int indto) const
{
  _activeCovList.clear();
  for (int i = inddeb; i < indto; i++) _activeCovList.push_back(i);
  _allActiveCov = false;
}

void CovList::setActiveCovList(const VectorInt& activeCovList, bool allActiveCov) const
{
  _activeCovList = activeCovList;
  _allActiveCov  = allActiveCov;
}

void CovList::_setContext(const CovContext& ctxt)
{
  for (auto& e: _covs)
  {
    e->setContext(ctxt);
  }
}

void CovList::copyCovContext(const CovContext& ctxt)
{
  int number = (int)_covs.size();
  for (int i = 0; i < number; i++) _covs[i]->copyCovContext(ctxt);
}

void CovList::normalize(double sill, int ivar, int jvar)
{
  double covval = 0.;
  for (int i = 0, n = getNCov(); i < n; i++) covval += _covs[i]->eval0(ivar, jvar);

  if (covval <= 0. || isEqual(covval, sill)) return;
  double ratio = sill / covval;

  for (int i = 0, n = getNCov(); i < n; i++)
    _covs[i]->setSill(_covs[i]->getSill(ivar, jvar) * ratio);
}

void CovList::setCovFiltered(int icov, bool filtered)
{
  if (!_isCovarianceIndexValid(icov)) return;
  _filtered[icov] = filtered;
  _updateLists();
}
