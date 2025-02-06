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

CovList::~CovList()
{
  delAllCov();
}

void CovList::addCovList(const CovList* covs)
{
  for (int icov = 0, ncov = covs->getNCov(); icov < ncov; icov++)
    addCov(covs->getCova(icov));
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

    if (! cov->getContext().isEqual(_covs[0]->getContext()))
    {
      messerr("Error: Covariances in the same CovList should share the same Context");
      messerr("Operation is cancelled");
      return;
    }
  }
  _covs.push_back(cov);
  _filtered.push_back(false);
  _updateLists();
}

void CovList::_updateLists()
{
  int ncov = getNCov();
  _allActiveCovList = VH::sequence(ncov);

  _activeCovList.clear();
  for (int icov = 0; icov < ncov; icov++)
    if (!_filtered[icov]) _activeCovList.push_back(icov);

  _allActiveCov = _activeCovList.size() == _allActiveCovList.size();
}

void CovList::delCov(int icov)
{
  if (! _isCovarianceIndexValid(icov)) return;
  delete _covs[icov];
  _covs.erase(_covs.begin() + icov);
  _filtered.erase(_filtered.begin() + icov);
  _delCov(icov);
  _updateLists();
}

void CovList::delAllCov()
{
  for (auto &e: _covs)
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
  for (const auto &e :_covs)
  {
    nostat = nostat || e->isNoStat();
  }
  return nostat;
}
void CovList::setFiltered(int icov, bool filtered)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _filtered[icov] = filtered;
  _updateLists();
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

double CovList::eval0(int ivar, int jvar, const CovCalcMode* mode) const
{
  double cov      = 0.;
  const VectorInt& list = _getListActiveCovariances(mode);
  for (int i = 0, n = (int)list.size(); i < n; i++)
    cov += _covs[list[i]]->eval0(ivar, jvar, mode);
  return cov;
}

/**
 * Calculate the Matrix of covariance for zero distance
 * @param mat   Covariance matrix (Dimension: nvar * nvar)
 * @param mode  Calculation Options
 *
 * @remarks: Matrix 'mat' should be dimensioned and initialized beforehand
 */
void CovList::addEval0CovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                    const CovCalcMode *mode) const
{
  const VectorInt& list = _getListActiveCovariances(mode);
  for (int i = 0, n = (int)list.size(); i < n; i++)
    _covs[list[i]]->addEval0CovMatBiPointInPlace(mat, mode);
}

void CovList::_optimizationSetTarget(const SpacePoint& pt) const
{
  for (int is = 0, ns = getNCov(); is < ns; is++)
    _covs[is]->optimizationSetTarget(pt);
}

void CovList::optimizationSetTargetByIndex(int iech) const
{
  for (int is = 0, ns = getNCov(); is < ns; is++)
    _covs[is]->optimizationSetTargetByIndex(iech);
}

double CovList::eval(const SpacePoint& p1,
                     const SpacePoint& p2,
                     int ivar,
                     int jvar,
                     const CovCalcMode* mode) const
{
  double cov      = 0.;
  const VectorInt& list = _getListActiveCovariances(mode);
  for (int i = 0, n = (int)list.size(); i < n; i++)
    cov += _covs[list[i]]->eval(p1, p2, ivar, jvar, mode);
  return cov;
}

double CovList::_loadAndEval(const SpacePoint& p1,
                          const SpacePoint&p2,
                          int ivar,
                          int jvar,
                          const CovCalcMode *mode) const
{
  double res      = 0.;
  const VectorInt& list = _getListActiveCovariances(mode);
  for (int i = 0, n = (int)list.size(); i < n; i++)
    res += _covs[list[i]]->loadAndEval(p1, p2, ivar, jvar, mode);
  return res;
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
void CovList::_addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                                const SpacePoint &p1,
                                                const SpacePoint &p2,
                                                const CovCalcMode *mode) const
{
  const VectorInt& list = _getListActiveCovariances(mode);
  for (int i = 0, n = (int)list.size(); i < n; i++)
    _covs[list[i]]->addEvalCovMatBiPointInPlace(mat, p1, p2, mode);
}

void CovList::_loadAndAddEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,const SpacePoint& p1,const SpacePoint&p2,
                                              const CovCalcMode *mode) const
{
  const VectorInt& list = _getListActiveCovariances(mode);
  for (int i = 0, n = (int)list.size(); i < n; i++)
    _covs[list[i]]->loadAndAddEvalCovMatBiPointInPlace(mat, p1, p2, mode);
}

void CovList::_load(const SpacePoint& p, bool case1) const
{
  const VectorInt& list = _getListActiveCovariances(nullptr);
  for (int i = 0, n = (int)list.size(); i < n; i++)
    _covs[list[i]]->load(p, case1);
}

String CovList::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (getNCov() <= 0) return sstr.str();

  for (int icov = 0, ncov = getNCov(); icov < ncov; icov++)
  {
    sstr << getCova(icov)->toString();
    if (isFiltered(icov))
      sstr << "  (This component is Filtered)" << std::endl;
  }
  sstr << std::endl;
 
  return sstr.str();
}


int  CovList::getNCov() const
{
  int ncov = (int) _covs.size();
  return ncov;
}

bool CovList::isFiltered(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return false;
  return _filtered[icov];
}

bool CovList::isAllActiveCovList() const
{
  for (int i=0, n=getNCov(); i<n; i++)
  {
    if (_filtered[i]) return false;
  }
  return true;
}

const CovBase* CovList::getCova(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return nullptr;
  return _covs[icov];
}

void CovList::setCova(int icov,const CovBase* covs)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covs[icov] = covs;
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
  if(! _isCovarianceIndexValid(icov)) return 0.;
  return _covs[icov]->getSill(ivar, jvar);
}
void CovList::setSill(int icov, int ivar, int jvar, double value)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covs[icov]->setSill(ivar, jvar, value);
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
    const CovBase* cova = getCova(icov);
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
      mat.setValue(ivar,jvar,getTotalSill(ivar,jvar));
  return mat;
}

bool CovList::_isCovarianceIndexValid(int icov) const
{
  return checkArg("Covariance Index", icov, getNCov());
}

void CovList::_optimizationPreProcess(const std::vector<SpacePoint>& p) const
{
  for (const auto &e :_covs)
  {
    e->optimizationPreProcess(p);
  }
}

void CovList::_optimizationPostProcess() const
{
	for (int is = 0, ns = getNCov(); is < ns; is++)
		_covs[is]->optimizationPostProcess();
}

void CovList::_manage(const Db* db1,const Db* db2)  const
{
  for (const auto &e : _covs)
  {
    e->manage(db1,db2);
  }
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
  for (const auto &e : _covs)
  {
    e->updateCovByPoints(icas1,iech1,icas2,iech2);
  }
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
