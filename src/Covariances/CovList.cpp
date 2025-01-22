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
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Space/ASpace.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "Db/Db.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"

#include <math.h>
#include <vector>

CovList::CovList(const ASpaceSharedPtr& space)
: ACov(space),
  _covs(),
  _filtered()
{
}


CovList::~CovList()
{
  delAllCov();
}

void CovList::addCovList(const CovList* covs)
{
  for (int icov = 0, ncov = covs->getCovaNumber(); icov < ncov; icov++)
    addCov(covs->getCova(icov));
}

void CovList::addCov(const CovBase* cov)
{
  if (getCovaNumber() == 0)
  {
    setNVar(cov->getNVariables());
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
}

void CovList::delCov(int icov)
{
  if (! _isCovarianceIndexValid(icov)) return;
  delete _covs[icov];
  _covs.erase(_covs.begin() + icov);
  _filtered.erase(_filtered.begin() + icov);
  _delCov(icov);
}

void CovList::delAllCov()
{
  for (auto &e: _covs)
  {
    delete e;
  }
  _covs.clear();
  _filtered.clear();
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
}

bool CovList::isConsistent(const ASpace* /*space*/) const
{
  /// TODO : CovList::isConsistent
  return true;
}

int CovList::getNVariables() const
{
  if (getCovaNumber() > 0)
    return _covs[0]->getNVariables();
  return 0;
}

bool CovList::_considerAllCovariances(const CovCalcMode* mode)
{
  if (mode == nullptr) return true;
  if (mode->isAllActiveCov()) return true;
  return false;
}

double CovList::eval0(int ivar, int jvar, const CovCalcMode* mode) const
{
  double cov = 0.;

  if (_considerAllCovariances(mode))
  {
    for (int i=0, n=getCovaNumber(); i<n; i++)
      cov += _covs[i]->eval0(ivar, jvar, mode);
  }
  else
  {
    for (int i=0, n=(int) mode->getActiveCovList().size(); i<n; i++)
      cov += _covs[mode->getActiveCovList(i)]->eval0(ivar, jvar, mode);
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
void CovList::addEval0CovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                    const CovCalcMode *mode) const
{
  if (_considerAllCovariances(mode))
  {
    for (int i=0, n=getCovaNumber(); i<n; i++)
    {
      _covs[i]->addEval0CovMatBiPointInPlace(mat, mode);
    }
  }
  else
  {
    for (int i=0, n=(int) mode->getActiveCovList().size(); i<n; i++)
    {
      _covs[mode->getActiveCovList(i)]->addEval0CovMatBiPointInPlace(mat, mode);
    }
  }
}

void CovList::_optimizationSetTarget(const SpacePoint& pt) const
{
  for (int is = 0, ns = getCovaNumber(); is < ns; is++)
    _covs[is]->optimizationSetTarget(pt);
}

void CovList::optimizationSetTargetByIndex(int iech) const
{
  for (int is = 0, ns = getCovaNumber(); is < ns; is++)
    _covs[is]->optimizationSetTargetByIndex(iech);
}

double CovList::eval(const SpacePoint& p1,
                           const SpacePoint& p2,
                           int ivar,
                           int jvar,
                           const CovCalcMode* mode) const
{
  double cov = 0.;

  if (_considerAllCovariances(mode))
  {
    for (int i=0, n=getCovaNumber(); i<n; i++)
      cov += _covs[i]->eval(p1, p2, ivar, jvar, mode);
  }
  else
  {
    for (int i=0, n=(int) mode->getActiveCovList().size(); i<n; i++)
      cov += _covs[mode->getActiveCovList(i)]->eval(p1, p2, ivar, jvar, mode);
  }
  return cov;
}

double CovList::_loadAndEval(const SpacePoint& p1,
                          const SpacePoint&p2,
                          int ivar,
                          int jvar,
                          const CovCalcMode *mode) const
{ 
  double res = 0.;
  for (const auto &e : _covs)
  {
    res += e->loadAndEval(p1, p2, ivar, jvar, mode);
  }
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
  if (_considerAllCovariances(mode))
  {
    for (int i=0, n=getCovaNumber(); i<n; i++)
    {
      _covs[i]->addEvalCovMatBiPointInPlace(mat,p1, p2, mode);
    }
  }
  else
  {
    for (int i=0, n=(int) mode->getActiveCovList().size(); i<n; i++)
    {
      _covs[mode->getActiveCovList(i)]->addEvalCovMatBiPointInPlace(mat,p1, p2, mode);
    }
  }
}

void CovList::_loadAndAddEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,const SpacePoint& p1,const SpacePoint&p2,
                                              const CovCalcMode *mode) const
{
  for (const auto &e : _covs)
  {
    e->loadAndAddEvalCovMatBiPointInPlace(mat,p1,p2,mode);
  }
}

String CovList::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (getCovaNumber() <= 0) return sstr.str();

  for (int icov = 0, ncov = getCovaNumber(); icov < ncov; icov++)
  {
    sstr << getCova(icov)->toString();
    if (isFiltered(icov))
      sstr << "  (This component is Filtered)" << std::endl;
  }
  sstr << std::endl;
 
  return sstr.str();
}

int  CovList::getCovaNumber() const
{
  int ncov = (int) _covs.size();
  return ncov;
}

bool CovList::isFiltered(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return false;
  return _filtered[icov];
}

VectorInt CovList::getActiveCovList() const
{
  VectorInt actives;
  for (int i=0, n=getCovaNumber(); i<n; i++)
  {
    if (_filtered[i]) continue;
    actives.push_back(i);
  }
  return actives;
}

bool CovList::isAllActiveCovList() const
{
  for (int i=0, n=getCovaNumber(); i<n; i++)
  {
    if (_filtered[i]) return false;
  }
  return true;
}

VectorInt CovList::getAllActiveCovList() const
{
  VectorInt actives;
  for (int i=0, n=getCovaNumber(); i<n; i++)
  {
    actives.push_back(i);
  }
  return actives;
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
const ECov& CovList::getType(int icov) const
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
  
const MatrixSquareSymmetric& CovList::getSill(int icov) const
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
  for (int icov = 0, ncov = getCovaNumber(); icov < ncov; icov++)
  {
    const CovBase* cova = getCova(icov);
    sill_total += cova->getSill(ivar, jvar);
  }
  return sill_total;
}

MatrixSquareSymmetric CovList::getTotalSills() const
{
  int nvar = getNVariables();
  MatrixSquareSymmetric mat(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
      mat.setValue(ivar,jvar,getTotalSill(ivar,jvar));
  return mat;
}

bool CovList::_isCovarianceIndexValid(int icov) const
{
  return checkArg("Covariance Index", icov, getCovaNumber());
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
	for (int is = 0, ns = getCovaNumber(); is < ns; is++)
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

