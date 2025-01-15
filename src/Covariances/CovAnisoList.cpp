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

#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovList.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Space/ASpace.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "Db/Db.hpp"
#include "Space/SpacePoint.hpp"

#include <math.h>
#include <vector>

CovAnisoList::CovAnisoList(const ASpace* space)
: CovList(space),
  _covAnisos()
{
}

CovAnisoList::CovAnisoList(const CovAnisoList &r)
: CovList(r._space),
  _covAnisos()
{
  _filtered = r._filtered;

  for (auto* e: r._covAnisos)
  {
    _pushCov(e->clone());
  }
}

CovAnisoList& CovAnisoList::operator=(const CovAnisoList &r)
{
  if (this != &r)
  {
    for (auto *e: r._covAnisos)
    {
     _pushCov(e->clone());
    }
    _filtered = r._filtered;
  }
  return *this;
}

CovAnisoList::~CovAnisoList()
{
  delAllCov();
}

void CovAnisoList::addCovList(const CovAnisoList* covs)
{
  for (int icov = 0, ncov = covs->getCovaNumber(); icov < ncov; icov++)
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
      messerr("Error: Covariances in the same CovAnisoList should share the same Context");
      messerr("Operation is cancelled");
      return;
    }
  }
  _pushCov(cov);
  _filtered.push_back(false);
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
void CovAnisoList::setFiltered(int icov, bool filtered)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _filtered[icov] = filtered;
}

bool CovAnisoList::isConsistent(const ASpace* /*space*/) const
{
  /// TODO : CovAnisoList::isConsistent
  return true;
}

int CovAnisoList::getNVariables() const
{
  if (getCovaNumber() > 0)
    return _covs[0]->getNVariables();
  return 0;
}


double CovAnisoList::eval0(int ivar, int jvar, const CovCalcMode* mode) const
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
void CovAnisoList::addEval0CovMatBiPointInPlace(MatrixSquareGeneral &mat,
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

/**
 * Evaluate the covariance rectangular matrix between samples of input 'db1' and
'db2'
 * @param db1 Input Db
 * @param db2 Output db
 * @param ivar0 Rank of the first variable (-1 for all variables)
 * @param jvar0 Rank of the second variable (-1 for all variables)
 * @param nbgh1 Vector of indices of active samples in db1 (optional)
 * @param nbgh2 Vector of indices of active samples in db2 (optional)
 * @param mode CovCalcMode structure
 * @param cleanOptim True if Optimization internal arrays must be cleaned at end
 * @return
 */
MatrixRectangular CovAnisoList::evalCovMatrixOptim(const Db* db1,
                                                   const Db* db2,
                                                   int ivar0,
                                                   int jvar0,
                                                   const VectorInt& nbgh1,
                                                   const VectorInt& nbgh2,
                                                   const CovCalcMode* mode,
                                                   bool cleanOptim) const
{
  MatrixRectangular mat;
  SpacePoint p2;
  if (db2 == nullptr) db2 = db1;
  VectorInt ivars = _getActiveVariables(ivar0);
  if (ivars.empty()) return mat;
  VectorInt jvars = _getActiveVariables(jvar0);
  if (jvars.empty()) return mat;

  // Prepare the Optimization for covariance calculation
  optimizationPreProcess(db1);

  // Create the sets of Vector of valid sample indices per variable (not masked and defined)
  VectorVectorInt index1 = db1->getMultipleRanksActive(ivars, nbgh1);
  VectorVectorInt index2 = db2->getMultipleRanksActive(jvars, nbgh2);

  // Creating the matrix
  int neq1 = VH::count(index1);
  int neq2 = VH::count(index2);
  if (neq1 <= 0 || neq2 <= 0)
  {
    messerr("The returned matrix does not have any valid sample for any valid variable");
    return mat;
  }
  mat.resize(neq1, neq2);

  // Loop on the second variable
  int icol = 0;
  for (int rvar2 = 0, nvar2 = (int) jvars.size(); rvar2 < nvar2; rvar2++)
  {
    int ivar2 = jvars[rvar2];

    // Loop on the second sample
    int nech2s = (int) index2[rvar2].size();
    for (int rech2 = 0; rech2 < nech2s; rech2++)
    {
      int iech2 = index2[rvar2][rech2];
      db2->getSampleAsSPInPlace(p2, iech2);
      optimizationSetTarget(p2);

      // Loop on the basic structures
      for (int i = 0, n = getCovaNumber(); i < n; i++)
         _covAnisos[i]->evalOptimInPlace(mat, ivars, index1, ivar2, icol, mode, false);
      icol++;
    }
  }

  if (cleanOptim) optimizationPostProcess();
  return mat;
}

void CovAnisoList::_optimizationSetTarget(const SpacePoint& pt) const
{
  for (int is = 0, ns = getCovaNumber(); is < ns; is++)
    _covs[is]->optimizationSetTarget(pt);
}

void CovAnisoList::optimizationSetTargetByIndex(int iech) const
{
  for (int is = 0, ns = getCovaNumber(); is < ns; is++)
    _covs[is]->optimizationSetTargetByIndex(iech);
}

/**
 * Evaluate the covariance matrix between samples of input 'db1'
 * @param db1 Input Db
 * @param ivar0 Rank of the first variable (-1 for all variables)
 * @param nbgh1 Vector of indices of active samples in db1 (optional)
 * @param mode CovCalcMode structure
 * @param cleanOptim When True, clean Optimization internal arrays at end
 * @return
 */
MatrixSquareSymmetric
CovAnisoList::evalCovMatrixSymmetricOptim(const Db* db1,
                                          int ivar0,
                                          const VectorInt& nbgh1,
                                          const CovCalcMode* mode,
                                          bool cleanOptim) const
{
  MatrixSquareSymmetric mat;
  SpacePoint p2;

  VectorInt ivars = _getActiveVariables(ivar0);
  if (ivars.empty()) return mat;

  // Prepare the Optimization for covariance calculation
  optimizationPreProcess(db1);

  // Create the sets of Vector of valid sample indices per variable (not masked and defined)
  VectorVectorInt index1 = db1->getMultipleRanksActive(ivars, nbgh1, true, true);

  // Creating the matrix
  int neq1 = VH::count(index1);
  if (neq1 <= 0)
  {
    messerr("The returned matrix does not have any valid sample for any valid variable");
    return mat;
  }
  mat.resize(neq1, neq1);

  // Loop on the second variable
  int icol = 0;
  for (int rvar2 = 0, nvar2 = (int) ivars.size(); rvar2 < nvar2; rvar2++)
  {
    int ivar2 = ivars[rvar2];

    // Loop on the second sample
    int nech2s = (int) index1[rvar2].size();
    for (int rech2 = 0; rech2 < nech2s; rech2++)
    {
      int iech2 = index1[rvar2][rech2];

      optimizationSetTargetByIndex(iech2);

      // Loop on the basic structures
      for (int i = 0, n = getCovaNumber(); i < n; i++)
         _covAnisos[i]->evalOptimInPlace(mat, ivars, index1, ivar2, icol, mode, true);

      icol++;
    }
  }

  // Update the matrix due to presence of Variance of Measurement Error
  _updateCovMatrixSymmetricVerr(db1, &mat, ivars, index1);

  if (cleanOptim) optimizationPostProcess();
  return mat;
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
void CovAnisoList::_addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
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

String CovAnisoList::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (getCovaNumber() <= 0) return sstr.str();

  for (int icov = 0, ncov = getCovaNumber(); icov < ncov; icov++)
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

int  CovAnisoList::getCovaNumber(bool skipNugget) const
{
  int ncov = (int) _covs.size();
  if (! skipNugget) return ncov;

  int nstruc = 0;
  for (int icov = 0; icov < ncov; icov++)
  {
    if (getCova(icov)->getType() != ECov::NUGGET) nstruc++;
  }
  return nstruc;
}

bool CovAnisoList::isFiltered(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return false;
  return _filtered[icov];
}

bool CovAnisoList::hasRange() const
{
  for (int i=0, n=getCovaNumber(); i<n; i++)
  {
    if (!getCova(i)->hasRange())
      return false;
  }
  return true;
}

bool CovAnisoList::isStationary() const
{
  for (int i=0, n=getCovaNumber(); i<n; i++)
  {
    if (getCova(i)->getMinOrder() >= 0)
      return false;
  }
  return true;
}

VectorInt CovAnisoList::getActiveCovList() const
{
  VectorInt actives;
  for (int i=0, n=getCovaNumber(); i<n; i++)
  {
    if (_filtered[i]) continue;
    actives.push_back(i);
  }
  return actives;
}

bool CovAnisoList::isAllActiveCovList() const
{
  for (int i=0, n=getCovaNumber(); i<n; i++)
  {
    if (_filtered[i]) return false;
  }
  return true;
}

VectorInt CovAnisoList::getAllActiveCovList() const
{
  VectorInt actives;
  for (int i=0, n=getCovaNumber(); i<n; i++)
  {
    actives.push_back(i);
  }
  return actives;
}

CovAniso CovAnisoList::extractCova(int icov) const
{
  return *(_covAnisos[icov]);
}

/**
 * @return The Minimum IRF-order induced by the covariances
 */
int CovAnisoList::getCovaMinIRFOrder() const
{
  int nmini = -1;
  for (unsigned i = 0, n = getCovaNumber(); i<n; i++)
  {
    int locmini = _covAnisos[i]->getMinOrder();
    if (locmini > nmini) nmini = locmini;
  }
  return nmini;
}

const CovAniso* CovAnisoList::getCova(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return nullptr;
  return _covAnisos[icov];
}
CovAniso* CovAnisoList::getCova(int icov)
{
  if (! _isCovarianceIndexValid(icov)) return nullptr;
  return _covAnisos[icov];
}
void CovAnisoList::setCovAniso(int icov, CovAniso* covs)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covAnisos[icov] = covs;
  _covs[icov] = covs;
}
const ECov& CovAnisoList::getType(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return ECov::UNKNOWN;
  return _covAnisos[icov]->getType();
}
String CovAnisoList::getCovName(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return String();
  return _covAnisos[icov]->getCovName();
}
double CovAnisoList::getParam(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return 0.;
  return _covAnisos[icov]->getParam();
}
double CovAnisoList::getRange(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return 0.;
  return _covAnisos[icov]->getRange();
}
VectorDouble CovAnisoList::getRanges(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return 0.;
  return _covAnisos[icov]->getRanges();
}

void CovAnisoList::setSill(int icov, int ivar, int jvar, double value)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covs[icov]->setSill(ivar, jvar, value);
}
void CovAnisoList::setRangeIsotropic(int icov, double range)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covAnisos[icov]->setRangeIsotropic(range);
}
void CovAnisoList::setParam(int icov, double value)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covAnisos[icov]->setParam(value);
}
void CovAnisoList::setMarkovCoeffs(int icov, const VectorDouble& coeffs)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covAnisos[icov]->setMarkovCoeffs(coeffs);
}
void CovAnisoList::setType(int icov, const ECov& type)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covAnisos[icov]->setType(type);
}

int CovAnisoList::getGradParamNumber(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return 0;
  return _covAnisos[icov]->getGradParamNumber();
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
  for (int icov = 0, ncov = getCovaNumber(); icov < ncov; icov++)
  {
    const CovAniso* cova = getCova(icov);
    if (cova->getMinOrder() >= 0) return TEST;
    sill_total += cova->getSill(ivar, jvar);
  }
  return sill_total;
}

MatrixSquareSymmetric CovAnisoList::getTotalSill() const
{
  int nvar = getNVariables();
  MatrixSquareSymmetric mat(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
      mat.setValue(ivar,jvar,getTotalSill(ivar,jvar));
  return mat;
}

bool CovAnisoList::_isCovarianceIndexValid(int icov) const
{
  return checkArg("Covariance Index", icov, getCovaNumber());
}

/**
 * Returns the largest range (in any direction in the heterotopic case)
 * @return
 */
double CovAnisoList::getMaximumDistance() const

{
  double maxdist = 0.;
  for (int icov = 0, ncov = getCovaNumber(); icov < ncov; icov++)
  {
    const CovAniso* cova = getCova(icov);
    if (! cova->hasRange()) continue;
    double range = cova->getRange();
    if (range > maxdist) maxdist = range;
  }
  return maxdist;
}

void CovAnisoList::copyCovContext(const CovContext& ctxt)
{
  int number = (int) _covAnisos.size();
  for (int i = 0; i < number; i++)
    _covAnisos[i]->copyCovContext(ctxt);
}

void CovAnisoList::normalize(double sill, int ivar, int jvar)
{
  double covval = 0.;
  for (int i=0, n=getCovaNumber(); i<n; i++)
    covval += _covs[i]->eval0(ivar, jvar);

  if (covval <= 0. || isEqual(covval, sill)) return;
  double ratio = sill / covval;

  for (int i=0, n=getCovaNumber(); i<n; i++)
  {
    CovAniso* cov = _covAnisos[i];
    cov->setSill(cov->getSill(ivar, jvar) * ratio);
  }
}

bool CovAnisoList::hasNugget() const
{
  for (int is = 0, ns = getCovaNumber(); is < ns; is++)
  {
    if (getType(is) == ECov::NUGGET) return true;
  }
  return false;
}

int CovAnisoList::getRankNugget() const
{
  for (int is = 0, ns = getCovaNumber(); is < ns; is++)
  {
    if (getType(is) == ECov::NUGGET) return is;
  }
  return -1;
}

void CovAnisoList::_optimizationPreProcess(const std::vector<SpacePoint>& p) const
{
  for (const auto &e :_covs)
  {
    e->optimizationPreProcess(p);
  }
}

void CovAnisoList::_optimizationPostProcess() const
{
	for (int is = 0, ns = getCovaNumber(); is < ns; is++)
		_covs[is]->optimizationPostProcess();
}

const CovAnisoList* CovAnisoList::createReduce(const VectorInt &validVars) const
{
  CovAnisoList* newcovlist = this->clone();

  for (int is = 0, ns = getCovaNumber(); is < ns; is++)
  {
    CovAniso* covs = newcovlist->getCova(is);
    newcovlist->setCovAniso(is,covs->createReduce(validVars));
  }
  return newcovlist;
}




