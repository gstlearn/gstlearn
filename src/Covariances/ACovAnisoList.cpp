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
#include "Covariances/ACovAnisoList.hpp"

#include "Covariances/CovCalcMode.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Space/ASpace.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Db/Db.hpp"
#include "Space/SpacePoint.hpp"

#include <algorithm>
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
  for (auto* e: r._covs)
    _covs.push_back(e->clone());
}

ACovAnisoList& ACovAnisoList::operator=(const ACovAnisoList &r)
{
  if (this != &r)
  {
    ACov::operator=(r);
    for (auto *e: r._covs)
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
  for (int icov = 0, ncov = covs->getCovaNumber(); icov < ncov; icov++)
    addCov(covs->getCova(icov));
}

void ACovAnisoList::addCov(const CovAniso* cov)
{
  if (getCovaNumber() > 0)
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

void ACovAnisoList::delCov(int icov)
{
  if (! _isCovarianceIndexValid(icov)) return;
  delete _covs[icov];
  _covs.erase(_covs.begin() + icov);
  _filtered.erase(_filtered.begin() + icov);
}

void ACovAnisoList::delAllCov()
{
  for (auto &e: _covs)
  {
    delete e;
  }
  _covs.clear();
  _filtered.clear();
}

bool ACovAnisoList::isNoStat() const
{
  bool nostat = false;
  for (const auto &e :_covs)
  {
    nostat = nostat || e->isNoStat();
  }
  return nostat;
}
void ACovAnisoList::setFiltered(int icov, bool filtered)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _filtered[icov] = filtered;
}

bool ACovAnisoList::isConsistent(const ASpace* /*space*/) const
{
  /// TODO : ACovAnisoList::isConsistent
  return true;
}

int ACovAnisoList::getNVariables() const
{
  if (getCovaNumber() > 0)
    return _covs[0]->getNVariables();
  return 0;
}

bool ACovAnisoList::_considerAllCovariances(const CovCalcMode* mode) const
{
  if (mode == nullptr) return true;
  if (mode->isAllActiveCov()) return true;
  return false;
}

double ACovAnisoList::eval0(int ivar, int jvar, const CovCalcMode* mode) const
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
void ACovAnisoList::addEval0CovMatBiPointInPlace(MatrixSquareSymmetric &mat,
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
 * Evaluate the covariance rectangular matrix between samples of input 'db1' and 'db2'
 * @param db1 Input Db
 * @param db2 Output db
 * @param ivar0 Rank of the first variable (-1 for all variables)
 * @param jvar0 Rank of the second variable (-1 for all variables)
 * @param nbgh1 Vector of indices of active samples in db1 (optional)
 * @param nbgh2 Vector of indices of active samples in db2 (optional)
 * @param mode CovCalcMode structure
 * @return
 */
MatrixRectangular ACovAnisoList::evalCovMatrixOptim(const Db *db1,
                                                    const Db *db2,
                                                    int ivar0,
                                                    int jvar0,
                                                    const VectorInt& nbgh1,
                                                    const VectorInt& nbgh2,
                                                    const CovCalcMode *mode) const
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
      p2.setIech(iech2);
      db2->getSampleAsSPInPlace(p2);
      optimizationSetTarget(p2);

      // Loop on the basic structures
      for (int i = 0, n = getCovaNumber(); i < n; i++)
         _covs[i]->evalOptimInPlace(mat, ivars, index1, ivar2, icol, mode, false);
      icol++;
    }
  }

  optimizationPostProcess();
  return mat;
}

void ACovAnisoList::_optimizationSetTarget(const SpacePoint& pt) const
{
  for (int is = 0, ns = getCovaNumber(); is < ns; is++)
    _covs[is]->optimizationSetTarget(pt);
}

void ACovAnisoList::optimizationSetTargetByIndex(int iech) const
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
 * @return
 */
MatrixSquareSymmetric ACovAnisoList::evalCovMatrixSymmetricOptim(const Db *db1,
                                                                 int ivar0,
                                                                 const VectorInt &nbgh1,
                                                                 const CovCalcMode *mode) const
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
         _covs[i]->evalOptimInPlace(mat, ivars, index1, ivar2, icol, mode, true);

      icol++;
    }
  }

  // Update the matrix due to presence of Variance of Measurement Error
  _updateCovMatrixSymmetricVerr(db1, &mat, ivars, index1);

  optimizationPostProcess();
  return mat;
}

double ACovAnisoList::eval(const SpacePoint& p1,
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

double ACovAnisoList::_loadAndEval(const SpacePoint& p1,
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
void ACovAnisoList::_addEvalCovMatBiPointInPlace(MatrixSquareSymmetric &mat,
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

void ACovAnisoList::_loadAndAddEvalCovMatBiPointInPlace(MatrixSquareSymmetric &mat,const SpacePoint& p1,const SpacePoint&p2,
                                              const CovCalcMode *mode) const
{
  for (const auto &e : _covs)
  {
    e->loadAndAddEvalCovMatBiPointInPlace(mat,p1,p2,mode);
  }
}

String ACovAnisoList::toString(const AStringFormat* /*strfmt*/) const
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

int  ACovAnisoList::getCovaNumber(bool skipNugget) const
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

bool ACovAnisoList::isFiltered(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return false;
  return _filtered[icov];
}

bool ACovAnisoList::hasRange() const
{
  for (int i=0, n=getCovaNumber(); i<n; i++)
  {
    if (!getCova(i)->hasRange())
      return false;
  }
  return true;
}

bool ACovAnisoList::isStationary() const
{
  for (int i=0, n=getCovaNumber(); i<n; i++)
  {
    if (getCova(i)->getMinOrder() >= 0)
      return false;
  }
  return true;
}

VectorInt ACovAnisoList::getActiveCovList() const
{
  VectorInt actives;
  for (int i=0, n=getCovaNumber(); i<n; i++)
  {
    if (_filtered[i]) continue;
    actives.push_back(i);
  }
  return actives;
}

bool ACovAnisoList::isAllActiveCovList() const
{
  for (int i=0, n=getCovaNumber(); i<n; i++)
  {
    if (_filtered[i]) return false;
  }
  return true;
}

VectorInt ACovAnisoList::getAllActiveCovList() const
{
  VectorInt actives;
  for (int i=0, n=getCovaNumber(); i<n; i++)
  {
    actives.push_back(i);
  }
  return actives;
}

CovAniso ACovAnisoList::extractCova(int icov) const
{
  return *(_covs[icov]);
}

/**
 * @return The Minimum IRF-order induced by the covariances
 */
int ACovAnisoList::getCovaMinIRFOrder() const
{
  int nmini = -1;
  for (unsigned i = 0, n = getCovaNumber(); i<n; i++)
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
double ACovAnisoList::getParam(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return 0.;
  return _covs[icov]->getParam();
}
double ACovAnisoList::getRange(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return 0.;
  return _covs[icov]->getRange();
}
VectorDouble ACovAnisoList::getRanges(int icov) const
{
  if (! _isCovarianceIndexValid(icov)) return 0.;
  return _covs[icov]->getRanges();
}
const MatrixSquareSymmetric& ACovAnisoList::getSill(int icov) const
{
  return _covs[icov]->getSill();
}
double ACovAnisoList::getSill(int icov, int ivar, int jvar) const
{
  if(! _isCovarianceIndexValid(icov)) return 0.;
  return _covs[icov]->getSill(ivar, jvar);
}
void ACovAnisoList::setSill(int icov, int ivar, int jvar, double value)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covs[icov]->setSill(ivar, jvar, value);
}
void ACovAnisoList::setRangeIsotropic(int icov, double range)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covs[icov]->setRangeIsotropic(range);
}
void ACovAnisoList::setParam(int icov, double value)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covs[icov]->setParam(value);
}
void ACovAnisoList::setMarkovCoeffs(int icov, VectorDouble coeffs)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covs[icov]->setMarkovCoeffs(coeffs);
}
void ACovAnisoList::setType(int icov, const ECov& type)
{
  if (! _isCovarianceIndexValid(icov)) return;
  _covs[icov]->setType(type);
}

int ACovAnisoList::getGradParamNumber(int icov) const
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
  for (int icov = 0, ncov = getCovaNumber(); icov < ncov; icov++)
  {
    const CovAniso* cova = getCova(icov);
    if (cova->getMinOrder() >= 0) return TEST;
    sill_total += cova->getSill(ivar, jvar);
  }
  return sill_total;
}

MatrixSquareSymmetric ACovAnisoList::getTotalSill() const
{
  int nvar = getNVariables();
  MatrixSquareSymmetric mat(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar <= ivar; jvar++)
      mat.setValue(ivar,jvar,getTotalSill(ivar,jvar));
  return mat;
}

bool ACovAnisoList::_isCovarianceIndexValid(int icov) const
{
  return checkArg("Covariance Index", icov, getCovaNumber());
}

/**
 * Returns the largest range (in any direction in the heterotopic case)
 * @return
 */
double ACovAnisoList::getMaximumDistance() const

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

void ACovAnisoList::copyCovContext(const CovContext& ctxt)
{
  int number = (int) _covs.size();
  for (int i = 0; i < number; i++)
    _covs[i]->copyCovContext(ctxt);
}

void ACovAnisoList::normalize(double sill, int ivar, int jvar)
{
  double covval = 0.;
  for (int i=0, n=getCovaNumber(); i<n; i++)
    covval += _covs[i]->eval0(ivar, jvar);

  if (covval <= 0. || areEqual(covval, sill)) return;
  double ratio = sill / covval;

  for (int i=0, n=getCovaNumber(); i<n; i++)
  {
    CovAniso* cov = _covs[i];
    cov->setSill(cov->getSill(ivar, jvar) * ratio);
  }
}

bool ACovAnisoList::hasNugget() const
{
  for (int is = 0, ns = getCovaNumber(); is < ns; is++)
  {
    if (getType(is) == ECov::NUGGET) return true;
  }
  return false;
}

int ACovAnisoList::getRankNugget() const
{
  for (int is = 0, ns = getCovaNumber(); is < ns; is++)
  {
    if (getType(is) == ECov::NUGGET) return is;
  }
  return -1;
}

// void ACovAnisoList::evalCovLHS(MatrixSquareSymmetric &mat,
//                           SpacePoint &pwork1,
//                           SpacePoint &pwork2,
//                           const Db* db, 
//                           const CovCalcMode *mode) const
// {
//   for (const auto &e : _covs)
//   {
//     e->evalCovLHS(mat, pwork1, pwork2, db, mode);
//   }

// }

// void ACovAnisoList::evalCovRHS(MatrixSquareSymmetric &mat,
//                           SpacePoint &pwork1,
//                           const Db* db,  SpacePoint& pout,  
//                           const CovCalcMode *mode) const
// {
//   for (const auto &e : _covs)
//   {
//     e->evalCovRHS(mat, pwork1, db, pout, mode);
//   }
  
// }

void ACovAnisoList::_optimizationPreProcess(const std::vector<SpacePoint>& p) const
{
  for (const auto &e :_covs)
  {
    e->optimizationPreProcess(p);
  }
}


void ACovAnisoList::_optimizationPostProcess() const
{
	for (int is = 0, ns = getCovaNumber(); is < ns; is++)
		_covs[is]->optimizationPostProcess();
}

const ACovAnisoList* ACovAnisoList::createReduce(const VectorInt &validVars) const
{
  ACovAnisoList* newcovlist = this->clone();

  for (int is = 0, ns = getCovaNumber(); is < ns; is++)
  {
    CovAniso* covs = newcovlist->getCova(is);
    newcovlist->setCova(is,covs->createReduce(validVars));
  }
  return newcovlist;
}

void ACovAnisoList::_manage(const Db* db1,const Db* db2)  const
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

void ACovAnisoList::updateCovByPoints(int icas1, int iech1, int icas2, int iech2) 
{
  for (const auto &e : _covs)
  {
    e->updateCovByPoints(icas1,iech1,icas2,iech2);
  }
}

