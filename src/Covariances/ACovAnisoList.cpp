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

#include "Space/ASpace.hpp"
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Db/Db.hpp"

#include <math.h>
#include <vector>

ACovAnisoList::ACovAnisoList(const ASpace* space)
: ACov(space),
  _covs(),
  _filtered(),
  _noStat()
{
}

ACovAnisoList::ACovAnisoList(const ACovAnisoList &r)
: ACov(r),
  _covs(),
  _filtered(r._filtered),
  _noStat(nullptr)
{
  for (auto e: r._covs)
    _covs.push_back(e->clone());
  if (r._noStat != nullptr)
    _noStat = dynamic_cast<ANoStat*>(r._noStat->clone());
}

ACovAnisoList& ACovAnisoList::operator=(const ACovAnisoList &r)
{
  if (this != &r)
  {
    ACov::operator=(r);
    for (auto e: r._covs)
      _covs.push_back(e->clone());
    _filtered = r._filtered;
    if (r._noStat != nullptr)
      _noStat = dynamic_cast<ANoStat*>(r._noStat->clone());
  }
  return *this;
}

ACovAnisoList::~ACovAnisoList()
{
  delAllCov();
  delNoStat();
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
  for (auto e: _covs)
  {
    delete e;
  }
  _covs.clear();
  _filtered.clear();
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
void ACovAnisoList::eval0MatInPlace(MatrixSquareGeneral &mat,
                                    const CovCalcMode *mode) const
{
  if (_considerAllCovariances(mode))
  {
    for (int i=0, n=getCovaNumber(); i<n; i++)
    {
      _covs[i]->eval0MatInPlace(mat, mode);
    }
  }
  else
  {
    for (int i=0, n=(int) mode->getActiveCovList().size(); i<n; i++)
    {
      _covs[mode->getActiveCovList(i)]->eval0MatInPlace(mat, mode);
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
      db2->getSampleCoordinatesAsSPInPlace(iech2, p2);
      optimizationSetTarget(p2);

      // Loop on the basic structures
      for (int i = 0, n = getCovaNumber(); i < n; i++)
         _covs[i]->evalOptimInPlace(mat, ivars, index1, ivar2, icol, mode);
      icol++;
    }
  }

  optimizationPostProcess();
  return mat;
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
  int icolVerr = -1;
  double verr = 0.;

  VectorInt ivars = _getActiveVariables(ivar0);
  if (ivars.empty()) return mat;

  // Prepare the Optimization for covariance calculation
  optimizationPreProcess(db1);

  // Create the sets of Vector of valid sample indices per variable (not masked and defined)
  VectorVectorInt index1 = db1->getMultipleRanksActive(ivars, nbgh1);

  // Update the Covariance matrix in case of presence of Variance of Measurement error
  bool useVerr = db1->hasLocVariable(ELoc::V);

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
    if (useVerr) icolVerr = db1->getColIdxByLocator(ELoc::V, ivar2);

    // Loop on the second sample
    int nech2s = (int) index1[rvar2].size();
    for (int rech2 = 0; rech2 < nech2s; rech2++)
    {
      int iech2 = index1[rvar2][rech2];
      if (icolVerr >= 0) verr = db1->getValueByColIdx(iech2, icolVerr);
      optimizationSetTarget(iech2);

      // Loop on the basic structures
      for (int i = 0, n = getCovaNumber(); i < n; i++)
         _covs[i]->evalOptimInPlace(mat, ivars, index1, ivar2, icol, mode);

      // Update the Diagonal due to the presence of Variance of Measurement Error
      if (icolVerr >= 0)
      {
        // Loop on the first variable
        int irow = 0;
        for (int rvar1 = 0, nvar1 = (int) ivars.size(); rvar1 < nvar1; rvar1++)
        {
          // Loop on the first sample
          int nech1s = (int) index1[rvar1].size();
          for (int rech1 = 0; rech1 < nech1s; rech1++)
          {
            if (irow == icol)
              mat.updValue(irow, icol, EOperator::ADD, verr);
            irow++;
          }
        }
      }
      icol++;
    }
  }

  optimizationPostProcess();
  return mat;
}


/**
 * Calculate the Matrix of covariance between two elements of two Dbs (defined beforehand)
 * @param icas1 Origin of the Db containing the first point
 * @param iech1 Rank of the first point
 * @param icas2 Origin of the Db containing the second point
 * @param iech2 Rank of the second point
 * @param mat   Covariance matrix (Dimension: nvar * nvar)
 * @param mode  Calculation Options
 *
 * @remarks: Matrix 'mat' should be dimensioned and initialized beforehand
 */
void ACovAnisoList::evalMatOptimInPlace(int icas1,
                                        int iech1,
                                        int icas2,
                                        int iech2,
                                        MatrixSquareGeneral &mat,
                                        const CovCalcMode *mode) const
{
  if (_considerAllCovariances(mode))
  {
    for (int i=0, n=getCovaNumber(); i<n; i++)
    {
      _covs[i]->evalMatOptimInPlace(icas1, iech1, icas2, iech2, mat, mode);
    }
  }
  else
  {
    for (int i=0, n=(int) mode->getActiveCovList().size(); i<n; i++)
    {
      _covs[mode->getActiveCovList(i)]->evalMatOptimInPlace(icas1, iech1, icas2, iech2, mat, mode);
    }
  }
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

/**
 * Calculate the Matrix of covariance between two space points
 * @param p1 Reference of the first space point
 * @param p2 Reference of the second space point
 * @param mat   Covariance matrix (Dimension: nvar * nvar)
 * @param mode  Calculation Options
 *
 * @remarks: Matrix 'mat' should be dimensioned and initialized beforehand
 */
void ACovAnisoList::evalMatInPlace(const SpacePoint &p1,
                                   const SpacePoint &p2,
                                   MatrixSquareGeneral &mat,
                                   const CovCalcMode *mode) const
{
  if (_considerAllCovariances(mode))
  {
    for (int i=0, n=getCovaNumber(); i<n; i++)
    {
      _covs[i]->evalMatInPlace(p1, p2, mat, mode);
    }
  }
  else
  {
    for (int i=0, n=(int) mode->getActiveCovList().size(); i<n; i++)
    {
      _covs[mode->getActiveCovList(i)]->evalMatInPlace(p1, p2, mat, mode);
    }
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

  // Non-stationary parameters
  if (_noStat != nullptr)
  {
    sstr << _noStat->toString();
  }
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

MatrixSquareGeneral ACovAnisoList::getTotalSill() const
{
  int nvar = getNVariables();
  MatrixSquareGeneral mat(nvar);
  for (int ivar = 0; ivar < nvar; ivar++)
    for (int jvar = 0; jvar < nvar; jvar++)
      mat.setValue(ivar,jvar,getTotalSill(ivar,jvar));
  return mat;
}

bool ACovAnisoList::_isCovarianceIndexValid(int icov) const
{
  if (icov >= (int) getCovaNumber())
  {
    mesArg("Covariance Index",icov,getCovaNumber());
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

bool ACovAnisoList::isOptimizationInitialized(const Db* db) const
{
  for (int is = 0, ns = getCovaNumber(); is < ns; is++)
    if (! _covs[is]->isOptimizationInitialized(db)) return false;
  return true;
}

void ACovAnisoList::optimizationPreProcess(const Db* db) const
{
	for (int is = 0, ns = getCovaNumber(); is < ns; is++)
		_covs[is]->optimizationPreProcess(db);
}

void ACovAnisoList::optimizationSetTarget(const SpacePoint& pt) const
{
  for (int is = 0, ns = getCovaNumber(); is < ns; is++)
    _covs[is]->optimizationSetTarget(pt);
}

void ACovAnisoList::optimizationSetTarget(int iech) const
{
  for (int is = 0, ns = getCovaNumber(); is < ns; is++)
    _covs[is]->optimizationSetTarget(iech);
}

void ACovAnisoList::optimizationPostProcess() const
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

void ACovAnisoList::delNoStat()
{
  delete _noStat;
  _noStat = nullptr;
}

/**
 * Define Non-stationary parameters
 * @param anostat ANoStat pointer will be duplicated
 * @return Error return code
 */
int ACovAnisoList::addNoStat(const ANoStat *anostat)
{
  if (anostat == nullptr) return 0;
  if (getNDim() > 3)
  {
    messerr("Non stationary model is restricted to Space Dimension <= 3");
    return 1;
  }

  for (int ipar = 0; ipar < (int) getNoStatElemNumber(); ipar++)
  {
    int icov = getNoStatElemIcov(ipar);
    EConsElem type = getNoStatElemType(ipar);

    // Check that the Non-stationary parameter is valid with respect
    // to the Model definition

    if (icov < 0 || icov >= getCovaNumber())
    {
      messerr("Invalid Covariance rank (%d) for the Non-Stationary Parameter (%d)",
              icov, ipar);
      return 1;
    }
    if (type == EConsElem::PARAM)
    {
      messerr("The current methodology does not handle constraint on third parameter");
      return 1;
    }
  }

  if (_noStat != nullptr) delete _noStat;
  _noStat = dynamic_cast<ANoStat*>(anostat->clone());
  return 0;
}

int ACovAnisoList::getNoStatElemNumber() const
{
  if (_noStat == nullptr) return 0;
  return _noStat->getNoStatElemNumber();
}

int ACovAnisoList::getNoStatElemIcov(int ipar) const
{
  if (_noStat == nullptr) return ITEST;
  return _noStat->getICov(ipar);
}
const EConsElem& ACovAnisoList::getNoStatElemType(int ipar) const
{
  if (_noStat == nullptr) return EConsElem::UNKNOWN;
  return _noStat->getType(ipar);
}

int ACovAnisoList::addNoStatElem(int igrf,
                                 int icov,
                                 const EConsElem &type,
                                 int iv1,
                                 int iv2)
{
  if (_noStat == nullptr) return 0;
  return _noStat->addNoStatElem(igrf, icov, type, iv1, iv2);
}

int ACovAnisoList::addNoStatElems(const VectorString &codes)
{
  if (_noStat == nullptr) return 0;
  return _noStat->addNoStatElems(codes);
}

CovParamId ACovAnisoList::getCovParamId(int ipar) const
{
  if (_noStat == nullptr) return CovParamId();
  return _noStat->getItems(ipar);
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
  double val1, val2;

  // If no non-stationary parameter is defined, simply skip
  if (! isNoStat()) return;
  int ndim = getNDim();

  // Loop on the elements that can be updated one-by-one

  for (int ipar = 0, npar = _noStat->getNoStatElemNumber(); ipar < npar; ipar++)
  {
    int icov = _noStat->getICov(ipar);
    EConsElem type = _noStat->getType(ipar);

    if (type == EConsElem::SILL)
    {
      if (_noStat->getInfoFromDb(ipar, icas1, iech1, icas2, iech2, &val1, &val2))
      {
        int iv1  = _noStat->getIV1(ipar);
        int iv2  = _noStat->getIV2(ipar);
        setSill(icov, iv1, iv2, sqrt(val1 * val2));
      }
    }
    else if (type == EConsElem::PARAM)

    {
      if (_noStat->getInfoFromDb(ipar, icas1, iech1, icas2, iech2, &val1, &val2))
        setParam(icov, 0.5 * (val1 + val2));
    }
  }

  // Loop on the other parameters (Anisotropy) that must be processed globally

  for (int icov = 0, ncov = getCovaNumber(); icov < ncov; icov++)
  {
    if (! _noStat->isDefinedforAnisotropy(icov)) continue;
    CovAniso* cova = getCova(icov);

    VectorDouble angle1;
    VectorDouble angle2;

    VectorDouble scale1;
    VectorDouble scale2;

    VectorDouble range1;
    VectorDouble range2;

    // Define the angles (for all space dimensions)
    bool flagRotTwo = false;
    bool flagRotOne = false;
    if (_noStat->isDefined(EConsElem::ANGLE, icov))
    {
      angle1 = cova->getAnisoAngles();
      angle2 = angle1;
      for (int idim = 0; idim < ndim; idim++)
      {
        if (_noStat->isDefined(EConsElem::ANGLE, icov, idim, 0))
        {
          int ipar = _noStat->getRank(EConsElem::ANGLE, icov, idim);
          if (ipar < 0) continue;
          flagRotOne = true;
          if (_noStat->getInfoFromDb(ipar, icas1, iech1, icas2, iech2, &angle1[idim], &angle2[idim]))
            flagRotTwo = true;
        }
      }
    }

    // Define the Theoretical ranges (for all space dimensions)

    bool flagScaleTwo = false;
    bool flagScaleOne = false;
    if (_noStat->isDefined(EConsElem::SCALE, icov))
    {
      scale1 = cova->getScales();
      scale2 = scale1;
      for (int idim = 0; idim < ndim; idim++)
      {
        if (_noStat->isDefined(EConsElem::SCALE, icov, idim, 0))
        {
          int ipar = _noStat->getRank(EConsElem::SCALE, icov, idim);
          if (ipar < 0) continue;
          flagScaleOne = true;
          if (_noStat->getInfoFromDb(ipar, icas1, iech1, icas2, iech2, &scale1[idim], &scale2[idim]))
            flagScaleTwo = true;
        }
      }
    }

    // Define the Practical ranges (for all space dimensions)

    bool flagRangeTwo = false;
    bool flagRangeOne = false;
    if (_noStat->isDefined(EConsElem::RANGE, icov))
    {
      range1 = cova->getRanges();
      range2 = range1;
      for (int idim = 0; idim < ndim; idim++)
      {
        if (_noStat->isDefined(EConsElem::RANGE, icov, idim))
        {
          int ipar = _noStat->getRank(EConsElem::RANGE, icov, idim);
          if (ipar < 0) continue;
          flagRangeOne = true;
          if (_noStat->getInfoFromDb(ipar, icas1, iech1, icas2, iech2, &range1[idim], &range2[idim]))
            flagRangeTwo = true;
        }
      }
    }

    // Update the Model
    double ratio = 1.;
    if (flagRotTwo || flagRangeTwo || flagScaleTwo)
    {
      // Extract the direct tensor at first point and square it
      cova->setRotationAnglesAndRadius(angle1, range1, scale1);
      MatrixSquareSymmetric direct1 = cova->getAniso().getTensorDirect2();
      double det1 = pow(direct1.determinant(), 0.25);

      // Extract the direct tensor at second point and square it
      cova->setRotationAnglesAndRadius(angle2, range2, scale2);
      MatrixSquareSymmetric direct2 = cova->getAniso().getTensorDirect2();
      double det2 = pow(direct2.determinant(), 0.25);

      // Calculate average squared tensor
      direct2.addMatInPlace(direct1, 0.5, 0.5);
      double detM = sqrt(direct2.determinant());

      // Update the tensor (squared version)
      Tensor tensor = cova->getAniso();
      tensor.setTensorDirect2(direct2);
      cova->setAniso(tensor);
      ratio = det1 * det2 / detM;
    }
    else if (flagRotOne || flagRangeOne || flagScaleOne)
    {
      // Simply update the model with one set of parameters
      cova->setRotationAnglesAndRadius(angle1, range1, scale1);
    }
    cova->setNoStatFactor(ratio);
  }
}

/**
 * Update the Model according to the Non-stationary parameters
 * @param imesh Rank of the target mesh
 */
void ACovAnisoList::updateCovByMesh(int imesh)
{
  // If no non-stationary parameter is defined, simply skip
  if (! isNoStat()) return;
  int ndim = getNDim();

  // Loop on the elements that can be updated one-by-one

  for (int ipar = 0; ipar < _noStat->getNoStatElemNumber(); ipar++)
  {
    int icov = _noStat->getICov(ipar);
    EConsElem type = _noStat->getType(ipar);

    if (type == EConsElem::SILL)
    {
      double sill = _noStat->getValueByParam(ipar, 0, imesh);
      int iv1  = _noStat->getIV1(ipar);
      int iv2  = _noStat->getIV2(ipar);
      setSill(icov, iv1, iv2, sill);
    }
  }

  // Loop on the other parameters (Anisotropy) that must be processed globally

  for (int icov = 0; icov < getCovaNumber(); icov++)
  {
    if (! _noStat->isDefinedforAnisotropy(icov)) continue;
    CovAniso* cova = getCova(icov);

    VectorDouble angles(cova->getAnisoAngles());
    VectorDouble scales(cova->getScales());
    VectorDouble ranges(cova->getRanges());

    // Define the angles (for all space dimensions)
    if (_noStat->isDefined(EConsElem::ANGLE, icov))
    {
      for (int idim = 0; idim < ndim; idim++)
      {
        if (_noStat->isDefined(EConsElem::ANGLE, icov, idim, 0))
        {
          int ipar = _noStat->getRank(EConsElem::ANGLE, icov, idim);
          if (ipar < 0) continue;
          angles[idim] = _noStat->getValueByParam(ipar, 0, imesh);
        }
      }
    }

    // Define the Theoretical ranges (for all space dimensions)
    if (_noStat->isDefined(EConsElem::SCALE, icov))
    {
      for (int idim = 0; idim < ndim; idim++)
      {
        if (_noStat->isDefined(EConsElem::SCALE, icov, idim))
        {
          int ipar = _noStat->getRank(EConsElem::SCALE, icov, idim);
          if (ipar < 0) continue;
          scales[idim] = _noStat->getValueByParam(ipar, 0, imesh);
        }
      }
    }

    // Define the Practical ranges (for all space dimensions)
    if (_noStat->isDefined(EConsElem::RANGE, icov))
    {
      for (int idim = 0; idim < ndim; idim++)
      {
        if (_noStat->isDefined(EConsElem::RANGE, icov, idim))
        {
          int ipar = _noStat->getRank(EConsElem::RANGE, icov, idim);
          if (ipar < 0) continue;
          ranges[idim] = _noStat->getValueByParam(ipar, 0, imesh);
        }
      }
    }

    // Exploit the Anisotropy
    cova->setRotationAnglesAndRadius(angles, ranges, scales);
  }
}
