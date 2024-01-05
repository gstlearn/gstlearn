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
  _noStat(),
  _matC()
{
}

ACovAnisoList::ACovAnisoList(const ACovAnisoList &r)
: ACov(r),
  _covs(),
  _filtered(r._filtered),
  _noStat(nullptr),
  _matC(r._matC)
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
    _matC = r._matC;
  }
  return *this;
}

ACovAnisoList::~ACovAnisoList()
{
  delAllCov();
  delete _noStat;
  _noStat = nullptr;
}

void ACovAnisoList::addCovList(const ACovAnisoList* covs)
{
  for (int icov = 0, ncov = covs->getCovNumber(); icov < ncov; icov++)
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
    for (int i=0, n=getCovNumber(); i<n; i++)
      cov += _covs[i]->eval0(ivar, jvar, mode);
  }
  else
  {
    for (int i=0, n=mode->getActiveCovList().size(); i<n; i++)
      cov += _covs[mode->getActiveCovList(i)]->eval0(ivar, jvar, mode);
  }
  return cov;
}

/**
 * Calculate the Matrix of covariance for zero distance
 * @param mat   Covariance matrix (Dimension: nvar * nvar)
 * @param mode  Calculation Options
 */
void ACovAnisoList::eval0MatInPlace(MatrixSquareGeneral &mat,
                                    const CovCalcMode *mode) const
{
  int nvar = mat.getNRows();
  if (_matC.getNRows() != nvar) _matC.reset(nvar,  nvar);

  mat.fill(0.);
  if (_considerAllCovariances(mode))
  {
    for (int i=0, n=getCovNumber(); i<n; i++)
    {
      _covs[i]->eval0MatInPlace(_matC, mode);
      mat.addMatrix(_matC);
    }
  }
  else
  {
    for (int i=0, n=mode->getActiveCovList().size(); i<n; i++)
    {
      _covs[mode->getActiveCovList(i)]->eval0MatInPlace(_matC, mode);
      mat.addMatrix(_matC);
    }
  }
}

/**
 * Evaluate the set of covariance vectors between samples of input 'db1' and
 * samples of output 'db2'
 * @param db1 Input Db
 * @param db2 Output db
 * @param ivar Rank of the first variable
 * @param jvar Rank of the second variable
 * @param mode CovCalcMode structure
 * @return
 */
VectorVectorDouble ACovAnisoList::evalCovMatrixOptim(const Db *db1,
                                                     const Db *db2,
                                                     int ivar,
                                                     int jvar,
                                                     const CovCalcMode *mode) const
{
  if (db2 == nullptr) db2 = db1;
  int nechtot2 = db2->getSampleNumber(false);
  int nech2 = db2->getSampleNumber(true);
  int nech1 = db1->getSampleNumber(true);
  VectorVectorDouble mat(nech2);

  for (auto &e : mat)
  {
    e = VectorDouble(nech1);
  }

  // Constitute the list of ALL samples contained in 'db1' (masked or active)
  std::vector<SpacePoint> p1s = db1->getSamplesAsSP();
  optimizationPreProcess(p1s);

  SpacePoint p2;

  int jech2 = 0;
  for (int iech2 = 0; iech2 < nechtot2; iech2++)
  {
    if (!db2->isActive(iech2)) continue;
    db2->getSampleCoordinatesAsSPInPlace(iech2, p2);
    optimizationSetTarget(p2);
    evalOptimInPlace(mat[jech2], ivar, jvar, mode);
    jech2++;
  }

  optimizationPostProcess();
  return mat;
}

void ACovAnisoList::evalOptimInPlace(VectorDouble &res,
                                     int ivar,
                                     int jvar,
                                     const CovCalcMode *mode) const
{
  for (auto &e : res)
    e = 0;
  for (int i = 0, n = getCovNumber(); i < n; i++)
    _covs[i]->evalOptimInPlace(res, ivar, jvar, mode);
}

void ACovAnisoList::evalMatOptimInPlace(int icas1,
                                        int iech1,
                                        int icas2,
                                        int iech2,
                                        MatrixSquareGeneral &mat,
                                        const CovCalcMode *mode) const
{
  int nvar = mat.getNRows();
  if (_matC.getNRows() != nvar) _matC.reset(nvar,  nvar);

  mat.fill(0.);
  if (_considerAllCovariances(mode))
  {
    for (unsigned int i=0, n=getCovNumber(); i<n; i++)
    {
      _covs[i]->evalMatOptimInPlace(icas1, iech1, icas2, iech2, _matC, mode);
      mat.addMatrix(_matC);
    }
  }
  else
  {
    for (int i=0, n=mode->getActiveCovList().size(); i<n; i++)
    {
      _covs[mode->getActiveCovList(i)]->evalMatOptimInPlace(icas1, iech1, icas2, iech2, _matC, mode);
      mat.addMatrix(_matC);
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
    for (unsigned int i=0, n=getCovNumber(); i<n; i++)
      cov += _covs[i]->eval(p1, p2, ivar, jvar, mode);
  }
  else
  {
    for (int i=0, n=mode->getActiveCovList().size(); i<n; i++)
      cov += _covs[mode->getActiveCovList(i)]->eval(p1, p2, ivar, jvar, mode);
  }
  return cov;
}

void ACovAnisoList::evalMatInPlace(const SpacePoint &p1,
                                   const SpacePoint &p2,
                                   MatrixSquareGeneral &mat,
                                   const CovCalcMode *mode) const
{
  int nvar = mat.getNRows();
  if (_matC.getNRows() != nvar) _matC.reset(nvar,  nvar);

  mat.fill(0.);
  if (_considerAllCovariances(mode))
  {
    for (unsigned int i=0, n=getCovNumber(); i<n; i++)
    {
      _covs[i]->evalMatInPlace(p1, p2, _matC, mode);
      mat.addMatrix(_matC);
    }
  }
  else
  {
    for (int i=0, n=mode->getActiveCovList().size(); i<n; i++)
    {
      _covs[mode->getActiveCovList(i)]->evalMatInPlace(p1, p2, _matC, mode);
      mat.addMatrix(_matC);
    }
  }
}

String ACovAnisoList::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (getCovNumber() <= 0) return sstr.str();

  for (int icov = 0, ncov = getCovNumber(); icov < ncov; icov++)
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

VectorInt ACovAnisoList::getActiveCovList() const
{
  VectorInt actives;
  for (unsigned int i=0, n=getCovNumber(); i<n; i++)
  {
    if (_filtered[i]) continue;
    actives.push_back(i);
  }
  return actives;
}

bool ACovAnisoList::isAllActiveCovList() const
{
  for (unsigned int i=0, n=getCovNumber(); i<n; i++)
  {
    if (_filtered[i]) return false;
  }
  return true;
}

VectorInt ACovAnisoList::getAllActiveCovList() const
{
  VectorInt actives;
  for (unsigned int i=0, n=getCovNumber(); i<n; i++)
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
  for (int icov = 0, ncov = getCovNumber(); icov < ncov; icov++)
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
  for (int icov = 0, ncov = getCovNumber(); icov < ncov; icov++)
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
  for (int is = 0, ns = getCovNumber(); is < ns; is++)
  {
    if (getType(is) == ECov::NUGGET) return true;
  }
  return false;
}

void ACovAnisoList::optimizationPreProcess(const std::vector<SpacePoint>& vec) const
{
	for (int is = 0, ns = getCovNumber(); is < ns; is++)
		_covs[is]->optimizationPreProcess(vec);
}

void ACovAnisoList::optimizationSetTarget(const SpacePoint& pt) const
{
  for (int is = 0, ns = getCovNumber(); is < ns; is++)
    _covs[is]->optimizationSetTarget(pt);
}

void ACovAnisoList::optimizationPostProcess() const
{
	for (int is = 0, ns = getCovNumber(); is < ns; is++)
		_covs[is]->optimizationPostProcess();
}

const ACovAnisoList* ACovAnisoList::reduce(const VectorInt &validVars) const
{
  ACovAnisoList* newcovlist = this->clone();

  for (int is = 0, ns = getCovNumber(); is < ns; is++)
  {
    CovAniso* covs = newcovlist->getCova(is);
    newcovlist->setCova(is,covs->reduce(validVars));
  }
  return newcovlist;
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

    if (icov < 0 || icov >= getCovNumber())
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
