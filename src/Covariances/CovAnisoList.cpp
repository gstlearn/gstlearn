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

#include "Basic/ListParams.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/CovList.hpp"
#include "Space/ASpace.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovLMGradient.hpp"
#include "Covariances/CovLMCConvolution.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Db/Db.hpp"
#include "Anamorphosis/AnamHermite.hpp"

#include <math.h>
#include <vector>

CovAnisoList::CovAnisoList(const CovContext& ctxt)
  : CovList(ctxt)
{
  setOptimEnabled(true);
}

CovAnisoList::CovAnisoList(const CovAnisoList& r)
  : CovList(r)
{
  
}

CovAnisoList& CovAnisoList::operator=(const CovAnisoList &r)
{

  

  if (this != &r)
  {
    CovList::operator=(r);
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
    addCov(covs->getCovAniso(icov));
}

void CovAnisoList::addCov(const CovBase* cov) 
{
  if (dynamic_cast<const CovAniso*>(cov) == nullptr)
  {
    messerr("The argument should be of type 'CovAniso*'");
    return;
  }
  CovList::addCov(cov);
}

const CovAniso* CovAnisoList::_getCovAniso(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return nullptr;
  const CovAniso* covaniso = dynamic_cast<const CovAniso*>(_covs[icov]);
  if (covaniso == nullptr)
  {
    messerr("The element 'icov' is not a CovAniso");
  }
  return covaniso;
}

CovAniso* CovAnisoList::_getCovAnisoModify(int icov)
{
  if (!_isCovarianceIndexValid(icov)) return nullptr;
  CovAniso* covaniso = dynamic_cast<CovAniso*>(_covs[icov]);
  if (covaniso == nullptr)
  {
    messerr("The element 'icov' is not a CovAniso");
  }
  return covaniso;
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
  for (const auto& j: list.getVector())
  {
    cov += _covs[j]->eval0(ivar, jvar, mode);
  }
  return cov;
}

String CovAnisoList::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  if (getNCov() <= 0) return sstr.str();

  for (int icov = 0, ncov = getNCov(); icov < ncov; icov++)
  {
    sstr << getCovAniso(icov)->toString();
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
    if (getCovAniso(icov)->getType() != ECov::NUGGET) nstruc++;
  }
  return nstruc;
}

bool CovAnisoList::hasRange() const
{
  for (int i = 0, n = getNCov(); i < n; i++)
  {
    if (!getCovAniso(i)->hasRange()) return false;
  }
  return true;
}

bool CovAnisoList::isStationary() const
{
  for (int i = 0, n = getNCov(); i < n; i++)
  {
    if (getCovAniso(i)->getMinOrder() >= 0) return false;
  }
  return true;
}

CovAniso CovAnisoList::extractCova(int icov) const
{
  const CovAniso* covaniso = _getCovAniso(icov);
  if (covaniso == nullptr)
    return CovAniso(ECov::NUGGET, CovContext());
  return *_getCovAniso(icov);
}

/**
 * @return The Minimum IRF-order induced by the covariances
 */
int CovAnisoList::getCovMinIRFOrder() const
{
  int nmini = -1;
  for (unsigned i = 0, n = getNCov(); i < n; i++)
  {
    const CovAniso* covaniso = _getCovAniso(i);
    if (covaniso == nullptr) continue; 
    int locmini = covaniso->getMinOrder();
    if (locmini > nmini) nmini = locmini;
  }
  return nmini;
}

CovAniso* CovAnisoList::getCovAniso(int icov) 
{
  if (!_isCovarianceIndexValid(icov)) return nullptr;
  return _getCovAnisoModify(icov);
}
const CovAniso* CovAnisoList::getCovAniso(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return nullptr;
  return _getCovAniso(icov);
}
void CovAnisoList::setCov(int icov, const CovBase* covs) 
// TODO rename into setOneCov
// to be different from the one in ModelGeneric
{
  if (dynamic_cast<const CovAniso*>(covs) == nullptr)
  {
    messerr("The argument should be of type 'CovAniso*'");
    return;
  }
  CovList::setCov(icov, covs);
}
const ECov& CovAnisoList::getCovType(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return ECov::UNKNOWN;
  const CovAniso* covaniso = _getCovAniso(icov);
  if (covaniso == nullptr)
  {
    messerr("The argument should be of type 'CovAniso*'");
    return CovList::getCovType(icov);
  }
  return covaniso->getType();
}

String CovAnisoList::getCovName(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return String();
  const CovAniso* covaniso = _getCovAniso(icov);
  if (covaniso == nullptr)
  {
    messerr("The argument should be of type 'CovAniso*'");
    return CovList::getCovName(icov);
  }
  return covaniso->getCovName();
}
double CovAnisoList::getParam(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return 0.;
  const CovAniso* covaniso = _getCovAniso(icov);
  if (covaniso == nullptr)
  {
    messerr("The argument should be of type 'CovAniso*'");
    return 1.;
  }
  return covaniso->getParam();
}
double CovAnisoList::getRange(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return 0.;
  const CovAniso* covaniso = _getCovAniso(icov);
  if (covaniso == nullptr)
  {
    messerr("The argument should be of type 'CovAniso*'");
    return 0.;
  }
  return covaniso->getRange();
}
VectorDouble CovAnisoList::getRanges(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return VectorDouble();
  const CovAniso* covaniso = _getCovAniso(icov);
  if (covaniso == nullptr)
  {
    messerr("The argument should be of type 'CovAniso*'");
    return VectorDouble();
  }
  return covaniso->getRanges();
}
VectorDouble CovAnisoList::getAngles(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return VectorDouble();
  const CovAniso* covaniso = _getCovAniso(icov);
  if (covaniso == nullptr)
  {
    messerr("The argument should be of type 'CovAniso*'");
    return VectorDouble();
  }
  return covaniso->getAnisoAngles();
}
void CovAnisoList::setRangeIsotropic(int icov, double range)
{
  if (!_isCovarianceIndexValid(icov)) return;
  CovAniso* covaniso = _getCovAnisoModify(icov);
  if (covaniso == nullptr)
  {
    messerr("The argument should be of type 'CovAniso*'");
    return;
  }
  covaniso->setRangeIsotropic(range);
}
void CovAnisoList::setParam(int icov, double value)
{
  if (!_isCovarianceIndexValid(icov)) return;
  CovAniso* covaniso = _getCovAnisoModify(icov);
  if (covaniso == nullptr)
  {
    messerr("The argument should be of type 'CovAniso*'");
    return;
  }
  covaniso->setParam(value);
}
void CovAnisoList::setMarkovCoeffs(int icov, const VectorDouble& coeffs)
{
  if (!_isCovarianceIndexValid(icov)) return;
  CovAniso* covaniso = _getCovAnisoModify(icov);
  if (covaniso == nullptr)
  {
    messerr("The argument should be of type 'CovAniso*'");
    return;
  }
  covaniso->setMarkovCoeffs(coeffs);
}
void CovAnisoList::setType(int icov, const ECov& type)
{
  if (!_isCovarianceIndexValid(icov)) return;
  CovAniso* covaniso = _getCovAnisoModify(icov);
  if (covaniso == nullptr)
  {
    messerr("The argument should be of type 'CovAniso*'");
    return;
  }
  covaniso->setType(type);
}

int CovAnisoList::getNGradParam(int icov) const
{
  if (!_isCovarianceIndexValid(icov)) return 0;
  const CovAniso* covaniso = _getCovAniso(icov);
  if (covaniso == nullptr)
  {
    messerr("The argument should be of type 'CovAniso*'");
    return 0;
  }
  return covaniso->getNGradParam();
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
    const CovAniso* cova = getCovAniso(icov);
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
    const CovAniso* cova = getCovAniso(icov);
    if (!cova->hasRange()) continue;
    double range = cova->getRange();
    if (range > maxdist) maxdist = range;
  }
  return maxdist;
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
    CovAniso* covs = newcovlist->getCovAniso(is);
    newcovlist->setCov(is, covs->createReduce(validVars));
  }
  newcovlist->setNVar((int) validVars.size());
  return newcovlist;
}

/**
 * Returns the Ball radius (from the first covariance of _covaList)
 * @return Value of the Ball Radius (if defined, i.e. for Numerical Gradient calculation)
 */
double CovAnisoList::getBallRadius() const
{
  // Check is performed on the first covariance
  const CovAniso* cova = getCovAniso(0);
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


void CovAnisoList::makeRangeNoStatDb(int icov, const String& namecol, int idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeRangeNoStatDb(namecol, idim);
}

void CovAnisoList::makeScaleNoStatDb(int icov, const String& namecol, int idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeScaleNoStatDb(namecol, idim, nullptr);

}
void CovAnisoList::makeAngleNoStatDb(int icov, const String& namecol, int idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeAngleNoStatDb(namecol, idim);
}

void CovAnisoList::makeTensorNoStatDb(int icov, const String& namecol, int idim, int jdim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeTensorNoStatDb(namecol, idim, jdim);

}
void CovAnisoList::makeParamNoStatDb(int icov, const String& namecol)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeParamNoStatDb(namecol);
}
void CovAnisoList::makeRangeNoStatFunctional(int icov, const AFunctional* func, int idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeRangeNoStatFunctional(func, idim);
}
void CovAnisoList::makeScaleNoStatFunctional(int icov, const AFunctional* func, int idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeScaleNoStatFunctional(func, idim);
}
void CovAnisoList::makeAngleNoStatFunctional(int icov, const AFunctional* func, int idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeAngleNoStatFunctional(func, idim);
}
void CovAnisoList::makeTensorNoStatFunctional(int icov, const AFunctional* func, int idim, int jdim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeTensorNoStatFunctional(func, idim, jdim);
}
void CovAnisoList::makeParamNoStatFunctional(int icov, const AFunctional* func)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeParamNoStatFunctional(func);
}
void CovAnisoList::makeRangeStationary(int icov, int idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeRangeStationary(idim);
}
void CovAnisoList::makeScaleStationary(int icov, int idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeScaleStationary(idim);
}
void CovAnisoList::makeAngleStationary(int icov, int idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeAngleStationary(idim);
}

void CovAnisoList::makeTensorStationary(int icov, int idim, int jdim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeTensorStationary(idim, jdim);
}
void CovAnisoList::makeParamStationary(int icov)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeParamStationary();
}

void CovAnisoList::appendParams(ListParams& listParams)
{
  if (!_sameRotation)
    CovList::appendParams(listParams);
  else
  {
    //Find the first structure with a rotation
    int icov = ITEST;
    int ncov = getNCov();
    int icur = 0;
    for (int jcov = 0; jcov < ncov; jcov++)
    {
      
      CovAniso* cova = getCovAniso(jcov);
      cova->appendParams(listParams);
      int icurnew = listParams.getNParams();
      int nadded = icurnew - icur;
      icur = icurnew;
      if (cova->getFlagRotation() && icov == ITEST)
      {
        icov = jcov;
      }
    }
  }
}