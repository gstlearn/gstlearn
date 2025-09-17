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
#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/ListParams.hpp"
#include "Basic/ParamInfo.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovContext.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovLMCAnamorphosis.hpp"
#include "Covariances/CovLMCConvolution.hpp"
#include "Covariances/CovLMCTapering.hpp"
#include "Covariances/CovList.hpp"
#include "Db/Db.hpp"
#include "Enum/EModelProperty.hpp"
#include "Space/ASpace.hpp"
#include "geoslib_define.h"

#include <cmath>
#include <cstddef>
#include <vector>

namespace gstlrn
{

CovAnisoList::CovAnisoList(const CovContext& ctxt)
  : CovList(ctxt)
{
  setOptimEnabled(true);
}

CovAnisoList::CovAnisoList(const CovAnisoList& r)
  : CovList(r)
{
}

CovAnisoList& CovAnisoList::operator=(const CovAnisoList& r)
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

void CovAnisoList::addCovList(const CovAnisoList& covs)
{
  for (Id icov = 0, ncov = covs.getNCov(); icov < ncov; icov++)
    addCov(*covs.getCovAniso(icov));
}

void CovAnisoList::addCov(const CovBase& cov)
{
  if (dynamic_cast<const CovAniso*>(&cov) == nullptr)
  {
    messerr("The argument should be of type 'CovAniso*'");
    return;
  }
  CovList::addCov(cov);
}

const CovAniso* CovAnisoList::_getCovAniso(Id icov) const
{
  if (!_isCovarianceIndexValid(icov)) return nullptr;
  const auto* covaniso = dynamic_cast<const CovAniso*>(_covs[icov].get());
  if (covaniso == nullptr)
  {
    messerr("The element 'icov' is not a CovAniso");
  }
  return covaniso;
}

CovAniso* CovAnisoList::_getCovAnisoModify(Id icov)
{
  if (!_isCovarianceIndexValid(icov)) return nullptr;
  auto covaniso = std::dynamic_pointer_cast<CovAniso>(_covs[icov]);
  if (covaniso == nullptr)
  {
    messerr("The element 'icov' is not a CovAniso");
  }
  return covaniso.get();
}

bool CovAnisoList::isConsistent(const ASpace* /*space*/) const
{
  /// TODO : CovAnisoList::isConsistent
  return true;
}

Id CovAnisoList::getNVar() const
{
  if (getNCov() > 0)
    return _covs[0]->getNVar();
  return 0;
}

double CovAnisoList::eval0(Id ivar, Id jvar, const CovCalcMode* mode) const
{
  double cov            = 0.;
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

  for (Id icov = 0, ncov = getNCov(); icov < ncov; icov++)
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

Id CovAnisoList::getNCov(bool skipNugget) const
{
  Id ncov = static_cast<Id>(_covs.size());
  if (!skipNugget) return ncov;

  Id nstruc = 0;
  for (Id icov = 0; icov < ncov; icov++)
  {
    if (getCovAniso(icov)->getType() != ECov::NUGGET) nstruc++;
  }
  return nstruc;
}

bool CovAnisoList::hasRange() const
{
  for (Id i = 0, n = getNCov(); i < n; i++)
  {
    if (!getCovAniso(i)->hasRange()) return false;
  }
  return true;
}

bool CovAnisoList::isStationary() const
{
  for (Id i = 0, n = getNCov(); i < n; i++)
  {
    if (getCovAniso(i)->getMinOrder() >= 0) return false;
  }
  return true;
}

CovAniso CovAnisoList::extractCova(Id icov) const
{
  const CovAniso* covaniso = _getCovAniso(icov);
  if (covaniso == nullptr)
    return CovAniso(ECov::NUGGET, CovContext());
  return *_getCovAniso(icov);
}

/**
 * @return The Minimum IRF-order induced by the covariances
 */
Id CovAnisoList::getCovMinIRFOrder() const
{
  Id nmini = -1;
  for (Id i = 0, n = getNCov(); i < n; i++)
  {
    const CovAniso* covaniso = _getCovAniso(i);
    if (covaniso == nullptr) continue;
    Id locmini = covaniso->getMinOrder();
    if (locmini > nmini) nmini = locmini;
  }
  return nmini;
}

CovAniso* CovAnisoList::getCovAniso(Id icov)
{
  if (!_isCovarianceIndexValid(icov)) return nullptr;
  return _getCovAnisoModify(icov);
}
const CovAniso* CovAnisoList::getCovAniso(Id icov) const
{
  if (!_isCovarianceIndexValid(icov)) return nullptr;
  return _getCovAniso(icov);
}
void CovAnisoList::setCov(Id icov, const CovBase* covs)
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
const ECov& CovAnisoList::getCovType(Id icov) const
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

String CovAnisoList::getCovName(Id icov) const
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
double CovAnisoList::getParam(Id icov) const
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
double CovAnisoList::getRange(Id icov) const
{
  if (!_isCovarianceIndexValid(icov)) return 0.;
  const CovAniso* covaniso = _getCovAniso(icov);
  if (covaniso == nullptr)
  {
    messerr("The argument should be of type 'CovAniso*'");
    return 0.;
  }
  return covaniso->getRangeIso();
}
VectorDouble CovAnisoList::getRanges(Id icov) const
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
VectorDouble CovAnisoList::getAngles(Id icov) const
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
void CovAnisoList::setRangeIsotropic(Id icov, double range)
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
void CovAnisoList::setParam(Id icov, double value)
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
void CovAnisoList::setMarkovCoeffs(Id icov, const VectorDouble& coeffs)
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
void CovAnisoList::setType(Id icov, const ECov& type)
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

Id CovAnisoList::getNGradParam(Id icov) const
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
double CovAnisoList::getTotalSill(Id ivar, Id jvar) const
{
  double sill_total = 0.;
  for (Id icov = 0, ncov = getNCov(); icov < ncov; icov++)
  {
    const CovAniso* cova = getCovAniso(icov);
    if (cova->getMinOrder() >= 0) return TEST;
    sill_total += cova->getSill(ivar, jvar);
  }
  return sill_total;
}

bool CovAnisoList::_isCovarianceIndexValid(Id icov) const
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
  for (Id icov = 0, ncov = getNCov(); icov < ncov; icov++)
  {
    const CovAniso* cova = getCovAniso(icov);
    if (!cova->hasRange()) continue;
    double range = cova->getRangeIso();
    if (range > maxdist) maxdist = range;
  }
  return maxdist;
}

bool CovAnisoList::hasNugget() const
{
  for (Id is = 0, ns = getNCov(); is < ns; is++)
  {
    if (getCovType(is) == ECov::NUGGET) return true;
  }
  return false;
}

Id CovAnisoList::getRankNugget() const
{
  for (Id is = 0, ns = getNCov(); is < ns; is++)
  {
    if (getCovType(is) == ECov::NUGGET) return is;
  }
  return -1;
}

const CovAnisoList* CovAnisoList::createReduce(const VectorInt& validVars) const
{
  CovAnisoList* newcovlist = this->clone();

  for (Id is = 0, ns = getNCov(); is < ns; is++)
  {
    CovAniso* covs = newcovlist->getCovAniso(is);
    newcovlist->setCov(is, covs->createReduce(validVars));
  }
  newcovlist->setNVar(static_cast<Id>(validVars.size()));
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

Id CovAnisoList::hasExternalCov() const
{
  for (Id icov = 0; icov < static_cast<Id>(getNCov()); icov++)
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

const gstlrn::AnamHermite* CovAnisoList::getAnamHermite() const
{
  const gstlrn::AAnam* anam = getAnam();
  if (anam == nullptr) return nullptr;
  const auto* anamH = dynamic_cast<const gstlrn::AnamHermite*>(anam);
  return anamH;
}

const EModelProperty& CovAnisoList::getCovMode() const
{
  if (dynamic_cast<const CovLMCTapering*>(this) != nullptr) return EModelProperty::TAPE;

  if (dynamic_cast<const CovLMCConvolution*>(this) != nullptr)
    return EModelProperty::CONV;

  if (dynamic_cast<const CovLMCAnamorphosis*>(this) != nullptr)
    return EModelProperty::ANAM;

  return EModelProperty::NONE;
}

void CovAnisoList::makeRangeNoStatDb(Id icov, const String& namecol, Id idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeRangeNoStatDb(namecol, idim);
}

void CovAnisoList::makeScaleNoStatDb(Id icov, const String& namecol, Id idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeScaleNoStatDb(namecol, idim, nullptr);
}
void CovAnisoList::makeAngleNoStatDb(Id icov, const String& namecol, Id idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeAngleNoStatDb(namecol, idim);
}

void CovAnisoList::makeTensorNoStatDb(Id icov, const String& namecol, Id idim, Id jdim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeTensorNoStatDb(namecol, idim, jdim);
}
void CovAnisoList::makeParamNoStatDb(Id icov, const String& namecol)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeParamNoStatDb(namecol);
}
void CovAnisoList::makeRangeNoStatFunctional(Id icov, const AFunctional* func, Id idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeRangeNoStatFunctional(func, idim);
}
void CovAnisoList::makeScaleNoStatFunctional(Id icov, const AFunctional* func, Id idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeScaleNoStatFunctional(func, idim);
}
void CovAnisoList::makeAngleNoStatFunctional(Id icov, const AFunctional* func, Id idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeAngleNoStatFunctional(func, idim);
}
void CovAnisoList::makeTensorNoStatFunctional(Id icov, const AFunctional* func, Id idim, Id jdim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeTensorNoStatFunctional(func, idim, jdim);
}
void CovAnisoList::makeParamNoStatFunctional(Id icov, const AFunctional* func)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeParamNoStatFunctional(func);
}
void CovAnisoList::makeRangeStationary(Id icov, Id idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeRangeStationary(idim);
}
void CovAnisoList::makeScaleStationary(Id icov, Id idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeScaleStationary(idim);
}
void CovAnisoList::makeAngleStationary(Id icov, Id idim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeAngleStationary(idim);
}

void CovAnisoList::makeTensorStationary(Id icov, Id idim, Id jdim)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeTensorStationary(idim, jdim);
}
void CovAnisoList::makeParamStationary(Id icov)
{
  if (!_isCovarianceIndexValid(icov)) return;
  getCovAniso(icov)->makeParamStationary();
}

void CovAnisoList::appendParams(ListParams& listParams,
                                std::vector<covmaptype>* gradFuncs)
{
  DECLARE_UNUSED(gradFuncs);
  if (!_sameRotation)
  {
    CovList::appendParams(listParams, gradFuncs);
    return;
  }

  // Find the first structure with a rotation
  bool found = false;
  auto ncov  = getNCov();
  std::vector<ParamInfo>* paramscur;
  std::vector<ParamInfo>* paramsref;
  std::vector<size_t> anglesrefLoc;
  for (Id jcov = 0; jcov < ncov; jcov++)
  {

    CovAniso* cova = getCovAniso(jcov);
    cova->appendParams(listParams, gradFuncs);
    if (cova->getFlagRotation()) // If rotation
    {
      paramscur = &cova->getAnglesParam(); // get the current angles

      if (!found) // If first time
      {
        found     = true;
        paramsref = paramscur;    // The ref becomes the current
        for (auto& p: *paramsref) // Create the vector of current locations
        {
          anglesrefLoc.push_back(p.getAddress());
        }
      }
      else
      {
        for (Id i = 0; i < static_cast<Id>(anglesrefLoc.size()); i++)
        {
          paramscur->at(i).setAddress(static_cast<Id>(anglesrefLoc[i]));
        }
      }
    }
  }
}
} // namespace gstlrn
