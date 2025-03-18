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
#include "Arrays/Array.hpp"
#include "Basic/AFunctional.hpp"
#include "Basic/AStringFormat.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/TabNoStat.hpp"
#include "Covariances/TabNoStatCovAniso.hpp"
#include "Db/Db.hpp"
#include "Covariances/NoStatArray.hpp"
#include "Covariances/CorAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Enum/EConsElem.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/FFT.hpp"
#include "Basic/Utilities.hpp"
#include "Space/ASpace.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/SpaceSN.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "geoslib_define.h"
#include <math.h>
#include <functional>
#include <vector>

static int NWGT[4]      = {2, 3, 4, 5};
static int NORWGT[4]    = {2, 6, 20, 70};
static int COVWGT[4][5] = {{2, -2, 0, 0, 0},
                           {6, -8, 2, 0, 0},
                           {20, -30, 12, -2, 0},
                           {70, -112, 56, -16, 2}};

CorAniso::CorAniso(const ECov& type, const CovContext& ctxt)
  : ACov(ctxt)
  , /// TODO : shared pointer
  _corfunc(CovFactory::createCovFunc(type, ctxt))
  , _aniso(ctxt.getSpace()->getNDim())
  , _tabNoStatCovAniso(nullptr)
  , _noStatFactor(1.)
{
  initFromContext();
}

CorAniso::CorAniso(const String& symbol, const CovContext& ctxt)
  : ACov(ctxt)
  , /// TODO : shared pointer
  _corfunc()
  , _aniso(ctxt.getSpace()->getNDim())
  , _tabNoStatCovAniso(nullptr)
  , _noStatFactor(1.)
{
  ECov covtype = CovFactory::identifyCovariance(symbol, ctxt);
  _corfunc     = CovFactory::createCovFunc(covtype, ctxt);
  initFromContext();
}

CorAniso::CorAniso(const ECov& type,
                   double range,
                   double param,
                   const CovContext& ctxt,
                   bool flagRange)
  : ACov(ctxt)
  , _corfunc(CovFactory::createCovFunc(type, ctxt))
  , _aniso(ctxt.getSpace()->getNDim())
  , _tabNoStatCovAniso(nullptr)
  , _noStatFactor(1.)
{
  initFromContext();

  // Param
  setParam(param);

  // Range
  if (flagRange)
    setRangeIsotropic(range);
  else
    setScale(range);
}

CorAniso::CorAniso(const CorAniso& r)
  : ACov(r)
  , _corfunc(CovFactory::duplicateCovFunc(*r._corfunc))
  , _aniso(r._aniso)
  , _tabNoStatCovAniso(new TabNoStatCovAniso(*r._tabNoStatCovAniso))
  , _noStatFactor(r._noStatFactor)
{
  _tabNoStat = _tabNoStatCovAniso;
}

CorAniso& CorAniso::operator=(const CorAniso& r)
{
  if (this != &r)
  {
    ACov::operator=(r);
    _corfunc           = CovFactory::duplicateCovFunc(*r._corfunc);
    _aniso             = r._aniso;
    _tabNoStatCovAniso = new TabNoStatCovAniso(*_tabNoStatCovAniso);
    _noStatFactor      = r._noStatFactor;
    _tabNoStat         = _tabNoStatCovAniso;
  }
  return *this;
}

CorAniso::~CorAniso()
{
  delete _corfunc;
}

TabNoStat* CorAniso::_createNoStatTab()
{
  _tabNoStatCovAniso = new TabNoStatCovAniso();
  return _tabNoStatCovAniso;
}
void CorAniso::computeCorrec()
{
  _corfunc->computeCorrec(getNDim());
}

void CorAniso::computeMarkovCoeffs()
{
  _corfunc->computeMarkovCoeffs(getNDim());
}

void CorAniso::_setContext(const CovContext& ctxt)
{
  _corfunc->setContext(ctxt);
  updateFromContext();
}

void CorAniso::setParam(double param)
{
  if (!_corfunc->hasParam()) return;
  _corfunc->setParam(param);
  updateFromContext();
}

void CorAniso::setRangeIsotropic(double range)
{
  if (!hasRange()) return;
  if (range <= EPSILON10)
  {
    messerr("Range is too small (%lf). It has been replaced by 1.", range);
    range = 1;
  }
  double scadef = _corfunc->getScadef();
  setScale(range / scadef);
}

void CorAniso::setRanges(const VectorDouble& ranges)
{
  if (!hasRange()) return;
  if ((int)ranges.size() != getNDim())
  {
    messerr("Inconsistency on Space Dimension");
    return;
  }
  for (unsigned int i = 0; i < ranges.size(); i++)
  {
    if (ranges[i] <= EPSILON10)
    {
      messerr("The range in Space dimension (%d) should not be too small", i);
    }
  }
  VectorDouble scales = ranges;
  double scadef       = _corfunc->getScadef();
  VH::divideConstant(scales, scadef);
  setScales(scales);
}

void CorAniso::setRange(int idim, double range)
{
  if (!hasRange()) return;
  if (range <= EPSILON10)
  {
    messerr("The range should not be too small");
    return;
  }
  double scadef = _corfunc->getScadef();
  setScale(idim, range / scadef);
}

void CorAniso::setScale(double scale)
{
  if (!hasRange()) return;
  if (scale <= EPSILON20) // should be less selective than setRange
  {
    messerr("A scale should not be too small");
    return;
  }
  _aniso.setRadiusIsotropic(scale);
  double scadef = _corfunc->getScadef();
  _corfunc->setField(scadef * scale);
}

void CorAniso::setScales(const VectorDouble& scales)
{
  if (!hasRange()) return;
  for (unsigned int i = 0; i < scales.size(); i++)
  {
    if (scales[i] <= EPSILON20) // should be less strict than setRange
    {
      messerr("The scale along Dimension (%d) should not be too small", i);
      return;
    }
  }
  _aniso.setRadiusVec(scales);
  double scadef = _corfunc->getScadef();
  _corfunc->setField(scadef * VH::maximum(scales));
}

void CorAniso::setScale(int idim, double scale)
{
  if (scale <= EPSILON10)
  {
    messerr("A scale should not be too small");
    return;
  }
  _aniso.setRadiusDir(idim, scale);
  double scadef = _corfunc->getScadef();
  _corfunc->setField(scadef * VH::maximum(_aniso.getRadius()));
}

void CorAniso::setAnisoRotation(const Rotation& rot)
{
  if (!hasRange()) return;
  _aniso.setRotation(rot);
}

void CorAniso::setAnisoRotation(const VectorDouble& rot)
{
  if (!hasRange()) return;
  int ndim = getNDim();
  if ((int)rot.size() != ndim * ndim)
  {
    messerr(
      "Dimension of 'rot' (%d) is not compatible with Space Dimension (%d)",
      (int)rot.size(), ndim);
  }
  Rotation r(ndim);
  r.setMatrixDirectVec(rot);
  _aniso.setRotation(r);
}

void CorAniso::setAnisoAngles(const VectorDouble& angles)
{
  if (!hasRange()) return;
  _aniso.setRotationAngles(angles);
}

void CorAniso::setAnisoAngle(int idim, double angle)
{
  if (!hasRange()) return;
  _aniso.setRotationAngle(idim, angle);
}

void CorAniso::setRotationAnglesAndRadius(const VectorDouble& angles,
                                          const VectorDouble& ranges,
                                          const VectorDouble& scales) const
{
  if (!hasRange()) return;

  VectorDouble scales_local;

  if (!scales.empty())
  {
    if (!ranges.empty())
    {
      messerr("You cannot define simultaneously 'ranges' and 'scales'");
      return;
    }

    if ((int)scales.size() != getNDim())
    {
      messerr("Inconsistency on Space Dimension");
      return;
    }
    for (unsigned int i = 0; i < scales.size(); i++)
    {
      if (scales[i] <= EPSILON20) // should be less strict than setRange
      {
        messerr("The scale along Dimension (%d) should not be too small", i);
        return;
      }
    }
    scales_local = scales;
  }

  if (!ranges.empty())
  {
    if ((int)ranges.size() != getNDim())
    {
      messerr("Inconsistency on Space Dimension");
      return;
    }
    for (unsigned int i = 0; i < ranges.size(); i++)
    {
      if (ranges[i] <= EPSILON10)
      {
        messerr("The range in Space dimension (%d) should not be too small", i);
      }
    }
    scales_local  = ranges;
    double scadef = _corfunc->getScadef();
    VH::divideConstant(scales_local, scadef);
  }

  // Perform the assignment and update the tensor

  _aniso.setRotationAnglesAndRadius(angles, scales_local);
}

bool CorAniso::isValidForTurningBand() const
{
  return _corfunc->isValidForTurningBand();
}
double CorAniso::simulateTurningBand(double t0, TurningBandOperate& operTB) const
{
  return _corfunc->simulateTurningBand(t0, operTB);
}
bool CorAniso::isValidForSpectral() const
{
  return _corfunc->isValidForSpectral();
}
MatrixRectangular CorAniso::simulateSpectralOmega(int nb) const
{
  return _corfunc->simulateSpectralOmega(nb);
}
bool CorAniso::isConsistent(const ASpace* space) const
{
  // Check against the Space Type
  if (space->getType() == ESpaceType::RN && !_corfunc->getCompatibleSpaceR()) return false;
  if (space->getType() == ESpaceType::SN && !_corfunc->getCompatibleSpaceS()) return false;

  // Check against the space dimension
  unsigned int maxndim = _corfunc->getMaxNDim();
  return maxndim <= 0 || (maxndim >= space->getNDim());
}

/**
 * Calculate the value of the covariance from the distance between bi-points
 * This distance has been calculated beforehand (possibly using anisotropy)
 * @param h    Input distance
 * @param mode Pointer to CovCalcMode structure (optional)
 * @return The covariance value
 */
double CorAniso::evalCorFromH(double h, const CovCalcMode* mode) const
{
  double cov = 0.;
  if (mode != nullptr)
  {
    int norder = mode->getOrderVario();
    if (norder == 0)
    {

      // Traditional Covariance or Variogram
      cov = _corfunc->evalCorFunc(h) * _noStatFactor;

      // Convert into a variogram
      if (mode->getAsVario()) cov = _corfunc->evalCorFunc(0) - cov;
    }
    else
    {
      double covcum = 0.;

      // Calculate High-order Variogram (only valuable when h != 0)
      for (int iwgt = 1, nwgt = NWGT[norder]; iwgt < nwgt; iwgt++)
      {
        double hp = h * (1. + iwgt);
        covcum += COVWGT[norder][iwgt] * _corfunc->evalCorFunc(hp);
      }
      cov = covcum / NORWGT[norder];
    }
  }
  else
  {
    cov = _corfunc->evalCorFunc(h) * _noStatFactor;
  }
  return cov;
}

double CorAniso::evalCor(const SpacePoint& p1,
                         const SpacePoint& p2,
                         const CovCalcMode* mode,
                         int ivar,
                         int jvar) const
{
  DECLARE_UNUSED(ivar, jvar);
  double h;
  if (p1.isProjected() && p2.isProjected())
    h = _pw2->getDistance(*_pw1);
  else
    h = getSpace()->getDistance(p1, p2, _aniso);
  return evalCorFromH(h, mode);
}

int CorAniso::addEvalCovVecRHSInPlace(vect vect,
                                      const VectorInt& index1,
                                      SpacePoint& pin,
                                      SpacePoint& pout,
                                      const int iech2) const
{
  if (!isOptimEnabled()) 
    return ACov::addEvalCovVecRHSInPlace(vect, index1, pin, pout, iech2);
  for (int i = 0; i < (int)vect.size(); i++)
  { 
    optimizationTransformSP(pin, pout);
    double h = pout.getDistance(_p1As[i]);
    vect[i] += evalCorFromH(h,nullptr);
  }
  return 0;
}
double CorAniso::_eval(const SpacePoint& p1,
                       const SpacePoint& p2,
                       int ivar,
                       int jvar,
                       const CovCalcMode* mode) const
{
  DECLARE_UNUSED(ivar, jvar)
  double cov = evalCor(p1, p2, mode);
  return cov;
}

double CorAniso::evalCovOnSphere(double alpha,
                                 int degree,
                                 bool flagScaleDistance,
                                 const CovCalcMode* mode) const
{
  if (!_corfunc->hasCovOnSphere()) return TEST;
  const ASpace* space    = getDefaultSpaceSh().get();
  const SpaceSN* spaceSn = dynamic_cast<const SpaceSN*>(space);
  if (spaceSn == nullptr) return TEST;

  double scale = getScale();
  if (flagScaleDistance)
  {
    double radius = spaceSn->getRadius();
    scale         = scale / radius;
    alpha         = alpha / radius;
  }

  double value = _corfunc->evalCovOnSphere(alpha, scale, degree);

  if (mode != nullptr && mode->getAsVario())
    value = _corfunc->evalCovOnSphere(0., scale, degree) - value;

  return value;
}

VectorDouble CorAniso::evalSpectrumOnSphere(int n, bool flagNormDistance, bool flagCumul) const
{
  if (!_corfunc->hasSpectrumOnSphere()) return VectorDouble();
  const ASpace* space    = getDefaultSpaceSh().get();
  const SpaceSN* spaceSn = dynamic_cast<const SpaceSN*>(space);
  if (spaceSn == nullptr) return VectorDouble();

  double scale = getScale();
  if (flagNormDistance)
  {
    double radius = spaceSn->getRadius();
    scale /= radius;
  }
  VectorDouble vec = _corfunc->evalSpectrumOnSphere(n, scale);
  if (flagCumul) VH::cumulateInPlace(vec);
  return vec;
}

void CorAniso::setMarkovCoeffs(const VectorDouble& coeffs)
{
  _corfunc->setMarkovCoeffs(coeffs);
  computeCorrec();
}

/* This function computes a polynomial P from two polynomials P1 and P2 and a small constant eps
 * P(x) = P1(x)^2 + x * P2(x)^2 + eps
 * This formula characterizes all the positive polynomials on R+.
 */
void CorAniso::setMarkovCoeffsBySquaredPolynomials(VectorDouble coeffs1,
                                                   VectorDouble coeffs2,
                                                   double eps)
{
  int size1 = (int)coeffs1.size();
  int size2 = (int)coeffs2.size();

  int size = MAX(2 * size1 - 1, 2 * size2);
  VectorDouble coeffs;
  coeffs.resize(size, 0.);

  for (int i = 0; i < size1; i++)
    for (int j = 0; j < size1; j++)
    {
      coeffs[i + j] += coeffs1[i] * coeffs1[j];
    }
  for (int i = 0; i < size2; i++)
    for (int j = 0; j < size2; j++)
    {
      coeffs[i + j + 1] += coeffs2[i] * coeffs2[j];
    }

  coeffs[0] += eps;
  setMarkovCoeffs(coeffs);
}

double CorAniso::getCorrec() const
{
  return _corfunc->getCorrec();
}

double CorAniso::getFullCorrec() const
{
  return _corfunc->getCorrec() / getDetTensor();
}

double CorAniso::getDetTensor() const
{
  VectorDouble scales = getScales();
  double detTensor    = 1.;
  for (auto& e: scales)
  {
    detTensor *= e;
  }
  return detTensor;
}

double CorAniso::evalSpectrum(const VectorDouble& freq, int ivar, int jvar) const
{
  DECLARE_UNUSED(ivar, jvar)
  if (!_corfunc->hasSpectrumOnRn()) return TEST;

  SpacePoint p1;
  SpacePoint p2;
  p2.setCoords(freq);
  double freqnorm = getSpace()->getFrequentialDistance(p1, p2, _aniso);
  double val      = _corfunc->evaluateSpectrum(freqnorm * freqnorm);
  return val / getCorrec();
}

double CorAniso::normalizeOnSphere(int n) const
{
  const ASpace* space    = getDefaultSpaceSh().get();
  const SpaceSN* spaceSn = dynamic_cast<const SpaceSN*>(space);
  double scale           = getScale();
  double radius          = spaceSn->getRadius();
  scale                  = scale / radius;
  return _corfunc->normalizeOnSphere(n, scale);
}

VectorDouble CorAniso::getMarkovCoeffs() const
{
  if (!_corfunc->hasMarkovCoeffs()) return VectorDouble();

  return _corfunc->getMarkovCoeffs();
}

VectorDouble CorAniso::evalCovOnSphereVec(const VectorDouble& alpha,
                                          int degree,
                                          bool flagScaleDistance,
                                          const CovCalcMode* mode) const
{
  int n = (int)alpha.size();
  VectorDouble vec(n);
  for (int i = 0; i < n; i++)
    vec[i] = evalCovOnSphere(alpha[i], degree, flagScaleDistance, mode);
  return vec;
}

String CorAniso::toStringParams(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt)
  std::stringstream sstr;

  if (_corfunc->hasRange() > 0)
  {

    // A sill is defined

    // Isotropy vs anisotropy

    if (_aniso.isIsotropic())
    {
      sstr << "- Range        = " << toDouble(getRange(0)) << std::endl;
      if (isAsymptotic())
        sstr << "- Theo. Range  = " << toDouble(getScale(0)) << std::endl;
    }
    else
    {
      sstr << toVector("- Ranges       = ", getRanges());
      if (isAsymptotic())
        sstr << toVector("- Theo. Ranges = ", getScales());
      if (!_aniso.getRotation().isIdentity())
      {
        VectorDouble angles = GeometryHelper::formatAngles(getAnisoAngles(), 180.);
        sstr << toVector("- Angles       = ", angles);
        sstr << toMatrix("- Rotation Matrix", VectorString(), VectorString(),
                         true, getNDim(), getNDim(), getAnisoRotMat().getValues());
      }
    }
  }
  else if (_corfunc->hasRange() < 0)
  {
    if (!_aniso.isIsotropic())
    {
      sstr << toVector("- Aniso, Coeff = ", _aniso.getRadius());
      if (!_aniso.getRotation().isIdentity())
      {
        VectorDouble angles = GeometryHelper::formatAngles(getAnisoAngles(), 180.);
        sstr << toVector("- Angles       = ", angles);
        sstr << toMatrix("- Rotation Matrix", VectorString(), VectorString(),
                         true, getNDim(), getNDim(), getAnisoRotMat().getValues());
      }
    }
  }
  return sstr.str();
}

String CorAniso::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  // Covariance Name
  sstr << _corfunc->toString(strfmt);
  sstr << toStringParams(strfmt);

  // Non-stationary parameters
  return sstr.str();
}

/*****************************************************************************/
/*!
 **  Update the Model in the case of Non-stationary parameters
 **  This requires the knowledge of the two end-points
 **
 ** \param[in]  covint       Internal structure for non-stationarity
 **                          or NULL (for stationary case)
 **
 *****************************************************************************/
void CorAniso::nostatUpdate(CovInternal* covint)
{
  if (covint == NULL) return;
  updateCovByPoints(covint->getIcas1(), covint->getIech1(),
                    covint->getIcas2(), covint->getIech2());
}

VectorDouble CorAniso::getRanges() const
{
  VectorDouble range = getScales();
  double scadef      = _corfunc->getScadef();
  if (!hasRange()) scadef = 0.;
  VH::multiplyConstant(range, scadef);
  return range;
}

void CorAniso::setType(const ECov& type)
{
  delete _corfunc;
  _corfunc = CovFactory::createCovFunc(type, getContext());
}

/**
 * This function returns the range in the isotropic case
 * In the anisotropic case, it returns the largest range over all directions
 * @return
 */
double CorAniso::getRange() const
{
  if (!hasRange()) return 0.;
  if (isIsotropic())
    return getRange(0);
  return VH::maximum(getRanges());
}

double CorAniso::getScale() const
{
  if (!hasRange()) return 0.;
  if (isIsotropic())
    return getScale(0);
  return VH::maximum(getScales());
}

VectorDouble CorAniso::getAnisoCoeffs() const
{
  VectorDouble coef = getRanges();
  double max        = VH::maximum(coef);
  if (isZero(max))
  {
    messerr("Range is null");
    return VectorDouble();
  }
  VH::divideConstant(coef, max);
  return coef;
}

/**
 * For compatibility, this function returns 0 if the Covariance has no Third Parameter
 *
 * @return Third parameter
 */
double CorAniso::getParam() const
{
  if (!hasParam())
    return 0.;
  return _corfunc->getParam();
}

void CorAniso::_initFromContext()
{
  int ndim = getNDim();
  _aniso.init(ndim);
  updateFromContext();
  createNoStatTab();
  setOptimEnabled(true);
}

void CorAniso::_updateFromContext()
{
  computeMarkovCoeffs();
  computeCorrec();
}

/**
 * Calculate the Integral Range in various Space Dimension (1, 2 or 3)
 * @return
 */
double CorAniso::getIntegralRange(int ndisc, double hmax) const
{
  int ndim = getNDim();
  SpacePoint dd(VectorDouble(ndim), -1);
  double delta = hmax / ndisc;
  double total = 0.;
  switch (ndim)
  {
    case 1:
      for (int j1 = -ndisc; j1 <= ndisc; j1++)
      {
        dd.setCoord(0, delta * j1);
        total += delta * evalCov(dd, SpacePoint());
      }
      break;

    case 2:
      for (int j1 = -ndisc; j1 <= ndisc; j1++)
        for (int j2 = -ndisc; j2 <= ndisc; j2++)
        {
          dd.setCoord(0, delta * j1);
          dd.setCoord(1, delta * j2);
          total += delta * delta * evalCov(dd, SpacePoint());
        }
      break;

    case 3:
      for (int j1 = -ndisc; j1 <= ndisc; j1++)
        for (int j2 = -ndisc; j2 <= ndisc; j2++)
          for (int j3 = -ndisc; j3 <= ndisc; j3++)
          {
            dd.setCoord(0, delta * j1);
            dd.setCoord(1, delta * j2);
            dd.setCoord(2, delta * j3);
            total += delta * delta * delta * evalCov(dd, SpacePoint());
          }
      break;

    default:
      my_throw("Integral Range has only been programmed for Space Dimension 1 to 3");
  }
  return total;
}

bool CorAniso::_isVariableValid(int ivar) const
{
  return checkArg("Rank of the Variable", ivar, getNVar());
}

int CorAniso::getNGradParam() const
{
  int ndim   = getNDim();
  int number = 0;

  // Anisotropy ranges
  if (hasRange())
  {
    // Anisotropy ranges
    number += ndim;

    // Anisotropy Rotation angles
    if (ndim == 2)
      number++;
    else
      number += ndim;
  }
  return number;
}

CorAniso* CorAniso::createIsotropic(const CovContext& ctxt,
                                    const ECov& type,
                                    double range,
                                    double param,
                                    bool flagRange)
{
  if (ctxt.getNVar() != 1)
  {
    messerr("This function is dedicated to the Monovariate case");
    return nullptr;
  }
  return new CorAniso(type, range, param, ctxt, flagRange);
}

CorAniso* CorAniso::createAnisotropic(const CovContext& ctxt,
                                      const ECov& type,
                                      const VectorDouble& ranges,
                                      double param,
                                      const VectorDouble& angles,
                                      bool flagRange)
{
  if (ctxt.getNVar() != 1)
  {
    messerr("This function is dedicated to the Monovariate case");
    return nullptr;
  }
  int ndim = (int)ranges.size();
  if ((int)ctxt.getNDim() != ndim)
  {
    messerr("Mismatch in Space Dimension between 'ranges'(%d) and 'ctxt'(%d)",
            ndim, ctxt.getNDim());
    return nullptr;
  }

  CorAniso* cov = new CorAniso(type, ctxt);
  if (flagRange)
    cov->setRanges(ranges);
  else
    cov->setScales(ranges);
  cov->setParam(param);
  if (!angles.empty()) cov->setAnisoAngles(angles);
  return cov;
}

CorAniso* CorAniso::createIsotropicMulti(const CovContext& ctxt,
                                         const ECov& type,
                                         double range,
                                         double param,
                                         bool flagRange)
{
  CorAniso* cov = new CorAniso(type, ctxt);

  if (flagRange)
    cov->setRangeIsotropic(range);
  else
    cov->setScale(range);
  cov->setParam(param);
  return cov;
}

CorAniso* CorAniso::createAnisotropicMulti(const CovContext& ctxt,
                                           const ECov& type,
                                           const VectorDouble& ranges,
                                           double param,
                                           const VectorDouble& angles,
                                           bool flagRange)
{

  int ndim = (int)ranges.size();
  if ((int)ctxt.getNDim() != ndim)
  {
    messerr("Mismatch in Space Dimension between 'ranges'(%d) and 'ctxt'(%d)",
            ndim, ctxt.getNDim());
    return nullptr;
  }

  CorAniso* cov = new CorAniso(type, ctxt);
  if (flagRange)
    cov->setRanges(ranges);
  else
    cov->setScales(ranges);
  cov->setParam(param);
  if (!angles.empty()) cov->setAnisoAngles(angles);
  return cov;
}

void CorAniso::_copyCovContext(const CovContext& ctxt)
{
  if (_corfunc != nullptr) _corfunc->copyCovContext(ctxt);
}

Array CorAniso::evalCovFFT(const VectorDouble& hmax,
                           int N,
                           int ivar,
                           int jvar) const
{
  if (!hasSpectrumOnRn()) return Array();

  std::function<double(const VectorDouble&)> funcSpectrum;
  funcSpectrum = [this, ivar, jvar](const VectorDouble& freq)
  {
    return evalSpectrum(freq, ivar, jvar) * getDetTensor();
  };
  return evalCovFFTSpatial(hmax, N, funcSpectrum);
}

void CorAniso::_optimizationPostProcess() const
{
}

/**
 * Transform a space point using the anisotropy tensor
 * @param ptin  Input Space Point
 * @param ptout Output Space Point
 */
void CorAniso::optimizationTransformSP(const SpacePoint& ptin,
                                       SpacePoint& ptout) const
{
  ptout.setIech(ptin.getIech());

  if (isOptimEnabled())
  {
    _aniso.applyInverseInPlace(ptin.getCoords(), ptout.getCoordRef());
    ptout.setProjected(true);
  }
}

/**
 * Transform a set of Space Points using the anisotropy tensor
 * The set of resulting Space Points are stored as private member of this.
 * Note that ALL samples are processed, regardless from the selection
 * @param mode 1 for _p1As and 2 for _p2As
 * @param ps Vector of SpacePoints
 */

void CorAniso::_optimizationPreProcess(int mode, const std::vector<SpacePoint>& ps) const
{
  if (!isOptimEnabled())
  {
    ACov::_optimizationPreProcess(mode, ps);
    return;
  }

  if (mode == 1)
    _p1As.clear();
  else
    _p2As.clear();

  SpacePoint pt(getSpace());
  for (int i = 0, n = (int)ps.size(); i < n; i++)
  {
    // Assert that the Space Point has been transformed
    pt.setIech(ps[i].getIech());

    // Perform the projection
    if (!ps[i].isFFFF())
      optimizationTransformSP(ps[i], pt);
    else
      pt.setFFFF();

    // Stack the Space point in _p1As or _p2As
    if (mode == 1)
      _p1As.push_back(pt);
    else
      _p2As.push_back(pt);
  }
}

/**
 * Checks that the Optimization has already been initiated, by:
 * - checking that the storage (for Sample Points projected in the Covariance
 * rotation system) is already allocated
 * - checking that the dimension of this storage is correct (only if 'db' is provided):
 * in particular, this check is not necessary when freeing this storage.
 */
bool CorAniso::isOptimizationInitialized(const std::vector<SpacePoint>& p1As,
                                         const Db* db)
{
  if (p1As.empty()) return false;
  if (db == nullptr) return true;
  int n = (int)p1As.size();
  return n == db->getNSample();
}

bool CorAniso::isNoStat() const
{
  return isNoStatForAnisotropy() || isNoStatForParam();
}
String CorAniso::toStringNoStat(const AStringFormat* strfmt, int i) const
{
  String sstr;
  sstr = _tabNoStatCovAniso->toStringInside(strfmt, i);
  return sstr;
}
///////////////////// Range ////////////////////////
void CorAniso::makeRangeNoStatDb(const String& namecol, int idim, const Db* db)
{
  if (!_checkTensor()) return;
  makeElemNoStat(EConsElem::RANGE, idim, 0, nullptr, db, namecol);
}

void CorAniso::makeRangeNoStatFunctional(const AFunctional* func, int idim)
{
  if (!_checkTensor()) return;
  makeElemNoStat(EConsElem::RANGE, idim, 0, func);
}

void CorAniso::makeRangeStationary(int idim)
{
  if (_tabNoStatCovAniso->removeElem(EConsElem::RANGE, idim) == 0 &&
      _tabNoStatCovAniso->removeElem(EConsElem::SCALE, idim) == 0)
  {
    messerr("This parameter was already stationary!");
  }
}

///////////////////// Scale ////////////////////////

void CorAniso::makeScaleNoStatDb(const String& namecol, int idim, const Db* db)
{
  if (!_checkTensor()) return;
  makeElemNoStat(EConsElem::SCALE, idim, 0, nullptr, db, namecol);
}

void CorAniso::makeScaleNoStatFunctional(const AFunctional* func, int idim)
{
  if (!_checkTensor()) return;
  makeElemNoStat(EConsElem::SCALE, idim, 0, func);
}

void CorAniso::makeScaleStationary(int idim)
{
  makeRangeStationary(idim);
}

///////////////////// Angle ////////////////////////

void CorAniso::makeAngleNoStatDb(const String& namecol, int idim, const Db* db)
{
  if (!_checkTensor()) return;
  makeElemNoStat(EConsElem::ANGLE, idim, 0, nullptr, db, namecol);
}

void CorAniso::makeAngleNoStatFunctional(const AFunctional* func, int idim)
{
  if (!_checkTensor()) return;
  makeElemNoStat(EConsElem::ANGLE, idim, 0, func);
}

void CorAniso::makeAngleStationary(int idim)
{
  if (_tabNoStatCovAniso->removeElem(EConsElem::ANGLE, idim) == 0)
  {
    messerr("This parameter was already stationary!");
  }
}
///////////////////// Tensor ////////////////////////

void CorAniso::makeTensorNoStatDb(const String& namecol, int idim, int jdim, const Db* db)
{
  if (!_checkRotation()) return;
  if (!_checkDims(idim, jdim)) return;
  makeElemNoStat(EConsElem::TENSOR, idim, jdim, nullptr, db, namecol);
}

void CorAniso::makeTensorNoStatFunctional(const AFunctional* func, int idim, int jdim)
{
  if (!_checkRotation()) return;
  if (!_checkDims(idim, jdim)) return;
  makeElemNoStat(EConsElem::TENSOR, idim, jdim, func);
}

void CorAniso::makeTensorStationary(int idim, int jdim)
{
  if (!_checkDims(idim, jdim)) return;
  if (_tabNoStatCovAniso->removeElem(EConsElem::TENSOR, idim, jdim) == 0)
  {
    messerr("This parameter was already stationary!");
  }
}

///////////////////// Param ////////////////////////

void CorAniso::makeParamNoStatDb(const String& namecol, const Db* db)
{
  if (!_checkParam()) return;
  makeElemNoStat(EConsElem::PARAM, 0, 0, nullptr, db, namecol);
}

void CorAniso::makeParamNoStatFunctional(const AFunctional* func)
{
  if (!_checkParam()) return;
  makeElemNoStat(EConsElem::PARAM, 0, 0, func);
}

void CorAniso::makeParamStationary()
{
  if (!_checkParam()) return;
  if (_tabNoStatCovAniso->removeElem(EConsElem::PARAM) == 0)
  {
    messerr("This parameter was already stationary!");
  }
}

/////////////////////////// Check functions ////////////////////:

bool CorAniso::_checkTensor() const
{
  if (isNoStatForTensor())
  {
    messerr("You have already defined non stationarity by using Tensor specifications");
    messerr("Use makeTensorStationary before specifying other non stationary parameters");
    messerr("for anisotropy.");
    return false;
  }
  return true;
}

bool CorAniso::_checkRotation() const
{
  if (isNoStatForRotation())
  {
    messerr("You have already defined non stationarity by using rotation");
    messerr("specifications (range, scale or angle).");
    messerr("Make these parameters stationary (e.g by makeRangeStationary) before specifying");
    messerr("non stationary tensors");
    return false;
  }
  return true;
}

bool CorAniso::_checkParam() const
{
  if (getType() != ECov::MATERN)
  {
    messerr("This covariance function has no parameters of this type");
    return false;
  }
  return true;
}

double CorAniso::getValue(const EConsElem& econs, int iv1, int iv2) const
{
  DECLARE_UNUSED(iv2)
  if (econs == EConsElem::RANGE)
    return getRange(iv1);
  if (econs == EConsElem::SCALE)
    return getScale(iv1);
  if (econs == EConsElem::ANGLE)
    return getAnisoAngles()[iv1];
  if (econs == EConsElem::PARAM)
    return getParam();
  return TEST;
}

void CorAniso::informMeshByMeshForAnisotropy(const AMesh* amesh) const
{
  for (const auto& e: _listaniso)
  {
    _tabNoStatCovAniso->informMeshByMesh(amesh, e);
  }
}

void CorAniso::informMeshByApexForAnisotropy(const AMesh* amesh) const
{
  for (const auto& e: _listaniso)
    _tabNoStatCovAniso->informMeshByApex(amesh, e);
}

void CorAniso::informDbInForAnisotropy(const Db* dbin) const
{
  for (const auto& e: _listaniso)
    _tabNoStatCovAniso->informDbIn(dbin, e);
}

void CorAniso::informDbOutForAnisotropy(const Db* dbout) const
{
  for (const auto& e: _listaniso)
    _tabNoStatCovAniso->informDbOut(dbout, e);
}

/**
 * Update the Model according to the Non-stationary parameters
 * @param icas1 Type of first Db: 1 for Input; 2 for Output
 * @param iech1 Rank of the target within Db1 (or -1)
 * @param icas2 Type of first Db: 1 for Input; 2 for Output
 * @param iech2 Rank of the target within Dbout (or -2)
 */
void CorAniso::updateCovByPoints(int icas1, int iech1, int icas2, int iech2)
{
  // If no non-stationary parameter is defined, simply skip
  if (!_tabNoStatCovAniso->isNoStat()) return;
  double val1, val2;

  int ndim = getNDim();

  const auto paramsnostat = _tabNoStatCovAniso->getTable();
  // Loop on the elements that can be updated one-by-one

  for (const auto& e: paramsnostat)
  {
    EConsElem type = e.first.getType();
    e.second->getValuesOnDb(icas1, iech1, &val1, icas2, iech2, &val2);

    if (type == EConsElem::PARAM)
    {
      setParam(0.5 * (val1 + val2));
    }
  }

  // Loop on the other parameters (Anisotropy) that must be processed globally

  if (!isNoStatForAnisotropy()) return;

  VectorDouble angle1;
  VectorDouble angle2;

  VectorDouble scale1;
  VectorDouble scale2;

  VectorDouble range1;
  VectorDouble range2;

  // Define the angles (for all space dimensions)
  bool flagRotTwo = false;
  bool flagRotOne = false;

  if (getNAngles() > 0)
  {
    angle1 = getAnisoAngles();
    angle2 = angle1;
    for (int idim = 0; idim < ndim; idim++)
    {
      if (_tabNoStatCovAniso->isElemDefined(EConsElem::ANGLE, idim))
      {
        auto noStat = _tabNoStatCovAniso->getElem(EConsElem::ANGLE, idim);
        flagRotOne  = true;
        if (noStat->getValuesOnDb(icas1, iech1, &angle1[idim], icas2, iech2, &angle2[idim]))
          flagRotTwo = true;
      }
    }
  }

  // Define the Theoretical ranges (for all space dimensions)

  bool flagScaleTwo = false;
  bool flagScaleOne = false;
  if (getNScales() > 0)
  {
    scale1 = getScales();
    scale2 = scale1;
    for (int idim = 0; idim < ndim; idim++)
    {
      if (_tabNoStatCovAniso->isElemDefined(EConsElem::SCALE, idim))
      {
        auto noStat  = _tabNoStatCovAniso->getElem(EConsElem::SCALE, idim);
        flagScaleOne = true;
        if (noStat->getValuesOnDb(icas1, iech1, &scale1[idim], icas2, iech2, &scale2[idim]))
          flagScaleTwo = true;
      }
    }
  }

  // Define the Practical ranges (for all space dimensions)

  bool flagRangeTwo = false;
  bool flagRangeOne = false;
  if (getNRanges() > 0)
  {
    range1 = getRanges();
    range2 = range1;
    for (int idim = 0; idim < ndim; idim++)
    {
      if (_tabNoStatCovAniso->isElemDefined(EConsElem::RANGE, idim))
      {
        auto noStat  = _tabNoStatCovAniso->getElem(EConsElem::RANGE, idim);
        flagRangeOne = true;
        if (noStat->getValuesOnDb(icas1, iech1, &range1[idim], icas2, iech2, &range2[idim]))
          flagRangeTwo = true;
      }
    }
  }

  // Update the Model
  double ratio = 1.;
  if (flagRotTwo || flagRangeTwo || flagScaleTwo)
  {
    // Extract the direct tensor at first point and square it
    setRotationAnglesAndRadius(angle1, range1, scale1);
    MatrixSquareSymmetric direct1 = getAniso().getTensorDirect2();
    double det1                   = pow(direct1.determinant(), 0.25);

    // Extract the direct tensor at second point and square it
    setRotationAnglesAndRadius(angle2, range2, scale2);
    MatrixSquareSymmetric direct2 = getAniso().getTensorDirect2();
    double det2                   = pow(direct2.determinant(), 0.25);

    // Calculate average squared tensor
    direct2.addMatInPlace(direct1, 0.5, 0.5);
    double detM = sqrt(direct2.determinant());

    // Update the tensor (squared version)
    Tensor tensor = getAniso();
    tensor.setTensorDirect2(direct2);
    setAniso(tensor);
    ratio = det1 * det2 / detM;
  }
  else if (flagRotOne || flagRangeOne || flagScaleOne)
  {
    // Simply update the model with one set of parameters
    setRotationAnglesAndRadius(angle1, range1, scale1);
  }
  setNoStatFactor(ratio);
}

void CorAniso::updateCovByMesh(int imesh, bool aniso) const
{
  // If no non-stationary parameter is defined, simply skip
  if (!_tabNoStatCovAniso->isNoStat()) return;

  int ndim = getNDim();

  // Loop on the elements that can be updated one-by-one
  if (!aniso)
  {
    const auto paramsnostat = _tabNoStatCovAniso->getTable();
    for (const auto& e: paramsnostat)
    {
      EConsElem type = e.first.getType();
    }
    return;
  }
  // Loop on the other parameters (Anisotropy) that must be processed globally

  if (!isNoStatForAnisotropy()) return;

  VectorDouble angles;
  VectorDouble scales;
  VectorDouble ranges;

  // Define the angles (for all space dimensions)
  if (getNAngles() > 0)
  {
    angles = getAnisoAngles();

    for (int idim = 0; idim < ndim; idim++)
    {
      if (_tabNoStatCovAniso->isElemDefined(EConsElem::ANGLE, idim))
      {
        auto noStat  = _tabNoStatCovAniso->getElem(EConsElem::ANGLE, idim);
        angles[idim] = noStat->getValueOnMeshByMesh(imesh);
      }
    }
  }

  // Define the Theoretical ranges (for all space dimensions)
  if (getNScales() > 0)
  {
    scales = getScales();
    for (int idim = 0; idim < ndim; idim++)
    {
      if (_tabNoStatCovAniso->isElemDefined(EConsElem::SCALE, idim))
      {
        auto noStat  = _tabNoStatCovAniso->getElem(EConsElem::SCALE, idim);
        scales[idim] = noStat->getValueOnMeshByMesh(imesh);
      }
    }
  }

  if (getNRanges() > 0)
  {
    ranges = getRanges();

    for (int idim = 0; idim < ndim; idim++)
    {
      if (_tabNoStatCovAniso->isElemDefined(EConsElem::RANGE, idim))
      {
        auto noStat  = _tabNoStatCovAniso->getElem(EConsElem::RANGE, idim);
        ranges[idim] = noStat->getValueOnMeshByMesh(imesh);
      }
    }
  }

  setRotationAnglesAndRadius(angles, ranges, scales);
  // TODO : This part is not finished
  if (isNoStatForTensor())
  {
    for (int idim = 0; idim < ndim; idim++)
    {
      for (int jdim = 0; jdim < ndim; jdim++)
        if (_tabNoStatCovAniso->isElemDefined(EConsElem::TENSOR, idim, jdim))
        {
          auto noStat = _tabNoStatCovAniso->getElem(EConsElem::TENSOR, idim, jdim);
        }
    }
  }
}

void CorAniso::_manage(const Db* db1, const Db* db2) const
{
  if (db1 != nullptr)
    informDbIn(db1);
  if (db2 != nullptr)
    informDbOut(db2);
}

/**
 * @brief Define the argument 'p' as the current target and project it if needed
 *
 * @param p Current space point
 *
 * @note: The target is stored in first position of '_p2As'
 * @note: Its pointer is stored in '_pw2' for quick covariance evaluation.
 */
void CorAniso::_optimizationSetTarget(SpacePoint& p) const
{
  if (!isOptimEnabled()) return;
  int iech = 0;
  optimizationTransformSP(p, _p2As[0]);
  p.setProjected(true);
  _pw2 = &_p2As[iech];
}
