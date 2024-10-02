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
#include "Covariances/TabNoStatCovAniso.hpp"
#include "Db/Db.hpp"
#include "Covariances/NoStatArray.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/NoStatFunctional.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Enum/EConsElem.hpp"
#include "Matrix/MatrixSquareGeneral.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/AException.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/FFT.hpp"
#include "Basic/Utilities.hpp"
#include "Space/ASpace.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpaceSN.hpp"
#include "Geometry/GeometryHelper.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "geoslib_define.h"
#include <math.h>
#include <functional>
#include <memory>

static int NWGT[4] = { 2, 3, 4, 5 };
static int NORWGT[4] = { 2, 6, 20, 70 };
static int COVWGT[4][5] = { { 2, -2, 0, 0, 0 },
                            { 6, -8, 2, 0, 0 },
                            { 20, -30, 12, -2, 0 },
                            { 70, -112, 56, -16, 2 } };

CovAniso::CovAniso(const ECov &type, const CovContext &ctxt)
    : ACov(ctxt.getSpace()), /// TODO : shared pointer
      _ctxt(ctxt),
      _cova(CovFactory::createCovFunc(type, ctxt)),
      _sill(),
      _aniso(ctxt.getSpace()->getNDim()),
      _tabNoStat(),
      _noStatFactor(1.),
      _optimEnabled(true),
      _p1As(),
      _p2A(ctxt.getSpace())
{
  _initFromContext();
}

CovAniso::CovAniso(const String &symbol, const CovContext &ctxt)
    : ACov(ctxt.getSpace()), /// TODO : shared pointer
      _ctxt(ctxt),
      _cova(),
      _sill(),
      _aniso(ctxt.getSpace()->getNDim()),
      _tabNoStat(),
      _noStatFactor(1.),
      _optimEnabled(true),
      _p1As(),
      _p2A(ctxt.getSpace())
{
  ECov covtype = CovFactory::identifyCovariance(symbol, ctxt);
  _cova = CovFactory::createCovFunc(covtype, ctxt);
  _initFromContext();
}

CovAniso::CovAniso(const ECov &type,
                   double range,
                   double param,
                   double sill,
                   const CovContext &ctxt,
                   bool flagRange)
    : ACov(ctxt.getSpace()), /// TODO : shared pointer
      _ctxt(ctxt),
      _cova(CovFactory::createCovFunc(type, ctxt)),
      _sill(),
      _aniso(ctxt.getSpace()->getNDim()),
      _tabNoStat(),
      _noStatFactor(1.),
      _optimEnabled(true),
      _p1As(),
      _p2A(ctxt.getSpace())
{
  _initFromContext();

  // Sill
  if (ctxt.getNVar() == 1)
    _sill.setValue(0, 0, sill);
  else
  {
    int nvar = ctxt.getNVar();
    _sill.fill(0);
    for (int ivar = 0; ivar < nvar; ivar++)
      _sill.setValue(ivar, ivar, sill);
  }

  // Param
  setParam(param);

  // Range
  if (flagRange)
    setRangeIsotropic(range);
  else
    setScale(range);
}

CovAniso::CovAniso(const CovAniso &r)
    : ACov(r),
      _ctxt(r._ctxt),
      _cova(CovFactory::duplicateCovFunc(*r._cova)),
      _sill(r._sill),
      _aniso(r._aniso),
      _tabNoStat(r._tabNoStat),
      _noStatFactor(r._noStatFactor),
      _optimEnabled(r._optimEnabled),
      _p1As(r._p1As),
      _p2A(r._p2A)
{
}

CovAniso& CovAniso::operator=(const CovAniso &r)
{
  if (this != &r)
  {
    ACov::operator =(r);
    _ctxt = r._ctxt;
    _cova = CovFactory::duplicateCovFunc(*r._cova);
    _sill = r._sill;
    _aniso = r._aniso;
    _tabNoStat = r._tabNoStat;
    _noStatFactor = r._noStatFactor;
    _optimEnabled = r._optimEnabled;
    _p1As = r._p1As;
    _p2A = r._p2A;
  }
  return *this;
}

CovAniso::~CovAniso()
{
  delete _cova;
}

void CovAniso::_computeCorrec()
{
  _cova->computeCorrec(getNDim());
}

void CovAniso::computeMarkovCoeffs()
{
  _cova->computeMarkovCoeffs(getNDim());
}

void CovAniso::setContext(const CovContext &ctxt)
{
  _ctxt = ctxt;
  _updateFromContext();
}

void CovAniso::setParam(double param)
{
  if (!_cova->hasParam()) return;
  _cova->setParam(param);
  _updateFromContext();
}

void CovAniso::setSill(double sill)
{
  if (getNVariables() != 1)
  {
    messerr("Number of provided sill doesn't match number of variables");
    return;
  }
  _sill.resetFromValue(1, 1, sill);
}

void CovAniso::setSill(const MatrixSquareSymmetric &sill)
{
  if (getNVariables() != sill.getNSize())
  {
    messerr("Number of provided sills doesn't match number of variables");
    return;
  }
  _sill = sill;
}

void CovAniso::setSill(const VectorDouble &sill)
{
  int size = static_cast<int>(sill.size());
  int nvar = getNVariables();
  if (size != nvar * nvar)
  {
    messerr("Number of provided sills doesn't match number of variables");
    return;
  }
  _sill.setValues(sill);
}

void CovAniso::setSill(int ivar, int jvar, double sill)
{
  if (!_isVariableValid(ivar)) return;
  if (!_isVariableValid(jvar)) return;
  /// TODO : Test if sill matrix is positive definite (if not, generate a warning)
  if (!_sill.isValid(ivar, jvar)) return;
  _sill.setValue(ivar, jvar, sill);
}

void CovAniso::initSill(double value)
{
  _sill.fill(value);
}

void CovAniso::setRangeIsotropic(double range)
{
  if (!hasRange()) return;
  if (range <= EPSILON10)
  {
    messerr("Range is too small (%lf). It has been replaced by 1.", range);
    range = 1;
  }
  double scadef = _cova->getScadef();
  setScale(range / scadef);
}

void CovAniso::setRanges(const VectorDouble &ranges)
{
  if (!hasRange()) return;
  if (ranges.size() != getNDim())
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
  double scadef = _cova->getScadef();
  VH::divideConstant(scales, scadef);
  setScales(scales);
}

void CovAniso::setRange(int idim, double range)
{
  if (!hasRange()) return;
  if (range <= EPSILON10)
  {
    messerr("The range should not be too small");
    return;
  }
  double scadef = _cova->getScadef();
  setScale(idim, range / scadef);
}

void CovAniso::setScale(double scale)
{
  if (!hasRange()) return;
  if (scale <= EPSILON20) // should be less selective than setRange
  {
    messerr("A scale should not be too small");
    return;
  }
  _aniso.setRadiusIsotropic(scale);
  double scadef = _cova->getScadef();
  _cova->setField(scadef * scale);
}

void CovAniso::setScales(const VectorDouble &scales)
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
  double scadef = _cova->getScadef();
  _cova->setField(scadef * VH::maximum(scales));
}

void CovAniso::setScale(int idim, double scale)
{
  if (scale <= EPSILON10)
  {
    messerr("A scale should not be too small");
    return;
  }
  _aniso.setRadiusDir(idim, scale);
  double scadef = _cova->getScadef();
  _cova->setField(scadef * VH::maximum(_aniso.getRadius()));
}

void CovAniso::setAnisoRotation(const Rotation &rot)
{
  if (!hasRange()) return;
  _aniso.setRotation(rot);
}

void CovAniso::setAnisoRotation(const VectorDouble &rot)
{
  if (!hasRange()) return;
  int ndim = getNDim();
  if ((int) rot.size() != ndim * ndim)
  {
    messerr(
        "Dimension of 'rot' (%d) is not compatible with Space Dimension (%d)",
        (int) rot.size(), ndim);
  }
  Rotation r(ndim);
  r.setMatrixDirectVec(rot);
  _aniso.setRotation(r);
}

void CovAniso::setAnisoAngles(const VectorDouble &angles)
{
  if (!hasRange()) return;
  _aniso.setRotationAngles(angles);
}

void CovAniso::setAnisoAngle(int idim, double angle)
{
  if (!hasRange()) return;
  _aniso.setRotationAngle(idim, angle);
}

void CovAniso::setRotationAnglesAndRadius(const VectorDouble &angles,
                                          const VectorDouble &ranges,
                                          const VectorDouble &scales)
{
  if (!hasRange()) return;

  VectorDouble scales_local;

  if (! scales.empty())
  {
    if (! ranges.empty())
    {
      messerr("You cannot define simultaneously 'ranges' and 'scales'");
      return;
    }

    if (scales.size() != getNDim())
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

  if (! ranges.empty())
  {
    if (ranges.size() != getNDim())
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
    scales_local = ranges;
    double scadef = _cova->getScadef();
    VH::divideConstant(scales_local, scadef);
  }

  // Perform the assignment and update the tensor

  _aniso.setRotationAnglesAndRadius(angles, scales_local);
}

bool CovAniso::isValidForTurningBand() const
{
  return _cova->isValidForTurningBand();
}
double CovAniso::simulateTurningBand(double t0, TurningBandOperate &operTB) const
{
  return _cova->simulateTurningBand(t0, operTB);
}
bool CovAniso::isValidForSpectral() const
{
  return _cova->isValidForSpectral();
}
MatrixRectangular CovAniso::simulateSpectralOmega(int nb) const
{
  return _cova->simulateSpectralOmega(nb);
}
bool CovAniso::isConsistent(const ASpace* space) const
{
  // Check against the Space Type
  if (space->getType() == ESpaceType::RN && ! _cova->getCompatibleSpaceR()) return false;
  if (space->getType() == ESpaceType::SN && ! _cova->getCompatibleSpaceS()) return false;

  // Check against the space dimension
  unsigned int maxndim = _cova->getMaxNDim();
  return maxndim <= 0 || (maxndim >= space->getNDim());
}

/**
 * Calculate the value of the covariance from the distance between bi-points
 * This distance has been calculated beforehand (possibly using anisotropy)
 * @param h    Input distance
 * @param mode Pointer to CovCalcMode structure (optional)
 * @return The covariance value
 */
double CovAniso::_evalCovFromH(double h, const CovCalcMode *mode) const
{
  double cov = 0.;
  if (mode != nullptr)
  {
    int norder = mode->getOrderVario();
    if (norder == 0)
    {

      // Traditional Covariance or Variogram
      cov = _cova->evalCov(h);

      // Convert into a variogram
      if (mode->getAsVario()) cov = _cova->evalCov(0) - cov;
    }
    else
    {
      double covcum = 0.;

      // Calculate High-order Variogram (only valuable when h != 0)
      for (int iwgt = 1, nwgt = NWGT[norder]; iwgt < nwgt; iwgt++)
      {
        double hp = h * (1. + iwgt);
        covcum += COVWGT[norder][iwgt] * _cova->evalCov(hp);
      }
      cov = covcum / NORWGT[norder];
    }
  }
  else
  {
    cov =  _cova->evalCov(h);
  }
  return cov;
}

double CovAniso::eval0(int ivar, int jvar, const CovCalcMode* mode) const
{
  double cov = _evalCovFromH(0, mode);

  if (mode == nullptr || ! mode->getUnitary())
    cov *= getSill(ivar, jvar);
  return (cov);
}

double CovAniso::eval(const SpacePoint &p1,
                      const SpacePoint &p2,
                      int ivar,
                      int jvar,
                      const CovCalcMode* mode) const
{
  double h = getSpace()->getDistance(p1, p2, _aniso);
  double cov = _evalCovFromH(h, mode);

  if (mode == nullptr || ! mode->getUnitary())
    cov *= getSill(ivar, jvar);
  return (cov);
}
/**
 * Calculate the Matrix of covariance for zero distance
 * @param mat   Covariance matrix (Dimension: nvar * nvar)
 * @param mode  Calculation Options
 *
 * @remarks: Matrix 'mat' should be dimensioned and initialized beforehand
 */
void CovAniso::eval0MatInPlace(MatrixSquareSymmetric &mat,
                               const CovCalcMode *mode) const
{
  double cov = _evalCovFromH(0, mode);

  if (mode == nullptr || ! mode->getUnitary())
    mat.addMatInPlace(_sill, 1., cov * _noStatFactor);
  else
  {
    MatrixSquareSymmetric identity = _sill;
    identity.setIdentity();
    mat.addMatInPlace(identity, 1., cov);
  }
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
void CovAniso::evalMatInPlace(const SpacePoint &p1,
                              const SpacePoint &p2,
                              MatrixSquareSymmetric &mat,
                              const CovCalcMode *mode) const
{
  double h = getSpace()->getDistance(p1, p2, _aniso);
  double cov = _evalCovFromH(h, mode);

  if (mode == nullptr || ! mode->getUnitary())
    mat.addMatInPlace(_sill, 1., cov * _noStatFactor);
  else
  {
    MatrixSquareSymmetric identity = _sill;
    identity.setIdentity();
    mat.addMatInPlace(identity, 1., cov);
  }
}

/**
 * Fill the vector of covariances between each valid SpacePoint (recorded in _p1As)
 * and the target (recorded in _p2A)
 * @param res  Vector of covariances
 * @param ivars Arrays of ranks for the first point
 * @param index1 Arrays of sample indices for the first point
 * @param ivar2 Rank of the variable for the second point
 * @param icol  Rank of the column (variable + sample) for the second point
 * @param mode CovCalcMode structure
 * @param flagSym True if used for a Symmetric matrix (should only fill upper triangle)
 *
 * @remark: The optimized version is not compatible with Franck's non-stationarity.
 * Then no correction must be applied to cov(h)
 */
void CovAniso::evalOptimInPlace(MatrixRectangular& res,
                                const VectorInt& ivars,
                                const VectorVectorInt& index1,
                                int ivar2,
                                int icol,
                                const CovCalcMode *mode,
                                bool flagSym) const
{
  double cov, hoptim;
  double sill = 1.;

  // Loop on the first variable
  int irow = 0;
  for (int rvar1 = 0, nvar1 = (int) ivars.size(); rvar1 < nvar1; rvar1++)
  {
    int ivar1 = ivars[rvar1];
    if (mode == nullptr || ! mode->getUnitary())
      sill = _sill.getValue(ivar1, ivar2);

    // Loop on the first sample
    int nech1s = (int) index1[rvar1].size();
    for (int rech1 = 0; rech1 < nech1s; rech1++)
    {
      if (!flagSym || irow <= icol)
      {
        int iech1 = index1[rvar1][rech1];
        hoptim = _p2A.getDistance(_p1As[iech1]);
        cov = _evalCovFromH(hoptim, mode);
        res.updValue(irow, icol, EOperator::ADD, sill * cov);
      }
      irow++;
    }
  }
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


void CovAniso::evalCovLHS(MatrixSquareSymmetric &mat,
                          SpacePoint &pwork1,
                          SpacePoint &pwork2,
                          const Db* db, 
                          const CovCalcMode *mode) const
{
  if (!_isOptimEnabled())
    ACov::evalCovLHS(mat, pwork1, pwork2, db, mode);
  else
  {
    SpacePoint* p1A = &_p1As[pwork1.getIech()]; 
    SpacePoint* p2A = &_p1As[pwork2.getIech()];
    _evalOptim(p1A,p2A,mat,mode);
  // Calculate covariance between two points
  }
}

void CovAniso::_evalOptim(SpacePoint* p1A, SpacePoint* p2A,
                          MatrixSquareSymmetric &mat,
                          const CovCalcMode *mode) const
{
  // Calculate covariance between two points
  double hoptim = p2A->getDistance(*p1A);
  double cov = _evalCovFromH(hoptim, mode);

  if (mode == nullptr || ! mode->getUnitary())
    mat.addMatInPlace(_sill, 1., cov * _noStatFactor);
  else
  {
    MatrixSquareSymmetric identity = _sill;
    identity.setIdentity();
    mat.addMatInPlace(identity, 1., cov);
  }
}

void CovAniso::evalCovRHS(MatrixSquareSymmetric &mat,
                          SpacePoint &pwork1,
                          const Db* db, SpacePoint& pout,  
                          const CovCalcMode *mode) const
{
  if (!_isOptimEnabled())
    ACov::evalCovRHS(mat, pwork1, db, pout, mode);
  else
  {
    SpacePoint* p1A = &_p1As[pwork1.getIech()]; 
    SpacePoint* p2A = &_p2A;
    _evalOptim(p1A,p2A,mat,mode);
  }
}

void CovAniso::evalMatOptimInPlace(int icas1,
                                   int iech1,
                                   int icas2,
                                   int iech2,
                                   MatrixSquareSymmetric &mat,
                                   const CovCalcMode *mode) const
{
  SpacePoint* p1A = (icas1 == 1) ? &_p1As[iech1] : &_p2A;
  SpacePoint* p2A = (icas2 == 1) ? &_p1As[iech2] : &_p2A;

  // Calculate covariance between two points
  double hoptim = p2A->getDistance(*p1A);
  double cov = _evalCovFromH(hoptim, mode);

  if (mode == nullptr || ! mode->getUnitary())
    mat.addMatInPlace(_sill, 1., cov * _noStatFactor);
  else
  {
    MatrixSquareSymmetric identity = _sill;
    identity.setIdentity();
    mat.addMatInPlace(identity, 1., cov);
  }
}

double CovAniso::evalCovOnSphere(double alpha,
                                 int degree,
                                 bool flagScaleDistance,
                                 const CovCalcMode* mode) const
{
  if (!_cova->hasCovOnSphere()) return TEST;
  const ASpace* space = getDefaultSpace();
  const SpaceSN* spaceSn = dynamic_cast<const SpaceSN*>(space);
  if (spaceSn == nullptr) return TEST;

  double scale = getScale();
  if (flagScaleDistance)
  {
    double radius = spaceSn->getRadius();
    scale = scale / radius;
    alpha = alpha / radius;
  }

  double value = _cova->evalCovOnSphere(alpha, scale, degree);

  if (mode != nullptr && mode->getAsVario())
    value = _cova->evalCovOnSphere(0., scale, degree) - value;

  if (mode == nullptr || ! mode->getUnitary())
    value *= getSill(0,0);

  return value;
}

VectorDouble CovAniso::evalSpectrumOnSphere(int n, bool flagNormDistance, bool flagCumul) const
{
  if (!_cova->hasSpectrumOnSphere()) return VectorDouble();
  const ASpace* space = getDefaultSpace();
  const SpaceSN* spaceSn = dynamic_cast<const SpaceSN*>(space);
  if (spaceSn == nullptr) return VectorDouble();

  double scale = getScale();
  if (flagNormDistance)
  {
    double radius = spaceSn->getRadius();
    scale /= radius;
  }
  VectorDouble vec = _cova->evalSpectrumOnSphere(n, scale);
  if (flagCumul) VH::cumulateInPlace(vec);
  return vec;
}

void CovAniso::setMarkovCoeffs(const VectorDouble& coeffs)
{
  _cova->setMarkovCoeffs(coeffs);
  _computeCorrec();
}

/* This function computes a polynomial P from two polynomials P1 and P2 and a small constant eps
 * P(x) = P1(x)^2 + x * P2(x)^2 + eps
 */
void CovAniso::setMarkovCoeffsBySquaredPolynomials(VectorDouble coeffs1,
                                                   VectorDouble coeffs2,
                                                   double eps)
{
  int size1 = (int) coeffs1.size();
  int size2 = (int) coeffs2.size();

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

double CovAniso::getCorrec() const
{
  return _cova->getCorrec();
}

double CovAniso::getFullCorrec() const
{
  return  _cova->getCorrec() / _getDetTensor();
}

double CovAniso::_getDetTensor() const
{
  VectorDouble scales = getScales();
  double detTensor = 1.;
  for (auto &e : scales)
  {
    detTensor *= e;
  }
  return detTensor;
}

double CovAniso::evalSpectrum(const VectorDouble& freq, int ivar, int jvar) const
{
  if (!_cova->hasSpectrumOnRn()) return TEST;

  double sill = getSill(ivar, jvar);

  SpacePoint p1;
  SpacePoint p2;
  p2.setCoords(freq);
  double freqnorm = getSpace()->getFrequentialDistance(p1, p2, _aniso);
  double val = _cova->evaluateSpectrum(freqnorm * freqnorm);
  return  sill * val / getCorrec();
}

double CovAniso::normalizeOnSphere(int n) const
{ 
  const ASpace* space = getDefaultSpace();
  const SpaceSN* spaceSn = dynamic_cast<const SpaceSN*>(space);
  double scale = getScale();
  double radius = spaceSn->getRadius();
  scale = scale / radius;
  return _cova->normalizeOnSphere(n,scale);
}

VectorDouble CovAniso::getMarkovCoeffs() const
{
  if (!_cova->hasMarkovCoeffs()) return VectorDouble();

  return _cova->getMarkovCoeffs();
}

VectorDouble CovAniso::evalCovOnSphereVec(const VectorDouble &alpha,
                                          int degree,
                                          bool flagScaleDistance,
                                          const CovCalcMode* mode) const
{
  int n = (int) alpha.size();
  VectorDouble vec(n);
  for (int i = 0; i < n; i++)
    vec[i] = evalCovOnSphere(alpha[i], degree, flagScaleDistance, mode);
  return vec;
}

String CovAniso::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  // Covariance Name
  sstr << _cova->toString();

  // Sill - Factor / Slope information
  if (_cova->hasRange() > 0)
  {

    // A sill is defined

    if (getNVariables() > 1)
    {
      sstr << toMatrix("- Sill matrix:", VectorString(), VectorString(), 0,
                      getNVariables(), getNVariables(), _sill.getValues());
    }
    else
    {
      sstr << "- Sill         = " << toDouble(_sill.getValue(0, 0)) << std::endl;
    }

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
  else if (_cova->hasRange() < 0)
  {
    // The sill is not defined: use slope instead

    if (getNVariables() > 1)
    {
      MatrixSquareSymmetric slopes = _sill;
      double range = getRange(0);
      for (int ivar = 0; ivar < getNVariables(); ivar++)
        for (int jvar = 0; jvar < getNVariables(); jvar++)
          slopes.setValue(ivar, jvar, _sill.getValue(ivar, jvar) / range);
      sstr << toMatrix("- Slope matrix:", VectorString(), VectorString(), 0,
                      getNVariables(), getNVariables(), slopes.getValues());
    }
    else
    {
      sstr << "- Slope        = " << toDouble(getSlope(0, 0)) << std::endl;
    }

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
  else
  {
    // Only sill is defined

    if (getNVariables() > 1)
    {
      sstr << toMatrix("- Sill matrix:", VectorString(), VectorString(), 0,
                      getNVariables(), getNVariables(), _sill.getValues());
    }
    else
    {
      sstr << "- Sill         = " << toDouble(_sill.getValue(0, 0)) << std::endl;
    }
  }
  // Non-stationary parameters
  if (isNoStat())
  {
    sstr << _tabNoStat.toString() << std::endl;
  }
  return sstr.str();
}

double CovAniso::getSill(int ivar, int jvar) const
{
  return _sill.getValue(ivar, jvar) * _noStatFactor;
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
void CovAniso::nostatUpdate(CovInternal *covint)
{
  if (covint == NULL) return;
  updateCovByPoints(covint->getIcas1(), covint->getIech1(),
                    covint->getIcas2(), covint->getIech2());
}

/**
 * Return the Slope calculated as the sill / range(idim=0)
 * @param ivar Rank of the first variable
 * @param jvar Rank of the second variable
 * @return
 */
double CovAniso::getSlope(int ivar, int jvar) const
{
  if (hasRange() == 0) return TEST;
  double range = getRange(0);
  return _sill.getValue(ivar, jvar) * _noStatFactor / range;
}

VectorDouble CovAniso::getRanges() const
{
  VectorDouble range = getScales();
  double scadef = _cova->getScadef();
  if (!hasRange()) scadef = 0.;
  VH::multiplyConstant(range, scadef);
  return range;
}

void CovAniso::setType(const ECov &type)
{
  delete _cova;
  _cova = CovFactory::createCovFunc(type, _ctxt);
}

/**
 * This function returns the range in the isotropic case
 * In the anisotropic case, it returns the largest range over all directions
 * @return
 */
double CovAniso::getRange() const
{
  if (!hasRange()) return 0.;
  if (isIsotropic())
    return getRange(0);
  return VH::maximum(getRanges());
}

double CovAniso::getScale() const
{
  if (!hasRange()) return 0.;
  if (isIsotropic())
    return getScale(0);
  return VH::maximum(getScales());
}

VectorDouble CovAniso::getAnisoCoeffs() const
{
  VectorDouble coef = getRanges();
  double max = VH::maximum(coef);
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
double CovAniso::getParam() const
{
  if (!hasParam())
    return 0.;
  return _cova->getParam();
}

void CovAniso::_initFromContext()
{
  int ndim = getNDim();
  int nvar = getNVariables();
  _sill.resetFromValue(nvar, nvar, 1.);
  _aniso.init(ndim);
  _updateFromContext();
}

void CovAniso::_updateFromContext()
{
  computeMarkovCoeffs();
  _computeCorrec();
}

/**
 * Calculate the Integral Range in various Space Dimension (1, 2 or 3)
 * @return
 */
double CovAniso::getIntegralRange(int ndisc, double hmax) const
{
  int ndim = getNDim();
  SpacePoint dd(VectorDouble(ndim),-1);
  double delta = hmax / ndisc;
  double total = 0.;
  switch (ndim)
  {
    case 1:
      for (int j1 = -ndisc; j1 <= ndisc; j1++)
      {
        dd.setCoord(0, delta * j1);
        total += delta * eval(dd, SpacePoint());
      }
      break;

    case 2:
      for (int j1 = -ndisc; j1 <= ndisc; j1++)
        for (int j2 = -ndisc; j2 <= ndisc; j2++)
        {
          dd.setCoord(0 , delta * j1);
          dd.setCoord(1, delta * j2);
          total += delta * delta * eval(dd, SpacePoint());
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
            total += delta * delta * delta * eval(dd, SpacePoint());
          }
      break;

    default:
      my_throw("Integral Range has only been programmed for Space Dimension 1 to 3");
  }
  return total;
}

bool CovAniso::_isVariableValid(int ivar) const
{
  return checkArg("Rank of the Variable", ivar, getNVariables());
}

int CovAniso::getGradParamNumber() const
{
  int ndim = getNDim();
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

CovAniso* CovAniso::createIsotropic(const CovContext &ctxt,
                                    const ECov &type,
                                    double range,
                                    double sill,
                                    double param,
                                    bool flagRange)
{
  if (ctxt.getNVar() != 1)
  {
    messerr("This function is dedicated to the Monovariate case");
    return nullptr;
  }
  return new CovAniso(type, range, param, sill, ctxt, flagRange);
}

CovAniso* CovAniso::createAnisotropic(const CovContext &ctxt,
                                      const ECov &type,
                                      const VectorDouble &ranges,
                                      double sill,
                                      double param,
                                      const VectorDouble &angles,
                                      bool flagRange)
{
  if (ctxt.getNVar() != 1)
  {
    messerr("This function is dedicated to the Monovariate case");
    return nullptr;
  }
  int ndim = (int) ranges.size();
  if ((int) ctxt.getNDim() != ndim)
  {
    messerr("Mismatch in Space Dimension between 'ranges'(%d) and 'ctxt'(%d)",
            ndim, ctxt.getNDim());
    return nullptr;
  }

  CovAniso *cov = new CovAniso(type, ctxt);
  if (flagRange)
    cov->setRanges(ranges);
  else
    cov->setScales(ranges);
  cov->setSill(sill);
  cov->setParam(param);
  if (!angles.empty()) cov->setAnisoAngles(angles);
  return cov;
}

CovAniso* CovAniso::createIsotropicMulti(const CovContext &ctxt,
                                         const ECov &type,
                                         double range,
                                         const MatrixSquareSymmetric &sills,
                                         double param,
                                         bool flagRange)
{
  CovAniso *cov = new CovAniso(type, ctxt);
  int nvar = sills.getNSize();
  if (ctxt.getNVar() != nvar)
  {
    messerr(
        "Mismatch in the number of variables between 'sills'(%d) and 'ctxt'(%d)",
        nvar, ctxt.getNVar());
    return nullptr;
  }
  if (flagRange)
    cov->setRangeIsotropic(range);
  else
    cov->setScale(range);
  cov->setSill(sills);
  cov->setParam(param);
  return cov;
}

CovAniso* CovAniso::createAnisotropicMulti(const CovContext &ctxt,
                                           const ECov &type,
                                           const VectorDouble &ranges,
                                           const MatrixSquareSymmetric& sills,
                                           double param,
                                           const VectorDouble &angles,
                                           bool flagRange)
{

  int nvar = sills.getNSize();
  if (ctxt.getNVar() != nvar)
  {
    messerr(
        "Mismatch in the number of variables between 'sills'(%d) and 'ctxt'(%d)",
        nvar, ctxt.getNVar());
    return nullptr;
  }
  int ndim = (int) ranges.size();
  if ((int) ctxt.getNDim() != ndim)
  {
    messerr("Mismatch in Space Dimension between 'ranges'(%d) and 'ctxt'(%d)",
            ndim, ctxt.getNDim());
    return nullptr;
  }

  CovAniso *cov = new CovAniso(type, ctxt);
  if (flagRange)
    cov->setRanges(ranges);
  else
    cov->setScales(ranges);
  cov->setSill(sills);
  cov->setParam(param);
  if (!angles.empty()) cov->setAnisoAngles(angles);
  return cov;
}

void CovAniso::copyCovContext(const CovContext &ctxt)
{
  _ctxt.copyCovContext(ctxt);
  if (_cova != nullptr) _cova->copyCovContext(ctxt);
}

Array CovAniso::evalCovFFT(const VectorDouble& hmax,
                           int N,
                           int ivar,
                           int jvar) const
{
  if (! hasSpectrumOnRn()) return Array();

  std::function<double(const VectorDouble&)> funcSpectrum;
  funcSpectrum = [this, ivar, jvar](const VectorDouble &freq)
  {
    return evalSpectrum(freq, ivar, jvar) * _getDetTensor();
  };
  return evalCovFFTSpatial(hmax, N, funcSpectrum);
}

CovAniso* CovAniso::createReduce(const VectorInt &validVars) const
{
  CovAniso* newCovAniso = this->clone();

  // Modify the CovContext
  int nvar = (int) validVars.size();
  CovContext ctxt = CovContext(nvar);

  // Modify the Matrix of sills
  newCovAniso->setContext(ctxt);
  MatrixSquareSymmetric* newsill = dynamic_cast<MatrixSquareSymmetric*>(MatrixFactory::createReduce(&_sill, validVars, validVars));
  newCovAniso->setSill(*newsill);
  return newCovAniso;
}

/**
 * Define the second Space Point by transforming the input Space Point 'pt'
 * on the basis of the current covariance
 *
 * @param pt Target sample provided as a Space Point
 */
void CovAniso::optimizationSetTarget(const SpacePoint& pt) const
{
  if (_isOptimEnabled())
  {
     _optimizationTransformSP(pt, _p2A);
  }
}

/**
 * Define the Second Space Point as coinciding with the Input Space Point 'iech'.
 * Note that, as the Input Space Points are already transformed in the basis
 * of the current structure, it is just an assignment.
 *
 * @param iech Rank of the sample among the recorded Space Points
 */
void CovAniso::optimizationSetTarget(int iech) const
{
  _p2A = _p1As[iech];
}

/**
 * Transform a space point using the anisotropy tensor
 * @param ptin  Input Space Point
 * @param ptout Output Space Point
 */
void CovAniso::_optimizationTransformSP(const SpacePoint& ptin, SpacePoint& ptout) const
{
  if (_isOptimEnabled())
	  _aniso.applyInverseInPlace(ptin.getCoords(), ptout.getCoordRef());
}

/**
 * Transform a set of Space Points using the anisotropy tensor
 * The set of resulting Space Points are stored as private member of this.
 * Note that ALL samples are processed, independently from the presence of a selection
 * or checking for heterotopy.
 * @param db Input Db
 */
void CovAniso::optimizationPreProcess(const Db* db) const
{
  if (!_isOptimEnabled()) return;
  if (isOptimizationInitialized(db)) return;
  std::vector<SpacePoint> p1s;

  db->getSamplesAsSP(p1s,_space);
  int n = (int) p1s.size();
	_p1As.resize(n);
	for(int i = 0; i < n ; i++)
	{
		_p1As[i] = SpacePoint(_space);
		if (! p1s[i].isFFFF())
		  _optimizationTransformSP(p1s[i], _p1As[i]);
		else
		  _p1As[i].setFFFF();
	}
  _p2A = SpacePoint(_space);
}

void CovAniso::optimizationPostProcess() const
{
  if (! isOptimizationInitialized()) return;
	_p1As.clear();
}

/**
 * Checks that the Optimization has already been initiated, by:
 * - checking that the storage (for Sample Points projected in the Covariance
 * rotation system) is already allocated
 * - checking that the dimension of this storage is correct (only if 'db' is provided):
 * in particular, this check is not necessary when freeing this storage.
 */
bool CovAniso::isOptimizationInitialized(const Db* db) const
{
  if (_p1As.empty()) return false;
  if (db == nullptr) return true;
  int n = (int) _p1As.size();
  return n == db->getSampleNumber();
}

double scale2range(const ECov &type, double scale, double param)
{
  CovContext ctxt = CovContext(1, 1);
  ACovFunc *cova = CovFactory::createCovFunc(type, ctxt);
  cova->setParam(param);
  double scadef = cova->getScadef();
  return scale * scadef;
}

double range2scale(const ECov &type, double range, double param)
{
  CovContext ctxt = CovContext(1, 1);
  ACovFunc *cova = CovFactory::createCovFunc(type, ctxt);
  cova->setParam(param);
  double scadef = cova->getScadef();
  return range / scadef;
}

  
// Set of functions to make parameters no stationary (or to make them back stationary).
// There is to types of non stationarities : NoStatDb in which the parameters are read in a
// DbGrid or NoStatFunctional for which you have to provide a function of the coordinates.
// Each parameter can have its own type of No stationarity and its own DbGrid in case
// of NoStatDb. 
// For specifying the NoStat DbGrid, you can first attach it by using attachNoStatDb.
// If not, you have to specify the DbGrid when you make the first parameter non stationary.

void CovAniso::attachNoStatDb(const Db* db)
{
  _tabNoStat.setDbNoStatRef(db);
}

bool CovAniso::_checkAndManageNoStatDb(const Db*&  db, const String& namecol)
{
 if (_tabNoStat.getDbNoStatRef() == nullptr && db == nullptr)
 {
  messerr("You have to define a Db (with attachNoStatDb or by specifying a Db here)");  
  return false;
 }
  _setNoStatDbIfNecessary(db);

 if (db->getUID(namecol)< 0)
 {
    messerr("You have to specified a name of a column of the reference Db");
    return false;
 }
 return true;
}

void CovAniso::_setNoStatDbIfNecessary(const Db*& db)
{
  if (_tabNoStat.getDbNoStatRef() == nullptr)
    attachNoStatDb(db);
  if (db == nullptr)
    db = _tabNoStat.getDbNoStatRef();
}

void CovAniso::_makeElemNoStat(const EConsElem &econs, int iv1, int iv2,const AFunctional* func, const Db* db, const String& namecol)
{
  std::shared_ptr<ANoStat> ns;
  if (func == nullptr)
  {
    if(!_checkAndManageNoStatDb(db,namecol)) return;
    ns = std::shared_ptr<ANoStat>(new NoStatArray(db,namecol));
  }
  else 
  {
    ns = std::unique_ptr<ANoStat>(new NoStatFunctional(func));
  }
   _tabNoStat.addElem(ns, econs,iv1,iv2);
  
}
///////////////////// Range ////////////////////////
void CovAniso::makeRangeNoStatDb(const String &namecol, int idim, const Db* db)
{   
  if(!_checkTensor()) return;
  _makeElemNoStat(EConsElem::RANGE, idim, 0,nullptr,db, namecol);
 
}

void CovAniso::makeRangeNoStatFunctional(const AFunctional *func, int idim)
{
  if(!_checkTensor()) return;
  _makeElemNoStat(EConsElem::RANGE, idim, 0,func);
}


void CovAniso::makeRangeStationary(int idim)
{
  if(_tabNoStat.removeElem(EConsElem::RANGE, idim) == 0 && 
     _tabNoStat.removeElem(EConsElem::SCALE, idim) == 0)
  {
    messerr("This parameter was already stationary!");
  }
}

///////////////////// Scale ////////////////////////

void CovAniso::makeScaleNoStatDb(const String &namecol, int idim, const Db* db)
{   
  if(!_checkTensor()) return;
  _makeElemNoStat(EConsElem::SCALE, idim, 0,nullptr,db, namecol);
}


void CovAniso::makeScaleNoStatFunctional(const AFunctional *func, int idim)
{
  if(!_checkTensor()) return;
  _makeElemNoStat(EConsElem::SCALE, idim, 0,func);
}

void CovAniso::makeScaleStationary(int idim)
{
  makeRangeStationary(idim);
}

///////////////////// Angle ////////////////////////


void CovAniso::makeAngleNoStatDb(const String &namecol, int idim, const Db* db)
{
  if(!_checkTensor()) return;
  _makeElemNoStat(EConsElem::ANGLE, idim, 0,nullptr,db, namecol);

}

void CovAniso::makeAngleNoStatFunctional(const AFunctional *func, int idim)
{
  if(!_checkTensor()) return;
  _makeElemNoStat(EConsElem::ANGLE, idim, 0,func);

}

void CovAniso::makeAngleStationary(int idim)
{
  if(_tabNoStat.removeElem(EConsElem::ANGLE, idim) == 0)
  {
    messerr("This parameter was already stationary!");
  }
}
///////////////////// Tensor ////////////////////////


void CovAniso::makeTensorNoStatDb(const String &namecol, int idim, int jdim,const Db* db)
{
  if(!_checkRotation()) return;
  if (!_checkDims(idim,jdim)) return;
  _makeElemNoStat(EConsElem::TENSOR, idim, jdim,nullptr,db, namecol);

}

void CovAniso::makeTensorNoStatFunctional(const AFunctional  *func, int idim, int jdim)
{
    if(!_checkRotation()) return;
    if (!_checkDims(idim,jdim)) return;
     _makeElemNoStat(EConsElem::TENSOR, idim, jdim,func);

}

void CovAniso::makeTensorStationary(int idim, int jdim)
{
  if (!_checkDims(idim,jdim)) return;
  if(_tabNoStat.removeElem(EConsElem::TENSOR, idim,jdim) == 0)
  {
    messerr("This parameter was already stationary!");
  }
}
///////////////////// Sill ////////////////////////

void CovAniso::makeSillNoStatDb(const String &namecol, int ivar, int jvar,const Db* db)
{
  if (!_checkSill(ivar,jvar)) return;
  _makeElemNoStat(EConsElem::SILL, ivar, jvar,nullptr,db, namecol);

}

void CovAniso::makeSillNoStatFunctional(const AFunctional  *func, int ivar, int jvar)
{
  if (!_checkSill(ivar,jvar)) return;
  _makeElemNoStat(EConsElem::SILL, ivar, jvar,func);

}
  
void CovAniso::makeSillStationary(int ivar, int jvar)
{
  if (!_checkSill(ivar,jvar)) return;
  if(_tabNoStat.removeElem(EConsElem::SILL, ivar,jvar) == 0)
  {
    messerr("This parameter was already stationary!");
  }
}

///////////////////// Param ////////////////////////

void CovAniso::makeParamNoStatDb(const String &namecol, const Db* db)
{
  if(!_checkParam()) return;
   _makeElemNoStat(EConsElem::PARAM, 0, 0,nullptr,db, namecol);

}

void CovAniso::makeParamNoStatFunctional(const AFunctional *func)
{
  if(!_checkParam()) return;
  _makeElemNoStat(EConsElem::PARAM, 0, 0,func);

}

void CovAniso::makeParamStationary()
{
  if(!_checkParam()) return;
  if(_tabNoStat.removeElem(EConsElem::PARAM) == 0)
  {
    messerr("This parameter was already stationary!");
  }
}

/////////////////////////// Check functions ////////////////////:

bool CovAniso::_checkTensor() const
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

bool CovAniso::_checkRotation() const
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

bool CovAniso::_checkSill(int ivar, int jvar) const
{
  int nvar = getNVariables();
  if ((ivar > nvar) || (jvar > nvar))
  {
    messerr("Your model has only %d variables.",nvar);
    return false;
  }
  return true;
}

bool CovAniso::_checkDims(int idim, int jdim) const
{
  int ndim = getNDim();
  if ((idim > ndim) || (jdim > ndim))
  {
    messerr("Your model is only in dimension %d.",ndim);
    return false;
  }
  return true;
}

bool CovAniso::_checkParam() const
{
  if (getType()!= ECov::MATERN)
  {
    messerr("This covariance function has no parameters of this type");
    return false;
  }
  return true;
}

/////////////  Functions to attach no stat information on various supports ////////
void CovAniso::informMeshByMesh(const AMesh* amesh) const
{
  _tabNoStat.informMeshByMesh(amesh);
}
void CovAniso::informMeshByApex(const AMesh* amesh) const
{
  _tabNoStat.informMeshByMesh(amesh);
}
void CovAniso::informDbIn(const Db* dbin) const
{
  _tabNoStat.informDbIn(dbin);
}
void CovAniso::informDbOut(const Db* dbout) const
{
  _tabNoStat.informDbOut(dbout);
}

double CovAniso::getValue(const EConsElem &econs,int iv1,int iv2) const
{
  if (econs == EConsElem::RANGE)
    return getRange(iv1);
  if (econs == EConsElem::SCALE)
    return getScale(iv1);
  if (econs == EConsElem::ANGLE)
    return getAnisoAngles()[iv1];
  if (econs == EConsElem::PARAM)
    return getParam();
  if (econs == EConsElem::SILL)
    return getSill(iv1,iv2);
  return TEST;
}

VectorDouble CovAniso::informCoords(const VectorVectorDouble& coords, 
                                    const EConsElem& econs,
                                    int iv1,
                                    int iv2) const
{
  VectorDouble result(coords[0].size(),getValue(econs,iv1,iv2));
  _tabNoStat.informCoords(coords,econs,iv1,iv2,result);
  return result;
}


void CovAniso::informMeshByMeshForAnisotropy(const AMesh* amesh) const
{
  for (const auto &e : _listaniso)
  {
    _tabNoStat.informMeshByMesh(amesh,e);
  }
}

void CovAniso::informMeshByApexForAnisotropy(const AMesh* amesh) const
{
   for (const auto &e : _listaniso)
      _tabNoStat.informMeshByMesh(amesh,e);
}

void CovAniso::informDbInForAnisotropy(const Db* dbin) const
{
   for (const auto &e :_listaniso)
      _tabNoStat.informDbIn(dbin,e);

}
void CovAniso::informDbOutForAnisotropy(const Db* dbout) const
{
   for (const auto &e: _listaniso)
      _tabNoStat.informDbOut(dbout,e);
}

void CovAniso::informMeshByMeshForSills(const AMesh* amesh) const
{
   _tabNoStat.informMeshByMesh(amesh,EConsElem::SILL);
}

void CovAniso::informMeshByApexForSills(const AMesh* amesh) const
{
   _tabNoStat.informMeshByApex(amesh,EConsElem::SILL);
}

void CovAniso::informDbInForSills(const Db* dbin) const
{
   _tabNoStat.informDbIn(dbin,EConsElem::SILL);
}

void CovAniso::informDbOutForSills(const Db* dbout) const
{
  _tabNoStat.informDbOut(dbout,EConsElem::SILL);
}



/**
 * Update the Model according to the Non-stationary parameters
 * @param icas1 Type of first Db: 1 for Input; 2 for Output
 * @param iech1 Rank of the target within Db1 (or -1)
 * @param icas2 Type of first Db: 1 for Input; 2 for Output
 * @param iech2 Rank of the target within Dbout (or -2)
 */
void CovAniso::updateCovByPoints(int icas1, int iech1, int icas2, int iech2) 
{
  // If no non-stationary parameter is defined, simply skip
  if (! isNoStat()) return;
  double val1, val2;

  
  int ndim = getNDim();

  const auto paramsnostat = _tabNoStat.getTable();
  // Loop on the elements that can be updated one-by-one

  for (const auto &e : paramsnostat)
  {
    EConsElem type = e.first.getType();
    e.second->getValuesOnDb( icas1, iech1, &val1, icas2, iech2, &val2);

    if (type == EConsElem::SILL)
    {
      int iv1 = e.first.getIV1();
      int iv2 = e.first.getIV2();
      setSill(iv1, iv2, sqrt(val1 * val2));
    }
    else if (type == EConsElem::PARAM)
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
      if (_tabNoStat.isElemDefined(EConsElem::ANGLE, idim))
      {
        auto noStat = _tabNoStat.getElem(EConsElem::ANGLE, idim);
        flagRotOne = true;                       
        if (noStat->getValuesOnDb(icas1, iech1,&angle1[idim], icas2, iech2,&angle2[idim]))
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
      if (_tabNoStat.isElemDefined(EConsElem::SCALE, idim))
      {
        auto noStat = _tabNoStat.getElem(EConsElem::SCALE, idim);
        flagScaleOne = true;
        if (noStat->getValuesOnDb( icas1, iech1, &scale1[idim], icas2, iech2, &scale2[idim]))
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
      if (_tabNoStat.isElemDefined(EConsElem::RANGE, idim))
      {
        auto noStat = _tabNoStat.getElem(EConsElem::RANGE, idim);
        flagRangeOne = true;
        if (noStat->getValuesOnDb(icas1, iech1,&range1[idim], icas2, iech2,  &range2[idim]))
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
    double det1 = pow(direct1.determinant(), 0.25);

  // Extract the direct tensor at second point and square it
    setRotationAnglesAndRadius(angle2, range2, scale2);
    MatrixSquareSymmetric direct2 = getAniso().getTensorDirect2();
    double det2 = pow(direct2.determinant(), 0.25);

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


void CovAniso::updateCovByMesh(int imesh,bool aniso)
{
  // If no non-stationary parameter is defined, simply skip
  if (! isNoStat()) return;
  int ndim = getNDim();

  // Loop on the elements that can be updated one-by-one
  if (!aniso)
  {
    const auto paramsnostat = _tabNoStat.getTable();
    for (const auto &e : paramsnostat)
    {
      EConsElem type = e.first.getType();
      if (type == EConsElem::SILL)
      {
        double sill = e.second->getValueOnMeshByApex(imesh);
        int iv1 = e.first.getIV1();
        int iv2 = e.first.getIV2();
        setSill(iv1, iv2, sill);
      }
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
      if (_tabNoStat.isElemDefined(EConsElem::ANGLE, idim))
      {
        auto noStat = _tabNoStat.getElem(EConsElem::ANGLE, idim);
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
      if (_tabNoStat.isElemDefined(EConsElem::SCALE, idim))
      {
        auto noStat = _tabNoStat.getElem(EConsElem::SCALE, idim);
        scales[idim] = noStat->getValueOnMeshByMesh(imesh);
      }
    }
  }
  
  if (getNRanges() > 0)
  {
    ranges = getRanges();

    for (int idim = 0; idim < ndim; idim++)
    {
      if (_tabNoStat.isElemDefined(EConsElem::RANGE, idim))
      {
        auto noStat = _tabNoStat.getElem(EConsElem::RANGE, idim);
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
        if (_tabNoStat.isElemDefined(EConsElem::TENSOR, idim,jdim))
        {
          auto noStat = _tabNoStat.getElem(EConsElem::TENSOR, idim, jdim);
          
      }
    }
  }
}

void CovAniso::makeStationary()
{
  _tabNoStat = TabNoStatCovAniso();
}

void CovAniso::_manage(const Db* db1,const Db* db2) const
{
  if (db1!=nullptr)
    informDbIn(db1);
  if (db2!=nullptr)
    informDbOut(db2);
}