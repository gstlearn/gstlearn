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
#include "Covariances/ACov.hpp"
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
#include "Basic/VectorNumT.hpp"
#include "Basic/FFT.hpp"
#include "Space/ASpace.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/SpaceSN.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "geoslib_define.h"
#include <math.h>
#include <functional>
#include <memory>
#include <ostream>
#include <vector>


CovAniso::CovAniso(const ECov &type, const CovContext &ctxt)
    : ACov(ctxt.getSpace()), /// TODO : shared pointer
      _cor(type, ctxt),
      _ctxt(ctxt),
      _sill(),
      _tabNoStat(),
      _optimEnabled(true)
{
  _initFromContext();
}

CovAniso::CovAniso(const String &symbol, const CovContext &ctxt)
    : ACov(ctxt.getSpace()), /// TODO : shared pointer
      _cor(symbol, ctxt),
      _ctxt(ctxt),
      _sill(),
      _tabNoStat(),
      _optimEnabled(true)
{
  ECov covtype = CovFactory::identifyCovariance(symbol, ctxt);
  _initFromContext();
}

CovAniso::CovAniso(const ECov &type,
                   double range,
                   double param,
                   double sill,
                   const CovContext &ctxt,
                   bool flagRange)
    : ACov(ctxt.getSpace()), /// TODO : shared pointer
      _cor(type, range,param, ctxt, flagRange),
      _ctxt(ctxt),
      _sill(),
      _tabNoStat(),
      _optimEnabled(true)
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
      _cor(r._cor),
      _ctxt(r._ctxt),
      _sill(r._sill),
      _tabNoStat(r._tabNoStat),
      _optimEnabled(r._optimEnabled)
{
}

CovAniso& CovAniso::operator=(const CovAniso &r)
{
  if (this != &r)
  {
    ACov::operator =(r);
    _cor = r._cor;
    _ctxt = r._ctxt;
    _sill = r._sill;
    _tabNoStat = r._tabNoStat;
    _optimEnabled = r._optimEnabled;
  }
  return *this;
}

CovAniso::~CovAniso()
{

}

void CovAniso::_computeCorrec()
{
  _cor.computeCorrec();
}

void CovAniso::computeMarkovCoeffs()
{
  _cor.computeMarkovCoeffs();
}

void CovAniso::setContext(const CovContext &ctxt)
{
  _ctxt = ctxt;
  _updateFromContext();
}

void CovAniso::setParam(double param)
{
  _cor.setParam(param);
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
  _cor.setRangeIsotropic(range);
}

void CovAniso::setRanges(const VectorDouble &ranges)
{
  _cor.setRanges(ranges);
}

void CovAniso::setRange(int idim, double range)
{
  _cor.setRange(idim, range);
}

void CovAniso::setScale(double scale)
{
  _cor.setScale(scale);
}

void CovAniso::setScales(const VectorDouble &scales)
{
  _cor.setScales(scales);
}

void CovAniso::setScale(int idim, double scale)
{
  _cor.setScale(idim, scale);
}

void CovAniso::setAnisoRotation(const Rotation &rot)
{
  _cor.setAnisoRotation(rot);
}

void CovAniso::setAnisoRotation(const VectorDouble &rot)
{
  _cor.setAnisoRotation(rot);
}

void CovAniso::setAnisoAngles(const VectorDouble &angles)
{
  _cor.setAnisoAngles(angles);
}

void CovAniso::setAnisoAngle(int idim, double angle)
{
  _cor.setAnisoAngle(idim, angle);
}

void CovAniso::setRotationAnglesAndRadius(const VectorDouble &angles,
                                          const VectorDouble &ranges,
                                          const VectorDouble &scales)
{
  _cor.setRotationAnglesAndRadius(angles, ranges, scales);
}

bool CovAniso::isValidForTurningBand() const
{
  return _cor.isValidForTurningBand();
}
double CovAniso::simulateTurningBand(double t0, TurningBandOperate &operTB) const
{
  return _cor.simulateTurningBand(t0, operTB);
}
bool CovAniso::isValidForSpectral() const
{
  return _cor.isValidForSpectral();
}
MatrixRectangular CovAniso::simulateSpectralOmega(int nb) const
{
  return _cor.simulateSpectralOmega(nb);
}
bool CovAniso::isConsistent(const ASpace* space) const
{
  return _cor.isConsistent(space);
}



double CovAniso::eval0(int ivar, int jvar, const CovCalcMode* mode) const
{
  double cov = _cor.evalCorFromH(0, mode);

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
  double cov = _cor.evalCor(p1,p2,mode);
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
void CovAniso::addEval0CovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                            const CovCalcMode *mode) const
{
  double cov = _cor.evalCorFromH(0, mode); 

  if (mode == nullptr || ! mode->getUnitary())
    mat.addMatInPlace(_sill, 1., cov);
  else
  {
    MatrixSquareGeneral identity = _sill;
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
void CovAniso::_addEvalCovMatBiPointInPlace(MatrixSquareGeneral &mat,
                                          const SpacePoint &p1,
                                          const SpacePoint &p2,
                                          const CovCalcMode *mode) const
{
  
  double cor = _cor.evalCor(p1,p2,mode);

  if (mode == nullptr || ! mode->getUnitary())
    mat.addMatInPlace(_sill, 1., cor);
  else
  {
    MatrixSquareGeneral identity = _sill;
    identity.setIdentity();
    mat.addMatInPlace(identity, 1., cor);
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
        cov = _cor.evalCorFromH(hoptim, mode);
        res.updValue(irow, icol, EOperator::ADD, sill * cov);
      }
      irow++;
    }
  }
}

void CovAniso::_evalOptim(SpacePoint* p1A, SpacePoint* p2A,
                          MatrixSquareGeneral &mat,
                          const CovCalcMode *mode) const
{
  // Calculate covariance between two points
  double hoptim = p2A->getDistance(*p1A);
  double cov = _cor.evalCorFromH(hoptim, mode);

  if (mode == nullptr || ! mode->getUnitary())
    mat.addMatInPlace(_sill, 1., cov);
  else
  {
    MatrixSquareGeneral identity = _sill;
    identity.setIdentity();
    mat.addMatInPlace(identity, 1., cov);
  }
}


double CovAniso::evalCovOnSphere(double alpha,
                                 int degree,
                                 bool flagScaleDistance,
                                 const CovCalcMode* mode) const
{
  double value = _cor.evalCovOnSphere(alpha, degree, flagScaleDistance,mode);

  if (mode == nullptr || ! mode->getUnitary())
    value *= getSill(0,0);

  return value;
}

VectorDouble CovAniso::evalSpectrumOnSphere(int n, bool flagNormDistance, bool flagCumul) const
{
  return _cor.evalSpectrumOnSphere(n, flagNormDistance, flagCumul);
}

void CovAniso::setMarkovCoeffs(const VectorDouble& coeffs)
{
  _cor.setMarkovCoeffs(coeffs);
}

/* This function computes a polynomial P from two polynomials P1 and P2 and a small constant eps
 * P(x) = P1(x)^2 + x * P2(x)^2 + eps
 */
void CovAniso::setMarkovCoeffsBySquaredPolynomials(VectorDouble coeffs1,
                                                   VectorDouble coeffs2,
                                                   double eps)
{
  _cor.setMarkovCoeffsBySquaredPolynomials(coeffs1, coeffs2, eps);
}

double CovAniso::getCorrec() const
{
  return _cor.getCorrec();
}

double CovAniso::getFullCorrec() const
{
  return  _cor.getFullCorrec();
}

double CovAniso::_getDetTensor() const
{
  return _cor.getDetTensor();
}

double CovAniso::evalSpectrum(const VectorDouble& freq, int ivar, int jvar) const
{
  if (!_cor.hasSpectrumOnRn()) return TEST;
  return _sill.getValue(ivar, jvar) * _cor.evalSpectrum(freq, ivar, jvar);
}

double CovAniso::normalizeOnSphere(int n) const
{ 
  return _cor.normalizeOnSphere(n);
}

VectorDouble CovAniso::getMarkovCoeffs() const
{
  return _cor.getMarkovCoeffs();
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

String CovAniso::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
 
  sstr << _cor.getCova()->toString();

  // Sill - Factor / Slope information
  if (_cor.hasRange() > 0)
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

  }
  else if (_cor.hasRange() < 0)
  {
    // The sill is not defined: use slope instead

    if (getNVariables() > 1)
    {
      MatrixSquareGeneral slopes = _sill;
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

   // Covariance Parameters
  sstr << _cor.toStringParams(strfmt);

  // Non-stationary parameters

  if (isNoStat())
  { 
    sstr << toTitle(1, "Non-Stationary Parameters");
    sstr << _tabNoStat.toString(strfmt);
    int i = _tabNoStat.getNSills();
    sstr << _cor.toStringNoStat(strfmt,i);
  }
  return sstr.str();
}

double CovAniso::getSill(int ivar, int jvar) const
{
  return _sill.getValue(ivar, jvar);
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
  return _sill.getValue(ivar, jvar) / range;
}

VectorDouble CovAniso::getRanges() const
{
  return _cor.getRanges();
}

void CovAniso::setType(const ECov &type)
{
 _cor.setType(type);
}

/**
 * This function returns the range in the isotropic case
 * In the anisotropic case, it returns the largest range over all directions
 * @return
 */
double CovAniso::getRange() const
{
 return _cor.getRange();
}

double CovAniso::getScale() const
{
 return _cor.getScale();
}

VectorDouble CovAniso::getAnisoCoeffs() const
{
  return _cor.getAnisoCoeffs();
}

/**
 * For compatibility, this function returns 0 if the Covariance has no Third Parameter
 *
 * @return Third parameter
 */
double CovAniso::getParam() const
{
  return _cor.getParam();
}

void CovAniso::_initFromContext()
{
  _cor.initFromContext();
  _sill.reset(_ctxt.getNVar(), _ctxt.getNVar());
}

void CovAniso::_updateFromContext()
{
  _cor.updateFromContext();
}

/**
 * Calculate the Integral Range in various Space Dimension (1, 2 or 3)
 * @return
 */
double CovAniso::getIntegralRange(int ndisc, double hmax) const
{
  return _sill.getValue(0, 0) * _cor.getIntegralRange(ndisc, hmax);
}

bool CovAniso::_isVariableValid(int ivar) const
{
  return checkArg("Rank of the Variable", ivar, getNVariables());
}

int CovAniso::getGradParamNumber() const
{
  return _cor.getGradParamNumber();
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
  _cor.copyCovContext(ctxt);
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
void CovAniso::_optimizationSetTarget(const SpacePoint& pt) const
{
  if (_isOptimEnabled())
  {  
    _optimizationTransformSP(pt, _p2A);
  }
  else 
  {
    _p2A = pt;
  }  
}

/**
 * Define the Second Space Point as coinciding with the Input Space Point 'iech'.
 * Note that, as the Input Space Points are already transformed in the basis
 * of the current structure, it is just an assignment.
 *
 * @param iech Rank of the sample among the recorded Space Points
 */
void CovAniso::optimizationSetTargetByIndex(int iech) const
{
  if (_isOptimPreProcessed)
  {
    _p2A = _p1As[iech];
    _p2A.setTarget(true);
  }
}

/**
 * Transform a space point using the anisotropy tensor
 * @param ptin  Input Space Point
 * @param ptout Output Space Point
 */
void CovAniso::_optimizationTransformSP(const SpacePoint& ptin, SpacePoint& ptout) const
{
  _cor.optimizationTransformSP(ptin, ptout);
}

void CovAniso::_optimizationPostProcess() const 
{
  _cor.optimizationPostProcess();
}

/**
 * Transform a set of Space Points using the anisotropy tensor
 * The set of resulting Space Points are stored as private member of this.
 * Note that ALL samples are processed, independently from the presence of a selection
 * or checking for heterotopy.
 * @param p vector of SpacePoints
 */
void CovAniso::_optimizationPreProcess(const std::vector<SpacePoint>& p) const
{
  if (!isOptimEnabled())
  {
     ACov::_optimizationPreProcess(p);
     return;
  }
  _cor.optimizationPreProcess(p,_p1As);
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
  _cor.attachNoStatDb(db);
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
  if (func == nullptr)
  {
    if(!_checkAndManageNoStatDb(db,namecol)) return;
  }

  if (econs != EConsElem::SILL)
  {
    _cor.makeElemNoStat(econs,iv1,iv2,func,db,namecol);
    return;
  }
  

  std::shared_ptr<ANoStat> ns;
  if (func == nullptr)
  {
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
  if(_checkAndManageNoStatDb(db,namecol))
  {
    _cor.makeRangeNoStatDb(namecol,idim,db);
  }
 
}

void CovAniso::makeRangeNoStatFunctional(const AFunctional *func, int idim)
{
  _cor.makeRangeNoStatFunctional(func,idim);
}


void CovAniso::makeRangeStationary(int idim)
{
  _cor.makeRangeStationary(idim);
}

///////////////////// Scale ////////////////////////

void CovAniso::makeScaleNoStatDb(const String &namecol, int idim, const Db* db)
{   
  if(_checkAndManageNoStatDb(db,namecol))
  {
    _cor.makeScaleNoStatDb(namecol,idim,db);
  }
}


void CovAniso::makeScaleNoStatFunctional(const AFunctional *func, int idim)
{
  _cor.makeScaleNoStatFunctional(func,idim);
}

void CovAniso::makeScaleStationary(int idim)
{
  _cor.makeScaleStationary(idim);
}

///////////////////// Angle ////////////////////////


void CovAniso::makeAngleNoStatDb(const String &namecol, int idim, const Db* db)
{
  if(_checkAndManageNoStatDb(db,namecol))
  {
   _cor.makeAngleNoStatDb(namecol,idim,db);
  }
}

void CovAniso::makeAngleNoStatFunctional(const AFunctional *func, int idim)
{
  _cor.makeAngleNoStatFunctional(func,idim);
}

void CovAniso::makeAngleStationary(int idim)
{
  _cor.makeAngleStationary(idim);
}
///////////////////// Tensor ////////////////////////


void CovAniso::makeTensorNoStatDb(const String &namecol, int idim, int jdim,const Db* db)
{
 if(_checkAndManageNoStatDb(db,namecol))
 {
    _cor.makeTensorNoStatDb(namecol,idim,jdim,db);
 } 
}

void CovAniso::makeTensorNoStatFunctional(const AFunctional  *func, int idim, int jdim)
{
   _cor.makeTensorNoStatFunctional(func,idim,jdim);
}

void CovAniso::makeTensorStationary(int idim, int jdim)
{
  _cor.makeTensorStationary(idim,jdim);
}
///////////////////// Sill ////////////////////////

void CovAniso::makeSillNoStatDb(const String &namecol, int ivar, int jvar,const Db* db)
{
  if (!_checkSill(ivar,jvar)) return;
  _makeElemNoStat(EConsElem::SILL, ivar, jvar,nullptr,db, namecol);
  _cor.checkAndManageNoStatDb(db,namecol);
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
 if(_checkAndManageNoStatDb(db,namecol))
 {
    _cor.makeParamNoStatDb(namecol,db);
 } 

}

void CovAniso::makeParamNoStatFunctional(const AFunctional *func)
{
  _cor.makeParamNoStatFunctional(func);

}

void CovAniso::makeParamStationary()
{
  _cor.makeParamStationary();
}

/////////////////////////// Check functions ////////////////////:

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


/////////////  Functions to attach no stat information on various supports ////////
void CovAniso::informMeshByMesh(const AMesh* amesh) const
{
  _tabNoStat.informMeshByMesh(amesh);
  _cor.informMeshByMesh(amesh);
}
void CovAniso::informMeshByApex(const AMesh* amesh) const
{
  _tabNoStat.informMeshByMesh(amesh);
  _cor.informMeshByMesh(amesh);
}
void CovAniso::informDbIn(const Db* dbin) const
{
  _tabNoStat.informDbIn(dbin);
  _cor.informDbIn(dbin);
}
void CovAniso::informDbOut(const Db* dbout) const
{
  _tabNoStat.informDbOut(dbout);
  _cor.informDbOut(dbout);
}

double CovAniso::getValue(const EConsElem &econs,int iv1,int iv2) const
{
  double val = _cor.getValue(econs,iv1,iv2);
  if (val == TEST)
  {
    if (econs == EConsElem::SILL)
    return getSill(iv1,iv2);
  }
  return val;
}

VectorDouble CovAniso::informCoords(const VectorVectorDouble& coords, 
                                    const EConsElem& econs,
                                    int iv1,
                                    int iv2) const
{
  if (econs == EConsElem::SILL)
  {
    VectorDouble result(coords[0].size(),getValue(econs,iv1,iv2));
    _tabNoStat.informCoords(coords,econs,iv1,iv2,result);
    return result;
  }
 
  return _cor.informCoords(coords,econs,iv1,iv2);
 
}


void CovAniso::informMeshByMeshForAnisotropy(const AMesh* amesh) const
{
  _cor.informMeshByMeshForAnisotropy(amesh);

}

void CovAniso::informMeshByApexForAnisotropy(const AMesh* amesh) const
{
   _cor.informMeshByApexForAnisotropy(amesh);
}

void CovAniso::informDbInForAnisotropy(const Db* dbin) const
{
  _cor.informDbInForAnisotropy(dbin);

}
void CovAniso::informDbOutForAnisotropy(const Db* dbout) const
{
   _cor.informDbOutForAnisotropy(dbout);  
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
  }
  _cor.updateCovByPoints(icas1, iech1 , icas2, iech2); 
}


void CovAniso::updateCovByMesh(int imesh,bool aniso)
{
  // If no non-stationary parameter is defined, simply skip
  if (! isNoStat()) return;

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
  }
 _cor.updateCovByMesh(imesh,aniso);
}

void CovAniso::makeStationary()
{
  _tabNoStat = TabNoStatCovAniso();
  _cor.makeStationary();
}

void CovAniso::_manage(const Db* db1,const Db* db2) const
{
  if (db1!=nullptr)
    informDbIn(db1);
  if (db2!=nullptr)
    informDbOut(db2);
  _cor.manage(db1,db2);
}