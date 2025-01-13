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
#include "Covariances/CorAniso.hpp"
#include "Covariances/CovProportional.hpp"
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
#include <ostream>
#include <vector>


CovAniso::CovAniso(const ECov &type, const CovContext &ctxt)
    : CovProportional(nullptr, MatrixSquareSymmetric(ctxt.getNVar())), /// TODO : shared pointer
      _corAniso(new CorAniso(type, ctxt)),
      _optimEnabled(true)
{
  CovProportional::setCor(_corAniso);
  _initFromContext();
}

CovAniso::CovAniso(const String &symbol, const CovContext &ctxt)
    :CovProportional(new CorAniso(symbol, ctxt)), /// TODO : shared pointer
      _corAniso((CorAniso*)getCor()),
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
    : CovProportional(new CorAniso(type, range,param, ctxt, flagRange), MatrixSquareSymmetric(ctxt.getNVar())), /// TODO : shared pointer
      _corAniso((CorAniso*)getCor()),
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
    : CovProportional(new CorAniso(*r._corAniso),r._sill), /// TODO : shared pointer
      _corAniso((CorAniso*)getCor()),
      _optimEnabled(r._optimEnabled)
{
  _ctxt.setNVar(r.getNVariables());
}

CovAniso& CovAniso::operator=(const CovAniso &r)
{
  if (this != &r)
  { 
     setCor(new CorAniso(*r._corAniso));
    _corAniso = (CorAniso*)getCor();
    _ctxt = r._ctxt;
    _sill = r._sill;
    _optimEnabled = r._optimEnabled;
  }
  return *this;
}

CovAniso::~CovAniso()
{
  delete _corAniso;
}

void CovAniso::_computeCorrec()
{
  _corAniso->computeCorrec();
}

void CovAniso::computeMarkovCoeffs()
{
  _corAniso->computeMarkovCoeffs();
}

void CovAniso::setParam(double param)
{
  _corAniso->setParam(param);
}


void CovAniso::setRangeIsotropic(double range)
{
  _corAniso->setRangeIsotropic(range);
}

void CovAniso::setRanges(const VectorDouble &ranges)
{
  _corAniso->setRanges(ranges);
}

void CovAniso::setRange(int idim, double range)
{
  _corAniso->setRange(idim, range);
}

void CovAniso::setScale(double scale)
{
  _corAniso->setScale(scale);
}

void CovAniso::setScales(const VectorDouble &scales)
{
  _corAniso->setScales(scales);
}

void CovAniso::setScale(int idim, double scale)
{
  _corAniso->setScale(idim, scale);
}

void CovAniso::setAnisoRotation(const Rotation &rot)
{
  _corAniso->setAnisoRotation(rot);
}

void CovAniso::setAnisoRotation(const VectorDouble &rot)
{
  _corAniso->setAnisoRotation(rot);
}

void CovAniso::setAnisoAngles(const VectorDouble &angles)
{
  _corAniso->setAnisoAngles(angles);
}

void CovAniso::setAnisoAngle(int idim, double angle)
{
  _corAniso->setAnisoAngle(idim, angle);
}

void CovAniso::setRotationAnglesAndRadius(const VectorDouble &angles,
                                          const VectorDouble &ranges,
                                          const VectorDouble &scales)
{
  _corAniso->setRotationAnglesAndRadius(angles, ranges, scales);
}

bool CovAniso::isValidForTurningBand() const
{
  return _corAniso->isValidForTurningBand();
}
double CovAniso::simulateTurningBand(double t0, TurningBandOperate &operTB) const
{
  return _corAniso->simulateTurningBand(t0, operTB);
}
bool CovAniso::isValidForSpectral() const
{
  return _corAniso->isValidForSpectral();
}
MatrixRectangular CovAniso::simulateSpectralOmega(int nb) const
{
  return _corAniso->simulateSpectralOmega(nb);
}


double CovAniso::eval0(int ivar, int jvar, const CovCalcMode* mode) const
{
  double cov = _corAniso->evalCorFromH(0, mode);

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
  double cov = _corAniso->evalCor(p1,p2,mode);
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
  double cov = _corAniso->evalCorFromH(0, mode); 

  if (mode == nullptr || ! mode->getUnitary())
    mat.addMatInPlace(_sill, 1., cov);
  else
  {
    mat.addMatInPlace(_workMat, 1., cov);
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
 * @remark: The optimized version is not compatible with non-stationarity (except sills).
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
        cov = _corAniso->evalCorFromH(hoptim, mode);
        res.updValue(irow, icol, EOperator::ADD, sill * cov);
      }
      irow++;
    }
  }
}



double CovAniso::evalCovOnSphere(double alpha,
                                 int degree,
                                 bool flagScaleDistance,
                                 const CovCalcMode* mode) const
{
  double value = _corAniso->evalCovOnSphere(alpha, degree, flagScaleDistance,mode);

  if (mode == nullptr || ! mode->getUnitary())
    value *= getSill(0,0);

  return value;
}

VectorDouble CovAniso::evalSpectrumOnSphere(int n, bool flagNormDistance, bool flagCumul) const
{
  return _corAniso->evalSpectrumOnSphere(n, flagNormDistance, flagCumul);
}

void CovAniso::setMarkovCoeffs(const VectorDouble& coeffs)
{
  _corAniso->setMarkovCoeffs(coeffs);
}

/* This function computes a polynomial P from two polynomials P1 and P2 and a small constant eps
 * P(x) = P1(x)^2 + x * P2(x)^2 + eps
 */
void CovAniso::setMarkovCoeffsBySquaredPolynomials(VectorDouble coeffs1,
                                                   VectorDouble coeffs2,
                                                   double eps)
{
  _corAniso->setMarkovCoeffsBySquaredPolynomials(coeffs1, coeffs2, eps);
}

double CovAniso::getCorrec() const
{
  return _corAniso->getCorrec();
}

double CovAniso::getFullCorrec() const
{
  return  _corAniso->getFullCorrec();
}

double CovAniso::_getDetTensor() const
{
  return _corAniso->getDetTensor();
}

double CovAniso::evalSpectrum(const VectorDouble& freq, int ivar, int jvar) const
{
  if (!_corAniso->hasSpectrumOnRn()) return TEST;
  return _sill.getValue(ivar, jvar) * _corAniso->evalSpectrum(freq, ivar, jvar);
}

double CovAniso::normalizeOnSphere(int n) const
{ 
  return _corAniso->normalizeOnSphere(n);
}

VectorDouble CovAniso::getMarkovCoeffs() const
{
  return _corAniso->getMarkovCoeffs();
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
 
  sstr << _corAniso->getCova()->toString();

  // Sill - Factor / Slope information
  if (_corAniso->hasRange() > 0)
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
  else if (_corAniso->hasRange() < 0)
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
  sstr << _corAniso->toStringParams(strfmt);

  // Non-stationary parameters

  if (isNoStat())
  { 
    sstr << toTitle(1, "Non-Stationary Parameters");
    sstr << _tabNoStat.toString(strfmt);
    int i = _tabNoStat.getNSills();
    sstr << _corAniso->toStringNoStat(strfmt,i);
  }
  return sstr.str();
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
  return _corAniso->getRanges();
}

void CovAniso::setType(const ECov &type)
{
 _corAniso->setType(type);
}

/**
 * This function returns the range in the isotropic case
 * In the anisotropic case, it returns the largest range over all directions
 * @return
 */
double CovAniso::getRange() const
{
 return _corAniso->getRange();
}

double CovAniso::getScale() const
{
 return _corAniso->getScale();
}

VectorDouble CovAniso::getAnisoCoeffs() const
{
  return _corAniso->getAnisoCoeffs();
}

/**
 * For compatibility, this function returns 0 if the Covariance has no Third Parameter
 *
 * @return Third parameter
 */
double CovAniso::getParam() const
{
  return _corAniso->getParam();
}

/**
 * Calculate the Integral Range in various Space Dimension (1, 2 or 3)
 * @return
 */
double CovAniso::getIntegralRange(int ndisc, double hmax) const
{
  return _sill.getValue(0, 0) * _corAniso->getIntegralRange(ndisc, hmax);
}


int CovAniso::getGradParamNumber() const
{
  return _corAniso->getGradParamNumber();
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

 
///////////////////// Range ////////////////////////
void CovAniso::makeRangeNoStatDb(const String &namecol, int idim, const Db* db)
{   
  if(_checkAndManageNoStatDb(db,namecol))
  {
    _corAniso->makeRangeNoStatDb(namecol,idim,db);
  }
 
}

void CovAniso::makeRangeNoStatFunctional(const AFunctional *func, int idim)
{
  _corAniso->makeRangeNoStatFunctional(func,idim);
}


void CovAniso::makeRangeStationary(int idim)
{
  _corAniso->makeRangeStationary(idim);
}

///////////////////// Scale ////////////////////////

void CovAniso::makeScaleNoStatDb(const String &namecol, int idim, const Db* db)
{   
  if(_checkAndManageNoStatDb(db,namecol))
  {
    _corAniso->makeScaleNoStatDb(namecol,idim,db);
  }
}


void CovAniso::makeScaleNoStatFunctional(const AFunctional *func, int idim)
{
  _corAniso->makeScaleNoStatFunctional(func,idim);
}

void CovAniso::makeScaleStationary(int idim)
{
  _corAniso->makeScaleStationary(idim);
}

///////////////////// Angle ////////////////////////


void CovAniso::makeAngleNoStatDb(const String &namecol, int idim, const Db* db)
{
  if(_checkAndManageNoStatDb(db,namecol))
  {
   _corAniso->makeAngleNoStatDb(namecol,idim,db);
  }
}

void CovAniso::makeAngleNoStatFunctional(const AFunctional *func, int idim)
{
  _corAniso->makeAngleNoStatFunctional(func,idim);
}

void CovAniso::makeAngleStationary(int idim)
{
  _corAniso->makeAngleStationary(idim);
}
///////////////////// Tensor ////////////////////////


void CovAniso::makeTensorNoStatDb(const String &namecol, int idim, int jdim,const Db* db)
{
 if(_checkAndManageNoStatDb(db,namecol))
 {
    _corAniso->makeTensorNoStatDb(namecol,idim,jdim,db);
 } 
}

void CovAniso::makeTensorNoStatFunctional(const AFunctional  *func, int idim, int jdim)
{
   _corAniso->makeTensorNoStatFunctional(func,idim,jdim);
}

void CovAniso::makeTensorStationary(int idim, int jdim)
{
  _corAniso->makeTensorStationary(idim,jdim);
}

///////////////////// Param ////////////////////////

void CovAniso::makeParamNoStatDb(const String &namecol, const Db* db)
{
 if(_checkAndManageNoStatDb(db,namecol))
 {
    _corAniso->makeParamNoStatDb(namecol,db);
 } 

}

void CovAniso::makeParamNoStatFunctional(const AFunctional *func)
{
  _corAniso->makeParamNoStatFunctional(func);

}

void CovAniso::makeParamStationary()
{
  _corAniso->makeParamStationary();
}


void CovAniso::informMeshByMeshForAnisotropy(const AMesh* amesh) const
{
  _corAniso->informMeshByMeshForAnisotropy(amesh);

}

void CovAniso::informMeshByApexForAnisotropy(const AMesh* amesh) const
{
   _corAniso->informMeshByApexForAnisotropy(amesh);
}

void CovAniso::informDbInForAnisotropy(const Db* dbin) const
{
  _corAniso->informDbInForAnisotropy(dbin);

}
void CovAniso::informDbOutForAnisotropy(const Db* dbout) const
{
   _corAniso->informDbOutForAnisotropy(dbout);  
}

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

void CovAniso::_optimizationTransformSP(const SpacePoint& ptin, SpacePoint& ptout) const
{
  _corAniso->optimizationTransformSP(ptin, ptout);
}
