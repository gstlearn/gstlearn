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
#include "geoslib_define.h"

#include "Arrays/Array.hpp"
#include "Basic/AFunctional.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CorAniso.hpp"
#include "Covariances/CovProportional.hpp"
#include "Db/Db.hpp"
#include "Covariances/NoStatArray.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/FFT.hpp"
#include "Space/ASpace.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/SpaceRN.hpp"
#include "Space/SpaceSN.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/MatrixFactory.hpp"

#include <math.h>
#include <functional>
#include <ostream>

CovAniso::CovAniso(const ECov& type, const CovContext& ctxt)
  : CovProportional(new CorAniso(type, ctxt), MatrixSymmetric(ctxt.getNVar()))
{
  CovProportional::setCor(getCorAniso());
  _initFromContext();
}

CovAniso::CovAniso(const String& symbol, const CovContext& ctxt)
  : CovProportional(new CorAniso(symbol, ctxt))
{
  ECov covtype = CovFactory::identifyCovariance(symbol, ctxt);
  _initFromContext();

}

CovAniso::CovAniso(const ECov& type,
                   double range,
                   double param,
                   double sill,
                   const CovContext& ctxt,
                   bool flagRange)
  : CovProportional(new CorAniso(type, range, param, ctxt, flagRange), MatrixSymmetric(ctxt.getNVar()))
{
  _initFromContext();

  // Sill
  if (ctxt.getNVar() == 1)
    _sillCur.setValue(0, 0, sill);
  else
  {
    int nvar = ctxt.getNVar();
    _sillCur.fill(0);
    for (int ivar = 0; ivar < nvar; ivar++)
      _sillCur.setValue(ivar, ivar, sill);
  }

  // Param
  setParam(param);

  // Range
  if (flagRange)
    setRangeIsotropic(range);
  else
    setScale(range);
}

CovAniso::CovAniso(const CovAniso& r)
  : CovProportional(r)
{
  _ctxt.setNVar(r.getNVar());
}

CovAniso& CovAniso::operator=(const CovAniso& r)
{
  if (this != &r)
  {
    setCor(new CorAniso(*r.getCorAniso()));
    _ctxt         = r._ctxt;
    _sillCur      = r._sillCur;
  }
  return *this;
}

CovAniso::~CovAniso()
{
  delete getCorAniso();
}

CorAniso* CovAniso::getCorAniso()
{
  return (CorAniso*)getCor();
}

MatrixDense CovAniso::simulateSpectralOmega(int nb) const
{
  return getCorAniso()->simulateSpectralOmega(nb);
}

double CovAniso::_getSillValue(int ivar, int jvar, const CovCalcMode* mode) const
{
  if (mode != nullptr && mode->getUnitary()) return 1.;
  return getSill(ivar, jvar);
}

double CovAniso::eval0(int ivar, int jvar, const CovCalcMode* mode) const
{
  double cov = getCorAniso()->evalCorFromH(0, mode);
  return cov * _getSillValue(ivar, jvar, mode);
}

double CovAniso::_eval(const SpacePoint& p1,
                       const SpacePoint& p2,
                       int ivar,
                       int jvar,
                       const CovCalcMode* mode) const
{
  double cov = getCorAniso()->evalCor(p1, p2, mode);
  return cov * _getSillValue(ivar, jvar, mode);
}

double CovAniso::evalCovOnSphere(double alpha,
                                 int degree,
                                 bool flagScaleDistance,
                                 const CovCalcMode* mode) const
{
  double value = getCorAniso()->evalCovOnSphere(alpha, degree, flagScaleDistance, mode);
  return value * _getSillValue(0, 0, mode);
}

VectorDouble CovAniso::evalSpectrumOnSphere(int n, bool flagNormDistance, bool flagCumul) const
{
  return getCorAniso()->evalSpectrumOnSphere(n, flagNormDistance, flagCumul);
}

double CovAniso::getCorrec() const
{
  return getCorAniso()->getCorrec();
}

double CovAniso::getFullCorrec() const
{
  return getCorAniso()->getFullCorrec();
}

double CovAniso::_getDetTensor() const
{
  return getCorAniso()->getDetTensor();
}

double CovAniso::evalSpectrum(const VectorDouble& freq, int ivar, int jvar) const
{
  if (!getCorAniso()->hasSpectrumOnRn()) return TEST;
  return _sillCur.getValue(ivar, jvar) * getCorAniso()->evalSpectrum(freq, ivar, jvar);
}

double CovAniso::normalizeOnSphere(int n) const
{
  return getCorAniso()->normalizeOnSphere(n);
}

VectorDouble CovAniso::getMarkovCoeffs() const
{
  return getCorAniso()->getMarkovCoeffs();
}

VectorDouble CovAniso::evalCovOnSphereVec(const VectorDouble& alpha,
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

String CovAniso::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;

  sstr << getCorAniso()->getCorFunc()->toString();

  // Sill - Factor / Slope information
  if (getCorAniso()->hasRange() >= 0)
  {
    // A sill is defined

    if (getNVar() > 1)
    {
      sstr << toMatrix("- Sill matrix:", VectorString(), VectorString(), 0,
                       getNVar(), getNVar(), _sillCur.getValues());
    }
    else
    {
      sstr << "- Sill         = " << toDouble(_sillCur.getValue(0, 0)) << std::endl;
    }
  }
  else
  {
    // The sill is not defined: use slope instead

    if (getNVar() > 1)
    {
      MatrixSquare slopes = _sillCur;
      double range               = getRange(0);
      for (int ivar = 0; ivar < getNVar(); ivar++)
        for (int jvar = 0; jvar < getNVar(); jvar++)
          slopes.setValue(ivar, jvar, _sillCur.getValue(ivar, jvar) / range);
      sstr << toMatrix("- Slope matrix:", VectorString(), VectorString(), 0,
                       getNVar(), getNVar(), slopes.getValues());
    }
    else
    {
      sstr << "- Slope        = " << toDouble(getSlope(0, 0)) << std::endl;
    }
  }

  // Covariance Parameters
  sstr << getCorAniso()->toStringParams(strfmt);

  // Non-stationary parameters

  if (isNoStat())
  {
    sstr << toTitle(1, "Non-Stationary Parameters");
    sstr << _tabNoStat->toString(strfmt);
    int i = getTabNoStatSills()->getNSills();
    sstr << getCorAniso()->toStringNoStat(strfmt, i);
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
  return _sillCur.getValue(ivar, jvar) / range;
}

VectorDouble CovAniso::getRanges() const
{
  return getCorAniso()->getRanges();
}

VectorDouble CovAniso::getAnisoCoeffs() const
{
  return getCorAniso()->getAnisoCoeffs();
}

/**
 * For compatibility, this function returns 0 if the Covariance has no Third Parameter
 *
 * @return Third parameter
 */
double CovAniso::getParam() const
{
  return getCorAniso()->getParam();
}

/**
 * Calculate the Integral Range in various Space Dimension (1, 2 or 3)
 * @return
 */
double CovAniso::getIntegralRange(int ndisc, double hmax) const
{
  return _sillCur.getValue(0, 0) * getCorAniso()->getIntegralRange(ndisc, hmax);
}

int CovAniso::getNGradParam() const
{
  return getCorAniso()->getNGradParam();
}

CovAniso* CovAniso::createIsotropic(const CovContext& ctxt,
                                    const ECov& type,
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

CovAniso* CovAniso::createAnisotropic(const CovContext& ctxt,
                                      const ECov& type,
                                      const VectorDouble& ranges,
                                      double sill,
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

  CovAniso* cov = new CovAniso(type, ctxt);
  if (flagRange)
    cov->setRanges(ranges);
  else
    cov->setScales(ranges);
  cov->setSill(sill);
  cov->setParam(param);
  if (!angles.empty()) cov->setAnisoAngles(angles);
  return cov;
}

CovAniso* CovAniso::createIsotropicMulti(const CovContext& ctxt,
                                         const ECov& type,
                                         double range,
                                         const MatrixSymmetric& sills,
                                         double param,
                                         bool flagRange)
{
  CovAniso* cov = new CovAniso(type, ctxt);
  int nvar      = sills.getNSize();
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

CovAniso* CovAniso::createAnisotropicMulti(const CovContext& ctxt,
                                           const ECov& type,
                                           const VectorDouble& ranges,
                                           const MatrixSymmetric& sills,
                                           double param,
                                           const VectorDouble& angles,
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
  int ndim = (int)ranges.size();
  if ((int)ctxt.getNDim() != ndim)
  {
    messerr("Mismatch in Space Dimension between 'ranges'(%d) and 'ctxt'(%d)",
            ndim, ctxt.getNDim());
    return nullptr;
  }

  CovAniso* cov = new CovAniso(type, ctxt);
  if (flagRange)
    cov->setRanges(ranges);
  else
    cov->setScales(ranges);
  cov->setSill(sills);
  cov->setParam(param);
  if (!angles.empty()) cov->setAnisoAngles(angles);
  return cov;
}

CovAniso* CovAniso::createFromParam(const ECov& type,
                                    double range,
                                    double sill,
                                    double param,
                                    const VectorDouble& ranges,
                                    const MatrixSymmetric& sills,
                                    const VectorDouble& angles,
                                    const ASpaceSharedPtr& space,
                                    bool flagRange)
{
  // Check consistency with parameters of the model

  int ndim = 0;
  if (!ranges.empty())
  {
    if (ndim > 0 && (int)ranges.size() != ndim)
    {
      messerr("Mismatch between the dimension of 'ranges' (%d)",
              (int)ranges.size());
      messerr("and the Space dimension stored in the Model (%d)", ndim);
      messerr("Operation is cancelled");
      return nullptr;
    }
    ndim = (int)ranges.size();
  }
  if (!angles.empty())
  {
    if (ndim > 0 && (int)angles.size() != ndim)
    {
      messerr("Mismatch between the dimension of 'angles' (%d)",
              (int)angles.size());
      messerr("and the Space dimension stored in the Model (%d)", ndim);
      messerr("Operation is cancelled");
      return nullptr;
    }
    ndim = (int)angles.size();
  }
  if (space != nullptr)
  {
    if (ndim > 0 && (int)space->getNDim() != ndim)
    {
      messerr("Mismatch between the space dimension in 'space' (%d)",
              (int)space->getNDim());
      messerr("and the Space dimension stored in the Model (%d)", ndim);
      messerr("Operation is cancelled");
      return nullptr;
    }
    ndim = (int)space->getNDim();
  }
  if (ndim <= 0)
  {
    messerr("You must define the SPace dimension");
    return nullptr;
  }

  int nvar = 0;
  if (!sills.empty())
  {
    if (nvar > 0 && nvar != sills.getNCols())
    {
      messerr("Mismatch between the number of rows 'sills' (%d)", sills.getNRows());
      messerr("and the Number of variables stored in the Model (%d)", nvar);
      messerr("Operation is cancelled");
      return nullptr;
    }
    nvar = (int)sqrt((double)sills.size());
  }
  if (nvar <= 0) nvar = 1;

  // Define the covariance

  const CovContext& ctxt = CovContext(nvar, space);
  CovAniso* cov          = new CovAniso(type, ctxt);

  // Define the Third parameter
  double parmax = cov->getParMax();
  if (param > parmax) param = parmax;
  cov->setParam(param);

  // Define the range
  if (!ranges.empty())
  {
    if (flagRange)
      cov->setRanges(ranges);
    else
      cov->setScales(ranges);
  }
  else
  {
    if (flagRange)
      cov->setRangeIsotropic(range);
    else
      cov->setScale(range);
  }

  // Define the sill
  if (!sills.empty())
    cov->setSill(sills);
  else
  {
    if (nvar <= 1)
      cov->setSill(sill);
    else
    {
      MatrixSymmetric locsills(nvar);
      locsills.setIdentity(sill);
      cov->setSill(locsills);
    }
  }

  if (!angles.empty()) cov->setAnisoAngles(angles);
  return cov;
}

Array CovAniso::evalCovFFT(const VectorDouble& hmax,
                           int N,
                           int ivar,
                           int jvar) const
{
  if (!hasSpectrumOnRn()) return Array();

  std::function<double(const VectorDouble&)> funcSpectrum;
  funcSpectrum = [this, ivar, jvar](const VectorDouble& freq)
  {
    return evalSpectrum(freq, ivar, jvar) * _getDetTensor();
  };
  return evalCovFFTSpatial(hmax, N, funcSpectrum);
}

CovAniso* CovAniso::createReduce(const VectorInt& validVars) const
{
  CovAniso* newCovAniso = this->clone();

  // Modify the CovContext
  int nvar        = (int)validVars.size();
  CovContext ctxt = CovContext(nvar);

  // Modify the Matrix of sills
  newCovAniso->setContext(ctxt);
  MatrixSymmetric* newsill = dynamic_cast<MatrixSymmetric*>(MatrixFactory::createReduce(&_sillCur, validVars, validVars));
  newCovAniso->setSill(*newsill);
  return newCovAniso;
}

double scale2range(const ECov& type, double scale, double param)
{
  CovContext ctxt = CovContext(1, 1);
  ACovFunc* cova  = CovFactory::createCovFunc(type, ctxt);
  cova->setParam(param);
  double scadef = cova->getScadef();
  return scale * scadef;
}

double range2scale(const ECov& type, double range, double param)
{
  CovContext ctxt = CovContext(1, 1);
  ACovFunc* cova  = CovFactory::createCovFunc(type, ctxt);
  cova->setParam(param);
  double scadef = cova->getScadef();
  return range / scadef;
}
