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
#include "Basic/AStringable.hpp"
#include "Basic/FFT.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/ACov.hpp"
#include "Covariances/CorAniso.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovProportional.hpp"
#include "Db/Db.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Space/ASpace.hpp"
#include "Space/ASpaceObject.hpp"
#include "Space/SpacePoint.hpp"
#include "Space/SpaceRN.hpp"
#include "Space/SpaceSN.hpp"

#include <cmath>
#include <functional>
#include <ostream>

namespace gstlrn
{

CovAniso::CovAniso(const ECov& type, const CovContext& ctxt)
  : CovProportional(nullptr, MatrixSymmetric(ctxt.getNVar()))
{
  auto tempcor = CorAniso(type, ctxt);
  CovProportional::setCor(&tempcor);
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
  : CovProportional(nullptr, MatrixSymmetric(ctxt.getNVar()))
{
  auto temp = CorAniso(type, range, param, ctxt, flagRange);
  setCor(&temp);
  _initFromContext();

  // Sill
  if (ctxt.getNVar() == 1)
    _sillCur.setValue(0, 0, sill);
  else
  {
    auto nvar = ctxt.getNVar();
    _sillCur.fill(0);
    for (Id ivar = 0; ivar < nvar; ivar++)
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
    _ctxt    = r._ctxt;
    _sillCur = r._sillCur;
  }
  return *this;
}

CovAniso::~CovAniso()
{
}

const CorAniso* CovAniso::getCorAniso() const
{
  return static_cast<const CorAniso*>(getCor());
}

CorAniso* CovAniso::getCorAnisoModify()
{
  return static_cast<CorAniso*>(getCorModify());
}

double CovAniso::_getSillValue(Id ivar, Id jvar, const CovCalcMode* mode) const
{
  if (mode != nullptr && mode->getUnitary()) return 1.;
  return getSill(ivar, jvar);
}

double CovAniso::eval0(Id ivar, Id jvar, const CovCalcMode* mode) const
{
  double cov = getCorAniso()->evalCorFromH(0, mode);
  return cov * _getSillValue(ivar, jvar, mode);
}

double CovAniso::_eval(const SpacePoint& p1,
                       const SpacePoint& p2,
                       Id ivar,
                       Id jvar,
                       const CovCalcMode* mode) const
{
  double cov = getCorAniso()->evalCor(p1, p2, mode);
  return cov * _getSillValue(ivar, jvar, mode);
}

double CovAniso::evalCovOnSphere(double alpha,
                                 Id degree,
                                 bool flagScaleDistance,
                                 const CovCalcMode* mode) const
{
  double value = getCorAniso()->evalCovOnSphere(alpha, degree, flagScaleDistance, mode);
  return value * _getSillValue(0, 0, mode);
}

VectorDouble CovAniso::evalSpectrumOnSphere(Id n, bool flagNormDistance, bool flagCumul) const
{
  return getCorAniso()->evalSpectrumOnSphere(n, flagNormDistance, flagCumul);
}

double CovAniso::evalSpectrum(const VectorDouble& freq, Id ivar, Id jvar) const
{
  if (!getCorAniso()->hasSpectrumOnRn()) return TEST;
  return _sillCur.getValue(ivar, jvar) * getCorAniso()->evalSpectrum(freq, ivar, jvar);
}

VectorDouble CovAniso::evalCovOnSphereVec(const VectorDouble& alpha,
                                          Id degree,
                                          bool flagScaleDistance,
                                          const CovCalcMode* mode) const
{
  Id n = static_cast<Id>(alpha.size());
  VectorDouble vec(n);
  for (Id i = 0; i < n; i++)
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
      double range        = getRange(0);
      for (Id ivar = 0; ivar < getNVar(); ivar++)
        for (Id jvar = 0; jvar < getNVar(); jvar++)
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
    auto i = getTabNoStatSills()->getNSills();
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
double CovAniso::getSlope(Id ivar, Id jvar) const
{
  if (hasRange() == 0) return TEST;
  double range = getRange(0);
  return _sillCur.getValue(ivar, jvar) / range;
}

/**
 * Calculate the Integral Range in various Space Dimension (1, 2 or 3)
 * @return
 */
double CovAniso::getIntegralRange(Id ndisc, double hmax) const
{
  return _sillCur.getValue(0, 0) * getCorAniso()->getIntegralRange(ndisc, hmax);
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
  Id ndim = static_cast<Id>(ranges.size());
  if (static_cast<Id>(ctxt.getNDim()) != ndim)
  {
    messerr("Mismatch in Space Dimension between 'ranges'(%d) and 'ctxt'(%d)",
            ndim, ctxt.getNDim());
    return nullptr;
  }

  auto* cov = new CovAniso(type, ctxt);
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
  auto* cov = new CovAniso(type, ctxt);
  auto nvar = sills.getNSize();
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

  auto nvar = sills.getNSize();
  if (ctxt.getNVar() != nvar)
  {
    messerr(
      "Mismatch in the number of variables between 'sills'(%d) and 'ctxt'(%d)",
      nvar, ctxt.getNVar());
    return nullptr;
  }
  Id ndim = static_cast<Id>(ranges.size());
  if (static_cast<Id>(ctxt.getNDim()) != ndim)
  {
    messerr("Mismatch in Space Dimension between 'ranges'(%d) and 'ctxt'(%d)",
            ndim, ctxt.getNDim());
    return nullptr;
  }

  auto* cov = new CovAniso(type, ctxt);
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

  Id ndim = 0;
  if (!ranges.empty())
  {
    if (ndim > 0 && static_cast<Id>(ranges.size()) != ndim)
    {
      messerr("Mismatch between the dimension of 'ranges' (%d)",
              static_cast<Id>(ranges.size()));
      messerr("and the Space dimension stored in the Model (%d)", ndim);
      messerr("Operation is cancelled");
      return nullptr;
    }
    ndim = static_cast<Id>(ranges.size());
  }
  if (!angles.empty())
  {
    if (ndim > 0 && static_cast<Id>(angles.size()) != ndim)
    {
      messerr("Mismatch between the dimension of 'angles' (%d)",
              static_cast<Id>(angles.size()));
      messerr("and the Space dimension stored in the Model (%d)", ndim);
      messerr("Operation is cancelled");
      return nullptr;
    }
    ndim = static_cast<Id>(angles.size());
  }
  if (space != nullptr)
  {
    if (ndim > 0 && static_cast<Id>(space->getNDim()) != ndim)
    {
      messerr("Mismatch between the space dimension in 'space' (%d)",
              static_cast<Id>(space->getNDim()));
      messerr("and the Space dimension stored in the Model (%d)", ndim);
      messerr("Operation is cancelled");
      return nullptr;
    }
    ndim = static_cast<Id>(space->getNDim());
  }
  if (ndim <= 0)
  {
    messerr("You must define the Space dimension");
    return nullptr;
  }

  Id nvar = 0;
  if (!sills.empty())
  {
    if (nvar > 0 && nvar != sills.getNCols())
    {
      messerr("Mismatch between the number of rows 'sills' (%d)", sills.getNRows());
      messerr("and the Number of variables stored in the Model (%d)", nvar);
      messerr("Operation is cancelled");
      return nullptr;
    }
    nvar = static_cast<Id>(sqrt(static_cast<double>(sills.size())));
  }
  if (nvar <= 0) nvar = 1;

  // Define the covariance

  const CovContext& ctxt = CovContext(nvar, space);
  auto* cov              = new CovAniso(type, ctxt);

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
                           Id N,
                           Id ivar,
                           Id jvar) const
{
  if (!hasSpectrumOnRn()) return Array();

  std::function<double(const VectorDouble&)> funcSpectrum;
  funcSpectrum = [this, ivar, jvar](const VectorDouble& freq)
  {
    return evalSpectrum(freq, ivar, jvar) * getDetTensor();
  };
  return evalCovFFTSpatial(hmax, N, funcSpectrum);
}

CovAniso* CovAniso::createReduce(const VectorInt& validVars) const
{
  CovAniso* newCovAniso = this->clone();

  // Modify the CovContext
  Id nvar = static_cast<Id>(validVars.size());
  CovContext ctxt(nvar);

  // Modify the Matrix of sills
  newCovAniso->setContext(ctxt);
  auto* newsill = dynamic_cast<MatrixSymmetric*>(MatrixFactory::createReduce(&_sillCur, validVars, validVars));
  newCovAniso->setSill(*newsill);
  return newCovAniso;
}

double scale2range(const ECov& type, double scale, double param)
{
  CovContext ctxt(1, 1);
  ACovFunc* cova = CovFactory::createCovFunc(type, ctxt);
  cova->setParam(param);
  double scadef = cova->getScadef();
  return scale * scadef;
}

double range2scale(const ECov& type, double range, double param)
{
  CovContext ctxt(1, 1);
  ACovFunc* cova = CovFactory::createCovFunc(type, ctxt);
  cova->setParam(param);
  double scadef = cova->getScadef();
  return range / scadef;
}
} // namespace gstlrn