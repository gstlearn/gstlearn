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
#include "Covariances/CorMatern.hpp"
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CorAniso.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Simulation/SpectrumOnRN.hpp"
#include "Simulation/SpectrumRN.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceComposite.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"

#include <cmath>
#include <memory>
#include <vector>

namespace gstlrn
{
CorMatern::CorMatern(const CovContext& ctxt, const ECov& type, const VectorDouble& params, const VectorDouble& kappas, const VectorDouble& ranges, const VectorDouble& angles, bool flagRange)
  : ACov(ctxt)
  , _nVar(static_cast<Id>(ctxt.getNVar()))
  , _corRef(std::shared_ptr<const CorAniso>(
      CorAniso::createAnisotropic(
        // CovContext(1, static_cast<Id>(ctxt.getNDim())),
        CovContext(1, ctxt.getSpace()),
        ECov::MATERN,
        ranges,
        params[0],
        angles,
        flagRange)))
  , _cor(*_corRef)
  , _Nu(_nVar)
  , _Kappa(_nVar)
  , _C0(_nVar)

{
  int ierr = 0;
  if (type != ECov::MATERN)
  {
    messerr("CorMatern: inconsistent type");
    ierr = 1;
  }
  if (static_cast<Id>(params.size()) != _nVar)
  {
    messerr("CorMatern: inconsistent size of params");
    messerr("CorMatern: params size = %d, context size = %d", params.size(), _nVar);
    ierr = 1;
  }
  if (static_cast<Id>(kappas.size()) != _nVar)
  {
    messerr("CorMatern: inconsistent size of coeffScales");
    messerr("CorMatern: coeffScales size = %d, context size = %d", kappas.size(), _nVar);
    ierr = 2;
  }

  //_coeffScales.push_front(1.);
  setContext(_cor.getContext());
  _ctxt.setNVar(_nVar);
  for (Id ivar = 0; ivar < _nVar; ivar++)
  {
    double nui    = params[ivar];
    double scalei = kappas[ivar];
    _Nu.setValue(ivar, ivar, nui);
    _Kappa.setValue(ivar, ivar, scalei);
    _C0.setValue(ivar, ivar, 1.);
    double lgni = loggamma(nui);
    for (Id jvar = ivar + 1; jvar < _nVar; jvar++)
    {
      double nuj     = params[jvar];
      double scalej  = kappas[jvar];
      double nuij    = computeNu(nui, nuj);
      double scaleij = computeKappa(scalei, scalej);
      _Nu.setValue(ivar, jvar, nuij);
      _Kappa.setValue(ivar, jvar, scaleij);
      double lgnj  = loggamma(nuj);
      double lgnij = loggamma(nuij);
      double val   = exp(lgnij - 0.5 * (lgni + lgnj));
      double tau   = pow(scalei, nui) * pow(scalej, nuj) / pow(scaleij, 2. * nuij);
      val *= tau;
      _C0.setValue(ivar, jvar, val);
    }
  }
  // checking the consistancy of the covariance matrix
  if (!_C0.isDefinitePositive())
  {
    messerr("CorMatern: inconsistent covariance matrix");
    ierr = 3;
  }
  if (ierr > 0)
  {
    _nVar  = 0;
    _C0    = MatrixSymmetric(0);
    _Kappa = MatrixSymmetric(0);
    _Nu    = MatrixSymmetric(0);
    return;
  }
}

CorMatern::CorMatern(const CorMatern& r)
  : ACov(r)
  , _cor(r._cor)
{
  _nVar   = r._nVar;
  _corRef = r._corRef;
  _C0     = r._C0;
  _Kappa  = r._Kappa;
  _Nu     = r._Nu;
}

CorMatern& CorMatern::operator=(const CorMatern& r)
{
  if (this != &r)
  {
    ACov::operator=(r);
    _nVar  = r._nVar;
    _cor   = r._cor;
    _C0    = r._C0;
    _Kappa = r._Kappa;
    _Nu    = r._Nu;
  }
  return *this;
}

CorMatern* CorMatern::create(
  const CovContext& ctxt,
  const ECov& type,
  const VectorDouble& params,
  const VectorDouble& kappas,
  const VectorDouble& ranges,
  const VectorDouble& angles,
  bool flagRange)
{
  Id ndim = ctxt.getNDim();
  Id nvar = ctxt.getNVar();
  if (ranges.length() != ndim)
  {
    messerr("Mismatch in Space Dimension between 'ctxt'(%d) and 'ranges'(%d)",
            ndim, ranges.length());
    return nullptr;
  }
  if (angles.length() != ndim)
  {
    messerr("Mismatch in Space Dimension between 'ctxt'(%d) and 'angles'(%d)",
            ndim, angles.length());
    return nullptr;
  }
  if (params.length() != nvar)
  {
    messerr("Mismatch in number of variables between 'ctxt'(%d) and 'params'(%d)",
            nvar, params.length());
    return nullptr;
  }
  if (kappas.length() != nvar)
  {
    messerr("Mismatch in number of variables between 'ctxt'(%d) and 'kappas'(%d)",
            nvar, kappas.length());
    return nullptr;
  }

  auto* cov = new CorMatern(ctxt, type, params, kappas, ranges, angles, flagRange);
  return cov;
}

double CorMatern::computeKappa(double kappai, double kappaj)
{
  double ci2   = kappai * kappai;
  double cj2   = kappaj * kappaj;
  double scale = sqrt(0.5 * (ci2 + cj2));
  return scale;
}

double CorMatern::evalSpectrumOnRN(const VectorDouble& freq, Id ivar, Id jvar) const
{
  double kappa        = _Kappa.getValue(ivar, jvar);
  VectorDouble scales = _corRef->getScales();
  VectorDouble angles = _corRef->getAnisoAngles();
  for (size_t idim = 0; idim < getSpace()->getNDim(); idim++)
  {
    scales[idim] /= kappa;
  }
  _cor.setRotationAnglesAndRadius(angles, VectorDouble(), scales);
  _cor.setParam(_Nu.getValue(ivar, jvar));
  return _C0.getValue(ivar, jvar) * _cor.evalSpectrumOnRN(freq, ivar, jvar);
}

double CorMatern::evalSpectrumRatio(const VectorDouble& freq, Id ivar, Id jvar, const ACov* cov0) const
{
  DECLARE_UNUSED(cov0)
  return evalSpectrumOnRN(freq, ivar, jvar) / _corRef->evalSpectrumOnRN(freq);
}

MatrixDense CorMatern::simulateSpectralOmega(Id nb) const
{
  return _corRef->simulateSpectralOmega(nb);
}

/*
SpectrumRN CorMatern::simulateSpectrumRN(Id ns, const ACov* cov0) const
{
  DECLARE_UNUSED(cov0)
  Id nvar = getNVar();

  // simulation of the frequencies
  MatrixDense omega = _corRef->simulateSpectralOmega(ns);

  // simulation of the normalizing factors gamma
  MatrixDense gamma(ns, getNVar());
  VectorDouble values(nvar);
  MatrixSymmetric H(nvar);
  for (Id ib = 0; ib < ns; ib++)
  {
    double val = sqrt(-log(law_uniform()) * 2 * nvar / ns);
    if (nvar == 1)
    {
      values[0] = val;
    }
    else
    {
      VectorDouble freq = omega.getRow(ib);
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        for (Id jvar = 0; jvar <= ivar; jvar++)
        {
          double ratioIS = evalSpectrumOnRN(freq, ivar, jvar) / _corRef->evalSpectrumOnRN(freq);
          H.setValue(ivar, jvar, ratioIS);
        }
      }
      // square root of the symmetric matrix
      if (H.squareRootInPlace(H) != 0)
      {
        message("Error in computing square root matrix\n");
      }
      Id icol        = law_int_uniform(0, nvar - 1);
      VectorDouble v = H.getColumn(icol);
      for (Id ivar = 0; ivar < nvar; ivar++)
      {
        values[ivar] = val * v[ivar];
      }
    }
    gamma.setRow(ib, values);
  }
  return SpectrumRN(gamma, omega);
}
*/

SpectrumOnRN* CorMatern::simulateOnRN(Id ns) const
{
  Id nvar = getNVar();

  // simulation of the random frequencies using _corRef
  MatrixDense omega = simulateSpectralOmega(ns);

  // simulation of the normalizing factors gamma
  MatrixDense gamma(ns, getNVar());
  VectorDouble values(nvar);
  MatrixSymmetric H(nvar);
  for (Id ib = 0; ib < ns; ib++)
  {
    double val        = sqrt(-log(law_uniform()) * nvar / ns);
    VectorDouble freq = omega.getRow(ib);
    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      for (Id jvar = 0; jvar <= ivar; jvar++)
      {
        double ratioIS = evalSpectrumOnRN(freq, ivar, jvar) / _corRef->evalSpectrumOnRN(freq);
        H.setValue(ivar, jvar, ratioIS);
      }
    }
    // square root of the symmetric matrix
    if (H.squareRootInPlace(H) != 0)
    {
      message("Error in computing square root matrix\n");
    }
    Id icol        = law_int_uniform(0, nvar - 1);
    VectorDouble v = H.getColumn(icol);
    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      values[ivar] = val * v[ivar];
    }
    gamma.setRow(ib, values);
  }

  // simulation of the random phase
  VectorDouble phi(ns);
  for (Id ib = 0; ib < ns; ib++)
  {
    phi[ib] = law_uniform(0.0, 2 * GV_PI);
  }

  // creating the spectrum
  auto* res = new SpectrumOnRNSimple(getNVar(), getNDim(), ns);
  res->setGamma(gamma);
  res->addFactor(omega, phi);
  return res;
};

CorMatern::~CorMatern()
{
}

void CorMatern::_optimizationSetTarget(SpacePoint& pt) const
{
  DECLARE_UNUSED(pt)
}

void CorMatern::_optimizationPreProcess(Id mode, const std::vector<SpacePoint>& ps) const
{
  DECLARE_UNUSED(mode)
  DECLARE_UNUSED(ps)
  // _covS->_optimizationPreProcess(p);
  // _covTemp->_optimizationPreProcess(p);
}

void CorMatern::_optimizationPostProcess() const
{
  //_covS->optimizationPostProcess();
  //_covTemp->optimizationPostProcess();
}

double CorMatern::_eval(const SpacePoint& p1,
                        const SpacePoint& p2,
                        Id ivar,
                        Id jvar,
                        const CovCalcMode* mode) const
{
  VectorDouble scales = _corRef->getScales();
  VectorDouble angles = _corRef->getAnisoAngles();
  double correcScale  = _Kappa.getValue(ivar, jvar);
  for (size_t idim = 0; idim < getSpace()->getNDim(); idim++)
  {
    scales[idim] /= correcScale;
  }
  _cor.setParam(_Nu.getValue(ivar, jvar));
  _cor.setRotationAnglesAndRadius(angles, VectorDouble(), scales);
  return _C0.getValue(ivar, jvar) * _cor.evalCov(p1, p2, ivar, jvar, mode);
}

} // namespace gstlrn