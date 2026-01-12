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
#include "Covariances/CorGaussianMixture.hpp"
#include "Basic/Law.hpp"
#include "Basic/MathFunc.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CorAniso.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceComposite.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"

#include <cmath>
#include <memory>

namespace gstlrn
{
CorGaussianMixture::CorGaussianMixture(
  const CovContext& ctxt,
  const ECov& type,
  const VectorDouble& params,
  const VectorDouble& kappas,
  const VectorDouble& ranges,
  const VectorDouble& angles,
  bool flagRange)
  : ACov(ctxt)
  , _nVar(static_cast<Id>(ctxt.getNVar()))
  , _corRef(std::shared_ptr<const CorAniso>(
      CorAniso::createAnisotropic(
        CovContext(1, static_cast<Id>(ctxt.getNDim())),
        type,
        ranges,
        params[0],
        angles,
        flagRange)))
  , _cor(*_corRef)
  , _Nu(_nVar)
  , _Kappa(_nVar)
  , _C0(_nVar)
  , _scaleGneiting(1.0)
{
  int ierr = 0;
  if ((type != ECov::CAUCHY) && (type != ECov::MATERN) && (type != ECov::GAUSSIAN))
  {
    messerr("CorGaussianMixture: inconsistent type");
    ierr = 1;
  }
  if (static_cast<Id>(params.size()) != _nVar)
  {
    messerr("CorGaussianMixture: inconsistent size of params");
    messerr("CorGaussianMixture: params size = %d, context size = %d", params.size(), _nVar);
    ierr = 2;
  }
  if (static_cast<Id>(kappas.size()) != _nVar)
  {
    messerr("CorGaussianMixture: inconsistent size of coeffScales");
    messerr("CorGaussianMixture: coeffScales size = %d, context size = %d", kappas.size(), _nVar);
    ierr = 3;
  }
  if (ierr == 0)
  {
    //_coeffScales.push_front(1.);
    setContext(_corRef->getContext());
    _ctxt.setNVar(_nVar);
    bool isMatern = (_corRef->getType() == ECov::MATERN);
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
        if (isMatern)
        {
          val *= tau;
        }
        else
        {
          val /= tau;
        }
        _C0.setValue(ivar, jvar, val);
      }
    }
    // checking the consistancy of the covariance matrix
    if (!_C0.isDefinitePositive())
    {
      messerr("CorGaussianMixture: inconsistent covariance matrix");
      ierr = 3;
    }
  }
  if (ierr != 0)
  {
    _nVar  = 0;
    _Nu    = MatrixSymmetric(0);
    _Kappa = MatrixSymmetric(0);
    _C0    = MatrixSymmetric(0);
    return;
  }
}

CorGaussianMixture::CorGaussianMixture(const CorGaussianMixture& r)
  : ACov(r)
  , _cor(r._cor)
{
  _nVar          = r._nVar;
  _corRef        = r._corRef;
  _Nu            = r._Nu;
  _Kappa         = r._Kappa;
  _C0            = r._C0;
  _scaleGneiting = r._scaleGneiting;
}

CorGaussianMixture& CorGaussianMixture::operator=(const CorGaussianMixture& r)
{
  if (this != &r)
  {
    ACov::operator=(r);
    _nVar          = r._nVar;
    _cor           = r._cor;
    _Nu            = r._Nu;
    _Kappa         = r._Kappa;
    _C0            = r._C0;
    _scaleGneiting = r._scaleGneiting;
  }
  return *this;
}

CorGaussianMixture::~CorGaussianMixture()
{
}

CorGaussianMixture* CorGaussianMixture::create(
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

  auto* cov = new CorGaussianMixture(ctxt, type, params, kappas, ranges, angles, flagRange);
  return cov;
}

double CorGaussianMixture::getScaleDim(Id idim) const
{
  return _corRef->getScale(idim);
}

VectorDouble CorGaussianMixture::getScales() const
{
  return _corRef->getScales();
};

void CorGaussianMixture::setScaleCor(Id idim, double scale)
{
  _cor.setScaleDim(idim, scale);
}

double CorGaussianMixture::computeKappa(double kappai, double kappaj) const
{
  double ci2    = kappai * kappai;
  double cj2    = kappaj * kappaj;
  double scale  = 1.0;
  if (_corRef->getType() == ECov::MATERN)
  {
    scale = sqrt(0.5 * (ci2 + cj2));
  }
  else if (_corRef->getType() == ECov::CAUCHY)
  {
    scale = 1 / sqrt(0.5 * (1 / ci2 + 1 / cj2));
  }
  return scale;
}

SpectrumRN CorGaussianMixture::simulateSpectrumRN(Id ns, const ACov* cov0) const
{
  DECLARE_UNUSED(cov0)
  Id nvar       = getNVar();
  Id ndim       = getNDim();
  double nu0    = _Nu.getDiagonal().minimum();
  double kappa0 = _Kappa.getDiagonal().maximum();
  MatrixDense omega(ns, ndim);
  MatrixDense gamma(ns, nvar);
  MatrixDense omega0(ns, ndim);
  VectorDouble xi0(ns);
  double ldf_xi  = 0.0;
  double ldf0_xi = 0.0;
  for (Id ib = 0; ib < ns; ib++)
  {
    double ll = -2.0 * log(law_uniform()) / ns;
    // simulation of random scale xi and the spatial frequency omega
    double xi = 0.0;
    if (_corRef->getType() == ECov::MATERN)
    {
      xi = 1 / law_gamma(nu0, (kappa0 * kappa0 / 4));
    }
    else if (_corRef->getType() == ECov::CAUCHY)
    {
      xi = law_gamma(nu0, 1 / (kappa0 * kappa0));
    }
    else // GAUSSIAN and nvar == 1
    {
      xi = _Kappa.getValue(0,0) * _Kappa.getValue(0,0);
    }
    for (Id idim = 0; idim < ndim; idim++)
    {
      omega.setValue(ib, idim, sqrt(2 * xi) * law_gaussian());
    }

    // simulation of the variable factors lambda
    if (_corRef->getType() == ECov::MATERN)
    {
      ldf0_xi = law_df_IGamma(xi, nu0, kappa0 * kappa0 / 4, true);
    }
    else if (_corRef->getType() == ECov::CAUCHY)
    {
      ldf0_xi = law_df_gamma(xi, nu0, 1 / (kappa0 * kappa0), true);
    }
    else
    {
      ldf0_xi = 1.0;
    }
    for (Id ivar = 0; ivar < nvar; ivar++)
    {
      double nu    = _Nu.getValue(ivar, ivar);
      double kappa = _Kappa.getValue(ivar, ivar);
      if (_corRef->getType() == ECov::MATERN)
      {
        ldf_xi = law_df_IGamma(xi, nu, kappa * kappa / 4, true);
      }
      else if (_corRef->getType() == ECov::CAUCHY)
      {
        ldf_xi = law_df_gamma(xi, nu, 1 / (kappa * kappa), true);
      }
      else
      {
        ldf_xi = 1.0;
      }
      gamma.setValue(ib, ivar, sqrt(ll * exp(ldf_xi - ldf0_xi)));
    }
    // storing omega0 and xi for Gneiting
    xi0[ib] = xi;
    omega0.setRow(ib, omega.getRow(ib));
  }
  // apply the geometrical anisotropy
  const auto& tensor = _corRef->getAniso().getTensorInverse();
  omega.prodMat(&tensor);
  return SpectrumRN(gamma, omega, omega0, xi0);
}

double CorGaussianMixture::_eval(const SpacePoint& p1,
                                 const SpacePoint& p2,
                                 Id ivar,
                                 Id jvar,
                                 const CovCalcMode* mode) const
{
  VectorDouble scales = _corRef->getScales();
  VectorDouble angles = _corRef->getAnisoAngles();
  double correcScale  = _Kappa.getValue(ivar, jvar) * _scaleGneiting;
  for (size_t idim = 0; idim < getSpace()->getNDim(); idim++)
  {
    scales[idim] /= correcScale;
  }
  _cor.setParam(_Nu.getValue(ivar, jvar));
  _cor.setRotationAnglesAndRadius(angles, VectorDouble(), scales);
  return _C0.getValue(ivar, jvar) * _cor.evalCov(p1, p2, ivar, jvar, mode);
}

} // namespace gstlrn