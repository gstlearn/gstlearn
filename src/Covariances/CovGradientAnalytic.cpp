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
#include "Covariances/CovGradientAnalytic.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovContext.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceComposite.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"

#define TR(i, j) (_Tr[3 * (i) + (j)])

namespace gstlrn
{
CovGradientAnalytic::CovGradientAnalytic(const CovAniso& cova)
  : CovGradientGeneric(cova)
  , _flagCalculateGG(true)
  , _p1()
  , _p2()
  , _covpp(TEST)
  , _covGp()
  , _covGG()
{
  setContext(cova.getContext());
  if (!_isValid()) return;
  _ctxt.setNVar(1);

  // Set the coordinate of memory space point to unprobable values
  _p1.setCoord(0, TEST);
  _p2.setCoord(0, TEST);

  // Initialize the working quantities
  _covpp = TEST;
  _covGp.resize(3);
  _covGG.resize(9);
  _dF.resize(3);
  _uF.resize(3);
  _hF.resize(3);
  _Tr.resize(9);
  _trttr.resize(9);
}

CovGradientAnalytic::CovGradientAnalytic(const CovGradientAnalytic& r)
  : CovGradientGeneric(r)
  , _flagCalculateGG(r._flagCalculateGG)
  , _p1(r._p1)
  , _p2(r._p2)
  , _covpp(r._covpp)
  , _covGp(r._covGp)
  , _covGG(r._covGG)
  , _dF(r._dF)
  , _uF(r._uF)
  , _hF(r._hF)
  , _Tr(r._Tr)
  , _trttr(r._trttr)
{
}

CovGradientAnalytic::~CovGradientAnalytic()
{
}

bool CovGradientAnalytic::_isValid() const
{
  auto ndim = getCovRef().getNDim();
  auto nvar = getCovRef().getNVar();
  if (ndim > 3)
  {
    messerr("The Gradient class (Analytic) is limited to 3-D");
    return false;
  }
  if (nvar != 1)
  {
    messerr("This class is limited to Monovariate case");
    return false;
  }
  return true;
}

void CovGradientAnalytic::_optimizationSetTarget(SpacePoint& pt) const
{
  DECLARE_UNUSED(pt)
}

/**
 * @brief According to the variable rank, call covariance between the variable and its derivatives
 *
 * @param p1 First point for covariance calculation
 * @param p2 Second point for covariance calculation
 * @param ivar Rank of the first variable (see remarks)
 * @param jvar Rank for the second variable (see remarks)
 * @param mode CovCalcMode structure
 * @return double
 *
 * @remark The use of this function is limited to the Monovariate case: therefore
 * @remark the arguments 'ivar' and 'jvar' cannot reference the variable rank.
 * @remark Instead the following convention applies for 'ivar' and 'jvar':
 * @remark - when =0, it refers to the variable itself
 * @remark - when >10, it applies to the gradient component (10 alongX, 11 along Y, ...)
 * @remark Note that if the space dimension is exceeded (it happens voluntarily in the code),
 * nothing particular must be done: the returned value will not be considered.
 */
double CovGradientAnalytic::_eval(const SpacePoint& p1,
                                  const SpacePoint& p2,
                                  Id ivar,
                                  Id jvar,
                                  const CovCalcMode* mode) const
{
  DECLARE_UNUSED(mode);
  if (_onePointHasChanged(p1, p2))
    _evalZAndGradients(p1, p2);

  Id idim = ivar - 1;
  Id jdim = jvar - 1;

  if (ivar == 0)
  {
    if (jvar == 0)
      return _covpp;

    return -_covGp[jdim];
  }
  if (jvar == 0)
    return _covGp[idim];
  return _covGG[3 * idim + jdim];
}

bool CovGradientAnalytic::_onePointHasChanged(const SpacePoint& p1, const SpacePoint& p2) const
{
  bool flagNew = false;
  for (Id idim = 0, ndim = getNDim(); idim < ndim; idim++)
  {
    if (p1.getCoord(idim) != _p1.getCoord(idim))
    {
      flagNew = true;
      break;
    }
    if (p2.getCoord(idim) != _p2.getCoord(idim))
    {
      flagNew = true;
      break;
    }
  }

  if (!flagNew) return false;
  _p1 = p1;
  _p2 = p2;
  return true;
}

String CovGradientAnalytic::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;

  sstr << "Covariance for Potential Model" << std::endl;
  return sstr.str();
}

const CovAniso* CovGradientAnalytic::_getCovRefAniso() const
{
  const auto* cova = dynamic_cast<const CovAniso*>(&getCovRef());
  return cova;
}

/**
 * Calculate the square of the transformation matrix which transforms
 * a vector into its isotropic equivalent
 */
void CovGradientAnalytic::_calculateTrTtr() const
{
  auto ndim = getNDim();
  _uF.fill(0.);
  _hF.fill(0.);
  _Tr.fill(0.);
  _trttr.fill(0.);

  // Matrix Tr = diag(coeffs) . R
  const MatrixSquare& mat    = _getCovRefAniso()->getAniso().getRotation().getMatrixDirect();
  const VectorDouble& scales = _getCovRefAniso()->getCorAniso()->getScales();

  for (Id i = 0; i < 3; i++)
    for (Id j = 0; j < 3; j++)
    {
      if (i >= ndim || j >= ndim) continue;
      TR(i, j) = mat.getValue(i, j) / scales[i];
    }

  // Calculate the t(Tr) %*% Tr matrix

  Id ecr = 0;
  for (Id i = 0; i < 3; i++)
    for (Id j = 0; j < 3; j++)
    {
      double prod = 0.;
      for (Id k = 0; k < 3; k++) prod += TR(k, i) * TR(k, j);
      _trttr[ecr++] = prod;
    }

  // Rotate the Vector

  _hF[0] = _Tr[0] * _dF[0] + _Tr[1] * _dF[1] + _Tr[2] * _dF[2];
  _hF[1] = _Tr[3] * _dF[0] + _Tr[4] * _dF[1] + _Tr[5] * _dF[2];
  _hF[2] = _Tr[6] * _dF[0] + _Tr[7] * _dF[1] + _Tr[8] * _dF[2];

  _uF[0] = _Tr[0] * _hF[0] + _Tr[3] * _hF[1] + _Tr[6] * _hF[2];
  _uF[1] = _Tr[1] * _hF[0] + _Tr[4] * _hF[1] + _Tr[7] * _hF[2];
  _uF[2] = _Tr[2] * _hF[0] + _Tr[5] * _hF[1] + _Tr[8] * _hF[2];
}

/**
 * Evaluates the covariance and gradient components
 * This function is restricted to the monovariate case
 * This function is limited to the only functions for which the
 * covariance(Point-Gradient) and covariance(Gradient-Gradient)
 * have been coded: i.e. Cubic or Gaussian, in addition to the nugget effect.
 *
 * If 'flag_grad' == 0, then the output array 'covGG' is not filled.
 *
 * @param p1     First point of the Increment
 * @param p2     Second point of the increment
 */

void CovGradientAnalytic::_evalZAndGradients(const SpacePoint& p1,
                                             const SpacePoint& p2) const
{
  _covpp = 0.;
  _covGp.fill(0);
  _covGG.fill(0);

  // Calculate the isotropic distance
  double hh = getSpace()->getDistance(p1, p2, _getCovRefAniso()->getAniso());

  //  Calculate the covariance
  const CovAniso* covref  = _getCovRefAniso();
  const ACovFunc* covfunc = covref->getCorFunc();
  double covar            = covref->getSill(0, 0) * covfunc->evalCorFunc(hh);
  _covpp += covar;
  if (covfunc->getType() == ECov::NUGGET) return;

  // Calculate distance and plunge into a 3-D vector
  VectorDouble d1 = VH::subtract(p1.getCoords(), p2.getCoords());
  for (Id i = 0; i < 3; i++)
    _dF[i] = (i < static_cast<Id>(d1.size())) ? d1[i] : 0.;

  _calculateTrTtr();
  double dcovsr = covref->getSill(0, 0) * covfunc->evalCovFirstDerivativeOverH(hh);

  //  Case where distance is null
  if (hh < EPSGRAD)
  {
    if (_flagCalculateGG)
    {
      for (Id ecr = 0; ecr < 9; ecr++)
        _covGG[ecr] -= dcovsr * _trttr[ecr];
    }
  }
  else
  {
    //  Calculate covariance between point and gradient
    for (Id ecr = 0; ecr < 3; ecr++)
      _covGp[ecr] += _uF[ecr] * dcovsr;

    //  Calculate the covariance between gradient and gradient
    if (_flagCalculateGG)
    {
      double d2cov = covref->getSill(0, 0) * covfunc->evalCovDerivative(2, hh);
      double a     = (dcovsr - d2cov) / (hh * hh);

      if (covref->isIsotropic())
      {

        //  Isotropic case

        double b = dcovsr * _trttr[0];
        for (Id i = 0, ecr = 0; i < 3; i++)
          for (Id j = 0; j < 3; j++, ecr++)
          {
            _covGG[ecr] += a * _uF[i] * _uF[j];
            if (i == j) _covGG[ecr] -= b;
          }
      }
      else
      {

        //  Anisotropic case

        for (Id i = 0, ecr = 0; i < 3; i++)
          for (Id j = 0; j < 3; j++, ecr++)
          {
            _covGG[ecr] += a * _uF[i] * _uF[j] - dcovsr * _trttr[ecr];
          }
      }
    }
  }
}

} // namespace gstlrn