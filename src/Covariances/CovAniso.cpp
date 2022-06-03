/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovFactory.hpp"
#include "Covariances/CovGradientNumerical.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/AException.hpp"
#include "Basic/Vector.hpp"
#include "Basic/Array.hpp"

#include "Space/ASpace.hpp"
#include "geoslib_f.h"
#include "geoslib_f_private.h"

#include <math.h>

static int NWGT[4] = { 2, 3, 4, 5 };
static int NORWGT[4] = { 2, 6, 20, 70 };
static int COVWGT[4][5] = { { 2, -2, 0, 0, 0 },
                            { 6, -8, 2, 0, 0 },
                            { 20, -30, 12, -2, 0 },
                            { 70, -112, 56, -16, 2 } };

CovAniso::CovAniso(const ECov& type, const CovContext& ctxt)
    : ACov(ctxt.getSpace()), /// TODO : shared pointer
      _ctxt(ctxt),
      _cova(CovFactory::createCovFunc(type, ctxt)),
      _sill(),
      _aniso(ctxt.getSpace()->getNDim())
{
  _updateFromContext();
}

CovAniso::CovAniso(const String& symbol, const CovContext& ctxt)
    : ACov(ctxt.getSpace()), /// TODO : shared pointer
      _ctxt(ctxt),
      _cova(),
      _sill(),
      _aniso(ctxt.getSpace()->getNDim())
{
  ECov covtype = CovFactory::identifyCovariance(symbol, ctxt);
  _cova = CovFactory::createCovFunc(covtype, ctxt);
  _updateFromContext();
}

CovAniso::CovAniso(const ECov& type,
                   double range,
                   double param,
                   double sill,
                   const CovContext& ctxt,
                   bool flagRange)
    : ACov(ctxt.getSpace()), /// TODO : shared pointer
      _ctxt(ctxt),
      _cova(CovFactory::createCovFunc(type, ctxt)),
      _sill(),
      _aniso(ctxt.getSpace()->getNDim())
{
  if (ctxt.getNVar() != 1)
  {
    messerr("This entry is dedicated to the Monovariate case");
    messerr("Additional parameters are ignored and stand constructor is used");
  }
  else
  {
    _updateFromContext();
    _sill.setValue(0, 0, sill);
    setParam(param);
    if (flagRange)
      setRange(range);
    else
      setScale(range);
  }
}

CovAniso::CovAniso(const CovAniso &r)
    : ACov(r),
      _ctxt(r._ctxt),
      _cova(CovFactory::duplicateCovFunc(*r._cova)),
      _sill(r._sill),
      _aniso(r._aniso)
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
  }
  return *this;
}

CovAniso::~CovAniso()
{
  delete _cova;
}

void CovAniso::setContext(const CovContext& ctxt)
{
  _ctxt = ctxt;
  _updateFromContext();
}

void CovAniso::setParam(double param)
{
  if (!_cova->hasParam()) return;
  _cova->setParam(param);
}

void CovAniso::setSill(double sill)
{
  if (getNVariables() != 1)
  {
    messerr("Number of provided sill doesn't match number of variables");
    return;
  }
  _sill.reset(1, 1, sill);
}

void CovAniso::setSill(const MatrixSquareGeneral& sill)
{
  if (getNVariables() != sill.getNSize())
  {
    messerr("Number of provided sills doesn't match number of variables");
    return;
  }
  _sill = sill;
}

void CovAniso::setSill(const VectorDouble& sill)
{
  int size = static_cast<int> (sill.size());
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
  if (! _isVariableValid(ivar)) return;
  if (! _isVariableValid(jvar)) return;
  /// TODO : Test if sill matrix is positive definite (if not, generate a warning)
  if (!_sill.isValid(ivar, jvar)) return;
  _sill.setValue(ivar, jvar, sill);
}

void CovAniso::initSill(double value)
{
  _sill.fill(value);
}

void CovAniso::setRange(double range)
{
  if (!hasRange()) return;
  if (range <= EPSILON20)
  {
    messerr("Range is too small (%lf). It has been replaced by 1.",range);
    range = 1;
  }
  double scadef = _cova->getScadef();
  setScale(range / scadef);
}

void CovAniso::setRanges(const VectorDouble& range)
{
  if (!hasRange()) return;
  if (range.size() != getNDim())
  {
    message("Inconsistency on Space Dimension");
    return;
  }
  for (unsigned int i = 0; i < range.size(); i++)
  {
    if (range[i] <= EPSILON20)
    {
      message("The range in Space dimension (%d) should not be too small",i);
    }
  }
  double scadef = _cova->getScadef();
  VectorDouble scale = range;
  ut_vector_divide_inplace(scale, scadef);
  setScales(scale);
}

void CovAniso::setRange(int idim, double range)
{
  if (!hasRange()) return;
  if (range <= EPSILON20)
  {
    message("The range should not be too small");
    return;
  }
  double scadef = _cova->getScadef();
  setScale(idim, range / scadef);
}

void CovAniso::setScale(double scale)
{
  if (!hasRange()) return;
  if (scale <= EPSILON20)
  {
    message("A scale should not be too small");
    return;
  }
  _aniso.setRadius(scale);
  double scadef = _cova->getScadef();
  _cova->setField(scadef * scale);
}

void CovAniso::setScales(const VectorDouble& scale)
{
  if (!hasRange()) return;
  for (unsigned int i = 0; i < scale.size(); i++)
  {
    if (scale[i] <= EPSILON20)
    {
      message("The scale in Space Dimension (%d) should not be too small",i);
      return;
    }
  }
  _aniso.setRadiusVec(scale);
  double scadef = _cova->getScadef();
  _cova->setField(scadef * ut_vector_max(scale));
}

void CovAniso::setScale(int idim, double scale)
{
  if (scale <= EPSILON20)
  {
    message("A scale should not be too small");
    return;
  }
  _aniso.setRadiusDir(idim, scale);
  double scadef = _cova->getScadef();
  _cova->setField(scadef * ut_vector_max(_aniso.getRadius()));
}

void CovAniso::setAnisoRotation(const Rotation& rot)
{
  if (!hasRange()) return;
  _aniso.setRotation(rot);
}

void CovAniso::setAnisoRotation(const VectorDouble& rot)
{
  if (!hasRange()) return;
  int ndim = getNDim();
  if ((int) rot.size() != ndim * ndim)
  {
    message("Dimension of 'rot' (%d) is not compatible with Space Dimension (%d)",
            (int) rot.size(),ndim);
  }
  Rotation r(ndim);
  r.setMatrixDirectByVector(rot);
  _aniso.setRotation(r);
}

void CovAniso::setAnisoAngles(const VectorDouble& angles)
{
  if (!hasRange()) return;
  _aniso.setRotationAngles(angles);
}

void CovAniso::setAnisoAngle(int idim, double angle)
{
  if (!hasRange()) return;
  _aniso.setRotationAngle(idim, angle);
}

bool CovAniso::isConsistent(const ASpace* /*space*/) const
{
  /// TODO : check something in CovAniso::isConsistent?
  return _cova->isConsistent();
}

double CovAniso::eval0(int ivar, int jvar, const CovCalcMode& mode) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;

  // Calculate unit distance by applying anisotropy
  double cov = _cova->evalCov(0);

  // Converting into a variogram is ignored as the result is obvious
  // and this could be not significant most of the time

  // Scale by the sill
  if (! mode.getUnitary())
  {
    double sill;
    if (mode.getEnvelop() != 0)
    {
      double sign = (mode.getEnvelop() > 0) ? 1 : -1;
      double coef = sqrt(getSill(ivar, ivar) * getSill(jvar, jvar));
      sill = sign * coef;
    }
    else
    {
      sill = getSill(ivar, jvar);
    }
    cov *= sill;
  }

  return (cov);
}

// TODO Replace p1 and p2 by SpaceTarget
double CovAniso::eval(int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode& mode) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  double cov = 0;

  // Calculate unit distance by applying anisotropy
  /// TODO: if composite space : return h1, h2, ... (number of sub-space)
  double h = getSpace()->getDistance(p1, p2, _aniso);

   // Calculation of Standardized Covariance

  int norder = mode.getOrderVario();
  if (norder > 0)
  {

    // Calculate High-order Variogram (only valuable when h != 0)

    for (int iwgt = 1; iwgt < NWGT[norder]; iwgt++)
    {
      double hp = h * (1. + iwgt);
      cov += COVWGT[norder][iwgt] * _cova->evalCov(hp);
    }
    cov /= NORWGT[norder];
  }
  else
  {
    cov = _cova->evalCov(h);

    // Convert into a variogram
    if (mode.getAsVario())
    {
      cov = _cova->evalCov(0) - cov;
    }
  }

  // Scale by the sill
  if (!mode.getUnitary())
  {
    double sill;
    if (mode.getEnvelop() != 0)
    {
      double sign = (mode.getEnvelop() > 0) ? 1 : -1;
      double coef = sqrt(getSill(ivar, ivar) * getSill(jvar, jvar));
      sill = sign * coef;
    }
    else
    {
      sill = getSill(ivar, jvar);
    }
    cov *= sill;
  }

  return (cov);
}

double CovAniso::evalCovOnSphere(double alpha, int degree, bool normalize) const
{
  if (! _cova->hasCovOnSphere()) return TEST;
  double radius;
  variety_get_characteristics(&radius);
  double scale = getScale() / radius;
  double sill = getSill(0, 0);

  double cov = _cova->evalCovOnSphere(alpha/radius, scale, degree);
  if (normalize)
  {
    double cov0 = _cova->evalCovOnSphere(0., scale, degree);
    cov /= cov0;
  }
  return sill * cov;
}


void CovAniso::setMarkovCoeffs(VectorDouble coeffs)
{
  _cova->setMarkovCoeffs(coeffs);
}

double CovAniso::evalSpectrum(double freq) const
{

  if (! _cova->hasSpectrum()) return TEST;
  double scale = getScale();
  double sill = getSill(0,0);
  double val = _cova->evaluateSpectrum(freq,scale,getNDim());
  return sill * val;
}


VectorDouble CovAniso::getMarkovCoeffs() const
{
  if( !_cova->hasMarkovCoeffs()) return VectorDouble();

  return _cova->getMarkovCoeffs();
}

VectorDouble CovAniso::evalCovOnSphere(const VectorDouble& alpha, int degree) const
{
  int n = (int) alpha.size();
  VectorDouble vec(n);
  double c0 = evalCovOnSphere(0., degree);
  for (int i = 0; i < n; i++)
    vec[i] = evalCovOnSphere(alpha[i], degree, false) / c0;

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
      sstr << toMatrix("- Sill matrix:",VectorString(),VectorString(),0,
                       getNVariables(),getNVariables(),_sill.getValues());
    }
    else
    {
      sstr << "- Sill         = " << toDouble(getSill(0, 0)) << std::endl;
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
        sstr << toVector("- Angles       = ", getAnisoAngles());
        sstr << toMatrix("- Rotation Matrix",VectorString(),VectorString(),true,
                         getNDim(),getNDim(),getAnisoRotMatVec());
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
          slopes.setValue(ivar, jvar, _sill.getValue(ivar,jvar) / range);
      sstr << toMatrix("- Slope matrix:",VectorString(),VectorString(),0,
                       getNVariables(),getNVariables(),slopes.getValues());
    }
    else
    {
      sstr << "- Slope        = " << toDouble(getSlope(0, 0));
    }

    if (! _aniso.isIsotropic())
    {
      sstr << toVector("- Aniso, Coeff = ", _aniso.getRadius());
      if (!_aniso.getRotation().isIdentity())
      {
        sstr << toVector("- Angles       = ", getAnisoAngles());
        sstr << toMatrix("- Rotation Matrix",VectorString(),VectorString(),true,
                         getNDim(),getNDim(),getAnisoRotMatVec());
      }
    }
  }
  else
  {
    // Only sill is defined

    if (getNVariables() > 1)
     {
       sstr << toMatrix("- Sill matrix:",VectorString(),VectorString(),0,
                        getNVariables(),getNVariables(),_sill.getValues());
     }
     else
     {
       sstr << "- Sill         = " << toDouble(getSill(0, 0)) << std::endl;
     }
  }

  return sstr.str();
}

double CovAniso::getSill(int ivar, int jvar) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  if (!_sill.isValid(ivar, jvar)) return TEST;
  return _sill.getValue(ivar, jvar);
}

/**
 * Return the Slope calculated as the sill / range(idim=0)
 * @param ivar Rank of the first variable
 * @param jvar Rank of the second variable
 * @return
 */
double CovAniso::getSlope(int ivar, int jvar) const
{
  if (! _isVariableValid(ivar)) return TEST;
  if (! _isVariableValid(jvar)) return TEST;
  if (!_sill.isValid(ivar, jvar)) return TEST;
  if (hasRange() == 0) return TEST;
  double range = getRange(0);
  return _sill.getValue(ivar, jvar) / range;
}

VectorDouble CovAniso::getRanges() const
{
  VectorDouble range = getScales();
  double scadef = _cova->getScadef();
  if (!hasRange()) scadef = 0.;
  ut_vector_multiply_inplace(range, scadef);
  return range;
}

void CovAniso::setType(const ECov& type)
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
  else
    return ut_vector_max(getRanges());
}

double CovAniso::getScale() const
{
  if (!hasRange()) return 0.;
  if (isIsotropic())
    return getScale(0);
  else
    return ut_vector_max(getScales());
}

const VectorDouble CovAniso::getAnisoCoeffs() const
{
  VectorDouble coef = getRanges();
  double max = ut_vector_max(coef);
  if (max < EPSILON10)
  {
    message("Range is null");
    return VectorDouble();
  }
  ut_vector_divide_inplace(coef, max);
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
  else
    return _cova->getParam();
}

void CovAniso::_updateFromContext()
{
  int ndim = getNDim();
  int nvar = getNVariables();
  _sill.reset(nvar, nvar, 1.);
  _aniso.init(ndim);
}

/**
 * Calculate the Integral Range in various Space Dimension (1, 2 or 3)
 * @return
 */
double CovAniso::getIntegralRange(int ndisc, double hmax) const
{
  int ndim = getNDim();
  VectorDouble dd(ndim);
  double delta = hmax / ndisc;
  double total = 0.;
  switch (ndim)
  {
    case 1:
      for (int j1 = -ndisc; j1 <= ndisc; j1++)
      {
        dd[0] = delta * j1;
        total += delta * eval(0,0,dd,SpacePoint());
      }
      break;

    case 2:
      for (int j1 = -ndisc; j1 <= ndisc; j1++)
        for (int j2 = -ndisc; j2 <= ndisc; j2++)
        {
          dd[0] = delta * j1;
          dd[1] = delta * j2;
          total += delta * delta * eval(0,0,dd,SpacePoint());
        }
      break;

    case 3:
      for (int j1 = -ndisc; j1 <= ndisc; j1++)
        for (int j2 = -ndisc; j2 <= ndisc; j2++)
          for (int j3 = -ndisc; j3 <= ndisc; j3++)
          {
            dd[0] = delta * j1;
            dd[1] = delta * j2;
            dd[2] = delta * j3;
            total += delta * delta * delta * eval(0,0,dd,SpacePoint());
          }
      break;

    default:
      my_throw("Integral Range has been programmed for Space Dimension 1 to 3");
  }
  return total;
}

bool CovAniso::_isVariableValid(int ivar) const
{
  if (ivar < 0 || ivar >= getNVariables())
  {
    mesArg("Rank of the Variable",1,getNVariables());
    return false;
  }
  return true;
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
      number ++;
    else
      number += ndim;
  }
  return number;
}

double CovAniso::scale2range(const ECov &type, double scale, double param)
{
  CovContext ctxt = CovContext(1, 1);
  ACovFunc *cova = CovFactory::createCovFunc(type, ctxt);
  cova->setParam(param);
  double factor = cova->getScadef();
  return scale * factor;
}

double CovAniso::range2scale(const ECov &type, double range, double param)
{
  CovContext ctxt = CovContext(1, 1);
  ACovFunc *cova = CovFactory::createCovFunc(type, ctxt);
  cova->setParam(param);
  double factor = cova->getScadef();
  return range / factor;
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
  return new CovAniso(type, range, param, sill, ctxt,flagRange);
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
  int ndim = (int) ranges.size();
  if ((int) ctxt.getNDim() != ndim)
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
  if (! angles.empty()) cov->setAnisoAngles(angles);
  return cov;
}

CovAniso* CovAniso::createIsotropicMulti(const CovContext& ctxt,
                                         const ECov& type,
                                         double range,
                                         const MatrixSquareGeneral& sills,
                                         double param,
                                         bool flagRange)
{
  CovAniso* cov = new CovAniso(type, ctxt);
  int nvar = sills.getNSize();
  if (ctxt.getNVar() != nvar)
  {
    messerr("Mismatch in the number of variables between 'sills'(%d) and 'ctxt'(%d)",
            nvar,ctxt.getNVar());
    return nullptr;
  }
  if (flagRange)
    cov->setRange(range);
  else
    cov->setScale(range);
  cov->setSill(sills);
  cov->setParam(param);
  return cov;
}

CovAniso* CovAniso::createAnisotropicMulti(const CovContext& ctxt,
                                           const ECov& type,
                                           const VectorDouble& ranges,
                                           const MatrixSquareGeneral& sills,
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
  int ndim = (int) ranges.size();
  if ((int) ctxt.getNDim() != ndim)
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
  if (! angles.empty()) cov->setAnisoAngles(angles);
  return cov;
}

void CovAniso::copyCovContext(const CovContext& ctxt)
{
  _ctxt.copyCovContext(ctxt);
  if (_cova != nullptr) _cova->copyCovContext(ctxt);
}

Array CovAniso::evalSpectrum(const VectorDouble& /*ext*/, int N) const
{
  int ndim = getNDim();
  VectorInt nxs(ndim);
  for (int idim = 0; idim < ndim; idim++)
    nxs[idim] = N;
  Array array(nxs);

  int ntotal = pow(N, ndim);
  for (int iad = 0; iad < ntotal; iad++)
  {
    VectorInt indices = array.RankToIndice(iad);
    double s = 0.;
    for (auto &e: indices)
      s += e * e;
    double res = evalSpectrum(s);
    array.setValue(indices, res);
  }
//  for (int idim = 0; idim < ndim; idim++)
//    for (int jdim = 0; jdim < ndim; jdim ++ )
//  {
//    for (int ix = -nx)
//  }
//  d = 2
//  a= np.pi * (N-1) * dx[0] #3 x la portée
//  ind = np.arange(0,int(N/2),2)
//  v = np.linspace(-1.,1.,N)
//  u = a/2 * v
//  deltau=a/(N-1)
//  normxi = np.array([i**2 + j**2 for i in u for j in u])#.reshape((len(u),len(u)))
//
//  fourier = 1./(2*np.pi)**d * np.array([model.getCova(0).evalSpectrum(i) for i in normxi.reshape(-1)])
//  # Préparation pour gstlearn
//  # On met la matrice dans un vector double
//  f = fourier.reshape(-1)
//  ve = gl.VectorDouble(f.shape[0])
//  for i in range(ve.size()):
//      ve[i]=f[i]
//  im = gl.VectorDouble()
//  gl.FFTn(2,[N,N],ve,im)
//
//  # Reformatage des sorties
//  w=np.array(ve).reshape((len(u),len(u)))
//
//  A = w[0,:][ind]
//  covZ=A*deltau**d
//  X= np.pi * (v[ind]-v[0]) /deltau
//
  return array;
}
